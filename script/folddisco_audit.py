"""Local FoldDisco provenance, query preparation and coverage calculations."""
from pathlib import Path
import hashlib
import json
import re

import numpy as np
import pandas as pd

EXPORTS = Path('data/exports/abp_site_domain')
RESIDUES = Path('data/filtered/details/3.interface_residues.csv')
CATALOG_FILES = [EXPORTS/'abp_site_table.csv', RESIDUES]


def residue_token(value):
    """Preserve signed PDB positions and insertion codes; never truncate fractions."""
    text=str(value).strip()
    if re.fullmatch(r'-?\d+\.0+',text):text=str(int(float(text)))
    match=re.fullmatch(r'(-?\d+)([A-Za-z]?)',text)
    return f'{int(match[1])}{match[2]}' if match else None


def residue_order(token):
    match=re.fullmatch(r'(-?\d+)([A-Za-z]?)',token)
    return int(match[1]),match[2]


def motif_catalog(site_table,residues):
    """Reconstruct script14's best-resolution choice; original submissions are absent."""
    columns=['query_abp','query_cluster','source_pdb','source_chain','interaction_id',
             'resolution_A','contact_positions','reconstructed_size','unparsed_contact_positions']
    if site_table.empty:return pd.DataFrame(columns=columns)
    table=site_table.copy()
    table['resolution_A']=pd.to_numeric(table.resolution.astype(str).str.extract(r'(\d+(?:\.\d+)?)',expand=False),errors='coerce')
    contacts={key:values.residue_number_structure.tolist() for key,values in residues.groupby(['interaction_id','chain'])}
    rows=[]
    for (abp,cluster),group in table.groupby(['abp_title','actin_site_cluster']):
        representative=group.sort_values('resolution_A',kind='stable',na_position='last').iloc[0]
        chain=str(representative.abp_chain);pdb=str(representative.pdb).lower()
        source=contacts.get((int(representative.interaction_id),f'{pdb}_{chain}'),[])
        positions=sorted({token for value in source if (token:=residue_token(value)) is not None},key=residue_order)
        invalid=sorted({str(value) for value in source if residue_token(value) is None})
        rows.append(dict(query_abp=abp,query_cluster=cluster,source_pdb=pdb,source_chain=chain,
                         interaction_id=int(representative.interaction_id),resolution_A=representative.resolution_A,
                         contact_positions=', '.join(positions),reconstructed_size=len(positions),
                         unparsed_contact_positions='; '.join(invalid)))
    return pd.DataFrame(rows,columns=columns)


def structure_path(row,root=Path('.')):
    base=root/EXPORTS;pdb=row['source_pdb'];chain=row['source_chain']
    candidates=[base/'abp_chains_disco'/f'{pdb}_{chain}.pdb']
    candidates+=sorted((base/'abp_chains').glob(f'*__{pdb}_{chain}.pdb'))
    candidates+=[root/'data/filtered/details/structures_files/assembly'/f'{pdb}{ext}' for ext in ['.pdb','.cif']]
    return next((p for p in candidates if p.is_file()),None)


def read_chain_ca(path,chain):
    """Read resolved C-alpha atoms only, preserving chain case and insertion codes."""
    from Bio.PDB import PDBParser,MMCIFParser
    from Bio.PDB.Polypeptide import protein_letters_3to1_extended
    parser=MMCIFParser(QUIET=True) if Path(path).suffix.lower()=='.cif' else PDBParser(QUIET=True)
    structure=parser.get_structure('query',str(path));model=next(iter(structure))
    if chain not in model:raise ValueError(f'Chain {chain!r} is absent from the local structure.')
    rows=[]
    for residue in model[chain]:
        if 'CA' not in residue or residue.resname not in protein_letters_3to1_extended:continue
        x,y,z=residue['CA'].coord
        rows.append(dict(position=f'{residue.id[1]}{residue.id[2].strip()}',residue=residue.resname,
                         x=float(x),y=float(y),z=float(z)))
    return pd.DataFrame(rows,columns=['position','residue','x','y','z'])


def query_chain(chain):
    """FoldDisco's parser only recognizes ASCII letters as chain prefixes.

    A single extracted numeric/symbol chain is exported as A. Residue numbers
    and coordinates are unchanged, and the original chain is kept in provenance.
    """
    return chain if re.fullmatch(r'[A-Za-z]', str(chain)) else 'A'


def validate_motif(text,available,chain):
    """Validate local PDB positions without claiming external FoldDisco execution."""
    parts=[part for part in re.split(r'[,;\s]+',str(text).strip()) if part]
    parsed=[residue_token(part) for part in parts]
    errors=[]
    invalid=[part for part,token in zip(parts,parsed) if token is None]
    if invalid:errors.append('Invalid position syntax: '+', '.join(invalid))
    positions=sorted({token for token in parsed if token is not None},key=residue_order)
    absent=sorted(set(positions)-set(available),key=residue_order)
    if absent:errors.append('No resolved C-alpha atom in the selected chain: '+', '.join(absent))
    if len(positions)<3:errors.append('At least three distinct resolved residues are required to prepare a motif.')
    compatible=all(re.fullmatch(r'\d+',p) for p in positions)
    # Negative numbers and insertion codes are not supported by the upstream
    # unsigned-integer/range parser. Never silently strip or renumber them.
    export_chain=query_chain(chain)
    query=','.join(f'{export_chain}{p}' for p in positions) if compatible and not errors else None
    return {'positions':positions,'errors':errors,'query':query,
            'submitted_chain':export_chain,
            'text_export_supported':compatible and not errors,'duplicate_tokens_removed':len(parsed)-len(set(parsed))}


def normalize_discovery(frame):
    """Expose the legacy source-label ratio and its distinct best-hit fallback.

    Script14 labels any hit from the same PDB as a source, regardless of chain;
    this is not a verified same-structure self match or a score bounded by one.
    """
    df=frame.copy();keys=['query_abp','query_cluster','db']
    source=df.is_source.astype(str).str.lower().eq('true')
    df['is_source']=source
    score=pd.to_numeric(df.idfscore,errors='coerce')
    df['idfscore']=score.where(np.isfinite(score))
    positive=df[df.idfscore.gt(0)]
    self_score=positive[source.loc[positive.index]].groupby(keys).idfscore.max()
    best=positive.groupby(keys).idfscore.max()
    denominator=self_score.reindex(best.index).fillna(best).rename('normalization_denominator')
    basis=pd.Series('best available hit',index=best.index)
    basis.loc[basis.index.isin(self_score.index)]='source-labelled hit (not verified self-match)'
    df=df.merge(pd.concat([denominator,basis.rename('normalization_basis')],axis=1),on=keys,how='left',validate='many_to_one')
    df['normalization_basis']=df.normalization_basis.fillna('unavailable (no positive finite score)')
    df['score_norm']=df.idfscore/df.normalization_denominator
    return df


def coverage_tables(catalog,discovery,local_pairs,root=Path('.')):
    """Separate expected motifs, saved result rows and available local structures."""
    keys=['query_abp','query_cluster']
    counts=discovery.groupby(keys+['db']).size().unstack(fill_value=0) if not discovery.empty else pd.DataFrame()
    if not counts.empty:counts=counts.rename(columns={'pdb':'PDB_result_rows','afdb':'AlphaFold_result_rows'}).reset_index()
    else:counts=pd.DataFrame(columns=keys+['PDB_result_rows','AlphaFold_result_rows'])
    size=discovery.groupby(keys).query_size.agg(lambda values:'; '.join(map(str,sorted(set(values.dropna()))))).rename('saved_query_sizes').reset_index() if not discovery.empty else pd.DataFrame(columns=keys+['saved_query_sizes'])
    detail=catalog.merge(counts,on=keys,how='outer').merge(size,on=keys,how='left')
    for column in ['PDB_result_rows','AlphaFold_result_rows']:
        if column not in detail:detail[column]=0
        detail[column]=detail[column].fillna(0).astype(int)
    detail['current_motif_available']=detail.contact_positions.fillna('').ne('')
    detail['local_structure_available']=[structure_path(row,root) is not None if row.get('current_motif_available') else False for row in detail.to_dict('records')]
    detail['saved_results_present']=detail[['PDB_result_rows','AlphaFold_result_rows']].sum(axis=1).gt(0)
    def compare_size(row):
        if not row.saved_results_present:return 'No saved result rows'
        if not row.current_motif_available:return 'No reconstructed motif'
        sizes={float(value) for value in str(row.saved_query_sizes).split('; ') if value and value!='nan'}
        return 'Same size only; positions unverified' if sizes=={float(row.reconstructed_size)} else 'Different query size; historical motif not recovered'
    detail['query_size_check']=detail.apply(compare_size,axis=1)
    detail['query_provenance']='Exact historical query positions and structure hash were not stored'
    local=set(local_pairs.query_abp) if not local_pairs.empty else set()
    rows=[]
    for name,group in detail.groupby('query_abp'):
        rows.append(dict(ABP=name,current_site_motifs=int(group.current_motif_available.sum()),
                         motifs_with_local_structure=int(group.local_structure_available.sum()),
                         sites_with_saved_PDB_hits=int(group.PDB_result_rows.gt(0).sum()),
                         sites_with_saved_AlphaFold_hits=int(group.AlphaFold_result_rows.gt(0).sum()),
                         motifs_without_saved_results=int((group.current_motif_available&~group.saved_results_present).sum()),
                         local_pairwise_results_present=name in local))
    return pd.DataFrame(rows),detail


def prepared_query_json(row,validation,path):
    return json.dumps(dict(status='Prepared locally; not submitted or recalculated',
        query_abp=row['query_abp'],query_cluster=row['query_cluster'],source_pdb=row['source_pdb'],
        source_chain=row['source_chain'],interaction_id=int(row['interaction_id']),
        submitted_chain=validation['submitted_chain'],
        selected_sites=row.get('clusters',[row['query_cluster']]),
        interaction_ids=row.get('interaction_ids',[int(row['interaction_id'])]),
        numbering='PDB author residue identifiers in the selected chain; not P60709',
        original_reconstructed_positions=row['contact_positions'],edited_positions=validation['positions'],
        folddisco_query=validation['query'],structure_file=Path(path).name,
        structure_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        historical_results='Not recalculated; historical exact query provenance unavailable'),indent=2)
