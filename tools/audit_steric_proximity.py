"""Reproducible descriptive distances after rigid ABP–actin placement on filaments."""
from collections import Counter
from datetime import datetime,timezone
import hashlib,json,os,sys
from pathlib import Path
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser,MMCIFParser
from Bio.PDB.Polypeptide import protein_letters_3to1_extended
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];os.chdir(ROOT);sys.path.insert(0,str(ROOT/'script'))
from scientific_analysis import sequence,isoform_map
from folddisco_audit import CATALOG_FILES,motif_catalog
from steric_screen import rigid_proximity


def main():
    data=pd.read_csv('data/filtered/filtered_all_data.csv',low_memory=False)
    interactions=pd.read_csv('data/filtered/details/1.interactions.csv').set_index('interaction_id')
    catalog=motif_catalog(*(pd.read_csv(p,low_memory=False) for p in CATALOG_FILES))
    refseq=sequence(ROOT/'data/P60709_ref.fasta');cache={};hashes={}
    def structure(pdb):
        if pdb in cache:return cache[pdb]
        base=ROOT/'data/filtered/details/structures_files/assembly'
        file=next((base/f'{pdb}{ext}' for ext in ['.pdb','.cif'] if (base/f'{pdb}{ext}').exists()),None)
        if file is None:raise ValueError('Local assembly unavailable')
        hashes[pdb]=hashlib.sha256(file.read_bytes()).hexdigest()
        parser=MMCIFParser(QUIET=True) if file.suffix=='.cif' else PDBParser(QUIET=True)
        model=next(iter(parser.get_structure(pdb,str(file))));out={}
        for c in model:
            residues=[r for r in c if r.resname in protein_letters_3to1_extended]
            ca=[r for r in residues if 'CA' in r];atoms=[];labels=[]
            for r in residues:
                for a in r:
                    if a.element in {'H','D'}:continue
                    if a.is_disordered() and a.get_altloc() not in {' ','A'}:continue
                    atoms.append(a.coord.astype(float));labels.append(f'{r.id[1]}{r.id[2].strip()}')
            label=f'{pdb}_{c.id}'
            mapping={}
            is_actin=((data.subunit_1.eq(label)&data.s1_actine)|(data.subunit_2.eq(label)&data.s2_actine)).any()
            if is_actin and len(ca)>=300:
                seq=''.join(protein_letters_3to1_extended[r.resname] for r in ca)
                mapping={p:ca[i-1]['CA'].coord.astype(float) for i,p in isoform_map(seq,refseq).items()}
            out[label]=dict(atoms=np.array(atoms),residues=labels,ca=mapping,is_actin=bool(is_actin))
        cache[pdb]=out;return out
    references=[]
    for pdb in ['3j8a','5yu8']:
        sub=data[data.pdb_id.eq(pdb)&data.s1_actine&data.s2_actine]
        edges={tuple(sorted((r.subunit_1,r.subunit_2))) for r in sub.itertuples()}
        degree=Counter(c for edge in edges for c in edge)
        anchor=sorted(degree,key=lambda c:(-degree[c],c))[0]
        chains=structure(pdb)
        others=[v['atoms'] for k,v in chains.items() if k!=anchor and v['is_actin'] and len(v['atoms'])]
        references.append(dict(pdb=pdb,anchor=anchor,degree=degree[anchor],ca=chains[anchor]['ca'],atoms=np.vstack(others),
                               neighbor_chains=[k for k,v in chains.items() if k!=anchor and v['is_actin']]))
    rows=[]
    for row in catalog.to_dict('records'):
        partner=f"{row['source_pdb']}_{row['source_chain']}";interaction=interactions.loc[row['interaction_id']]
        chain_a,chain_b=interaction['chain_A_id'],interaction['chain_B_id']
        anchor=chain_b if chain_a==partner else chain_a if chain_b==partner else None
        for reference in references:
            item=dict(ABP=row['query_abp'],site=row['query_cluster'],interaction_id=row['interaction_id'],
                source_PDB=row['source_pdb'],source_actin=anchor,source_ABP=partner,
                reference_PDB=reference['pdb'].upper(),reference_anchor=reference['anchor'])
            try:
                chains=structure(row['source_pdb']);a=chains[anchor];b=chains[partner]
                native=float(cKDTree(a['atoms']).query(b['atoms'])[0].min())
                item['source_pair_minimum_distance_A']=native
                if native>8:raise ValueError('Source chains are not adjacent in this assembly coordinate frame (>8 Å).')
                result=rigid_proximity(a['ca'],reference['ca'],b['atoms'],b['residues'],reference['atoms'])
                item.update(result,state='measured' if result['anchor_CA_RMSD_A']<=3 else 'poor anchor fit (>3 Å): inspect')
            except (ValueError,KeyError) as exc:item.update(state=f'unavailable: {exc}')
            rows.append(item)
        print(row['query_abp'],row['query_cluster'],flush=True)
    out=ROOT/'reports/scientific_audit';frame=pd.DataFrame(rows);frame.to_csv(out/'rigid_filament_proximity.csv',index=False)
    files=[*CATALOG_FILES,Path('data/filtered/filtered_all_data.csv'),Path('data/filtered/details/1.interactions.csv')]
    manifest=dict(created_utc=datetime.now(timezone.utc).isoformat(),dataset=ROOT.name,
        scope='Best-resolution observed ABP/site representatives; finite reference assemblies; protein heavy atoms only, first model, altloc blank or A.',
        references=[{k:v for k,v in ref.items() if k not in {'ca','atoms'}} for ref in references],
        method='Reference anchor: highest actin–actin contact degree in retained source table; ties lexicographic. P60709 global sequence mapping (BLOSUM62, gap open -10, extend -0.5, consensus optimal mappings). At least 300 shared C-alpha positions. Rigid actin fit; transform partner with the same rotation/translation. Nearest heavy-atom distances to all other reference actins. Count distinct partner atoms and residues below 1.5, 2.0, 2.5 Å. No ASA threshold.',
        limits='Descriptive geometric screen, not a validated steric-clash/competition prediction or an automatic major/minor assignment. No atom-specific van der Waals radii, relaxation, conformational ensembles or infinite filament. Source assembly adjacency checked at 8 Å; fits above 3 Å flagged.',
        states=frame.state.value_counts().to_dict(),source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        structure_sha256=hashes,code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'script/steric_screen.py']})
    (out/'rigid_filament_proximity_manifest.json').write_text(json.dumps(manifest,indent=2))

if __name__=='__main__':main()
