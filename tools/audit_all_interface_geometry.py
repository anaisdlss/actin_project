"""Audit all locally retained homo chain pairs against two documented reference ensembles."""
from datetime import datetime,timezone
import hashlib,json,os,sys
from pathlib import Path
import pandas as pd
from Bio.PDB import PDBParser,MMCIFParser
from Bio.PDB.Polypeptide import protein_letters_3to1_extended
ROOT=Path(__file__).resolve().parents[1];os.chdir(ROOT)
sys.path.insert(0,str(ROOT/'script'))
from scientific_analysis import sequence,isoform_map
from interface_geometry import compare_pair


def main():
    reference=sequence(ROOT/'data/P60709_ref.fasta')
    data=pd.read_csv('data/filtered/filtered_all_data.csv',low_memory=False)
    homo=data[data.s1_actine & data.s2_actine].drop_duplicates(['subunit_1','subunit_2','cluster_data_70']).copy()
    homo['resolution_numeric']=pd.to_numeric(homo.resolution.astype(str).str.extract(r'(\d+(?:\.\d+)?)',expand=False),errors='coerce')
    templates=homo[homo.pdb_id.isin(['3j8a','5yu8'])].sort_values(['resolution_numeric','subunit_1','subunit_2']).drop_duplicates(['pdb_id','cluster_data_70'])
    source_hashes={}
    def coordinates(pdb,needed):
        base=ROOT/'data/filtered/details/structures_files/assembly'
        file=next((base/f'{pdb}{ext}' for ext in ['.pdb','.cif'] if (base/f'{pdb}{ext}').exists()),None)
        if file is None:return {}
        source_hashes[pdb]=hashlib.sha256(file.read_bytes()).hexdigest()
        parser=MMCIFParser(QUIET=True) if file.suffix=='.cif' else PDBParser(QUIET=True)
        model=next(iter(parser.get_structure(pdb,str(file))))
        out={}
        for chain in model:
            label=f'{pdb}_{chain.id}'
            if label not in needed:continue
            residues=[r for r in chain if 'CA' in r and r.resname in protein_letters_3to1_extended]
            if len(residues)<300:continue
            seq=''.join(protein_letters_3to1_extended[r.resname] for r in residues)
            mapping=isoform_map(seq,reference)
            out[label]={p:residues[i-1]['CA'].coord.astype(float) for i,p in mapping.items()}
        return out
    template_coords={}
    for pdb,group in templates.groupby('pdb_id'):
        template_coords.update(coordinates(pdb,set(group.subunit_1)|set(group.subunit_2)))
    rows=[];coverage=[]
    for pdb,group in homo.groupby('pdb_id'):
        try:coords=coordinates(pdb,set(group.subunit_1)|set(group.subunit_2))
        except (ValueError,KeyError,OSError) as exc:
            coords={};print(pdb,str(exc),flush=True)
        for r in group.itertuples():
            count=0
            for template in templates.itertuples():
                result=compare_pair(coords.get(r.subunit_1,{}),coords.get(r.subunit_2,{}),
                    template_coords.get(template.subunit_1,{}),template_coords.get(template.subunit_2,{}))
                if result is None:continue
                rows.append(dict(pdb_id=pdb,chain_A=r.subunit_1,chain_B=r.subunit_2,c70=r.cluster_data_70,
                    template_pdb=template.pdb_id,template_c70=template.cluster_data_70,
                    template_chain_A=template.subunit_1,template_chain_B=template.subunit_2,**result));count+=1
            coverage.append(dict(pdb_id=pdb,chain_A=r.subunit_1,chain_B=r.subunit_2,c70=r.cluster_data_70,
                template_comparisons=count,state='measured' if count else 'unavailable: missing structure or fewer than 300 unambiguously mapped C-alpha atoms'))
        print(pdb,len(group),'pairs',flush=True)
    out=ROOT/'reports/scientific_audit';out.mkdir(parents=True,exist_ok=True)
    measures=pd.DataFrame(rows);inventory=pd.DataFrame(coverage)
    measures.to_csv(out/'all_interface_geometry.csv',index=False);inventory.to_csv(out/'all_interface_geometry_coverage.csv',index=False)
    # Nearest reference geometry per pair and ensemble; then one median per PDB.
    nearest=measures.sort_values('neighbor_CA_RMSD_A').drop_duplicates(['pdb_id','chain_A','chain_B','c70','template_pdb'])
    pdb=nearest.groupby(['c70','template_pdb','pdb_id']).neighbor_CA_RMSD_A.median().reset_index()
    summary=pdb.groupby(['c70','template_pdb']).neighbor_CA_RMSD_A.agg(PDB_count='size',median_neighbor_RMSD_A='median',minimum_PDB_median_A='min',maximum_PDB_median_A='max').reset_index()
    summary.to_csv(out/'all_interface_geometry_summary.csv',index=False)
    manifest=dict(created_utc=datetime.now(timezone.utc).isoformat(),dataset=ROOT.name,
        expected_pairs=len(inventory),measured_pairs=int(inventory.state.eq('measured').sum()),
        expected_clusters=int(inventory.c70.nunique()),measured_clusters=int(inventory.loc[inventory.state.eq('measured'),'c70'].nunique()),
        references=['3J8A: F-actin with tropomyosin','5YU8: cofilin-decorated actin'],
        method='P60709 sequence mapping: BLOSUM62 global alignment, open -10, extend -0.5; consensus of optimal mappings. At least 300 C-alpha positions in each subunit. Fit anchor; measure neighbour under same transform. Minimum over chain orientations and reference pairs, separately for each reference PDB. Summaries use one median per PDB.',
        limits='Geometry similarity is descriptive, not a major/minor biological assignment, steric-clash test or validation of deposited assemblies. Coverage includes only locally retained PDB chain pairs. Templates themselves remain in the inventory and can yield zero RMSD.',
        structure_sha256=source_hashes,source_table_sha256=hashlib.sha256((ROOT/'data/filtered/filtered_all_data.csv').read_bytes()).hexdigest(),
        code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'script/interface_geometry.py']})
    (out/'all_interface_geometry_manifest.json').write_text(json.dumps(manifest,indent=2));print({k:manifest[k] for k in ['expected_pairs','measured_pairs','expected_clusters','measured_clusters']},flush=True)

if __name__=='__main__':main()
