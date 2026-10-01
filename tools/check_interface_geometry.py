"""Representative dimer-geometry check: F-actin/tropomyosin vs cofilactin.
Fits one actin subunit, measures displacement of the second without refitting it.
This is an illustrative structural check, not validation of every deposited assembly.
"""
import sys,os,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1_extended
from Bio.SVDSuperimposer import SVDSuperimposer
ROOT=Path(__file__).resolve().parents[1];os.chdir(ROOT);sys.path.insert(0,str(ROOT/'script'))
from scientific_analysis import sequence,isoform_map
DATA=ROOT/'data/filtered/details/structures_files/assembly'
OUT=ROOT/'reports/scientific_audit';OUT.mkdir(parents=True,exist_ok=True)
ref=sequence(ROOT/'data/P60709_ref.fasta');coords={}
for pdb in ['3j8a','5yu8','6vao']:
    model=PDBParser(QUIET=True).get_structure(pdb,DATA/f'{pdb}.pdb')[0]
    for chain in model:
        residues=[r for r in chain if 'CA' in r and r.resname in protein_letters_3to1_extended]
        seq=''.join(protein_letters_3to1_extended[r.resname] for r in residues)
        if len(seq)<300:continue
        mapping=isoform_map(seq,ref)
        coords[f'{pdb}_{chain.id}']={p:residues[i-1]['CA'].coord.astype(float) for i,p in mapping.items()}
x=pd.read_csv('data/filtered/filtered_all_data.csv',low_memory=False)
x=x[x.pdb_id.isin(['3j8a','5yu8','6vao']) & x.s1_actine & x.s2_actine].drop_duplicates(['subunit_1','subunit_2'])
templates=x[x.pdb_id.isin(['3j8a','5yu8'])].drop_duplicates(['pdb_id','cluster_data_70'])
rows=[]
for r in x.itertuples():
    for t in templates.itertuples():
        best=None
        for flip in [False,True]:
            ta,tb=(t.subunit_2,t.subunit_1) if flip else (t.subunit_1,t.subunit_2)
            ap=sorted(set(coords[r.subunit_1])&set(coords[ta]));bp=sorted(set(coords[r.subunit_2])&set(coords[tb]))
            if len(ap)<300 or len(bp)<300:raise ValueError('Insufficient mapped C-alpha coverage')
            target=np.array([coords[ta][p] for p in ap]);mobile=np.array([coords[r.subunit_1][p] for p in ap])
            fit=SVDSuperimposer();fit.set(target,mobile);fit.run();rot,trans=fit.get_rotran()
            expected=np.array([coords[tb][p] for p in bp]);neighbor=np.array([coords[r.subunit_2][p] for p in bp])@rot+trans
            rms=float(np.sqrt(np.mean(np.sum((neighbor-expected)**2,axis=1))))
            result=dict(pdb_id=r.pdb_id,chain_A=r.subunit_1,chain_B=r.subunit_2,c70=r.cluster_data_70,
                template_pdb=t.pdb_id,template_c70=t.cluster_data_70,template_A=ta,template_B=tb,
                n_anchor_CA=len(ap),n_neighbor_CA=len(bp),anchor_CA_RMSD_A=fit.get_rms(),neighbor_CA_RMSD_A=rms)
            if best is None or rms<best['neighbor_CA_RMSD_A']:best=result
        rows.append(best)
result=pd.DataFrame(rows);result.to_csv(OUT/'representative_geometry.csv',index=False)
summary=result.groupby(['pdb_id','c70','template_pdb','template_c70']).agg(n_pairs=('neighbor_CA_RMSD_A','size'),median_neighbor_RMSD_A=('neighbor_CA_RMSD_A','median'),median_anchor_RMSD_A=('anchor_CA_RMSD_A','median')).reset_index()
summary.to_csv(OUT/'representative_geometry_summary.csv',index=False)
manifest={'reference_sequence':'P60709','fit':'Global sequence alignment (BLOSUM62, gap open -10, extend -0.5), consensus across optimal alignments; C-alpha Kabsch fit of one actin. Neighbor RMSD under SAME transform. Best of two chain assignments per reference pair.',
 'scope':'Representative structures only; does not establish steric clashes or classification of all interfaces.',
 'sources':{pdb:{'url':f'https://www.rcsb.org/structure/{pdb.upper()}','sha256':hashlib.sha256((DATA/f'{pdb}.pdb').read_bytes()).hexdigest()} for pdb in ['3j8a','5yu8','6vao']},
 'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(OUT/'representative_geometry_manifest.json').write_text(json.dumps(manifest,indent=2))
print(summary.to_string(index=False))
