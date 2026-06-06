#!/usr/bin/env python3
"""
Écrit la conservation ProteoCast (moyenne de sensibilité mutationnelle par résidu)
dans le B-factor d'une structure d'actine humaine, + script PyMOL.

B-factor = conservation = -mean(Variant_score) par position (plus grand = plus sensible/conservé).
Numérotation du PDB = position P60709 (vérifié : 100 % identique, resnum 1..375).

Usage : python script/bfactor_conservation.py
Sorties : data/proteocast/actin_human_conservation.pdb + load_conservation.pml
"""
from pathlib import Path
import pandas as pd
from Bio.PDB import PDBParser, PDBIO
import warnings
warnings.filterwarnings('ignore')

ROOT = Path(__file__).resolve().parent.parent
PDB_IN = ROOT / 'actin_human_pdb.pdb'
CSV = ROOT / 'data/proteocast/conservation_vs_asa_per_position.csv'
OUT_PDB = ROOT / 'data/proteocast/actin_human_conservation.pdb'
OUT_PML = ROOT / 'data/proteocast/load_conservation.pml'

df = pd.read_csv(CSV)
cons = dict(zip(df['Residue'].astype(int), df['conservation']))

struct = PDBParser(QUIET=True).get_structure('actin', str(PDB_IN))
n_set, n_miss = 0, 0
for atom in struct.get_atoms():
    res = atom.get_parent()
    if res.id[0] != ' ':
        continue
    v = cons.get(res.id[1])
    if v is None:
        atom.set_bfactor(0.0); n_miss += 1
    else:
        atom.set_bfactor(round(float(v), 2)); n_set += 1

io = PDBIO(); io.set_structure(struct); io.save(str(OUT_PDB))
vmin, vmax = min(cons.values()), max(cons.values())
print(f'B-factor écrit : {OUT_PDB}')
print(f'  atomes renseignés : {n_set} | sans valeur : {n_miss}')
print(f'  conservation min={vmin:.2f}  max={vmax:.2f}')

pml = f"""# Conservation ProteoCast (sensibilité mutationnelle moyenne) en B-factor
load {OUT_PDB.name}, actin
bg_color white
hide everything
show cartoon, actin
# bleu = tolérant (peu conservé) -> rouge = sensible (très conservé)
spectrum b, blue_white_red, actin, minimum={vmin:.2f}, maximum={vmax:.2f}
set cartoon_fancy_helices, 1
set ray_shadows, 0
# surface optionnelle :
# show surface, actin ; set transparency, 0.2
"""
OUT_PML.write_text(pml)
print(f'Script PyMOL : {OUT_PML}')
