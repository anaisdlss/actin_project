#!/usr/bin/env python
"""
Etape 2 — Extraire la chaine ABP seule de chaque complexe representatif.

Lit data/exports/abp_site_domain/abp_representatives.csv, recupere le fichier
assembly {iid}_{pdb}.pdb|.cif, isole la chaine ABP -> un .pdb par ABP pour foldseek.

Sortie : data/exports/abp_site_domain/abp_chains/<slug>.pdb
"""
import re
from pathlib import Path
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/exports/abp_site_domain"
ASM = ROOT / "data/filtered/details/structures_files/assembly"
CHAINS = OUT / "abp_chains"
CHAINS.mkdir(parents=True, exist_ok=True)

rep = pd.read_csv(OUT / "abp_representatives.csv")


def slug(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(s)).strip("_")[:50]


class ChainSel(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        # garder seulement les acides amines (hetero flag vide)
        return residue.id[0] == " "


io = PDBIO()
written, errors = [], []
for _, r in rep.iterrows():
    base = f"{r.pdb}"
    f = ASM / (base + ".pdb")
    parser = PDBParser(QUIET=True)
    if not f.exists():
        f = ASM / (base + ".cif")
        parser = MMCIFParser(QUIET=True)
    try:
        struct = parser.get_structure(base, str(f))
        model = next(iter(struct))
        if r.abp_chain not in [c.id for c in model]:
            errors.append((r.abp_title, base, f"chaine {r.abp_chain} absente: {[c.id for c in model]}"))
            continue
        name = f"{slug(r.abp_title)}__{r.pdb}_{r.abp_chain}"
        io.set_structure(struct)
        io.save(str(CHAINS / (name + ".pdb")), ChainSel(r.abp_chain))
        written.append(name)
    except Exception as e:
        errors.append((r.abp_title, base, str(e)))

print(f"chaines ecrites : {len(written)} dans {CHAINS}")
if errors:
    print(f"\nerreurs ({len(errors)}):")
    for t, b, e in errors:
        print(f"  {t} [{b}] : {e}")
