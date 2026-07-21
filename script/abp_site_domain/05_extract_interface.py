#!/usr/bin/env python
"""
Etape 6 — Extraire UNIQUEMENT les residus de l'ABP qui touchent l'actine
(empreinte d'interface), au lieu de la proteine entiere.

Variante reglable par PAD : nb de residus de contexte ajoutes de part et d'autre
de chaque residu d'interface (PAD=0 = residus de contact stricts).

Sortie : data/exports/abp_site_domain/abp_interface[_padN]/<slug>.pdb
"""
import re
import sys
from pathlib import Path
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select

PAD = int(sys.argv[1]) if len(sys.argv) > 1 else 0

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/exports/abp_site_domain"
ASM = ROOT / "data/filtered/details/structures_files/assembly"
DST = OUT / (f"abp_interface_pad{PAD}" if PAD else "abp_interface")
DST.mkdir(parents=True, exist_ok=True)

rep = pd.read_csv(OUT / "abp_representatives.csv")
res = pd.read_csv(ROOT / "data/filtered/details/3.interface_residues.csv")


def slug(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(s)).strip("_")[:50]


class ResSel(Select):
    def __init__(self, chain_id, resids):
        self.chain_id = chain_id
        self.resids = resids

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return residue.id[0] == " " and residue.id[1] in self.resids


io = PDBIO()
written, skipped = [], []
for _, r in rep.iterrows():
    iid = int(r.interaction_id)
    abp_full = f"{r.pdb}_{r.abp_chain}"
    sub = res[(res.interaction_id == iid) & (res.chain == abp_full)]
    contacts = sorted(pd.to_numeric(sub.residue_number_structure, errors="coerce").dropna().astype(int))
    if not contacts:
        skipped.append((r.abp_title, iid, abp_full, "0 residu interface"))
        continue
    keep = set()
    for c in contacts:
        for d in range(-PAD, PAD + 1):
            keep.add(c + d)

    base = f"{r.pdb}"
    f = ASM / (base + ".pdb"); parser = PDBParser(QUIET=True)
    if not f.exists():
        f = ASM / (base + ".cif"); parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure(base, str(f))
    name = f"{slug(r.abp_title)}__{r.pdb}_{r.abp_chain}"
    io.set_structure(struct)
    io.save(str(DST / (name + ".pdb")), ResSel(r.abp_chain, keep))
    written.append((name, len(contacts)))

print(f"PAD={PAD} -> {DST}")
print(f"interfaces ecrites : {len(written)}")
print("taille des empreintes (nb residus de contact) :")
import numpy as np
sizes = [n for _, n in written]
print(f"  min={min(sizes)} median={int(np.median(sizes))} max={max(sizes)}")
if skipped:
    print(f"ignores ({len(skipped)}): {[s[0] for s in skipped]}")
