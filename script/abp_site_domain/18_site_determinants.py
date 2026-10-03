#!/usr/bin/env python
"""
Pourquoi des ABP différentes ciblent-elles les mêmes résidus d'actine ?
=> propriété de la SURFACE de l'actine, pas des ABP.

Pour chaque résidu canonique de l'actine : nb de familles d'ABP distinctes qui le
touchent, croisé avec exposition (RSA) et conservation (ProteoCast actine).

Sorties :
  data/exports/abp_site_domain/actin_residue_determinants.csv
  data/exports/abp_site_domain/figure_site_determinants.png
"""
import sys
import streamlit
from pathlib import Path
from collections import defaultdict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/exports/abp_site_domain"

df = pd.read_csv(ROOT / "data/filtered/filtered_all_data.csv", low_memory=False)
di = pd.read_csv(ROOT / "data/filtered/details/1.interactions.csv")[
    ["interaction_id", "chain_A_id", "chain_B_id"]]
res = pd.read_csv(ROOT / "data/filtered/details/3.interface_residues.csv")
res["canon"] = pd.to_numeric(res["residue_number_canon_mafft"], errors="coerce")
fam = pd.read_csv(OUT / "familles.csv")
fam_of = {a.strip(): r.famille for _, r in fam.iterrows() for a in str(r.membres).split(" · ")}
sys.path.insert(0, str(ROOT / "script"))
from scientific_analysis import canonical_conservation
from footprint_comparison import footprint_records
import numbering
cons = canonical_conservation(ROOT)

fam_by_canon = defaultdict(set)
records = footprint_records(df, di, res)
for r in records[records.kind.eq('abp') & records.asa.gt(0)].itertuples():
    fa = fam_of.get(r.group)
    if not fa:
        continue
    fam_by_canon[numbering.to_canon(int(r.position))].add(fa)

nfam = pd.Series({c: len(s) for c, s in fam_by_canon.items()}, name="n_familles")
t = cons.set_index("canon").join(nfam).copy()
t["n_familles"] = t["n_familles"].fillna(0).astype(int)
t["rsa"] = pd.to_numeric(t["rsa"], errors="coerce")
t.reset_index().to_csv(OUT / "actin_residue_determinants.csv", index=False)

sub = t.dropna(subset=["conservation", "rsa"])
rho_rsa, p_rsa = spearmanr(sub.n_familles, sub.rsa)
rho_cons, p_cons = spearmanr(sub.n_familles, sub.conservation)
print(f"n_familles vs RSA          : rho={rho_rsa:.2f} (p={p_rsa:.1e})")
print(f"n_familles vs conservation : rho={rho_cons:.2f} (p={p_cons:.1e})")

# figure
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
a1.scatter(sub.rsa, sub.n_familles, s=14, alpha=0.5, color="#E69F00")
a1.set_xlabel("RSA of human actin alone (8DNH chain B)")
a1.set_ylabel("Observed contacting ABP families")
a1.set_title(f"RSA association   (ρ={rho_rsa:.2f})")
a2.scatter(sub.conservation, sub.n_familles, s=14, alpha=0.5, color="#0072B2")
a2.set_xlabel("Mutational sensitivity (ProteoCast)")
a2.set_ylabel("Observed contacting ABP families")
a2.set_title(f"Sensitivity association   (ρ={rho_cons:.2f})")
fig.suptitle("Observed family contacts, reference RSA and mutational sensitivity",
             fontsize=13, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(OUT / "figure_site_determinants.png", dpi=150, bbox_inches="tight")
print("écrit : actin_residue_determinants.csv + figure_site_determinants.png")
