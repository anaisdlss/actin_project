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
cons = pd.read_csv(ROOT / "data/proteocast/conservation_vs_asa_per_position.csv")

m = df.merge(di, left_on=["subunit_1", "subunit_2"],
             right_on=["chain_A_id", "chain_B_id"], how="left")
m = m[(m.s1_actine) & (~m.s2_actine)]

fam_by_canon = defaultdict(set)
for _, r in m.iterrows():
    fa = fam_of.get(r.subunit_2_title)
    if not fa:
        continue
    rr = res[(res.interaction_id == r.interaction_id) & (res.chain == r.subunit_1)]
    for c in rr.canon.dropna().astype(int):
        fam_by_canon[c].add(fa)

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
a1.scatter(sub.rsa, sub.n_familles, s=14, alpha=0.5, color="#e76f51")
a1.set_xlabel("Exposition au solvant du résidu d'actine (RSA)")
a1.set_ylabel("Nb de familles d'ABP qui le touchent")
a1.set_title(f"Exposition → convergence   (ρ={rho_rsa:.2f})")
a2.scatter(sub.conservation, sub.n_familles, s=14, alpha=0.5, color="#2a9d8f")
a2.set_xlabel("Conservation du résidu d'actine (ProteoCast)")
a2.set_ylabel("Nb de familles d'ABP qui le touchent")
a2.set_title(f"Conservation → peu prédictif   (ρ={rho_cons:.2f})")
fig.suptitle("Pourquoi ces résidus d'actine ? — c'est l'EXPOSITION qui attire les ABP",
             fontsize=13, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(OUT / "figure_site_determinants.png", dpi=150, bbox_inches="tight")
print("écrit : actin_residue_determinants.csv + figure_site_determinants.png")
