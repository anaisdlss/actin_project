#!/usr/bin/env python
"""
Regroupement des ABP en FAMILLES, puis test des convergences INTER-FAMILLES.

1. Famille = composante connexe d'un graphe où deux ABP sont reliés s'ils
   partagent ≥1 domaine Pfam (hors tags de fusion) OU le même accession UniProt.
2. Pour chaque cluster de site : familles distinctes présentes.
   Cluster "multi-familles" = ≥2 familles → vraie convergence.
3. Pour chaque PAIRE de familles différentes d'un même site, meilleur TM-score
   entre leurs membres (protéine entière) + %id.

Sorties :
  data/exports/abp_site_domain/familles.csv
  data/exports/abp_site_domain/convergences_inter_familles.csv
  data/exports/abp_site_domain/figure_convergence_familles.png
"""
import re
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/exports/abp_site_domain"
TAG_DROP = {"PF01547", "PF01035", "PF02870"}
UNI_DROP = {"P0AEX9"}  # MBP (tag de cristallisation) — ne doit pas relier les fusions
# constructions de fusion/chimères : faussent le graphe de familles (ponts artificiels) → exclues
EXCLUDE_FUSION = {
    "Myosin-14,Alpha-actinin A",  # chimère myosine+actinine → pont myosines/CH
    "Maltose/maltodextrin-binding periplasmic protein,Adenylate cyclase ExoY",
    "Maltose/maltodextrin-binding periplasmic protein,TccC3",
    "Coronin-1B,Methylated-DNA--protein-cysteine methyltransferase",
}

ip = pd.read_csv(OUT / "abp_interpro.csv")
m = pd.read_csv(OUT / "abp_master.csv")
t = pd.read_csv(OUT / "abp_site_table.csv")
wp = pd.read_csv(OUT / "whole_pairs_all.tsv", sep="\t", header=None,
                 names=["q", "tt", "qtm", "ttm", "rmsd", "fident", "alnlen"])


def slug(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(s)).strip("_")[:50]


abps = sorted(a for a in ip.abp_title.unique() if a not in EXCLUDE_FUSION)
pfam = {r.abp_title: ((set() if pd.isna(r.pfam_acc) else set(str(r.pfam_acc).split(";"))) - TAG_DROP) for _, r in ip.iterrows()}
pfam_names = {r.abp_title: (str(r.pfam_domains) if pd.notna(r.pfam_domains) else "") for _, r in ip.iterrows()}
uni = {r.abp_title: ((set() if pd.isna(r.uniprot) else set(str(r.uniprot).split(";"))) - UNI_DROP) for _, r in ip.iterrows()}
stem = {r.abp_title: f"{slug(r.abp_title)}__{r.pdb}_{r.abp_chain}" for _, r in m.iterrows()}
TM = {(r.q, r.tt): (r.qtm, r.fident) for _, r in wp.iterrows()}
def pair(a, b):
    c = [TM.get((stem[a], stem[b])), TM.get((stem[b], stem[a]))]
    c = [x for x in c if x is not None]
    return max(c, key=lambda x: x[0]) if c else None


# --- union-find pour les familles ---
parent = {a: a for a in abps}
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
def union(a, b):
    parent[find(a)] = find(b)

for a, b in combinations(abps, 2):
    if (pfam[a] and pfam[b] and (pfam[a] & pfam[b])) or (uni[a] and uni[b] and (uni[a] & uni[b])):
        union(a, b)

fam = {}
for a in abps:
    fam.setdefault(find(a), []).append(a)

# nom de famille = domaine Pfam le plus partagé, sinon nom du membre
def fam_label(members):
    common = set.intersection(*[pfam[x] for x in members]) if all(pfam[x] for x in members) else set()
    if common:
        # nom lisible du 1er domaine commun
        names = pfam_names[members[0]].split(" | ")
        return names[0] if names else members[0]
    return members[0] + (f" (+{len(members)-1})" if len(members) > 1 else "")

fam_of = {}
fam_rows = []
for i, (root, members) in enumerate(sorted(fam.items(), key=lambda kv: -len(kv[1])), 1):
    lab = f"F{i:02d} · {fam_label(members)}"
    for x in members:
        fam_of[x] = lab
    fam_rows.append(dict(famille=lab, n_ABP=len(members), membres=" · ".join(members)))
pd.DataFrame(fam_rows).to_csv(OUT / "familles.csv", index=False)

print(f"{len(abps)} ABP regroupés en {len(fam)} familles")
for r in fam_rows:
    print(f"  {r['famille'][:45]:45s} ({r['n_ABP']}) : {r['membres'][:80]}")

# --- convergences inter-familles par cluster de site ---
pt = t.dropna(subset=["actin_site_cluster"]).drop_duplicates(["actin_site_cluster", "abp_title"])
conv = []
multi_fam_clusters = 0
for site, g in pt.groupby("actin_site_cluster"):
    members = sorted(a for a in g.abp_title.unique() if a in fam_of)  # exclut les fusions
    fams = sorted({fam_of[a] for a in members})
    if len(fams) < 2:
        continue
    multi_fam_clusters += 1
    # pour chaque paire de familles, meilleure paire d'ABP (TM max)
    by_fam = {}
    for a in members:
        by_fam.setdefault(fam_of[a], []).append(a)
    for fa, fb in combinations(fams, 2):
        best = None
        for a in by_fam[fa]:
            for b in by_fam[fb]:
                r = pair(a, b)
                if r and (best is None or r[0] > best[0]):
                    best = (r[0], r[1], a, b)
        if best:
            tm, fid, a, b = best
            conv.append(dict(site=site, familleA=fa, familleB=fb, TM=round(tm, 2),
                             pid=round(100*fid), repA=a, repB=b,
                             convergence="repliement proche" if tm >= 0.5 else "structures différentes"))
        else:
            conv.append(dict(site=site, familleA=fa, familleB=fb, TM=np.nan, pid=np.nan,
                             repA=by_fam[fa][0], repB=by_fam[fb][0], convergence="non alignable"))

cv = pd.DataFrame(conv)
cv.to_csv(OUT / "convergences_inter_familles.csv", index=False)

n_pairs = len(cv)
n_ge = (cv.TM >= 0.5).sum()
n_lt = ((cv.TM < 0.5)).sum()
n_na = cv.TM.isna().sum()
print(f"\nClusters multi-familles : {multi_fam_clusters}")
print(f"Paires de familles différentes au même site : {n_pairs}")
print(f"  TM ≥ 0.50 (repliement proche)   : {n_ge}  ({100*n_ge/n_pairs:.0f}%)")
print(f"  TM < 0.50 (structures différentes): {n_lt}  ({100*n_lt/n_pairs:.0f}%)")
print(f"  non alignable                    : {n_na}")
print(f"TM médian inter-familles : {cv.TM.median():.2f}")

# --- figure ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), gridspec_kw={"width_ratios": [1, 1.3]})
vals = cv.TM.dropna()
ax1.hist(vals, bins=np.arange(0, 1.01, 0.1), color="#e76f51", edgecolor="white")
ax1.axvline(0.5, ls="--", color="grey")
ax1.text(0.51, ax1.get_ylim()[1]*0.9, "seuil 0.5", color="grey")
ax1.set_xlabel("TM-score (meilleure paire entre 2 familles d'un même site)")
ax1.set_ylabel("Nb de paires de familles")
ax1.set_title(f"Convergences inter-familles (n={n_pairs} paires)")

cats = ["structures\ndifférentes\n(TM<0.5)", "repliement\nproche\n(TM≥0.5)", "non\nalignable"]
ax2.bar(cats, [n_lt, n_ge, n_na], color=["#e76f51", "#8ab17d", "#bbbbbb"], edgecolor="black")
for i, vv in enumerate([n_lt, n_ge, n_na]):
    ax2.text(i, vv + 0.3, str(vv), ha="center", fontweight="bold")
ax2.set_ylabel("Nb de paires de familles différentes")
ax2.set_title("Quand 2 familles partagent un site, se ressemblent-elles ?")
fig.suptitle(f"Après regroupement en familles : {len(fam)} familles, "
             f"{multi_fam_clusters} sites multi-familles", fontsize=13, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(OUT / "figure_convergence_familles.png", dpi=150, bbox_inches="tight")
print(f"\nfigures + csv écrits dans {OUT}")
