#!/usr/bin/env python
"""
Etape 5 — Croisement : les ABPs d'un meme cluster de site de liaison actine
partagent-ils le meme fold (Foldseek) et le meme domaine (InterPro) ?

Sorties :
  data/exports/abp_site_domain/site_cluster_summary.csv
  data/exports/abp_site_domain/abp_master.csv
  data/exports/abp_site_domain/figure_site_vs_domaine.png
"""
import re
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/exports/abp_site_domain"

t = pd.read_csv(OUT / "abp_site_table.csv")
rep = pd.read_csv(OUT / "abp_representatives.csv")
ip = pd.read_csv(OUT / "abp_interpro.csv")
fc = pd.read_csv(OUT / "fold_cluster_cluster.tsv", sep="\t", header=None, names=["rep", "mem"])


def slug(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(s)).strip("_")[:50]


# abp_title -> stem -> fold cluster (rep)
rep["stem"] = rep.apply(lambda r: f"{slug(r.abp_title)}__{r.pdb}_{r.abp_chain}", axis=1)
stem2fold = dict(zip(fc.mem, fc.rep))
rep["fold_cluster"] = rep.stem.map(stem2fold)

# table maitresse par ABP
master = rep[["abp_title", "pdb", "abp_chain", "fold_cluster"]].merge(
    ip[["abp_title", "uniprot", "pfam_domains", "interpro_domains"]], on="abp_title", how="left")
master["pfam_domains"] = master["pfam_domains"].fillna("(non annote)")
master.to_csv(OUT / "abp_master.csv", index=False)

# (site cluster, abp) uniques
pairs = t.dropna(subset=["actin_site_cluster"]).drop_duplicates(["actin_site_cluster", "abp_title"])
pairs = pairs.merge(master[["abp_title", "fold_cluster", "pfam_domains"]], on="abp_title", how="left")

rows = []
for site, g in pairs.groupby("actin_site_cluster"):
    abps = sorted(g.abp_title.unique())
    folds = set(g.fold_cluster.dropna())
    doms = set(g.pfam_domains)
    rows.append(dict(
        site_cluster=site, n_abp=len(abps),
        n_folds=len(folds), n_domaines=len(doms),
        mono_fold=(len(folds) == 1), mono_domaine=(len(doms) == 1),
        abps=" / ".join(a[:25] for a in abps),
    ))
summ = pd.DataFrame(rows).sort_values(["n_abp", "site_cluster"], ascending=[False, True])
summ.to_csv(OUT / "site_cluster_summary.csv", index=False)

multi = summ[summ.n_abp >= 2]
print(f"clusters de site total              : {len(summ)}")
print(f"  dont mono-ABP (non informatifs)   : {(summ.n_abp == 1).sum()}")
print(f"  dont multi-ABP (>=2)              : {len(multi)}")
print(f"\nParmi les {len(multi)} clusters multi-ABP :")
print(f"  meme FOLD (Foldseek)     : {multi.mono_fold.sum()}/{len(multi)} ({100*multi.mono_fold.mean():.0f}%)")
print(f"  meme DOMAINE (Pfam)      : {multi.mono_domaine.sum()}/{len(multi)} ({100*multi.mono_domaine.mean():.0f}%)")
print(f"\nDetail des clusters multi-ABP :")
print(multi[["site_cluster", "n_abp", "n_folds", "mono_fold", "abps"]].to_string(index=False))

# ---------- FIGURE ----------
fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# Panel A : mono vs mixte (fold et domaine)
ax = axes[0]
cats = ["Même fold\n(Foldseek)", "Même domaine\n(Pfam)"]
mono = [multi.mono_fold.sum(), multi.mono_domaine.sum()]
mixte = [(~multi.mono_fold).sum(), (~multi.mono_domaine).sum()]
x = np.arange(len(cats))
ax.bar(x, mono, 0.55, label="1 seul (convergence d'une famille)", color="#2a9d8f")
ax.bar(x, mixte, 0.55, bottom=mono, label="plusieurs (familles différentes)", color="#e76f51")
for i, (m, mx) in enumerate(zip(mono, mixte)):
    ax.text(i, m / 2, str(m), ha="center", va="center", color="white", fontweight="bold")
    if mx:
        ax.text(i, m + mx / 2, str(mx), ha="center", va="center", color="white", fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels(cats)
ax.set_ylabel("Nb de clusters de site multi-ABP")
ax.set_title(f"Clusters de site touchés par ≥2 ABP (n={len(multi)})")
ax.legend(fontsize=9)

# Panel B : detail par cluster, ABP colorés par fold
ax = axes[1]
mdet = multi.sort_values("n_abp")
folds_all = sorted(pairs.fold_cluster.dropna().unique())
cmap = plt.get_cmap("tab20")
fold_color = {f: cmap(i % 20) for i, f in enumerate(folds_all)}
for yi, (_, r) in enumerate(mdet.iterrows()):
    g = pairs[pairs.actin_site_cluster == r.site_cluster].drop_duplicates("abp_title")
    for xi, (_, a) in enumerate(g.iterrows()):
        ax.scatter(xi, yi, s=260, color=fold_color.get(a.fold_cluster, "grey"),
                   edgecolor="black", linewidth=0.5)
        ax.text(xi, yi, a.abp_title[:14], ha="center", va="center", fontsize=6)
ax.set_yticks(range(len(mdet)))
ax.set_yticklabels([f"{r.site_cluster} ({'mono' if r.mono_fold else 'MIXTE'})"
                    for _, r in mdet.iterrows()], fontsize=8)
ax.set_xlabel("ABPs du cluster (couleur = fold Foldseek)")
ax.set_title("Composition de chaque cluster multi-ABP")
ax.set_xlim(-0.6, mdet.n_abp.max() - 0.4)
fig.suptitle("Les ABPs d'un même site de l'actine partagent-ils le même domaine ?",
             fontsize=14, y=1.02)
fig.tight_layout()
fig.savefig(OUT / "figure_site_vs_domaine.png", dpi=150, bbox_inches="tight")
print(f"\nfigure : {OUT/'figure_site_vs_domaine.png'}")
