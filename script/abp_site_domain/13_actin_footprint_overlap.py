#!/usr/bin/env python
"""
Recouvrement des EMPREINTES SUR L'ACTINE entre ABP d'un même site.

Pour chaque ABP : ensemble des résidus canoniques (MAFFT) de l'actine qu'il contacte.
Puis, par cluster de site, recouvrement de Jaccard entre paires d'ABP :
  - même famille  vs  familles différentes (fusions exclues, familles de familles.csv)

Question : des familles différentes au même site touchent-elles LES MÊMES résidus
d'actine, ou seulement une région voisine ?

Sorties :
  data/exports/abp_site_domain/actin_footprint_overlap.csv
  data/exports/abp_site_domain/figure_footprint_overlap.png
"""
from pathlib import Path
import sys
import streamlit
from itertools import combinations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/exports/abp_site_domain"
EXCLUDE_FUSION = {
    "Myosin-14,Alpha-actinin A",
    "Maltose/maltodextrin-binding periplasmic protein,Adenylate cyclase ExoY",
    "Maltose/maltodextrin-binding periplasmic protein,TccC3",
    "Coronin-1B,Methylated-DNA--protein-cysteine methyltransferase",
}

df = pd.read_csv(ROOT / "data/filtered/filtered_all_data.csv", low_memory=False)
di = pd.read_csv(ROOT / "data/filtered/details/1.interactions.csv")[["interaction_id", "chain_A_id", "chain_B_id", "resolution"]].rename(columns={"resolution": "res"})
resdf = pd.read_csv(ROOT / "data/filtered/details/3.interface_residues.csv")
sys.path.insert(0, str(ROOT / "script"))
import numbering
resdf["canon"] = numbering.uniprot_series(resdf["residue_number_canon_mafft"])
fams = pd.read_csv(OUT / "familles.csv")

fam_of = {}
for _, r in fams.iterrows():
    for a in str(r.membres).split(" · "):
        fam_of[a.strip()] = r.famille

m = df.merge(di, left_on=["subunit_1", "subunit_2"], right_on=["chain_A_id", "chain_B_id"], how="left")
m = m[m.s1_actine.ne(m.s2_actine)].copy()
m["abp_title"] = np.where(m.s1_actine, m.subunit_2_title, m.subunit_1_title)
m["actin_chain"] = np.where(m.s1_actine, m.subunit_1, m.subunit_2)
m["res"] = pd.to_numeric(m.res, errors="coerce")

# empreinte actine = résidus canoniques contactés côté chaîne actine (subunit_1)
def footprint(iid, actin_chain):
    s = resdf[(resdf.interaction_id == iid) & (resdf.chain == actin_chain)]
    return set(s.canon.dropna().astype(int))

# 1 interaction représentante par (site, ABP) : meilleure résolution
m["site"] = np.where(m.s1_actine, m.s1_binding_site_cluster_data_70, m.s2_binding_site_cluster_data_70)
sites_abp = (m.dropna(subset=["site", "interaction_id"])
             .sort_values(["res", "interaction_id", "actin_chain"], kind="stable")
             .drop_duplicates(["site", "abp_title"]))

# pré-calcul empreintes
fp = {}
for _, r in sites_abp.iterrows():
    fp[(r.site, r.abp_title)] = footprint(int(r.interaction_id), r.actin_chain)

rows = []
for site, g in sites_abp.groupby("site"):
    abps = [a for a in g.abp_title.unique() if a in fam_of]
    if len(abps) < 2:
        continue
    for a, b in combinations(sorted(abps), 2):
        fa, fb = fp[(site, a)], fp[(site, b)]
        if not fa or not fb:
            continue
        inter = len(fa & fb); union = len(fa | fb)
        jac = inter / union if union else np.nan
        rows.append(dict(site=site, a=a, b=b,
                         fam_a=fam_of[a], fam_b=fam_of[b],
                         same_fam=(fam_of[a] == fam_of[b]),
                         n_a=len(fa), n_b=len(fb), partages=inter,
                         jaccard=round(jac, 2),
                         recouvrement_min=round(inter / min(len(fa), len(fb)), 2)))

ov = pd.DataFrame(rows)
# dédoublonner les couples de familles différentes (1 paire représentante = jaccard max)
ov["fam_couple"] = ov.apply(lambda r: tuple(sorted((r.fam_a, r.fam_b))), axis=1)
diff = ov[~ov.same_fam].sort_values("jaccard", ascending=False).groupby(["site", "fam_couple"], as_index=False).first()
same = ov[ov.same_fam]
ov.to_csv(OUT / "actin_footprint_overlap.csv", index=False)

print(f"paires même famille        : {len(same)}  (Jaccard médian {same.jaccard.median():.2f})")
print(f"couples familles diff.     : {len(diff)}  (Jaccard médian {diff.jaccard.median():.2f})")
print(f"  -> recouvrement_min médian familles diff. : {diff.recouvrement_min.median():.2f}")
print(f"  -> couples diff. avec ≥50% de résidus communs : {(diff.recouvrement_min>=0.5).sum()}/{len(diff)}")
print("\n=== top couples de familles différentes par recouvrement ===")
show = diff.sort_values("recouvrement_min", ascending=False)
show["fA"] = show.fam_a.str.replace(r"F\d+ . ", "", regex=True).str[:20]
show["fB"] = show.fam_b.str.replace(r"F\d+ . ", "", regex=True).str[:20]
print(show[["site", "fA", "fB", "partages", "n_a", "n_b", "jaccard", "recouvrement_min"]].head(12).to_string(index=False))

# figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
ax1.boxplot([same.jaccard.dropna(), diff.jaccard.dropna()], tick_labels=["Same\nfamily", "Different\nfamilies"],
            patch_artist=True, boxprops=dict(facecolor="#d9eaf5"), medianprops=dict(color="black"))
rng = np.random.default_rng(0)  # Reproducible display jitter; it never changes the measurements.
ax1.scatter(rng.normal(1, 0.05, len(same)), same.jaccard, color="#0072B2", alpha=0.5, s=18)
ax1.scatter(rng.normal(2, 0.05, len(diff)), diff.jaccard, color="#E69F00", alpha=0.6, s=22)
ax1.set_ylabel("Jaccard of contacted actin residues")
ax1.set_title("Observed actin footprint overlap")

ax2.hist(diff.recouvrement_min.dropna(), bins=np.arange(0, 1.01, 0.1), color="#E69F00", edgecolor="white")
ax2.axvline(0.5, ls="--", color="grey")
ax2.set_xlabel("Shared residues / smaller footprint")
ax2.set_ylabel("Number of pairs of different families")
ax2.set_title("Overlap between families sharing an actin site")
fig.suptitle("Actin footprints of representative interactions (P60709 positions)", fontsize=13, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(OUT / "figure_footprint_overlap.png", dpi=150, bbox_inches="tight")
print(f"\nfigure : {OUT/'figure_footprint_overlap.png'}")
