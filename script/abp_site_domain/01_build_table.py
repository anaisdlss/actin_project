#!/usr/bin/env python
"""
Etape 1 — Table ABP x cluster de site de liaison actine.

Question : les ABPs qui touchent le MEME patch de surface de l'actine
(meme binding_site_cluster) partagent-ils le meme domaine / fold ?

Sortie :
  data/exports/abp_site_domain/abp_site_table.csv  (toutes les interactions ABP-actine)
  data/exports/abp_site_domain/abp_representatives.csv (1 structure par ABP)
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/exports/abp_site_domain"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(ROOT / "data/filtered/filtered_all_data.csv", low_memory=False)
di = pd.read_csv(ROOT / "data/filtered/details/1.interactions.csv")

# interaction_id + resolution via les paires de chaines
di_key = di[["interaction_id", "chain_A_id", "chain_B_id", "resolution"]].rename(
    columns={"resolution": "res_struct"})
m = df.merge(di_key, left_on=["subunit_1", "subunit_2"],
             right_on=["chain_A_id", "chain_B_id"], how="left")

recs = []
for _, r in m.iterrows():
    # garder uniquement actine <-> ABP (exclure filament actine-actine)
    if r.s1_actine and not r.s2_actine:        # s1 = actine, s2 = ABP
        abp_sub, abp_title = r.subunit_2, r.subunit_2_title
        site = r.s1_binding_site_cluster_data_70
        sup = r.get("s1_supercluster")
    elif r.s2_actine and not r.s1_actine:      # s2 = actine, s1 = ABP
        abp_sub, abp_title = r.subunit_1, r.subunit_1_title
        site = r.s2_binding_site_cluster_data_70
        sup = r.get("s2_supercluster")
    else:
        continue
    if pd.isna(r.interaction_id) or pd.isna(abp_sub):
        continue
    pdb = str(abp_sub).split("_")[0]
    chain = str(abp_sub).split("_")[1]
    recs.append(dict(
        interaction_id=int(r.interaction_id), pdb=pdb, abp_chain=chain,
        abp_subunit=abp_sub, abp_title=abp_title,
        actin_site_cluster=site, supercluster=sup,
        resolution=r.res_struct,
    ))

t = pd.DataFrame(recs)
t.to_csv(OUT / "abp_site_table.csv", index=False)

# 1 representant structural par ABP title : meilleure resolution dispo
t_res = t.dropna(subset=["abp_title"]).copy()
t_res["resolution"] = pd.to_numeric(t_res["resolution"], errors="coerce")
rep = (t_res.sort_values("resolution")
            .groupby("abp_title", as_index=False).first())
rep.to_csv(OUT / "abp_representatives.csv", index=False)

print(f"interactions ABP-actine : {len(t)}")
print(f"ABP titles distincts    : {t.abp_title.nunique()}")
print(f"clusters de site (_70)  : {t.actin_site_cluster.nunique()}")
print(f"representants ecrits     : {len(rep)}")
print("\nExemple representants :")
print(rep[["abp_title", "pdb", "abp_chain", "resolution", "actin_site_cluster"]].head(10).to_string(index=False))
