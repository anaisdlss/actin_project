#!/usr/bin/env python
# coding: utf-8

import os as _os
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # cwd = racine projet (robuste, peu importe d'où on lance)


try:
    from IPython.display import display
except Exception:
    def display(*a, **k):
        for x in a: print(x)



# # Méga-clusters de sites de liaison actine (S1 et S2 confondus)
# 
# Regroupe **tous** les clusters de sites de liaison actine — **sans distinguer S1 et S2** : un même
# site est la même région de surface de l'actine qu'il soit observé côté S1 (actine = subunit_1) ou
# côté S2 (actine = subunit_2). On **unit donc les deux empreintes** par identifiant de cluster, puis
# on regroupe par **contenance à 100 %** : chaque méga-cluster est un cluster maximal qui recouvre
# entièrement les empreintes des plus petits qu'il contient.
# 
# Les 4 familles principales restent **6685_1, 6685_2, 6685_3, 6685_4**.

# In[1]:


import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict


# ## 1. Chargement des données

# In[2]:


ROOT = Path(".")  # depuis notebooks/

df_all  = pd.read_csv(ROOT / "data/filtered/filtered_all_data.csv", low_memory=False)
df_int  = pd.read_csv(ROOT / "data/filtered/details/1.interactions.csv")
df_res  = pd.read_csv(ROOT / "data/filtered/details/3.interface_residues.csv")

# Garder uniquement les résidus avec position canonique connue
df_res["canon"] = pd.to_numeric(df_res["residue_number_canon_mafft"], errors="coerce")
df_res = df_res[df_res["canon"].notna()].copy()

print(f"filtered_all_data : {len(df_all)} lignes, {df_all.shape[1]} colonnes")
print(f"interactions      : {len(df_int)} lignes")
print(f"interface_residues: {len(df_res)} lignes")


# ## 2. Mapping interaction_id → clusters S1 / S2

# In[3]:


# Joindre les clusters de sites de liaison via chain_A_id / chain_B_id
# df_all a subunit_1, subunit_2 ; df_int a chain_A_id, chain_B_id, interaction_id
imap = df_int[["interaction_id", "chain_A_id", "chain_B_id"]].merge(
    df_all[[
        "subunit_1", "subunit_2",
        "s1_binding_site_cluster_data_70",
        "s2_binding_site_cluster_data_70",
        "s1_actine", "s2_actine"
    ]],
    left_on=["chain_A_id", "chain_B_id"],
    right_on=["subunit_1", "subunit_2"],
    how="left"
).drop(columns=["subunit_1", "subunit_2"])

imap["s1_actine"] = imap["s1_actine"].fillna(False).astype(bool)
imap["s2_actine"] = imap["s2_actine"].fillna(False).astype(bool)

print(imap.head(3))


# ## 3. Calcul des empreintes de résidus canoniques actine

# In[4]:


def compute_footprints(imap_subset, cluster_col, actin_chain_col):
    """
    Pour chaque cluster de sites de liaison, calcule l'ensemble des résidus
    canoniques actine contactés (union sur toutes les interactions du cluster).

    actin_chain_col : 'chain_A_id' si actin est S1, 'chain_B_id' si actin est S2
    """
    footprints = {}
    for cluster, grp in imap_subset.groupby(cluster_col):
        if pd.isna(cluster):
            continue
        canon_positions = set()
        for _, row in grp.iterrows():
            iid        = row["interaction_id"]
            actin_ch   = row[actin_chain_col]
            res_subset = df_res[
                (df_res["interaction_id"] == iid) &
                (df_res["chain"] == actin_ch)
            ]
            canon_positions.update(res_subset["canon"].dropna().astype(int).tolist())
        if canon_positions:
            footprints[str(cluster)] = canon_positions
    return footprints

# S1 : actin est chain_A (s1_actine=True)
imap_s1 = imap[imap["s1_actine"] & imap["s1_binding_site_cluster_data_70"].notna()]
fp_s1   = compute_footprints(imap_s1, "s1_binding_site_cluster_data_70", "chain_A_id")

# S2 : actin est chain_B (s2_actine=True)
imap_s2 = imap[imap["s2_actine"] & imap["s2_binding_site_cluster_data_70"].notna()]
fp_s2   = compute_footprints(imap_s2, "s2_binding_site_cluster_data_70", "chain_B_id")

print(f"Empreintes S1 calculées : {len(fp_s1)} clusters")
print(f"Empreintes S2 calculées : {len(fp_s2)} clusters")


# ## 4. Construction des super-clusters par contenance à 100% (sans chaînage)
# 
# Un **super-cluster = un cluster maximal** (dont l'empreinte de résidus canoniques actine
# n'est strictement contenue dans aucune autre) **+ tous les clusters dont l'empreinte y est
# entièrement incluse**.
# 
# Deux règles pour éviter le faux regroupement observé avec l'ancien union-find :
# - **les clusters maximaux ne fusionnent jamais entre eux** : deux gros clusters ne sont pas
#   réunis sous prétexte qu'un même petit site est inclus dans les deux (c'est le chaînage qui
#   faisait atterrir Myosin-7 dans `6685_3` alors qu'elle contacte une région à 0 résidu commun) ;
# - un petit site contenu dans **plusieurs** maximaux **appartient à chacun** (multi-appartenance) —
#   chaque site est ainsi garanti d'être à 100 % dans tout super-cluster qui le porte.

# In[5]:


# Confondre S1 et S2 : un site de liaison actine est la MÊME région qu'il soit
# vu côté S1 (actine=subunit_1) ou côté S2 (actine=subunit_2). On unit donc les
# deux empreintes par identifiant de cluster avant de regrouper.
fp_all = {}
for _d in (fp_s1, fp_s2):
    for _c, _f in _d.items():
        fp_all.setdefault(_c, set()).update(_f)
print(f"Empreintes fusionnées S1 ∪ S2 : {len(fp_all)} clusters de sites")


def build_superclusters(footprints, main_clusters):
    """
    Méga-clusters par CONTENANCE 100% SANS chaînage.

    Un méga-cluster = un cluster MAXIMAL (empreinte non strictement incluse dans
    une autre) + tous les clusters dont l'empreinte y est entièrement contenue.
      - Les clusters maximaux ne fusionnent jamais entre eux (pas d'union-find).
      - Un site contenu dans plusieurs maximaux appartient à CHACUN (multi-appartenance).

    Retourne :
      cluster_to_supers : dict {cluster_id -> [supercluster_names triés]}
      components         : dict {supercluster_name -> [membres]}
    """
    clusters = list(footprints.keys())

    # Clusters maximaux : aucun autre cluster d'empreinte STRICTEMENT plus grande
    # ne les contient à 100%.
    maximal = [
        c for c in clusters
        if not any(
            o != c
            and footprints[c] <= footprints[o]
            and len(footprints[o]) > len(footprints[c])
            for o in clusters
        )
    ]

    # Dédupliquer les maximaux d'empreinte identique : un seul nom (id le + petit)
    rep = {}
    by_fp = {}
    for m in sorted(maximal):
        key = frozenset(footprints[m])
        rep[m] = by_fp.setdefault(key, m)

    # Multi-appartenance : chaque cluster rejoint TOUT maximal qui le contient à 100%
    cluster_to_supers = {}
    for c in clusters:
        supers = sorted({
            rep[m] for m in maximal
            if footprints[c] <= footprints[m]
        })
        cluster_to_supers[c] = supers

    components = defaultdict(list)
    for c, supers in cluster_to_supers.items():
        for s in supers:
            components[s].append(c)

    n_super = len(set(rep.values()))
    n_multi = sum(1 for v in cluster_to_supers.values() if len(v) > 1)
    print(f"  {n_super} méga-clusters (maximaux) · "
          f"{n_multi}/{len(clusters)} sites en multi-appartenance")

    return cluster_to_supers, components


MAIN = ["6685_1", "6685_2", "6685_3", "6685_4"]

print("=== Méga-clusters (S1 + S2 confondus) ===")
c2s, comp = build_superclusters(fp_all, MAIN)
print(f"  {len(comp)} méga-clusters au total")


# ## 5. Vue détaillée des super-clusters

# In[6]:


def show_superclusters(c2s, footprints):
    """c2s : {cluster -> [megaclusters]} (multi-appartenance, S1+S2 confondus)."""
    rows = []
    for c, supers in c2s.items():
        for s in supers:
            rows.append({
                "binding_site_cluster": c,
                "supercluster": s,
                "n_residues": len(footprints.get(c, set())),
                "is_main_cluster": c in ["6685_1", "6685_2", "6685_3", "6685_4"],
            })
    return pd.DataFrame(rows)

df_sc = show_superclusters(c2s, fp_all)

print("=== Méga-clusters (S1 + S2 confondus) ===")
display(df_sc.groupby("supercluster")[["binding_site_cluster", "n_residues"]].agg(list))


# ## 6. Ajout des colonnes dans filtered_all_data.csv

# In[7]:


df_out = df_all.copy()

# s1/s2_supercluster : même mapping fusionné c2s (S1+S2 confondus), liste des
# méga-clusters contenant le site, jointe par ';'. Rempli si le côté est actine.
def _map_super(side_actine, bs_col):
    def f(row):
        if not row.get(side_actine, False):
            return None
        c = str(row[bs_col]) if pd.notna(row[bs_col]) else None
        supers = c2s.get(c) if c else None
        return ";".join(supers) if supers else None
    return f

df_out["s1_supercluster"] = df_out.apply(
    _map_super("s1_actine", "s1_binding_site_cluster_data_70"), axis=1)
df_out["s2_supercluster"] = df_out.apply(
    _map_super("s2_actine", "s2_binding_site_cluster_data_70"), axis=1)

# Nombre de méga-clusters par site (repère les sites ambigus en multi-appartenance)
df_out["s1_n_superclusters"] = df_out["s1_supercluster"].apply(
    lambda v: len(v.split(";")) if isinstance(v, str) and v else 0)
df_out["s2_n_superclusters"] = df_out["s2_supercluster"].apply(
    lambda v: len(v.split(";")) if isinstance(v, str) and v else 0)

print("s1_supercluster — distribution :")
print(df_out["s1_supercluster"].value_counts(dropna=False).head(10))
print("\ns2_supercluster — distribution :")
print(df_out["s2_supercluster"].value_counts(dropna=False).head(10))


# In[8]:


# Sauvegarder filtered_all_data.csv avec les nouvelles colonnes
out_path = ROOT / "data/filtered/filtered_all_data.csv"
df_out.to_csv(out_path, index=False)
print(f"Sauvegardé : {out_path}")
print(f"Nouvelles colonnes : s1_supercluster, s2_supercluster")


# ## 7. CSV de traçabilité

# In[9]:


df_trace = df_sc.sort_values(
    ["supercluster", "n_residues"], ascending=[True, False]
).reset_index(drop=True)

trace_path = ROOT / "data/filtered/binding_site_superclusters.csv"
df_trace.to_csv(trace_path, index=False)
print(f"CSV de traçabilité : {trace_path}  ({len(df_trace)} appartenances)")
display(df_trace)


# ## 8. Résumé des empreintes des 4 familles principales

# In[10]:


print("Empreintes des 4 familles principales (résidus canoniques actine) :")
for m in MAIN:
    fp_s1_m = fp_s1.get(m, set())
    fp_s2_m = fp_s2.get(m, set())
    print(f"\n{m}")
    if fp_s1_m:
        print(f"  S1 : {len(fp_s1_m)} résidus — min={min(fp_s1_m)}, max={max(fp_s1_m)}")
        print(f"       {sorted(fp_s1_m)[:10]}{'...' if len(fp_s1_m)>10 else ''}")
    if fp_s2_m:
        print(f"  S2 : {len(fp_s2_m)} résidus — min={min(fp_s2_m)}, max={max(fp_s2_m)}")
        print(f"       {sorted(fp_s2_m)[:10]}{'...' if len(fp_s2_m)>10 else ''}")

