import streamlit as st
import pandas as pd
import numpy as np
import os
from pathlib import Path as _Path

_ABP_HM_FILES = [
    "data/filtered/details/3.interface_residues.csv",
    "data/filtered/details/1.interactions.csv",
    "data/filtered/filtered_all_data.csv",
    "data/filtered/proteins_per_pdb.csv",
]

_S1_SUPER_FILES = [
    "data/filtered/details/3.interface_residues.csv",
    "data/filtered/details/1.interactions.csv",
    "data/filtered/filtered_all_data.csv",
]


@st.cache_data(show_spinner="Calcul heatmap ABP…")
def _build_abp_heatmap_data(*_):
    if not all(os.path.exists(f) for f in _ABP_HM_FILES):
        return None
    df3 = pd.read_csv(_ABP_HM_FILES[0])
    df_int_h = pd.read_csv(_ABP_HM_FILES[1])
    df_all_h = pd.read_csv(_ABP_HM_FILES[2])
    df_pp_h = pd.read_csv(_ABP_HM_FILES[3])

    df3["buried_ASA_percent"] = pd.to_numeric(
        df3["buried_ASA_percent"].astype(
            str).str.replace("%", "", regex=False),
        errors="coerce",
    )
    df3["residue_number_canon_mafft"] = pd.to_numeric(
        df3["residue_number_canon_mafft"], errors="coerce"
    )
    df3 = df3[df3["residue_number_canon_mafft"].notna()
              & df3["buried_ASA_percent"].notna()].copy()

    _actin_ch = set(df_pp_h[df_pp_h["is_actin"]]["chain"].str.lower())
    homo_iids = set(
        df_int_h[df_int_h["chain_B_id"].str.lower().isin(_actin_ch)
                 ]["interaction_id"]
    )
    het_int = (
        df_int_h[~df_int_h["interaction_id"].isin(homo_iids)]
        .merge(
            df_all_h[["subunit_1", "subunit_2",
                      "subunit_2_title", "s2_actine", "cluster_data_70"]],
            left_on=["chain_A_id", "chain_B_id"],
            right_on=["subunit_1", "subunit_2"], how="left",
        )
        .drop_duplicates("interaction_id")
    )
    het_int = het_int[het_int["s2_actine"].fillna(False) == False].copy()
    het_int["abp_name"] = (
        het_int["subunit_2_title"].fillna("Unknown")
        .str.replace(r"\s*\(.*?\)", "", regex=True).str.strip().str[:50]
    )
    _s1ch = het_int.set_index("interaction_id")["chain_A_id"].str.lower()
    _abpn = het_int.set_index("interaction_id")["abp_name"]
    _c70n = het_int.set_index("interaction_id")["cluster_data_70"]

    het_ids = set(het_int["interaction_id"])
    df3_h = df3[df3["interaction_id"].isin(het_ids)].copy()
    df3_h["_s1c"] = df3_h["interaction_id"].map(_s1ch)
    df3_h = df3_h[df3_h["chain"].str.lower() == df3_h["_s1c"]].copy()
    df3_h["abp"] = df3_h["interaction_id"].map(_abpn)
    df3_h["canon"] = df3_h["residue_number_canon_mafft"].astype(int)
    df3_h["c70"] = df3_h["interaction_id"].map(_c70n)

    # Nombre d'interactions (pour label et tri)
    abp_freq = (
        df3_h.groupby("abp")["interaction_id"].nunique()
        .sort_values(ascending=False)
    )

    # ── Moyenne équitable par cluster C70 ────────────────────────────────────
    # Étape 1 : n interactions par (abp, c70)
    _abp_c70_n = df3_h.groupby(["abp", "c70"])["interaction_id"].nunique()

    # Étape 2 : somme ASA par (abp, c70, canon) — 0 implicite si absent
    _agg_c70 = (
        df3_h.groupby(["abp", "c70", "canon"])["buried_ASA_percent"]
        .sum().reset_index(name="asa_sum")
    )

    # Étape 3 : moyenne par c70 (dénominateur = toutes les interactions du c70)
    _agg_c70["asa_c70"] = _agg_c70.apply(
        lambda r: r["asa_sum"] / max(_abp_c70_n.get((r["abp"], r["c70"]), 1), 1),
        axis=1,
    )

    # Étape 4 : n clusters C70 par ABP
    _n_c70_abp = df3_h.groupby("abp")["c70"].nunique()

    # Étape 5 : moyenne équitable (chaque c70 pèse 1)
    _agg_eq = (
        _agg_c70.groupby(["abp", "canon"])["asa_c70"]
        .sum().reset_index()
    )
    _agg_eq["buried_ASA_percent"] = (
        _agg_eq["asa_c70"] / _agg_eq["abp"].map(_n_c70_abp)
    )
    agg = _agg_eq[["abp", "canon", "buried_ASA_percent"]]

    n_abp = max(abp_freq.shape[0], 1)
    res_freq = agg.groupby("canon")["abp"].nunique() / n_abp

    pivot = agg.pivot(index="abp", columns="canon",
                      values="buried_ASA_percent")
    pivot = pivot.loc[abp_freq.index.intersection(pivot.index)]
    return pivot, abp_freq, res_freq


@st.cache_data(show_spinner="Computing C70 competition (overlap)…")
def _build_c70_jaccard_edges(jaccard_threshold: float, *_mtimes):
    """
    Pour chaque cluster C70, construit l'empreinte de résidus canonical actin contactés.
    Deux ABPs sont en compétition si au moins une paire de leurs C70 a
    recouvrement ≥ seuil, où recouvrement = |A ∩ B| / min(|A|, |B|).
    Cela détecte le cas où une grande interface contient totalement une petite.
    Retourne : dict {(abp_a, abp_b): n_paires_au_dessus_seuil}
    """
    import numpy as np  # noqa: F401

    if not all(os.path.exists(f) for f in _S1_SUPER_FILES):
        return {}

    df3 = pd.read_csv(_S1_SUPER_FILES[0])
    df3 = df3[df3["residue_number_canon_mafft"].notna()][[
        "interaction_id", "chain", "residue_number_canon_mafft"
    ]]
    df1 = pd.read_csv(_S1_SUPER_FILES[1],
                      usecols=["interaction_id", "chain_A_id", "chain_B_id"])
    df_all_s = pd.read_csv(_S1_SUPER_FILES[2])
    df_all_s["s1_actine"] = df_all_s["s1_actine"].fillna(False).astype(bool)
    df_all_s["s2_actine"] = df_all_s["s2_actine"].fillna(False).astype(bool)

    pp_c70 = pd.read_csv("data/filtered/proteins_per_pdb.csv")
    abp_pp_c70 = pp_c70[~pp_c70["is_actin"]].copy()
    abp_pp_c70["chain_low"] = abp_pp_c70["chain"].str.lower()

    # Interactions hétéro (les deux directions) avec cluster C70 et chaîne actin
    h1 = df_all_s[df_all_s["s1_actine"] & ~df_all_s["s2_actine"]][[
        "subunit_1", "subunit_2", "cluster_data_70"
    ]].rename(columns={"subunit_1": "actin_chain", "subunit_2": "abp_chain"})
    h2 = df_all_s[~df_all_s["s1_actine"] & df_all_s["s2_actine"]][[
        "subunit_2", "subunit_1", "cluster_data_70"
    ]].rename(columns={"subunit_2": "actin_chain", "subunit_1": "abp_chain"})
    hetero_c70 = pd.concat([h1, h2], ignore_index=True).dropna(
        subset=["cluster_data_70"])
    hetero_c70["abp_chain_low"] = hetero_c70["abp_chain"].str.lower()

    # Récupérer interaction_id via 1.interactions.csv
    hetero_c70 = hetero_c70.merge(
        df1,
        left_on=["actin_chain", "abp_chain"],
        right_on=["chain_A_id", "chain_B_id"],
        how="inner",
    )

    # Empreinte par cluster C70 : ensemble de résidus canonical actin
    merged_c70_df3 = hetero_c70.merge(df3, on="interaction_id", how="inner")
    merged_c70_df3 = merged_c70_df3[
        merged_c70_df3["chain"] == merged_c70_df3["actin_chain"]
    ]
    c70_footprint: dict = {}
    for c70, grp in merged_c70_df3.groupby("cluster_data_70"):
        c70_footprint[str(c70)] = set(
            grp["residue_number_canon_mafft"].unique())

    # ABP → ensemble de C70 avec leur empreinte connue
    hetero_c70_abp = hetero_c70.merge(
        abp_pp_c70[["chain_low", "protein"]],
        left_on="abp_chain_low", right_on="chain_low", how="inner",
    )
    abp_to_c70s: dict = {}
    for _, _r in hetero_c70_abp.dropna(subset=["protein", "cluster_data_70"]).iterrows():
        _c = str(_r["cluster_data_70"])
        if _c in c70_footprint:
            abp_to_c70s.setdefault(_r["protein"], set()).add(_c)

    # Arêtes : au moins 1 paire de C70 avec recouvrement ≥ seuil
    abp_list_c70 = sorted(abp_to_c70s.keys())
    edge_wts_c70: dict = {}
    abp_c70_comp: dict = {}  # ABP → set de C70 impliqués dans au moins 1 compétition
    for _ii, _a in enumerate(abp_list_c70):
        for _b in abp_list_c70[_ii + 1:]:
            n_pairs = 0
            for _ca in abp_to_c70s[_a]:
                fa = c70_footprint[_ca]
                for _cb in abp_to_c70s[_b]:
                    fb = c70_footprint[_cb]
                    _inter = len(fa & fb)
                    _min_sz = min(len(fa), len(fb))
                    if _min_sz > 0 and _inter / _min_sz >= jaccard_threshold:
                        n_pairs += 1
                        abp_c70_comp.setdefault(_a, set()).add(_ca)
                        abp_c70_comp.setdefault(_b, set()).add(_cb)
            if n_pairs > 0:
                _k = (min(_a, _b), max(_a, _b))
                edge_wts_c70[_k] = n_pairs

    abp_c70_comp_count = {k: len(v) for k, v in abp_c70_comp.items()}
    return edge_wts_c70, abp_c70_comp_count, c70_footprint, abp_to_c70s


@st.cache_data(show_spinner="Calcul super-clusters S1…")
def _build_s1_superclusters(jaccard_threshold: float, *_mtimes):
    """
    Regroupe les sites S1 dont les résidus canonical actin se chevauchent
    (Jaccard ≥ jaccard_threshold) en super-clusters.

    Retourne :
      s1_to_super : dict  s1_site → supercluster_id (int)
      n_clusters  : int   nombre de super-clusters
    """
    import numpy as np
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform

    if not all(os.path.exists(f) for f in _S1_SUPER_FILES):
        return {}, 0

    df3 = pd.read_csv(_S1_SUPER_FILES[0])
    df1 = pd.read_csv(_S1_SUPER_FILES[1], usecols=[
                      "interaction_id", "chain_A_id", "chain_B_id"])
    df_all_s = pd.read_csv(_S1_SUPER_FILES[2])
    df_all_s["s1_actine"] = df_all_s["s1_actine"].fillna(False).astype(bool)
    df_all_s["s2_actine"] = df_all_s["s2_actine"].fillna(False).astype(bool)

    # Interactions hétéro (actin = S1) — filtered_all_data n'a pas d'interaction_id
    hetero_base = df_all_s[
        df_all_s["s1_actine"] & ~df_all_s["s2_actine"]
    ][["subunit_1", "subunit_2", "s1_binding_site_cluster_data_70"]].dropna(
        subset=["s1_binding_site_cluster_data_70"]
    )
    # Obtenir interaction_id via 1.interactions.csv
    hetero_s = hetero_base.merge(
        df1,
        left_on=["subunit_1", "subunit_2"],
        right_on=["chain_A_id", "chain_B_id"],
        how="inner",
    )[["interaction_id", "subunit_1", "s1_binding_site_cluster_data_70"]].rename(columns={
        "subunit_1": "actin_chain",
        "s1_binding_site_cluster_data_70": "s1_site",
    })

    # Joindre avec df3 : garder uniquement les résidus actin (chain == actin_chain)
    df3_canon = df3[df3["residue_number_canon_mafft"].notna()][
        ["interaction_id", "chain", "residue_number_canon_mafft"]
    ]
    merged_s = hetero_s.merge(df3_canon, on="interaction_id", how="inner")
    merged_s = merged_s[merged_s["chain"] == merged_s["actin_chain"]]

    # Profil de chaque S1 = ensemble de positions canonical contactées
    s1_profiles: dict[str, set] = {}
    for s1, grp in merged_s.groupby("s1_site"):
        s1_profiles[s1] = set(grp["residue_number_canon_mafft"].unique())

    s1_sites = sorted(s1_profiles.keys())
    n = len(s1_sites)
    if n == 0:
        return {}, 0
    if n == 1:
        return {s1_sites[0]: 1}, 1

    # Matrice de distance Jaccard (distance = 1 − similarité)
    dist_mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            a, b = s1_profiles[s1_sites[i]], s1_profiles[s1_sites[j]]
            union = len(a | b)
            jac = len(a & b) / union if union > 0 else 0.0
            dist_mat[i, j] = dist_mat[j, i] = 1.0 - jac

    Z = linkage(squareform(dist_mat), method="average")
    labels = fcluster(Z, t=1.0 - jaccard_threshold, criterion="distance")

    s1_to_super = {s1_sites[i]: int(labels[i]) for i in range(n)}

    # Sauvegarde CSV pour analyses extérieures
    _sc_rows = [
        {"s1_binding_site": s1, "supercluster": sc,
         "n_residues": len(s1_profiles[s1]),
         "jaccard_threshold": jaccard_threshold}
        for s1, sc in s1_to_super.items()
    ]
    pd.DataFrame(_sc_rows).sort_values("supercluster").to_csv(
        "data/filtered/s1_superclusters.csv", index=False
    )

    return s1_to_super, int(labels.max())
