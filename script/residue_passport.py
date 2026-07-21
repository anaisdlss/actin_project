"""Table « passeport résidu » de l'actin.

Fondation commune aux vues interactives (résidu cliquable, sélecteurs cluster /
ABP / paire d'ABP / séquence query). On joint, par position canonical (MAFFT) de
l'actin :

  * conservation évolutive + sensibilité mutationnelle ProteoCast
    (data/proteocast/conservation_vs_asa_per_position.csv)
  * les ABP qui contactent ce résidu, avec leur % ASA enfouie
    (3.interface_residues.csv, côté actin = chaîne A des interactions hétéro)
  * l'aa partenaire de l'ABP + type de contact (van der Waals / liaison H / pont
    salin) au niveau contact résidu-résidu (4.inter-residue_contacts.csv)
  * les clusters C70 et sites de liaison S1 d'appartenance
    (filtered_all_data.csv)

Convention (identique au reste de l'app, cf. data_analysis._build_abp_heatmap_data
et streamlit._abp_actin_footprint) : pour une interaction hétéro actin-ABP,
l'actin est le sous-unité 1 = chaîne A, l'ABP le sous-unité 2 = chaîne B.
"""

import os
import pandas as pd
import numpy as np
import streamlit as st

_PP_FILES = [
    "data/filtered/details/1.interactions.csv",
    "data/filtered/details/3.interface_residues.csv",
    "data/filtered/details/4.inter-residue_contacts.csv",
    "data/filtered/filtered_all_data.csv",
    "data/proteocast/conservation_vs_asa_per_position.csv",
]


def pp_mtimes():
    """Empreinte temporelle des fichiers sources (clé de cache)."""
    return tuple(os.path.getmtime(f) if os.path.exists(f) else 0.0
                 for f in _PP_FILES)


def _clean_abp_name(s: pd.Series) -> pd.Series:
    """Nom d'ABP propre : retire la parenthèse d'espèce, tronque à 50 car."""
    return (s.fillna("Unknown")
            .astype(str)
            .str.replace(r"\s*\(.*?\)", "", regex=True)
            .str.strip().str[:50])


def _join_semicol(vals) -> str:
    """Liste triée unique -> chaîne « a ; b ; c » (pour affichage / tri)."""
    u = sorted({str(v) for v in vals if pd.notna(v) and str(v) != "nan"})
    return " ; ".join(u)


@st.cache_data(show_spinner="Building the residue-passport table…")
def build_passport(_mtimes):
    """Renvoie un dict de DataFrames, ou None si les sources manquent.

    Clés :
      * ``pos``       : 1 ligne / position canonical (conservation + agrégats ABP)
      * ``res_abp``   : 1 ligne / (canon, ABP) agrégée (% ASA, n interactions, clusters)
      * ``res_long``  : 1 ligne / (interaction, résidu actin) — pour filtrer par cluster/ABP
      * ``contacts``  : 1 ligne / contact résidu-résidu (aa partenaire, type de liaison)
      * ``abp_list``  : ABP triés par nombre de résidus d'actin touchés
      * ``c70_list``  : clusters C70 actin-ABP présents
      * ``site_list`` : sites de liaison S1 présents
    """
    f_int, f_iface, f_con, f_all, f_cons = _PP_FILES
    if not all(os.path.exists(f) for f in (f_int, f_iface, f_all)):
        return None

    inter = pd.read_csv(f_int)[["interaction_id", "chain_A_id", "chain_B_id"]]
    alld = pd.read_csv(f_all, low_memory=False)

    # Interactions hétéro actin (S1) — ABP (S2, non-actin)
    het = alld[(alld["s1_actine"].fillna(False).astype(bool)) &
               (~alld["s2_actine"].fillna(False).astype(bool))].copy()
    het["abp"] = _clean_abp_name(het["subunit_2_title"])
    het = het[["subunit_1", "subunit_2", "abp",
               "cluster_data_70", "s1_binding_site_cluster_data_70"]]

    # interaction_id -> (abp, c70, site_s1, chaîne actin)
    m = inter.merge(
        het, left_on=["chain_A_id", "chain_B_id"],
        right_on=["subunit_1", "subunit_2"], how="inner",
    ).rename(columns={
        "cluster_data_70": "c70",
        "s1_binding_site_cluster_data_70": "s1_site",
    })
    m = m[["interaction_id", "chain_A_id", "abp", "c70", "s1_site"]]

    # ── res_abp : % ASA du résidu actin par ABP (source = 3.interface_residues) ─
    iface = pd.read_csv(f_iface)
    iface["asa_pct"] = pd.to_numeric(
        iface["buried_ASA_percent"].astype(str).str.replace("%", "", regex=False),
        errors="coerce")
    iface["canon"] = pd.to_numeric(
        iface["residue_number_canon_mafft"], errors="coerce")
    iface = iface[iface["canon"].notna() & iface["asa_pct"].notna()].copy()
    iface["canon"] = iface["canon"].astype(int)

    im = iface.merge(m, on="interaction_id", how="inner")
    # ne garder que le côté actin (chaîne A de l'interaction)
    im = im[im["chain"].str.lower() == im["chain_A_id"].str.lower()].copy()

    res_abp = (
        im.groupby(["canon", "abp"])
        .agg(asa_max=("asa_pct", "max"),
             asa_mean=("asa_pct", "mean"),
             n_int=("interaction_id", "nunique"),
             c70s=("c70", _join_semicol),
             sites=("s1_site", _join_semicol))
        .reset_index()
    )
    res_abp["asa_mean"] = res_abp["asa_mean"].round(1)
    res_abp["asa_max"] = res_abp["asa_max"].round(1)

    # ── res_long : 1 ligne par (interaction, résidu actin) ─────────────────────
    # Base commune pour filtrer par cluster C70 / ABP / site S1.
    res_long = im.rename(columns={"residue_name": "actin_aa"})[
        ["canon", "actin_aa", "abp", "c70", "s1_site", "asa_pct", "interaction_id"]
    ].copy()
    res_long["c70"] = res_long["c70"].astype(str)
    res_long["s1_site"] = res_long["s1_site"].astype(str)

    # ── contacts : détail résidu-résidu (aa partenaire ABP + type de liaison) ───
    contacts = None
    if os.path.exists(f_con):
        con = pd.read_csv(f_con)
        con = con.merge(m, on="interaction_id", how="inner")
        con = con[con["chain_A_id_x"].str.lower() ==
                  con["chain_A_id_y"].str.lower()] if "chain_A_id_x" in con else con
        con["canon"] = pd.to_numeric(con["residue_A_canon_mafft"], errors="coerce")
        con = con[con["canon"].notna()].copy()
        con["canon"] = con["canon"].astype(int)
        con["contact_type"] = con["contact_type"].fillna("van der Waals")
        con["asa_pct"] = pd.to_numeric(con["asa_pct_A"], errors="coerce")
        contacts = con.rename(columns={
            "residue_A_name": "actin_aa",
            "residue_A_sequence": "actin_resnum",
            "residue_B_name": "abp_aa",
            "residue_B_sequence": "abp_resnum",
            "contact_area": "area",
        })[["canon", "actin_aa", "actin_resnum", "abp", "abp_aa", "abp_resnum",
            "asa_pct", "area", "contact_type", "c70", "s1_site", "interaction_id"]]

    # ── pos : conservation par position + agrégats ABP ─────────────────────────
    if os.path.exists(f_cons):
        pos = pd.read_csv(f_cons)
    else:
        pos = pd.DataFrame({"canon": sorted(res_abp["canon"].unique())})
    pos["canon"] = pd.to_numeric(pos["canon"], errors="coerce")
    pos = pos[pos["canon"].notna()].copy()
    pos["canon"] = pos["canon"].astype(int)

    # aa canonical de l'actin par position (depuis les résidus d'interface)
    aa_by_canon = (im.dropna(subset=["residue_name"])
                   .drop_duplicates("canon")
                   .set_index("canon")["residue_name"])
    pos["actin_aa"] = pos["canon"].map(aa_by_canon)

    # agrégats ABP par position
    abp_by_canon = res_abp.groupby("canon")["abp"].apply(
        lambda s: sorted(set(s)))
    asa_max_by_canon = res_abp.groupby("canon")["asa_max"].max()
    c70_by_canon = im.groupby("canon")["c70"].apply(_join_semicol)
    site_by_canon = im.groupby("canon")["s1_site"].apply(_join_semicol)

    pos["abp_list"] = pos["canon"].map(abp_by_canon)
    pos["abp_list"] = pos["abp_list"].apply(
        lambda v: v if isinstance(v, list) else [])
    pos["n_abp"] = pos["abp_list"].apply(len)
    pos["asa_max"] = pos["canon"].map(asa_max_by_canon)
    pos["c70_list"] = pos["canon"].map(c70_by_canon).fillna("")
    pos["site_list"] = pos["canon"].map(site_by_canon).fillna("")

    # listes de référence pour les sélecteurs
    abp_list = (res_abp.groupby("abp")["canon"].nunique()
                .sort_values(ascending=False).index.tolist())
    c70_list = sorted(im["c70"].dropna().astype(str).unique())
    site_list = sorted(im["s1_site"].dropna().astype(str).unique())

    return {
        "pos": pos,
        "res_abp": res_abp,
        "res_long": res_long,
        "contacts": contacts,
        "abp_list": abp_list,
        "c70_list": c70_list,
        "site_list": site_list,
    }
