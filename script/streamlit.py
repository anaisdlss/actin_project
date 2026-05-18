# streamlit run script/streamlit.py
import streamlit as st
import pandas as pd
import networkx as nx
import py3Dmol
from Bio import AlignIO
from pyvis.network import Network
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import ast
import csv
import os
import re
import subprocess
import sys
import numpy as np
from collections import Counter as _Counter
from pathlib import Path as _Path
import matplotlib
matplotlib.use("Agg")

_GRAPH_COMPONENT = st.components.v1.declare_component(
    "actin_graph",
    path=str(_Path(__file__).parent / "graph_component"),
)

# ── Bipartite network helpers ─────────────────────────────────────────────────
# Incrémenter pour invalider le cache en mémoire après un correctif logique
_BIP_CACHE_VERSION = 3

_ACTIN_COLOR_NET = "#5B9BD5"
_ABP_PALETTE_NET = [
    "#E8735A", "#F0A500", "#7DBF6E", "#C580C3", "#70B8D4",
    "#A0522D", "#DB7093", "#8FBC8F", "#FF7F50", "#9370DB",
    "#20B2AA", "#DAA520", "#BA55D3", "#2E8B57", "#CD853F",
]
_BIPARTITE_FILES = [
    "data/filtered/details/3.interface_residues.csv",
    "data/filtered/details/1.interactions.csv",
    "data/filtered/filtered_all_data.csv",
    "data/filtered/patches_infos_s1_binding_site.csv",
    "data/filtered/patches_infos_cluster_data_70.csv",
    "data/filtered/proteins_per_pdb.csv",
    "data/filtered/details/4.inter-residue_contacts.csv",
    "data/filtered/binding_site_roles.csv",
]


def _bip_mtimes():
    return tuple(os.path.getmtime(f) if os.path.exists(f) else 0.0 for f in _BIPARTITE_FILES)


@st.cache_data(show_spinner="Chargement données réseau…")
def _load_bipartite_base(_v, *_mtimes):
    """Load & pre-process all data for S1 bipartite networks (cached per file state)."""
    import re as _re
    from collections import defaultdict as _dd

    if not all(os.path.exists(f) for f in _BIPARTITE_FILES):
        return None

    df_res3 = pd.read_csv(_BIPARTITE_FILES[0])
    df_int_b = pd.read_csv(_BIPARTITE_FILES[1])
    df_all_b = pd.read_csv(_BIPARTITE_FILES[2])
    df_s1_b = pd.read_csv(_BIPARTITE_FILES[3])
    df_c70_b = pd.read_csv(_BIPARTITE_FILES[4])
    df_pp_b = pd.read_csv(_BIPARTITE_FILES[5])

    df_res3["buried_ASA_percent"] = pd.to_numeric(
        df_res3["buried_ASA_percent"].astype(
            str).str.replace("%", "", regex=False),
        errors="coerce",
    )

    s1_chain = df_int_b.set_index("interaction_id")["chain_A_id"].str.lower()
    s2_chain = df_int_b.set_index("interaction_id")["chain_B_id"].str.lower()
    actin_ch = set(df_pp_b[df_pp_b["is_actin"]]["chain"].str.lower())
    homo_ids = set(
        df_int_b[df_int_b["chain_B_id"].str.lower().isin(actin_ch)]["interaction_id"])

    mc = (
        df_int_b.merge(
            df_all_b[["subunit_1", "subunit_2",
                      "s1_binding_site_cluster_data_70",
                      "s2_binding_site_cluster_data_70"]],
            left_on=["chain_A_id", "chain_B_id"],
            right_on=["subunit_1", "subunit_2"], how="inner",
        )[["interaction_id",
            "s1_binding_site_cluster_data_70",
            "s2_binding_site_cluster_data_70"]]
        .drop_duplicates("interaction_id")
        .set_index("interaction_id")
    )

    # ── Correction swap S1/S2 : même logique que le réseau (ancres étendues via CSV) ──
    _bsr = (pd.read_csv("data/filtered/binding_site_roles.csv")
            if os.path.exists("data/filtered/binding_site_roles.csv") else None)
    _FIXED_S1_B = (set(_bsr[_bsr["role"] == "S1"]["patch"].astype(str)) if _bsr is not None
                   else {"6685_1", "6685_3"})
    _FIXED_S2_B = (set(_bsr[_bsr["role"] == "S2"]["patch"].astype(str)) if _bsr is not None
                   else {"6685_2", "6685_4"})
    _s1p_b = mc["s1_binding_site_cluster_data_70"].astype(str)
    _s2p_b = mc["s2_binding_site_cluster_data_70"].astype(str)

    _homo_df2 = df_all_b[df_all_b["s2_actine"].fillna(False) == True]
    _gs1b = _homo_df2["s1_binding_site_cluster_data_70"].astype(str).value_counts()
    _gs2b = _homo_df2["s2_binding_site_cluster_data_70"].astype(str).value_counts()
    _can_s1b = {p: int(_gs1b.get(p, 0)) >= int(_gs2b.get(p, 0))
                for p in set(_gs1b.index) | set(_gs2b.index)}

    _prb: dict = {p: "S1" for p in _FIXED_S1_B}
    _prb.update({p: "S2" for p in _FIXED_S2_B})
    _chg = True
    while _chg:
        _chg = False
        for _iid in mc.index:
            _c1, _c2 = _s1p_b[_iid], _s2p_b[_iid]
            if _c1 in _prb and _c2 not in _prb:
                _prb[_c2] = "S2" if _prb[_c1] == "S1" else "S1"
                _chg = True
            elif _c2 in _prb and _c1 not in _prb:
                _prb[_c1] = "S1" if _prb[_c2] == "S2" else "S2"
                _chg = True
    for _iid in mc.index:
        for _pp in (_s1p_b[_iid], _s2p_b[_iid]):
            if _pp not in _prb:
                _prb[_pp] = "S1" if _can_s1b.get(_pp, True) else "S2"

    # Patch corrigé par rôle + s1_chain/s2_chain corrigés
    _s1pat, _s2pat = {}, {}
    _swbip = set()
    for _iid in mc.index:
        _c1, _c2 = _s1p_b[_iid], _s2p_b[_iid]
        if _prb.get(_c1) == "S1":
            _s1pat[_iid], _s2pat[_iid] = _c1, _c2
        elif _prb.get(_c2) == "S1":
            _s1pat[_iid], _s2pat[_iid] = _c2, _c1
            _swbip.add(_iid)
        else:
            _s1pat[_iid], _s2pat[_iid] = _c1, _c2

    _s1c, _s2c = s1_chain.copy(), s2_chain.copy()
    for _iid in _swbip:
        if _iid in _s1c.index and _iid in _s2c.index:
            _s1c[_iid], _s2c[_iid] = _s2c[_iid], _s1c[_iid]
    s1_chain, s2_chain = _s1c, _s2c
    # ──────────────────────────────────────────────────────────────────────────

    tmp = df_res3.copy()
    tmp["s1_chain"] = tmp["interaction_id"].map(s1_chain)
    tmp["s2_chain"] = tmp["interaction_id"].map(s2_chain)

    s1r = tmp[
        tmp["residue_number_canon_mafft"].notna() &
        (tmp["chain"].str.lower() == tmp["s1_chain"])
    ].copy()
    s1r["canon"] = s1r["residue_number_canon_mafft"].astype(int)
    s1r["patch"] = s1r["interaction_id"].map(_s1pat)

    s2r = tmp[
        tmp["interaction_id"].isin(homo_ids) &
        tmp["residue_number_canon_mafft"].notna() &
        (tmp["chain"].str.lower() == tmp["s2_chain"])
    ].copy()
    s2r["canon"] = s2r["residue_number_canon_mafft"].astype(int)
    s2r["patch"] = s2r["interaction_id"].map(_s2pat)

    df_res = pd.concat([s1r, s2r], ignore_index=True)
    df_res["buried_ASA_percent"] = pd.to_numeric(
        df_res["buried_ASA_percent"], errors="coerce")

    def _parse(s):
        s = _re.sub(r"np\.int64\((\d+)\)", r"\1", str(s))
        return {int(x) for x in _re.findall(r"\d+", s)}

    id_to_c70 = {}
    for _, r in df_c70_b.iterrows():
        for iid in _parse(r["ids_interactions"]):
            id_to_c70[iid] = r["patch"]

    df_int_meta = (
        df_int_b.merge(
            df_all_b[["subunit_1", "subunit_2",
                      "subunit_1_title", "subunit_2_title",
                      "s1_actine", "s2_actine"]],
            left_on=["chain_A_id", "chain_B_id"],
            right_on=["subunit_1", "subunit_2"], how="left",
        )
        .drop_duplicates("interaction_id")
        .set_index("interaction_id")
    )

    homo_c, hetero_c = [], []
    for _, row in df_s1_b.iterrows():
        ids = _parse(row["ids_interactions"])
        itype = str(row.get("interaction_type", "")).strip().lower()
        n = int(row["n_interactions"])
        entry = {"patch": row["patch"], "type": itype, "iids": ids, "n": n}
        if itype == "homo":
            homo_c.append(entry)
        elif itype == "mixed":
            hi = ids & homo_ids
            he = ids - homo_ids
            if hi:
                homo_c.append(
                    {**entry, "iids": hi, "n": len(hi), "type": "homo"})
            if he:
                hetero_c.append(
                    {**entry, "iids": he, "n": len(he), "type": "hetero"})
        else:
            hetero_c.append(entry)

    s2p2i = (
        df_all_b[df_all_b["s2_actine"] == True]
        .merge(df_int_b[["interaction_id", "chain_A_id", "chain_B_id"]],
               left_on=["subunit_1", "subunit_2"],
               right_on=["chain_A_id", "chain_B_id"], how="left")
        .dropna(subset=["interaction_id"])
        .groupby("s2_binding_site_cluster_data_70")["interaction_id"]
        .apply(set).to_dict()
    )
    for c in homo_c:
        c["iids"] |= s2p2i.get(c["patch"], set())
        c["n"] = len(c["iids"])

    p2c = _dd(list)
    for c in homo_c + hetero_c:
        p2c[c["patch"]].append({**c, "iids": list(c["iids"])})

    return df_res, df_int_meta, id_to_c70, dict(p2c)


@st.cache_data(show_spinner="Génération du réseau…")
def _build_bipartite_html(patch, _v, *_mtimes):
    """Build PyVis interactive network HTML for one S1 patch (physics layout, white bg).

    Returns (html_str, n_residues, n_proteins, n_total) or (None, 0, 0, 0).
    """
    import matplotlib.colors as _mc
    import matplotlib as _mpl_inner  # noqa: F841

    bip = _load_bipartite_base(_v, *_mtimes)
    if bip is None:
        return None, 0, 0, 0
    df_res, df_int_meta, id_to_c70, patch_clusters = bip

    clusters = patch_clusters.get(patch, [])
    if not clusters:
        return None, 0, 0, 0

    all_iids = set()
    for c in clusters:
        all_iids |= set(c["iids"])
    n_total = sum(c["n"] for c in clusters)

    sub = df_res[
        (df_res["interaction_id"].isin(all_iids)) &
        (df_res["patch"] == patch)
    ].copy()
    if sub.empty:
        return None, 0, 0, 0

    partner_map = {}
    for iid in all_iids:
        if iid not in df_int_meta.index:
            continue
        m = df_int_meta.loc[iid]
        s2t = m.get("subunit_2_title")
        s2a = bool(m.get("s2_actine", False))
        partner_map[iid] = "Actine" if s2a else (
            str(s2t) if pd.notna(s2t) else None)

    sub["is_s1_res"] = sub["chain"].str.lower() == sub["s1_chain"].str.lower()
    sub["partner"] = np.where(
        sub["is_s1_res"],
        sub["interaction_id"].map(partner_map),
        "Actine",
    )
    sub = sub[sub["partner"].notna()].copy()
    if sub.empty:
        return None, 0, 0, 0

    sub["c70"] = sub["interaction_id"].map(id_to_c70)
    # Nombre total d'iids par C70 (depuis all_iids — inclut les interactions sans résidu = 0%)
    c70_total: dict = {}
    for iid in all_iids:
        c70 = id_to_c70.get(iid)
        if c70:
            c70_total[c70] = c70_total.get(c70, 0) + 1
    # ASA équitable : sum/n_total par C70, puis moyenne sur TOUS les C70 du patch
    # (les C70 sans le résidu contribuent 0 → pénalisent les résidus absents)
    sub_c70 = sub[sub["c70"].notna()]
    if not sub_c70.empty and c70_total:
        c70_profiles = []
        for c70_val, n in c70_total.items():
            c70_sub = sub_c70[sub_c70["c70"] == c70_val]
            if not c70_sub.empty:
                area = (
                    c70_sub.groupby(["interaction_id", "canon"])[
                        "buried_ASA_percent"]
                    .max()
                    .groupby(level="canon").sum() / n
                )
                c70_profiles.append(area)
            # C70 sans aucun résidu → contribution 0 pour tous les résidus (pas ajouté)
        n_c70_total = len(c70_total)
        # sum() / n_c70_total inclut implicitement les 0 des C70 absents
        asa_c70 = pd.concat(c70_profiles).groupby(level=0).sum() / n_c70_total
    else:
        asa_c70 = pd.Series(dtype=float)
    # fallback : max par iid puis sum/n_total pour les résidus sans C70 mappé
    asa_raw = (
        sub.groupby(["interaction_id", "canon"])["buried_ASA_percent"]
        .max()
        .groupby(level="canon").sum() / len(all_iids)
    )

    edge_df = (
        sub.groupby(["canon", "partner"])["interaction_id"]
        .nunique()
        .reset_index(name="count")
    )
    all_residues = sorted(edge_df["canon"].unique().tolist())
    asa_mean = asa_c70.reindex(all_residues).fillna(asa_raw)
    prot_max_per_edge = edge_df.groupby("partner")["count"].transform("max")
    edge_df["norm"] = edge_df["count"] / prot_max_per_edge

    all_proteins = sorted(edge_df["partner"].unique().tolist())
    abp_colors = {}
    cidx = 0
    for p in all_proteins:
        if p != "Actine":
            abp_colors[p] = _ABP_PALETTE_NET[cidx % len(_ABP_PALETTE_NET)]
            cidx += 1

    cmap_asa = _mpl_inner.colormaps["YlOrRd"]
    norm_asa = _mc.Normalize(vmin=0, vmax=100)

    def _hexc(v):
        r, g, b, _ = cmap_asa(norm_asa(float(v)))
        return "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))

    n_res = len(all_residues)
    n_prot = len(all_proteins)
    h_px = 520

    net = Network(height=f"{h_px}px", width="100%",
                  bgcolor="#ffffff", font_color="#222")

    # Fréquences = nb d'iids distincts où le résidu/protéine apparaît au moins une fois
    res_freq = sub.groupby("canon")["interaction_id"].nunique()
    prot_freq = sub.groupby("partner")["interaction_id"].nunique()
    res_max = res_freq.max() if res_freq.max() > 0 else 1
    prot_max = prot_freq.max() if prot_freq.max() > 0 else 1

    # Distribution des acides aminés par position canonique
    # Pour les clusters homo, is_s1_res peut être False pour tous → on prend tout sub
    # (S1 et S2 sont les deux actin, même séquence, mêmes AA aux mêmes positions)
    _s1_sub_aa = sub.drop_duplicates(["interaction_id", "canon"])
    _aa_counts = (
        _s1_sub_aa.groupby(["canon", "residue_name"])["interaction_id"]
        .nunique().reset_index(name="cnt")
    )
    _aa_counts["pct"] = _aa_counts["cnt"] / n_total * 100
    _aa_counts["letter"] = _aa_counts["residue_name"].str.upper().fillna("?")
    _aa_by_pos = {
        pos: grp.sort_values("pct", ascending=False)
        for pos, grp in _aa_counts.groupby("canon")
    }

    # Residue nodes — couleur = ASA buried, taille = nb interactions
    # widthConstraint force la taille pour shape="circle" (size est ignoré par vis.js)
    for pos in all_residues:
        asa_v = float(asa_mean[pos]) if pos in asa_mean.index else 0.0
        bg = _hexc(asa_v)
        tc = "#222" if asa_v < 55 else "#fff"
        freq = int(res_freq.get(pos, 1))
        # 36–78 px (plus grand pour afficher pos + AA)
        sz = 36 + int(42 * freq / res_max)

        _pos_aa = _aa_by_pos.get(pos, pd.DataFrame())
        _top_letter = _pos_aa.iloc[0]["letter"] if not _pos_aa.empty else "?"
        _top_pct = _pos_aa.iloc[0]["pct"] if not _pos_aa.empty else 0.0
        _aa_tip = " · ".join(
            f"{r['letter']} : {r['pct']:.1f} %"
            for _, r in _pos_aa.iterrows()
        ) if not _pos_aa.empty else ""

        net.add_node(
            f"r{pos}", label=f"{_top_letter}{pos}",
            color={"background": bg, "border": "#888",
                   "highlight": {"background": bg, "border": "#E05000"},
                   "hover":     {"background": bg, "border": "#E05000"}},
            shape="circle",
            widthConstraint={"minimum": sz, "maximum": sz},
            font={"color": tc, "multi": False},
            title=(
                f"{_top_letter}{pos}\n"
                f"ASA buried : {asa_v:.1f} % · {freq} interactions"
                + (f"\nAA : {_aa_tip}" if _aa_tip and len(_pos_aa) > 1 else "")
            ),
            borderWidth=1.5, borderWidthSelected=2.5,
        )

    # Protein nodes — rectangle fixe
    PROT_W = 120
    PROT_H = 40
    prot_degree = edge_df.groupby("partner")["canon"].nunique()
    for prot in all_proteins:
        col = "#E03030" if prot == "Actine" else abp_colors.get(prot, "#888")
        freq = int(prot_freq.get(prot, 1))
        net.add_node(
            f"p_{prot}", label=prot,
            color={"background": col, "border": "#444",
                   "highlight": {"background": col, "border": "#E05000"},
                   "hover":     {"background": col, "border": "#E05000"}},
            shape="box",
            widthConstraint={"minimum": PROT_W},
            heightConstraint={"minimum": PROT_H},
            font={"color": "#ffffff"},
            title=f"{prot} — {freq} interactions · {int(prot_degree.get(prot, 0))} résidus connectés",
            borderWidth=1.5, borderWidthSelected=2.5,
        )

    # Edges — width per-protein normalized, smooth dynamic curves
    for _, row in edge_df.iterrows():
        nv = float(row["norm"])
        cnt = int(row["count"])
        net.add_edge(
            f"r{row['canon']}", f"p_{row['partner']}",
            width=0.5 + 4.5 * nv,
            color={"color":     "rgba(80,120,200,0.30)",
                   "highlight": "rgba(220,80,0,0.90)",
                   "hover":     "rgba(220,80,0,0.65)"},
            title=f"{cnt} interactions",
            smooth={"enabled": True, "type": "dynamic"},
        )

    net.set_options("""{
      "layout": {
        "randomSeed": 42
      },
      "physics": {
        "solver": "forceAtlas2Based",
        "forceAtlas2Based": {
          "gravitationalConstant": -60,
          "centralGravity": 0.01,
          "springLength": 100,
          "springConstant": 0.15,
          "damping": 0.5,
          "avoidOverlap": 1.0
        },
        "maxVelocity": 60,
        "timestep": 0.35,
        "stabilization": {
          "enabled": true,
          "iterations": 500,
          "updateInterval": 25,
          "fit": true
        },
        "minVelocity": 0.5
      },
      "interaction": {
        "hover": true, "tooltipDelay": 80,
        "zoomView": true, "dragView": true, "dragNodes": true,
        "navigationButtons": false, "keyboard": false
      },
      "edges": {
        "smooth": {"enabled": true, "type": "dynamic"},
        "selectionWidth": 2
      },
      "nodes": {
        "font": {"size": 18, "face": "monospace"},
        "shadow": {"enabled": true, "color": "rgba(0,0,0,0.12)", "size": 8, "x": 2, "y": 3}
      }
    }""")

    html = net.generate_html()

    # Légende — Actine uniquement si présente (interactions homo)
    legend_rows = []
    if "Actine" in all_proteins:
        legend_rows.append(
            '<div style="margin:4px 0;display:flex;align-items:center">'
            '<span style="display:inline-block;width:14px;height:14px;background:#E03030;'
            'border-radius:3px;margin-right:8px;flex-shrink:0"></span>Actine</div>'
        )
    for p in all_proteins:
        if p == "Actine":
            continue
        col = abp_colors.get(p, "#888")
        legend_rows.append(
            '<div style="margin:4px 0;display:flex;align-items:center">'
            f'<span style="display:inline-block;width:14px;height:14px;background:{col};'
            f'border-radius:3px;margin-right:8px;flex-shrink:0"></span>{p}</div>'
        )
    legend_html = (
        '<div style="position:fixed;top:10px;left:10px;'
        'background:rgba(255,255,255,0.96);border:1px solid #dde;border-radius:10px;'
        'padding:12px 14px;font-family:\'Segoe UI\',sans-serif;font-size:11px;color:#333;'
        'z-index:999;max-width:240px;max-height:82vh;overflow-y:auto;'
        'box-shadow:0 2px 12px rgba(0,0,0,0.10);">'
        '<div style="font-weight:700;color:#555;margin-bottom:8px;font-size:12px">Partenaires</div>'
        + "".join(legend_rows)
        + '<div style="margin-top:12px;font-weight:700;color:#555;font-size:11px">Résidus — ASA buried</div>'
        '<div style="background:linear-gradient(to right,#FFFFCC,#FD8D3C,#800026);'
        'height:10px;border-radius:4px;margin:5px 0 2px"></div>'
        '<div style="display:flex;justify-content:space-between;font-size:9px;color:#999">'
        '<span>0 %</span><span>50 %</span><span>100 %</span></div>'
        '</div>'
    )
    html = html.replace("</body>", legend_html + "\n</body>")

    # Inject: immediate fit + freeze+fit after stabilisation
    inject_js = (
        "setTimeout(function(){ network.fit(); }, 80);\n"
        "network.once('stabilizationIterationsDone', function(){\n"
        "  network.setOptions({physics:{enabled:false}});\n"
        "  setTimeout(function(){ network.fit({animation:{duration:300}}); }, 50);\n"
        "});\n"
    )
    for pat in [
        "network = new vis.Network(container, data, options);",
        "var network = new vis.Network(container, data, options);",
    ]:
        if pat in html:
            html = html.replace(pat, pat + "\n" + inject_js)
            break

    return html, n_res, n_prot, n_total


@st.cache_data(show_spinner="Chargement contacts inter-résidus…")
def _load_res4(_v, *_mtimes):
    """Charge et pre-traite 4.inter-residue_contacts.csv une seule fois."""
    p = _BIPARTITE_FILES[6]
    if not os.path.exists(p):
        return None
    df = pd.read_csv(p)
    df["contact_area"] = pd.to_numeric(
        df["contact_area"], errors="coerce").fillna(0.0)
    df["residue_B_canon_mafft"] = pd.to_numeric(
        df["residue_B_canon_mafft"], errors="coerce")
    df["residue_B_sequence"] = pd.to_numeric(
        df["residue_B_sequence"], errors="coerce")
    return df


@st.cache_data(show_spinner="Génération réseau C70…")
def _build_bipartite_c70_html(patch_c70, bipartite, _v, *_mtimes):
    """Réseau bipartite interactif pour un patch C70 : résidus actine (S1) ↔ résidus ABP (S2)."""
    import matplotlib.colors as _mc
    import matplotlib as _mpl_c70
    import re as _re_c70

    if not all(os.path.exists(f) for f in _BIPARTITE_FILES):
        return None, 0, 0, 0

    df_int = pd.read_csv(_BIPARTITE_FILES[1])
    df_all = pd.read_csv(_BIPARTITE_FILES[2])
    df_c70b = pd.read_csv(_BIPARTITE_FILES[4])

    def _p(s):
        s = _re_c70.sub(r"np\.int64\((\d+)\)", r"\1", str(s))
        return {int(x) for x in _re_c70.findall(r"\d+", s)}

    rows_c70 = df_c70b[df_c70b["patch"].astype(str) == str(patch_c70)]
    if rows_c70.empty:
        return None, 0, 0, 0
    all_iids: set = set()
    for _, r in rows_c70.iterrows():
        all_iids |= _p(r["ids_interactions"])
    if not all_iids:
        return None, 0, 0, 0
    n_total = len(all_iids)

    s1_chain = df_int.set_index("interaction_id")["chain_A_id"].str.lower()
    s2_chain = df_int.set_index("interaction_id")["chain_B_id"].str.lower()

    df_all["s2_actine"] = df_all["s2_actine"].fillna(False)

    # ── Orientation par propagation depuis ancres fixes ──
    # Priorité 1 : ancres fixes étendues chargées depuis binding_site_roles.csv
    # Priorité 2 : propagation locale au cluster
    # Fallback   : vote majoritaire global sur tous les homo
    _bsr2 = (pd.read_csv("data/filtered/binding_site_roles.csv")
             if os.path.exists("data/filtered/binding_site_roles.csv") else None)
    _FIXED_S1 = (set(_bsr2[_bsr2["role"] == "S1"]["patch"].astype(str)) if _bsr2 is not None
                 else {"6685_1", "6685_3"})
    _FIXED_S2 = (set(_bsr2[_bsr2["role"] == "S2"]["patch"].astype(str)) if _bsr2 is not None
                 else {"6685_2", "6685_4"})

    _homo_df = df_all[df_all["s2_actine"] == True]
    _gs1 = _homo_df["s1_binding_site_cluster_data_70"].astype(
        str).value_counts()
    _gs2 = _homo_df["s2_binding_site_cluster_data_70"].astype(
        str).value_counts()
    _is_canonical_s1 = {
        p: int(_gs1.get(p, 0)) >= int(_gs2.get(p, 0))
        for p in set(_gs1.index) | set(_gs2.index)
    }

    _cmap_raw = df_int.merge(
        df_all[["subunit_1", "subunit_2",
                "s1_binding_site_cluster_data_70",
                "s2_binding_site_cluster_data_70"]],
        left_on=["chain_A_id", "chain_B_id"],
        right_on=["subunit_1", "subunit_2"], how="left"
    ).drop_duplicates("interaction_id").set_index("interaction_id")
    _s1p_map = _cmap_raw["s1_binding_site_cluster_data_70"].astype(str)
    _s2p_map = _cmap_raw["s2_binding_site_cluster_data_70"].astype(str)

    # Initialiser avec les ancres fixes
    _patch_role: dict = {}
    for _p in _FIXED_S1:
        _patch_role[_p] = "S1"
    for _p in _FIXED_S2:
        _patch_role[_p] = "S2"

    # Propagation itérative locale au cluster
    _changed = True
    while _changed:
        _changed = False
        for _iid in all_iids:
            if _iid not in _cmap_raw.index:
                continue
            _cs1 = _s1p_map[_iid]
            _cs2 = _s2p_map[_iid]
            if _cs1 in _patch_role and _cs2 not in _patch_role:
                # cs1 connu : si c'est S1 → cs2 est S2, sinon cs2 est S1
                _patch_role[_cs2] = "S2" if _patch_role[_cs1] == "S1" else "S1"
                _changed = True
            elif _cs2 in _patch_role and _cs1 not in _patch_role:
                # cs2 connu : si c'est S2 → cs1 est S1, sinon cs1 est S2
                _patch_role[_cs1] = "S1" if _patch_role[_cs2] == "S2" else "S2"
                _changed = True

    # Fallback : vote majoritaire global pour patches non atteints par propagation
    for _iid in all_iids:
        if _iid not in _cmap_raw.index:
            continue
        for _pp in (_s1p_map[_iid], _s2p_map[_iid]):
            if _pp not in _patch_role:
                _patch_role[_pp] = "S1" if _is_canonical_s1.get(
                    _pp, True) else "S2"

    # Appliquer le swap pour les interactions dont le S1 actuel devrait être S2
    _s1_chain = s1_chain.copy()
    _s2_chain = s2_chain.copy()
    _swap_iids: set = set()
    for _iid in all_iids:
        if _iid not in _cmap_raw.index:
            continue
        _cs1 = _s1p_map[_iid]
        _cs2 = _s2p_map[_iid]
        if _patch_role.get(_cs1) == "S2" and _patch_role.get(_cs2) == "S1":
            if _iid in _s1_chain.index and _iid in _s2_chain.index:
                _s1_chain[_iid], _s2_chain[_iid] = _s2_chain[_iid], _s1_chain[_iid]
                _swap_iids.add(_iid)

    s1_chain, s2_chain = _s1_chain, _s2_chain

    # Patches S1 et S2 présents dans ce cluster (pour la légende)
    _cluster_s1_patches = sorted({
        p for _iid in all_iids if _iid in _cmap_raw.index
        for p in (_s1p_map[_iid], _s2p_map[_iid])
        if _patch_role.get(p) == "S1"
    })
    _cluster_s2_patches = sorted({
        p for _iid in all_iids if _iid in _cmap_raw.index
        for p in (_s1p_map[_iid], _s2p_map[_iid])
        if _patch_role.get(p) == "S2"
    })

    meta = (
        df_int.merge(
            df_all[["subunit_1", "subunit_2", "subunit_2_title", "s2_actine"]],
            left_on=["chain_A_id", "chain_B_id"],
            right_on=["subunit_1", "subunit_2"], how="left")
        .drop_duplicates("interaction_id").set_index("interaction_id")
    )

    def _partner(iid):
        if iid not in meta.index:
            return None
        r = meta.loc[iid]
        return "Actine" if r["s2_actine"] else (
            str(r["subunit_2_title"]) if pd.notna(r["subunit_2_title"]) else None)

    _AA3TO1 = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
        "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
        "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }

    def _aa1(aa3):
        return _AA3TO1.get(str(aa3).upper(), str(aa3)[:1])

    # Mapping interaction_id → couple (s1_patch × s2_patch)
    # Réutilise _cmap_raw ; couple inversé pour les interactions swappées
    _cm_s1p = _cmap_raw["s1_binding_site_cluster_data_70"].astype(str)
    _cm_s2p = _cmap_raw["s2_binding_site_cluster_data_70"].astype(str)
    _is_sw = _cmap_raw.index.isin(_swap_iids)
    iid_to_couple = pd.Series(
        np.where(_is_sw, _cm_s2p + "×" + _cm_s1p, _cm_s1p + "×" + _cm_s2p),
        index=_cmap_raw.index,
    )

    # Nombre d'interactions par couple dans ce cluster
    _couple_sizes = iid_to_couple[iid_to_couple.index.isin(
        all_iids)].value_counts()

    def _aa_majority_eq(grp):
        """AA le plus fréquent en moyenne équitable par couple."""
        majorities = [
            c_grp["residue_name"].value_counts().index[0]
            for _, c_grp in grp.groupby("couple")
        ]
        return pd.Series(majorities).value_counts().index[0]

    def _aa_dist_eq(grp):
        """Distribution AA en proportions moyennées équitablement par couple."""
        couple_dists = [
            c_grp["residue_name"].value_counts(normalize=True)
            for _, c_grp in grp.groupby("couple")
        ]
        if not couple_dists:
            return ""
        combined = pd.concat(couple_dists, axis=1).fillna(
            0).mean(axis=1).sort_values(ascending=False)
        total = combined.sum()
        if total == 0:
            return ""
        return " / ".join(
            f"{_aa1(aa)} ({100*v/total:.0f}%)" for aa, v in combined.items() if v > 0
        )

    # Charger table 4 (contient asa_pct_A et asa_pct_B)
    df_res4 = _load_res4(_v, *_mtimes)
    if df_res4 is None:
        return None, 0, 0, 0
    t4 = df_res4[df_res4["interaction_id"].isin(all_iids)].copy()
    t4["partner"] = t4["interaction_id"].map(_partner)
    t4["couple"] = t4["interaction_id"].map(iid_to_couple)

    # ── Correction du swap PPI3D dans table 4 ────────────────────────────────
    # PPI3D écrit parfois les chaînes dans l'ordre PDB (table 4) ≠ ordre du
    # résumé d'interaction (chain_A_id = S1).  On détecte le décalage en
    # comparant residue_A de la table 4 avec les résidus de chain_A dans la
    # table 3 (qui respecte toujours l'ordre du résumé).
    _df3_t4fix = pd.read_csv(_BIPARTITE_FILES[0])
    _df3_t4fix["residue_number_canon_mafft"] = pd.to_numeric(
        _df3_t4fix["residue_number_canon_mafft"], errors="coerce")
    _df3_t4fix = _df3_t4fix[_df3_t4fix["interaction_id"].isin(all_iids)].copy()
    _df3_t4fix["chain_lower"] = _df3_t4fix["chain"].str.lower()
    _int_ch = df_int[df_int["interaction_id"].isin(all_iids)].set_index("interaction_id")
    _t4_swapped_ppi3d: set = set()
    for _iid in all_iids:
        if _iid not in _int_ch.index:
            continue
        _cA = str(_int_ch.at[_iid, "chain_A_id"]).lower()
        _cB = str(_int_ch.at[_iid, "chain_B_id"]).lower()
        _sub3 = _df3_t4fix[_df3_t4fix["interaction_id"] == _iid]
        _t3A = set(_sub3[_sub3["chain_lower"] == _cA]["residue_number_canon_mafft"].dropna().astype(int))
        _t3B = set(_sub3[_sub3["chain_lower"] == _cB]["residue_number_canon_mafft"].dropna().astype(int))
        if not _t3A and not _t3B:
            continue
        _t4A = set(df_res4[df_res4["interaction_id"] == _iid]["residue_A_canon_mafft"].dropna().astype(int))
        if not _t4A:
            continue
        _ov_A = len(_t4A & _t3A) / max(len(_t4A | _t3A), 1)
        _ov_B = len(_t4A & _t3B) / max(len(_t4A | _t3B), 1)
        if _ov_B > _ov_A:
            _t4_swapped_ppi3d.add(_iid)
    if _t4_swapped_ppi3d:
        _sw_ppi3d = t4["interaction_id"].isin(_t4_swapped_ppi3d)
        t4.loc[_sw_ppi3d, ["residue_A_canon_mafft", "residue_B_canon_mafft"]] = \
            t4.loc[_sw_ppi3d, ["residue_B_canon_mafft", "residue_A_canon_mafft"]].values
        t4.loc[_sw_ppi3d, ["asa_pct_A", "asa_pct_B"]] = \
            t4.loc[_sw_ppi3d, ["asa_pct_B", "asa_pct_A"]].values
        t4.loc[_sw_ppi3d, ["residue_A_name", "residue_B_name"]] = \
            t4.loc[_sw_ppi3d, ["residue_B_name", "residue_A_name"]].values
    # ──────────────────────────────────────────────────────────────────────────

    # Pour les interactions swappées S1↔S2 : inverser A↔B dans table 4
    if _swap_iids:
        _sw = t4["interaction_id"].isin(_swap_iids)
        t4.loc[_sw, ["residue_A_canon_mafft", "residue_B_canon_mafft"]] = \
            t4.loc[_sw, ["residue_B_canon_mafft", "residue_A_canon_mafft"]].values
        t4.loc[_sw, ["asa_pct_A", "asa_pct_B"]] = \
            t4.loc[_sw, ["asa_pct_B", "asa_pct_A"]].values
        t4.loc[_sw, ["residue_A_name", "residue_B_name"]] = \
            t4.loc[_sw, ["residue_B_name", "residue_A_name"]].values

    # Stats S1 — moyenne équitable par couple : mean(sum/n_couple par canon)
    t4_s1 = t4[t4["residue_A_canon_mafft"].notna()].copy()
    if t4_s1.empty:
        return None, 0, 0, 0
    t4_s1["canon"] = t4_s1["residue_A_canon_mafft"].astype(int)
    t4_s1["residue_name"] = t4_s1["residue_A_name"]
    _s1_profiles = []
    for _c, _grp in t4_s1.groupby("couple"):
        _n = _grp["interaction_id"].nunique()
        _prof = (_grp.groupby(["interaction_id", "canon"])["asa_pct_A"].max()
                 .groupby(level="canon").sum() / float(_n))
        _s1_profiles.append(_prof)
    asa_s1 = (pd.concat(_s1_profiles, axis=1).fillna(0).mean(axis=1)
              if _s1_profiles else pd.Series(dtype=float))

    # Stats AA S1 depuis t4_s1
    s1_aa_majority = t4_s1.groupby("canon").apply(_aa_majority_eq)
    s1_aa_dist = t4_s1.groupby("canon").apply(_aa_dist_eq)

    # Fréquence équitable : fraction d'apparition moyennée sur TOUS les couples
    # (couple absent = 0, cohérent avec l'ASA qui fait fillna(0) avant mean)
    _n_couples_all = max(len(_couple_sizes), 1)

    def _equitable_freq(grp):
        total = sum(
            c_grp["interaction_id"].nunique(
            ) / max(_couple_sizes.get(couple, 1), 1)
            for couple, c_grp in grp.groupby("couple")
        )
        return total / _n_couples_all

    # Fréquence S1 depuis t4_s1 (même source que l'ASA) — cohérence taille + tooltip
    freq_s1 = t4_s1.groupby("canon").apply(_equitable_freq)
    s1_max = max(freq_s1.max(), 1e-9)
    n_couples_total = _couple_sizes.shape[0]
    s1_n_iids = t4_s1.groupby("canon")["interaction_id"].nunique()
    s1_n_couples = t4_s1.groupby("canon")["couple"].nunique()

    # Moyenne asa_pct_B — équitable par couple
    t4_s2 = t4[t4["partner"].notna() & t4["residue_B_canon_mafft"].notna()].copy()
    if t4_s2.empty:
        return None, 0, 0, 0
    t4_s2["s2_pos4"] = t4_s2["residue_B_canon_mafft"].astype(int)
    t4_s2["node_id4"] = t4_s2["partner"].astype(
        str).str[:15] + "_" + t4_s2["s2_pos4"].astype(str)
    t4_s2["residue_name"] = t4_s2["residue_B_name"]
    _s2_profiles = []
    for _c, _grp in t4_s2.groupby("couple"):
        _n = _grp["interaction_id"].nunique()
        _prof = (_grp.groupby(["interaction_id", "node_id4"])["asa_pct_B"].max()
                 .reset_index().groupby("node_id4")["asa_pct_B"].sum() / float(_n))
        _s2_profiles.append(_prof)
    avg_asa_s2 = (pd.concat(_s2_profiles, axis=1).fillna(0).mean(axis=1)
                  if _s2_profiles else pd.Series(dtype=float))
    # Fréquence S2 depuis t4_s2 (même source que l'ASA)
    freq_s2 = t4_s2.groupby("node_id4").apply(_equitable_freq)
    min_freq = 0.05
    top_s2 = freq_s2[freq_s2 >= min_freq].index.tolist()
    if not top_s2:
        top_s2 = freq_s2.nlargest(20).index.tolist()
    s2_ca_max = max(avg_asa_s2[avg_asa_s2.index.isin(
        top_s2)].max() if top_s2 else 1.0, 1.0)
    s2_n_iids = t4_s2.groupby("node_id4")["interaction_id"].nunique()
    s2_n_couples = t4_s2.groupby("node_id4")["couple"].nunique()

    # Stats AA S2 depuis t4_s2
    nd_to_pos = t4_s2.drop_duplicates(
        "node_id4").set_index("node_id4")["s2_pos4"]
    nd_to_partner = t4_s2.drop_duplicates(
        "node_id4").set_index("node_id4")["partner"]
    nd_to_majority = t4_s2.groupby("node_id4").apply(_aa_majority_eq)
    nd_to_aa_dist = t4_s2.groupby("node_id4").apply(_aa_dist_eq)
    nd_to_label = nd_to_majority.map(_aa1) + nd_to_pos.astype(str)

    # Arêtes : co-occurrence S1 canon × S2 node, épaisseur = nb couples, couleur = type contact
    s1_with_type = t4_s1[["interaction_id",
                          "canon", "contact_type", "couple"]].copy()
    s2_ids = t4_s2[t4_s2["node_id4"].isin(
        top_s2)][["interaction_id", "node_id4"]].copy()
    edges_raw = s1_with_type.merge(s2_ids, on="interaction_id")
    if edges_raw.empty:
        return None, 0, 0, 0
    edge_counts = (edges_raw.groupby(["canon", "node_id4"])["couple"]
                   .nunique().reset_index())
    edge_counts.columns = ["s1_canon", "s2_node", "n_couples"]
    edge_type = (edges_raw.groupby(["canon", "node_id4"])["contact_type"]
                 .agg(lambda x: x.dropna().value_counts().index[0]
                      if x.dropna().shape[0] > 0 else "").reset_index())
    edge_type.columns = ["s1_canon", "s2_node", "contact_type"]
    edge_df = edge_counts.merge(edge_type, on=["s1_canon", "s2_node"])

    connected_s1 = set(edge_df["s1_canon"])
    all_s1 = sorted(c for c in freq_s1.index if c in connected_s1)

    all_partners = sorted(set(
        t4_s2[t4_s2["node_id4"].isin(top_s2)]["partner"].dropna().unique()
    ))

    # Gradient S1 : YlOrRd (jaune → rouge) par ASA buried, normalisé au max du cluster
    s1_ca_max = max(float(asa_s1.max()) if not asa_s1.empty else 1.0, 1.0)
    cmap_s1 = _mpl_c70.colormaps["YlOrRd"]
    norm_s1 = _mc.Normalize(vmin=0, vmax=s1_ca_max)

    def _hex_s1(v):
        r, g, b, _ = cmap_s1(norm_s1(float(v)))
        return "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))

    # Gradient S2 hétéro : YlGn (jaune → vert) par aire de contact moyenne
    cmap_s2 = _mpl_c70.colormaps["YlGn"]
    norm_s2 = _mc.Normalize(vmin=0, vmax=s2_ca_max)

    def _hex_s2(ca_val):
        r, g, b, _ = cmap_s2(norm_s2(float(ca_val)))
        return "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))

    # Gradient S2 homo actine : rose (blanc rosé → rose vif → rose foncé)
    cmap_s2_homo = _mpl_c70.colors.LinearSegmentedColormap.from_list(
        "pink_grad", ["#FFF0F5", "#FF69B4", "#C71585"])
    norm_s2_homo = _mc.Normalize(vmin=0, vmax=s2_ca_max)

    def _hex_s2_homo(ca_val):
        r, g, b, _ = cmap_s2_homo(norm_s2_homo(float(ca_val)))
        return "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))

    # ── Positions bipartites (S1 gauche, S2 droite, tri pour min. croisements) ──
    n_s1 = len(all_s1)
    n_s2 = len(top_s2)

    if bipartite:
        # Trier S2 par position S1 moyenne connectée → minimise les croisements
        s2_mean_s1 = {
            nid: float(np.mean(edge_df[edge_df["s2_node"] == nid]["s1_canon"]))
            if not edge_df[edge_df["s2_node"] == nid].empty else 0.0
            for nid in top_s2
        }
        top_s2_ordered = sorted(top_s2, key=lambda nid: s2_mean_s1.get(nid, 0))
        h_net = max(700, max(n_s1, n_s2) * 55)
    else:
        top_s2_ordered = top_s2
        h_net = 600

    net = Network(height=f"{h_net}px", width="100%",
                  bgcolor="#ffffff", font_color="#222")

    # Taille commune S1 et S2 : même échelle de fréquence → même rayon visuel
    freq_s2_max = max(freq_s2[top_s2].max(), 1)
    _freq_max = max(s1_max, freq_s2_max)

    def _node_radius(freq):
        return 40 + int(30 * freq / _freq_max)

    # Nœuds S1
    for i, pos in enumerate(all_s1):
        asa_v = float(asa_s1[pos]) if pos in asa_s1.index else 0.0
        bg = _hex_s1(asa_v)
        tc = "#222" if asa_v < 55 else "#fff"
        _s1_freq = freq_s1.get(pos, 0.0)
        diam = _node_radius(_s1_freq) * 2
        _s1_aa = _aa1(s1_aa_majority.get(pos, "?"))
        _s1_dist = s1_aa_dist.get(pos, "")
        _s1_dist_str = f"AA : {_s1_dist}" if _s1_dist else ""
        _n_iids_s1 = int(s1_n_iids.get(pos, 0))
        _n_c_s1 = int(s1_n_couples.get(pos, 0))
        kwargs = {}
        if bipartite:
            y_pos = int((i - (n_s1 - 1) / 2) * (h_net / max(n_s1, 1)))
            kwargs = {"x": -450, "y": y_pos, "fixed": True, "physics": False}
        net.add_node(
            f"s1_{pos}", label=f"<b>{_s1_aa}{pos}</b>",
            color={"background": bg, "border": "#888",
                   "highlight": {"background": bg, "border": "#E05000"},
                   "hover":     {"background": bg, "border": "#E05000"}},
            shape="circle",
            widthConstraint={"minimum": diam, "maximum": diam},
            font={"color": tc},
            title=f"résidu {_s1_aa}{pos} · ASA buried : {asa_v:.1f} % · {_n_iids_s1} interactions · {_n_c_s1}/{n_couples_total} couples{_s1_dist_str}",
            borderWidth=1.5, borderWidthSelected=2.5,
            **kwargs,
        )

    # Nœuds S2
    for i, nid in enumerate(top_s2_ordered):
        partner = str(nd_to_partner.get(nid, "?"))
        label = str(nd_to_label.get(nid, nid))
        _s2_freq = freq_s2.get(nid, 0.0)
        ca_val = float(avg_asa_s2.get(nid, 0.0))
        aa_dist = str(nd_to_aa_dist.get(nid, ""))
        sz = _node_radius(_s2_freq)
        if partner == "Actine":
            col = _hex_s2_homo(ca_val)
            _bord_h = "#880044"
        else:
            col = _hex_s2(ca_val)
            _bord_h = "#007700"
        tc_s2 = "#222"
        _variants = f"AA : {aa_dist}" if aa_dist else ""
        _n_iids_s2 = int(s2_n_iids.get(nid, 0))
        _n_c_s2 = int(s2_n_couples.get(nid, 0))
        diam_s2 = sz * 2
        kwargs = {}
        if bipartite:
            y_pos = int((i - (n_s2 - 1) / 2) * (h_net / max(n_s2, 1)))
            kwargs = {"x": 450, "y": y_pos, "fixed": True, "physics": False}
        net.add_node(
            f"s2_{nid}", label=f"<b>{label}</b>",
            color={"background": col, "border": "#555",
                   "highlight": {"background": col, "border": _bord_h},
                   "hover":     {"background": col, "border": _bord_h}},
            shape="circle",
            widthConstraint={"minimum": diam_s2, "maximum": diam_s2},
            font={"color": tc_s2},
            title=f"{partner} — résidu {label} · % ASA moy : {ca_val:.1f} % · {_n_iids_s2} interactions · {_n_c_s2}/{n_couples_total} couples{_variants}",
            borderWidth=1.5,
            **kwargs,
        )

    def _edge_col(ct):
        ct_l = str(ct).lower()
        if "salt" in ct_l:
            return "rgba(210,175,0,0.55)"
        if "h-bond" in ct_l or "hydrogen" in ct_l:
            return "rgba(68,136,238,0.30)"
        return "rgba(50,50,50,0.20)"

    # Arêtes
    e_max = max(edge_df["n_couples"].max(), 1)
    for _, row in edge_df.iterrows():
        nc = int(row["n_couples"])
        ecol = _edge_col(row["contact_type"])
        net.add_edge(
            f"s1_{int(row['s1_canon'])}", f"s2_{row['s2_node']}",
            width=6.0 if (nc == n_couples_total and n_couples_total > 1) else 0.5,
            color={"color": ecol, "highlight": "#FF4400", "hover": "#FF4400"},
            title=f"{nc}/{n_couples_total} couples · {row['contact_type']}",
            smooth={"enabled": True, "type": "curvedCW" if bipartite else "dynamic",
                    "roundness": 0.15},
        )

    if bipartite:
        net.set_options("""{
          "layout": {"randomSeed": 42},
          "physics": {"enabled": false},
          "interaction": {"hover": true, "tooltipDelay": 80,
                          "zoomView": true, "dragView": true, "dragNodes": true},
          "edges": {"smooth": {"enabled": true, "type": "curvedCW", "roundness": 0.15},
                    "selectionWidth": 2},
          "nodes": {"font": {"face": "monospace", "size": 25, "multi": "html"},
                    "shadow": {"enabled": true, "color": "rgba(0,0,0,0.10)",
                               "size": 6, "x": 2, "y": 2}}
        }""")
    else:
        net.set_options("""{
          "layout": {"randomSeed": 42},
          "physics": {
            "solver": "forceAtlas2Based",
            "forceAtlas2Based": {
              "gravitationalConstant": -200,
              "centralGravity": 0.15,
              "springLength": 30,
              "springConstant": 0.4,
              "damping": 0.8,
              "avoidOverlap": 1.0
            },
            "maxVelocity": 100,
            "timestep": 0.25,
            "stabilization": {
              "enabled": true,
              "iterations": 1500,
              "updateInterval": 25,
              "fit": true
            },
            "minVelocity": 0.1
          },
          "interaction": {"hover": true, "tooltipDelay": 80,
                          "zoomView": true, "dragView": true, "dragNodes": true},
          "edges": {"smooth": {"enabled": true, "type": "dynamic"}, "selectionWidth": 2},
          "nodes": {"font": {"face": "monospace", "size": 25, "multi": "html"},
                    "shadow": {"enabled": true, "color": "rgba(0,0,0,0.10)",
                               "size": 6, "x": 2, "y": 2}}
        }""")

    html = net.generate_html()

    # Légende
    s2_partner_names = [p for p in all_partners if p != "Actine"]
    _has_homo = "Actine" in all_partners
    _has_hetero = bool(s2_partner_names)
    _leg = (
        '<div style="position:fixed;top:10px;left:10px;'
        'background:rgba(255,255,255,0.96);border:1px solid #dde;border-radius:10px;'
        'padding:12px 14px;font-family:\'Segoe UI\',sans-serif;font-size:11px;color:#333;'
        'z-index:999;max-width:230px;max-height:80vh;overflow-y:auto;'
        'box-shadow:0 2px 12px rgba(0,0,0,0.10);">'
        '<div style="font-weight:700;color:#555;margin-bottom:6px;font-size:12px">'
        'Résidus S1 (actine) — % ASA buried moy</div>'
        '<div style="background:linear-gradient(to right,#FFFFCC,#FD8D3C,#800026);'
        'height:10px;border-radius:4px;margin:4px 0 2px"></div>'
        f'<div style="display:flex;justify-content:space-between;font-size:9px;color:#999;'
        f'margin-bottom:10px"><span>0</span><span>max ({s1_ca_max:.1f} %)</span></div>'
    )
    if _has_homo:
        _leg += (
            '<div style="font-weight:700;color:#555;margin-bottom:6px;font-size:12px">'
            'Résidus S2 (actine homo) — % ASA buried moy</div>'
            '<div style="background:linear-gradient(to right,#FFF0F5,#FF69B4,#C71585);'
            'height:10px;border-radius:4px;margin:4px 0 2px"></div>'
            f'<div style="display:flex;justify-content:space-between;font-size:9px;color:#999;'
            f'margin-bottom:10px"><span>0</span><span>max ({s2_ca_max:.1f} %)</span></div>'
        )
    if _has_hetero:
        _leg += (
            '<div style="font-weight:700;color:#555;margin-bottom:6px;font-size:12px">'
            'Résidus S2 (ABP) — % ASA buried moy</div>'
            '<div style="background:linear-gradient(to right,#FFFFCC,#78C679,#006837);'
            'height:10px;border-radius:4px;margin:4px 0 2px"></div>'
            f'<div style="display:flex;justify-content:space-between;font-size:9px;color:#999;'
            f'margin-bottom:10px"><span>0</span><span>max ({s2_ca_max:.1f} %)</span></div>'
            + '<div style="font-weight:700;color:#555;margin-bottom:6px;font-size:11px">'
            'Partenaires</div>'
            + "".join(
                f'<div style="margin:3px 0;font-size:10px;color:#444">• {p[:35]}</div>'
                for p in s2_partner_names
            )
        )
    # Couples (s1_patch × s2_patch) présents dans ce cluster
    _couple_counts_leg = (
        iid_to_couple[iid_to_couple.index.isin(all_iids)]
        .value_counts()
        .sort_values(ascending=False)
    )
    _leg += (
        '<div style="border-top:1px solid #eee;margin:8px 0 6px"></div>'
        '<div style="font-weight:700;color:#555;margin-bottom:4px;font-size:12px">'
        'Combinaisons S1×S2</div>'
        + "".join(
            f'<div style="margin:2px 0;font-size:10px;color:#444">'
            f'{c} <span style="color:#999">({n})</span></div>'
            for c, n in _couple_counts_leg.items()
        )
    )
    # Binding sites S1 / S2
    _leg += (
        '<div style="border-top:1px solid #eee;margin:8px 0 6px"></div>'
        '<div style="font-weight:700;color:#555;margin-bottom:4px;font-size:12px">'
        'Binding sites S1</div>'
        + "".join(
            f'<div style="margin:2px 0;font-size:10px;color:#444">{p}</div>'
            for p in _cluster_s1_patches
        )
        + '<div style="font-weight:700;color:#555;margin:6px 0 4px;font-size:12px">'
        'Binding sites S2</div>'
        + "".join(
            f'<div style="margin:2px 0;font-size:10px;color:#444">{p}</div>'
            for p in _cluster_s2_patches
        )
    )
    legend_html = _leg + '</div>'
    html = html.replace("</body>", legend_html + "\n</body>")

    inject_js = (
        "setTimeout(function(){ network.fit(); }, 80);\n"
        "network.once('stabilizationIterationsDone', function(){\n"
        "  network.setOptions({physics:{enabled:false}});\n"
        "  setTimeout(function(){ network.fit({animation:{duration:300}}); }, 50);\n"
        "});\n"
    )
    for pat in [
        "network = new vis.Network(container, data, options);",
        "var network = new vis.Network(container, data, options);",
    ]:
        if pat in html:
            html = html.replace(pat, pat + "\n" + inject_js)
            break

    # ── Représentation 3D du couple représentatif ──────────────────────────────
    _html_3d = None
    try:
        from collections import defaultdict as _dd3
        _df8_3d = pd.read_csv("data/filtered/details/8.structures.csv")

        # Couple dominant dans ce cluster
        _dom_couple = _couple_sizes.idxmax()
        _dom_iids = [iid for iid in all_iids
                     if iid in iid_to_couple.index and iid_to_couple[iid] == _dom_couple]

        # Séquences pour ce couple (inverser si swappé)
        _seq_mrg = df_int[df_int["interaction_id"].isin(_dom_iids)].merge(
            df_all[["subunit_1", "subunit_2", "s1_sequence", "s2_sequence"]],
            left_on=["chain_A_id", "chain_B_id"],
            right_on=["subunit_1", "subunit_2"], how="left"
        ).drop_duplicates("interaction_id").set_index("interaction_id")
        for _iid in _dom_iids:
            if _iid in _swap_iids and _iid in _seq_mrg.index:
                _seq_mrg.loc[_iid, ["s1_sequence", "s2_sequence"]] = \
                    _seq_mrg.loc[_iid, ["s2_sequence", "s1_sequence"]].values
        _seq_mrg["_sp"] = (_seq_mrg["s1_sequence"].str.lower().fillna("") +
                           "|||" + _seq_mrg["s2_sequence"].str.lower().fillna(""))

        # Paire de séquences la plus fréquente → représentant
        _best_sp = _seq_mrg["_sp"].value_counts()
        _rep_iid3d, _pdb3d = None, None
        if not _best_sp.empty:
            for _riid in _seq_mrg[_seq_mrg["_sp"] == _best_sp.index[0]].index:
                _r8 = _df8_3d[_df8_3d["interaction_id"] == _riid]
                if _r8.empty:
                    continue
                _pp = str(_r8.iloc[0]["pairwise_pdb_file"])
                if _pp and _pp != "nan" and os.path.exists(_pp):
                    _rep_iid3d, _pdb3d = _riid, _pp
                    break

        if _rep_iid3d is not None:
            # Détecter le swap via s1_chain (plus fiable que _swap_iids pour les cas
            # où _patch_role n'a pas propagé le rôle d'un patch inconnu)
            _rep_s1_sub = s1_chain.get(_rep_iid3d, "").lower()
            _rep_cha_sub = (meta.loc[_rep_iid3d, "chain_A_id"].lower()
                            if _rep_iid3d in meta.index else "")
            if _rep_s1_sub and _rep_cha_sub:
                _is_sw3d = (_rep_s1_sub != _rep_cha_sub)
            else:
                _is_sw3d = _rep_iid3d in _swap_iids
            _chain_s1_3d = "B" if _is_sw3d else "A"
            _chain_s2_3d = "A" if _is_sw3d else "B"
            _sc_s1 = "residue_B_structure" if _is_sw3d else "residue_A_structure"
            _sc_s2 = "residue_A_structure" if _is_sw3d else "residue_B_structure"

            _is_homo3d = _partner(_rep_iid3d) == "Actine"

            _pdb_id_3d = _df8_3d.loc[_df8_3d["interaction_id"] == _rep_iid3d,
                                     "pdb_id"].iat[0]
            _same_pdb_iids = (
                set(_df8_3d.loc[_df8_3d["pdb_id"] == _pdb_id_3d, "interaction_id"])
                & all_iids
            )
            # Utiliser df_res4 (données brutes, non-swappées) pour mapper les
            # numéros de résidus physiques (chain A/B) directement à l'ASA.
            # Pas besoin de filtrer par _swap_iids : l'ASA raw est toujours
            # asa_pct_A = chaîne A physique, asa_pct_B = chaîne B physique.
            _t4r = df_res4[df_res4["interaction_id"].isin(_same_pdb_iids)]
            # Colonne ASA correspondant à la chaîne physique S1 / S2
            _asa_s1_col = "asa_pct_B" if _is_sw3d else "asa_pct_A"
            _asa_s2_col = "asa_pct_A" if _is_sw3d else "asa_pct_B"

            # B-factor = asa_pct réel de ce PDB (max par résidu si contacts multiples)
            _s1_bfac, _s2_bfac = {}, {}
            for _, _rw in _t4r.iterrows():
                _rs1 = _rw.get(_sc_s1)
                _v1 = _rw.get(_asa_s1_col)
                if pd.notna(_rs1) and pd.notna(_v1):
                    _k = int(_rs1)
                    if float(_v1) > _s1_bfac.get(_k, 0.0):
                        _s1_bfac[_k] = float(_v1)
                _rs2 = _rw.get(_sc_s2)
                _v2 = _rw.get(_asa_s2_col)
                if pd.notna(_rs2) and pd.notna(_v2):
                    _k2 = int(_rs2)
                    if float(_v2) > _s2_bfac.get(_k2, 0.0):
                        _s2_bfac[_k2] = float(_v2)

            # Réécrire les B-factors : 0 hors interface (→ 1er stop gradient),
            # valeur ASA pour interface — même convention que bfactor_cluster
            _pdb_lines = []
            for _line in open(_pdb3d):
                if _line.startswith(("ATOM", "HETATM")) and len(_line) > 26:
                    _ch = _line[21]
                    try:
                        _ri = int(_line[22:26].strip())
                    except ValueError:
                        _pdb_lines.append(_line)
                        continue
                    if _ch == _chain_s1_3d:
                        _bv = _s1_bfac.get(_ri, 0.0)
                    elif _ch == _chain_s2_3d:
                        _bv = _s2_bfac.get(_ri, 0.0)
                    else:
                        _bv = 0.0
                    _line = f"{_line[:60]}{_bv:6.2f}{_line[66:]}"
                _pdb_lines.append(_line)
            _pdb_mod = "".join(_pdb_lines)

            # Normalisation locale : max observé dans ce PDB (relatif à la structure)
            _s1_bfac_max = max(max(_s1_bfac.values()) if _s1_bfac else 0.0, 1.0)
            _s2_bfac_max = max(max(_s2_bfac.values()) if _s2_bfac else 0.0, 1.0)

            _ylord = ["#FFFFCC", "#FFF0A9", "#FEE186", "#FECA65", "#FDAA48",
                      "#FC8C3B", "#FC5A2D", "#EC2D21", "#D30F20", "#AF0026", "#800026"]
            _s1_sch = {"prop": "b", "gradient": "linear",
                       "colors": _ylord, "min": 0, "max": _s1_bfac_max}
            _s2_cols = (["#FFF0F5", "#FFD6E0", "#FFB6C1", "#FF85A1",
                         "#FF69B4", "#FF1493", "#C71585"]
                        if _is_homo3d else
                        ["#FFFFCC", "#D9F0A3", "#ADDD8E", "#78C679",
                         "#41AB5D", "#238443", "#006837"])
            _s2_sch = {"prop": "b", "gradient": "linear",
                       "colors": _s2_cols, "min": 0, "max": _s2_bfac_max}

            _v3 = py3Dmol.view(width=700, height=430)
            _v3.addModel(_pdb_mod, "pdb")
            _v3.setStyle({}, {})
            _v3.addSurface(py3Dmol.SES,
                           {"opacity": 1, "colorscheme": _s1_sch},
                           {"chain": _chain_s1_3d})
            _v3.addSurface(py3Dmol.SES,
                           {"opacity": 1, "colorscheme": _s2_sch},
                           {"chain": _chain_s2_3d})
            _v3.setBackgroundColor("white")
            _v3.zoomTo()
            _html_3d = _v3._make_html()
    except Exception:
        _html_3d = None

    return html, len(all_s1), len(top_s2), n_total, _html_3d


_BFACTOR_CLUSTER_DIR = "data/filtered/details/structures_files/bfactor_cluster"


@st.cache_data(show_spinner="Génération interface 3D S1…")
def _build_s1_3d_html(patch_s1, _v, *_mtimes):
    """Vue 3D pour un patch S1 binding site.

    Utilise le PDB pré-calculé bfactor_cluster/{patch}.pdb (chaîne A,
    b-factors = % ASA cluster-level). Retourne (html, s1_ca_max, 1.0, True) ou None.
    """
    pdb_path = os.path.join(_BFACTOR_CLUSTER_DIR, f"{patch_s1}.pdb")
    if not os.path.exists(pdb_path):
        return None
    try:
        pdb_content = open(pdb_path).read()
        # détecter la chaîne et les b-factors depuis le PDB
        chain_id = "A"
        bfacs = []
        for _ln in pdb_content.splitlines():
            if _ln.startswith("ATOM") and len(_ln) > 65:
                if chain_id == "A":
                    chain_id = _ln[21]  # colonne 22 PDB (0-indexed 21)
                try:
                    bfacs.append(float(_ln[60:66]))
                except ValueError:
                    pass
        s1_ca_max = max(bfacs) if bfacs else 1.0

        _ylord = ["#FFFFCC", "#FFF0A9", "#FEE186", "#FECA65", "#FDAA48",
                  "#FC8C3B", "#FC5A2D", "#EC2D21", "#D30F20", "#AF0026", "#800026"]
        _s1_sch = {"prop": "b", "gradient": "linear",
                   "colors": _ylord, "min": 0, "max": 100}

        _v3 = py3Dmol.view(width=500, height=450)
        _v3.addModel(pdb_content, "pdb")
        _v3.setStyle({}, {})
        _v3.addSurface(py3Dmol.SES,
                       {"opacity": 1, "colorscheme": _s1_sch},
                       {"chain": chain_id})
        _v3.setBackgroundColor("white")
        _v3.zoomTo({"chain": chain_id})
        _v3.zoom(0.7)
        import re as _re_fog
        _html3d = _v3._make_html()
        # disable fog — inject before first render() call
        _html3d_nofog = _re_fog.sub(
            r'(viewer_\w+)\.render\(\)',
            r'\1.setFogParameters({density:0});\1.render()',
            _html3d, count=1,
        )
        _html3d = _html3d_nofog if _html3d_nofog != _html3d else _html3d
        return _html3d, s1_ca_max, 1.0, True
    except Exception:
        return None


st.set_page_config(layout="wide", page_title="Analyse actine-ABP - PPI3D")

_CHAIN_ID_RE = re.compile(r'^[0-9][a-z0-9]{3}_[A-Za-z]')


def norm_chain_id(s):
    """PDB ID en minuscules, lettre de chaîne préservée. Ex: '6VEC_A' → '6vec_A'."""
    parts = str(s).split("_", 1)
    return parts[0].lower() + "_" + parts[1] if len(parts) == 2 else str(s).lower()


@st.cache_data
def _cluster_protein_index(aln_dir: str) -> dict:
    """Pour chaque cluster : liste de {name, ids} avec fallback sur le nom le plus commun."""
    result = {}
    for f in sorted(_Path(aln_dir).glob("*.aln")):
        try:
            aln = AlignIO.read(str(f), "fasta")
            entries = []
            for r in aln:
                subunits = r.id.split("|")
                desc = r.description
                if desc.startswith(r.id):
                    desc = desc[len(r.id):].strip()
                names = [p.strip() for p in desc.split(" | ") if p.strip()]
                entries.append({"ids": subunits, "names": names})

            # Noms "réels" (pas des IDs de chaîne PDB) présents dans l'alignement
            real_names = [
                n for e in entries for n in e["names"]
                if not _CHAIN_ID_RE.match(n)
            ]
            fallback = _Counter(real_names).most_common(1)[
                0][0] if real_names else ""

            result[f.stem] = [
                {
                    "name": (e["names"][0] if e["names"] else "") or fallback or e["ids"][0],
                    "ids": e["ids"],
                }
                for e in entries
            ]
        except Exception:
            result[f.stem] = []
    return result


st.title("Analyse des interactions actine-actine et actine-ABP - PPI3D")

with st.sidebar:
    st.markdown("## Sommaire")
    st.markdown("""
- [Téléchargement des données](#telechargement-des-donnees)
- [Données filtrées](#donnees-filtrees-s1-actine)
- [Structures PDB valides](#structures-pdb-valides)
- [Clusters d'interactions](#clusters-d-interactions)
- [ABP](#abp)
""")

# ---------------------------------------------------------------------------
# Section téléchargement
# ---------------------------------------------------------------------------

st.header("Téléchargement des données")

STEPS = {
    "1/13":  "Téléchargement du summary PPI3D (BLAST)",
    "2/13":  "Téléchargement des entrées PDB",
    "3/13":  "Téléchargement de toutes les données (cluster table)",
    "4/13":  "Filtrage des structures (≥ 4 actines connectées) - notebook",
    "5/13":  "Téléchargement des interactions d'interface",
    "6/13":  "Alignement MAFFT par cluster de séquences",
    "7/13":  "Analyse des clusters d'interaction C70 - notebook",
    "8/13":  "Calcul B-factors interface C70 par cluster",
    "9/13":  "Génération script PyMOL surface complète C70",
    "10/13": "Génération scripts PyMOL par site S1",
    "11/13": "Analyse interface par cluster C70 - notebook",
    "12/13": "Heatmap S1 binding site et références clusters - notebook",
    "13/13": "Calcul B-factors S1 par cluster",
}

# Fichier de sortie attendu pour chaque étape
STEP_OUTPUT_FILES = {
    "1/13":  "data/raw/ppi3d_actin_summary.csv",
    "2/13":  "data/raw/pdb_entry_results.csv",
    "3/13":  "data/raw/all_data.csv",
    "4/13":  "data/filtered/filtered_pdb_entry.csv",
    "5/13":  "data/filtered/details/1.interactions.csv",
    "6/13":  "data/alignments/.done",
    "7/13":  "data/filtered/patches_infos_cluster_data_70.csv",
    "8/13":  "data/filtered/details/structures_files/bfactor_c70_interface",
    "9/13":  "data/filtered/details/structures_files/bfactor_c70_interface/view_full_surface.pml",
    "10/13": "data/filtered/details/structures_files/bfactor_c70_interface/by_s1_cluster",
    "11/13": "visualisations/actin_c70_contacts",
    "12/13": "visualisations/actin_s1_all_equitable_heatmap.png",
    "13/13": "data/filtered/details/structures_files/bfactor_cluster",
}

STEP_KEYS = list(STEPS.keys())
TOTAL = len(STEPS)

SKIP_KEYWORDS = ["Nothing to do", "unchanged",
                 "Using existing", "already up", "Déjà à jour",
                 "No changes detected"]


def step_md(key, label, state):
    if state == "pending":
        return f"⬜ **{key}** — {label}"
    if state == "running":
        return f"🔄 **{key}** — {label} *(en cours...)*"
    if state == "done":
        return f"✅ **{key}** — {label} — *données téléchargées*"
    if state == "skipped":
        return f"✅ **{key}** — {label} — *déjà à jour*"
    if state == "error":
        return f"❌ **{key}** — {label} — *erreur*"
    return f"⬜ **{key}** — {label}"


def initial_state(key):
    path = STEP_OUTPUT_FILES.get(key)
    if path and os.path.exists(path):
        return "skipped"
    return "pending"


def parse_sub_progress(line):
    """
    Retourne (current, total, label) si la ligne contient une progression,
    sinon None.
    Patterns détectés :
      - Étape 2 : "1/148 4A7F"          → PDB X/total
      - Étape 5 : "Downloading 1/3"     → pass X/total
      - Étape 6 : "200,000 lignes traitées (471 Mo en RAM)..."
    """
    # Étape 2 : "X/Y PDBID"
    m = re.match(r"^(\d+)/(\d+)\s+\w+", line.strip())
    if m:
        return int(m.group(1)), int(m.group(2)), f"PDB {line.strip()}"

    # Étape 5 : "Downloading X/Y"
    m = re.search(r"Downloading\s+(\d+)/(\d+)", line)
    if m:
        return int(m.group(1)), int(m.group(2)), line.strip()

    return None


# Bouton en haut, étapes en dessous

clicked = st.button("Lancer le téléchargement")

progress_bar = st.empty()
log_placeholder = st.empty()

placeholders = {}
for key, label in STEPS.items():
    placeholders[key] = st.empty()
    placeholders[key].markdown(step_md(key, label, initial_state(key)))

sub_progress_bar = st.empty()
sub_progress_text = st.empty()

if "just_downloaded" not in st.session_state:
    st.session_state.just_downloaded = False

all_done = all(initial_state(k) == "skipped" for k in STEP_KEYS)
if all_done and not st.session_state.just_downloaded:
    st.write("✅ Données déjà téléchargées")


if clicked:
    st.session_state.just_downloaded = True

    for key, label in STEPS.items():
        init = initial_state(key)
        placeholders[key].markdown(
            step_md(key, label, init if init == "skipped" else "pending"))

    pb = progress_bar.progress(0, text="Démarrage...")

    current_idx = -1
    step_had_skip = {}
    log_lines = []
    step_log_lines = []

    cmd = [sys.executable, "-u", "-m", "script.data_extract.pipeline_data"]
    if sys.platform == "darwin":
        cmd = ["caffeinate", "-i"] + cmd
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    for line in proc.stdout:
        line = line.rstrip()

        # Détecter le début d'une nouvelle étape
        for i, key in enumerate(STEP_KEYS):
            if "ETAPE" in line and key in line:
                if current_idx >= 0:
                    prev_key = STEP_KEYS[current_idx]
                    state = "skipped" if step_had_skip.get(
                        prev_key) else "done"
                    placeholders[prev_key].markdown(
                        step_md(prev_key, STEPS[prev_key], state))

                current_idx = i
                step_log_lines = []
                placeholders[key].markdown(step_md(key, STEPS[key], "running"))
                pb.progress(i / TOTAL, text=f"Étape {key} — {STEPS[key]}")
                sub_progress_bar.empty()
                sub_progress_text.empty()
                log_placeholder.empty()
                break

        # Sous-progression de l'étape courante
        result = parse_sub_progress(line)
        if result is not None:
            current, total, label = result
            if current is not None and total is not None and total > 0:
                sub_progress_bar.progress(current / total)
                sub_progress_text.caption(label)
            else:
                sub_progress_text.caption(label)

        # Suivre les mots-clés "skip"
        if current_idx >= 0:
            current_key = STEP_KEYS[current_idx]
            if any(kw in line for kw in SKIP_KEYWORDS):
                step_had_skip[current_key] = True

        if line:
            log_lines.append(line)
            step_log_lines.append(line)
            log_placeholder.text("\n".join(step_log_lines[-15:]))

    proc.wait()

    sub_progress_bar.empty()
    sub_progress_text.empty()
    log_placeholder.empty()

    if current_idx >= 0:
        last_key = STEP_KEYS[current_idx]
        state = "skipped" if step_had_skip.get(last_key) else "done"
        if proc.returncode != 0:
            state = "error"
        placeholders[last_key].markdown(
            step_md(last_key, STEPS[last_key], state))

    if proc.returncode == 0:
        pb.progress(1.0, text="Téléchargement terminé ✓")
    else:
        pb.progress(1.0, text="Erreur — pipeline arrêté")
        st.error("Pipeline arrêté — dernières lignes de log :")
        st.code("\n".join(log_lines[-40:]), language=None)


# ---------------------------------------------------------------------------
# Section données filtrées
# ---------------------------------------------------------------------------

st.divider()
st.header("Données filtrées (s1 actine)")


@st.cache_data
def _load_pdb_file(path, mtime=None):
    with open(path) as f:
        return f.read()


def extract_chain(pdb_data: str, chain_letter: str) -> str:
    lines = [l for l in pdb_data.splitlines()
             if l.startswith('ATOM')
             and len(l) > 21 and l[21] == chain_letter]
    return "\n".join(lines)


@st.cache_data
def load_csv(path, mtime=None):
    with open(path, newline="") as f:
        sample = f.read(10000)
    try:
        sep = csv.Sniffer().sniff(sample, delimiters=";,\t").delimiter
    except csv.Error:
        sep = ","
    return pd.read_csv(path, sep=sep)


def read_csv(path):
    """Charge un CSV en invalidant le cache si le fichier a changé."""
    mtime = os.path.getmtime(path) if os.path.exists(path) else None
    return load_csv(path, mtime=mtime)


TABLES = {
    "filtered_summary.csv": "data/filtered/filtered_summary.csv",
    "filtered_pdb_entry.csv": "data/filtered/filtered_pdb_entry.csv",
    "filtered_all_data.csv": "data/filtered/filtered_all_data.csv",
    "1.interactions.csv": "data/filtered/details/1.interactions.csv",
    "2.proteins.csv": "data/filtered/details/2.proteins.csv",
    "3.interface_residues.csv": "data/filtered/details/3.interface_residues.csv",
    "4.inter-residue_contacts.csv": "data/filtered/details/4.inter-residue_contacts.csv",
    "5.ligands.csv": "data/filtered/details/5.ligands.csv",
    "6.meta_alignement.csv": "data/filtered/details/6.meta_alignement.csv",
    "7.alignment_sequences.csv": "data/filtered/details/7.alignment_sequences.csv",
    "8.structures.csv": "data/filtered/details/8.structures.csv",
}

available_tables = {t: p for t, p in TABLES.items() if os.path.exists(p)}

if available_tables:
    selected = st.selectbox("Choisir un tableau",
                            list(available_tables.keys()))
    df = read_csv(available_tables[selected])
    hide_constant = st.checkbox(
        "Masquer les colonnes sans variation", value=False)
    if hide_constant:
        cols_to_show = [c for c in df.columns if df[c].nunique() > 1]
    else:
        cols_to_show = list(df.columns)
    st.caption(
        f"{len(df):,} lignes · {len(cols_to_show)} colonnes affichées (sur {len(df.columns)} total)")
    st.dataframe(df[cols_to_show], width="stretch")

# Métrique PDB
METRICS_FILES = {
    "pdb": "data/filtered/filtered_pdb_entry.csv",
    "interactions": "data/filtered/details/1.interactions.csv",
    "proteins": "data/filtered/details/2.proteins.csv",
    "residues": "data/filtered/details/3.interface_residues.csv",
    "ligands": "data/filtered/details/5.ligands.csv",
}

if all(os.path.exists(p) for p in METRICS_FILES.values()):
    df_pdb = read_csv(METRICS_FILES["pdb"])
    df_int = read_csv(METRICS_FILES["interactions"])
    df_prot = read_csv(METRICS_FILES["proteins"])
    df_res = read_csv(METRICS_FILES["residues"])
    df_lig = read_csv(METRICS_FILES["ligands"])

    nb_homo = int((df_pdb["Interface type"] == "homo").sum()
                  ) if "Interface type" in df_pdb.columns else 0
    nb_hetero = int((df_pdb["Interface type"] == "hetero").sum(
    )) if "Interface type" in df_pdb.columns else 0

    col1, col2, col3 = st.columns(3)
    col1.metric("Structures PDB",
                f"{df_pdb['pdb_id'].nunique():,}".replace(',', ' '))
    col2.metric("Interactions protéine-protéine",
                f"{len(df_int):,}".replace(',', ' '))
    col3.metric("Protéines partenaires uniques (ABP)",
                f"{df_prot[~df_prot['protein_name'].str.lower().str.contains('actin', na=False)]['protein_name'].nunique():,}".replace(',', ' '))


# ---------------------------------------------------------------------------
# Section PDB valides — explorateur
# ---------------------------------------------------------------------------

pdb_filt_path = "data/filtered/filtered_pdb_entry.csv"
if os.path.exists(pdb_filt_path):
    st.divider()
    st.header("Structures PDB valides")

    df_entry = read_csv(pdb_filt_path)
    pdb_ids = sorted(df_entry["pdb_id"].str.upper().unique())

    summary_path = "data/filtered/filtered_summary.csv"
    if os.path.exists(summary_path):
        df_sum = read_csv(summary_path)
        unique = df_sum[["PDB ID", "Structure title"]
                        ].drop_duplicates("PDB ID")
        title_map = dict(
            zip(unique["PDB ID"].str.upper(), unique["Structure title"]))
    else:
        title_map = {}

    col_sel, _ = st.columns([1, 1])
    with col_sel:
        selected_pdb = st.selectbox(
            f"Choisir une structure ({len(pdb_ids)} PDB valides)",
            pdb_ids,
            format_func=lambda x: f"{x} — {title_map[x]}" if x in title_map else x,
            key="pdb_selector",
        )

    sub = df_entry[df_entry["pdb_id"].str.upper() == selected_pdb]

    # Reset de l'interaction sélectionnée quand on change de PDB
    if st.session_state.get("last_pdb") != selected_pdb:
        st.session_state["sel_inter"] = None
        st.session_state["sel_node"] = None
        st.session_state["last_pdb"] = selected_pdb

    if "sel_inter" not in st.session_state:
        st.session_state["sel_inter"] = None
    if "sel_node" not in st.session_state:
        st.session_state["sel_node"] = None

    int1_path = "data/filtered/details/1.interactions.csv"
    pp_path = "data/filtered/proteins_per_pdb.csv"
    struct_csv = "data/filtered/details/8.structures.csv"
    cont_csv = "data/filtered/details/4.inter-residue_contacts.csv"

    # Données communes aux deux colonnes
    df_int1_g = read_csv(int1_path) if os.path.exists(int1_path) else None
    sub_int1_g = df_int1_g[df_int1_g["pdb_id"].str.upper(
    ) == selected_pdb] if df_int1_g is not None else None

    raw_path = "data/raw/pdb_entry_results.csv"
    sub_raw_g = None
    if os.path.exists(raw_path):
        df_raw_g = read_csv(raw_path)
        sub_raw_g = df_raw_g[df_raw_g["pdb_id"].str.upper() == selected_pdb]

    col_net, col_3d = st.columns(2)

    # ── Colonne gauche : réseau interactif ────────────────────────────────
    with col_net:
        st.markdown(f"**Réseau d'interactions — {selected_pdb}**")
        if os.path.exists(pp_path) and sub_int1_g is not None:
            df_pp_g = read_csv(pp_path)
            sub_pp_g = df_pp_g[df_pp_g["pdb_id"].str.upper()
                               == selected_pdb]

            # Nœuds
            g_nodes = []
            for _, nr in sub_pp_g.drop_duplicates("chain").iterrows():
                is_act = bool(nr["is_actin"])
                chain_norm = norm_chain_id(str(nr["chain"]))
                g_nodes.append({
                    "id": chain_norm,
                    "label": chain_norm.split("_")[-1].upper(),
                    "title": str(nr["protein"]),
                    "color": "#E67E22" if is_act else "#2ECC71",
                    "size": 15,
                })

            # Mapping pair → interaction_id
            pair_to_iid = {}
            if sub_int1_g is not None:
                for _, er in sub_int1_g.iterrows():
                    ca = norm_chain_id(str(er["chain_A_id"]))
                    cb = norm_chain_id(str(er["chain_B_id"]))
                    pair_to_iid[frozenset([ca, cb])] = str(
                        er["interaction_id"])

            # Arêtes depuis pdb_entry_results
            g_edges = []
            edge_ids_set = set()
            sel_inter_str = str(st.session_state.get("sel_inter", ""))
            src = sub_raw_g if sub_raw_g is not None else (
                sub_int1_g.rename(
                    columns={"chain_A_id": "Interactor 1", "chain_B_id": "Interactor 2"})
                if sub_int1_g is not None else None
            )
            if src is not None:
                for _, er in src.iterrows():
                    ca = norm_chain_id(str(er["Interactor 1"]))
                    cb = norm_chain_id(str(er["Interactor 2"]))
                    iid = pair_to_iid.get(frozenset([ca, cb]))
                    eid = iid if iid else f"pair_{ca}_{cb}"
                    if iid:
                        edge_ids_set.add(eid)
                    is_sel = (eid == sel_inter_str)
                    g_edges.append({
                        "id": eid, "source": ca, "target": cb,
                        "width": 4 if is_sel else 1,
                        "color": "#FFD700" if is_sel else "#aaaaaa",
                    })

            if "sel_node" not in st.session_state:
                st.session_state["sel_node"] = None

            sel = _GRAPH_COMPONENT(
                nodes=g_nodes, edges=g_edges, height=450, key="actin_graph"
            )

            if sel is not None:
                if sel.get("type") == "edge" and sel.get("id") in edge_ids_set:
                    new_iid = int(sel["id"])
                    if new_iid != st.session_state.get("sel_inter"):
                        st.session_state["sel_inter"] = new_iid
                        st.session_state["sel_node"] = None
                        st.rerun()
                elif sel.get("type") == "node":
                    new_node = sel["id"]
                    if new_node != st.session_state.get("sel_node"):
                        st.session_state["sel_node"] = new_node
                        st.session_state["sel_inter"] = None
                        st.rerun()
            elif sel is None:
                if st.session_state.get("sel_node") or st.session_state.get("sel_inter"):
                    st.session_state["sel_node"] = None
                    st.session_state["sel_inter"] = None
                    st.rerun()

    # ── Colonne droite : viewer 3D ────────────────────────────────────────
    with col_3d:
        st.markdown("**Visualisation 3D — contacts d'interface**")
        if os.path.exists(struct_csv) and os.path.exists(cont_csv):
            df_struct8 = read_csv(struct_csv)
            sub_s8 = df_struct8[df_struct8["pdb_id"].str.upper()
                                == selected_pdb]

            pdb_file = None
            for _, r in sub_s8.iterrows():
                p = str(r.get("biological_assembly_pdb_file", "")).strip()
                if p and p.lower() != "nan" and os.path.exists(p):
                    pdb_file = p
                    break
            if pdb_file is None:
                for _, r in sub_s8.iterrows():
                    p = str(r.get("biological_assembly_cif_file", "")).strip()
                    if p and p.lower() != "nan" and os.path.exists(p):
                        pdb_file = p
                        break

            if pdb_file:
                int_ids = sub_s8["interaction_id"].tolist()
                df_cont4 = read_csv(cont_csv)

                sel_inter = st.session_state.get("sel_inter")
                sel_node = st.session_state.get("sel_node")

                # chain_colors : dict chain_letter → couleur (None = pas de sélection)
                chain_colors = None

                if sel_inter and sel_inter in int_ids:
                    # Cas 1 : une arête est sélectionnée — A cyan, B magenta
                    sub_c = df_cont4[df_cont4["interaction_id"]
                                     == sel_inter]
                    sel_row = df_int1_g[df_int1_g["interaction_id"]
                                        == sel_inter].iloc[0]
                    chain_colors = {
                        str(sel_row["chain_A_id"]).split("_")[-1]: "#FFD700",
                        str(sel_row["chain_B_id"]).split("_")[-1]: "#39FF14",
                    }

                elif sel_node:
                    # Cas 2 : un nœud est sélectionné
                    # contacts : depuis 1.interactions.csv (données structurales)
                    node_interactions = df_int1_g[
                        (df_int1_g["chain_A_id"] == sel_node) |
                        (df_int1_g["chain_B_id"] == sel_node)
                    ]["interaction_id"].tolist()
                    sub_c = df_cont4[df_cont4["interaction_id"].isin(
                        node_interactions)]

                    # couleurs : tous les partenaires depuis pdb_entry_results
                    sel_letter = sel_node.split("_")[-1].upper()
                    chain_colors = {sel_letter: "#FFD700"}
                    if sub_raw_g is not None:
                        for _, r in sub_raw_g.iterrows():
                            i1 = norm_chain_id(str(r["Interactor 1"]))
                            i2 = norm_chain_id(str(r["Interactor 2"]))
                            if i1 == sel_node or i2 == sel_node:
                                partner = i2 if i1 == sel_node else i1
                                cl = partner.split("_")[-1].upper()
                                chain_colors.setdefault(cl, "#00E5FF")

                else:
                    # Cas 3 : rien de sélectionné
                    sub_c = df_cont4[df_cont4["interaction_id"].isin(
                        int_ids)]

                contact_res: dict = {}
                for _, cr in sub_c.iterrows():
                    for cc, rc in [("chain_A_id", "residue_A_structure"),
                                   ("chain_B_id", "residue_B_structure")]:
                        chain = str(cr.get(cc, "")).split("_")[-1]
                        resi = cr.get(rc)
                        if pd.notna(resi) and chain:
                            try:
                                contact_res.setdefault(
                                    chain, set()).add(int(float(resi)))
                            except (ValueError, TypeError):
                                pass

                fmt = "cif" if pdb_file.endswith(".cif") else "pdb"
                pdb_data = _load_pdb_file(
                    pdb_file, mtime=os.path.getmtime(pdb_file))

                # Cache : ne recalcule le HTML que si sel_inter ou le PDB change
                viewer_key = (selected_pdb, sel_inter, sel_node)
                if st.session_state.get("viewer_key") != viewer_key:
                    view = py3Dmol.view(width=580, height=450)
                    view.addModel(pdb_data, fmt)
                    if chain_colors:
                        view.setStyle(
                            {}, {"cartoon": {"color": "#BDBDBD", "opacity": 0.8}})
                        for chain, col in chain_colors.items():
                            view.addStyle(
                                {"chain": chain},
                                {"cartoon": {"color": col, "opacity": 1.0}}
                            )
                    else:
                        view.setStyle(
                            {}, {"cartoon": {"color": "white", "opacity": 0.9}})
                    for chain, residues in contact_res.items():
                        resi_list = list(residues)
                        # Cartoon rouge pour le squelette des résidus de contact
                        view.addStyle(
                            {"chain": chain, "resi": resi_list},
                            {"cartoon": {"color": "red"}}
                        )
                        # Sticks rouges pour les chaînes latérales
                        view.addStyle(
                            {"chain": chain, "resi": resi_list},
                            {"stick": {"color": "red", "radius": 0.15}}
                        )
                    view.setBackgroundColor("#1e1e1e")
                    view.zoomTo()
                    st.session_state["viewer_html"] = view._make_html()
                    st.session_state["viewer_key"] = viewer_key

                st.components.v1.html(
                    st.session_state["viewer_html"], height=470, scrolling=False)
            else:
                st.info("Fichier PDB non disponible.")

    # ── Tableau interactions du PDB ───────────────────────────────────────
    all_data_path_pdb = "data/filtered/filtered_all_data.csv"
    if os.path.exists(all_data_path_pdb):
        df_all_pdb = read_csv(all_data_path_pdb)
        sub_all_pdb = df_all_pdb[df_all_pdb["pdb_id"].str.upper(
        ) == selected_pdb]
        cols_to_show = [c for c in [
            "subunit_1", "subunit_2",
            "cluster_data_70",
            "s1_binding_site_cluster_data_70",
            "s2_binding_site_cluster_data_70",
        ] if c in sub_all_pdb.columns]
        if not sub_all_pdb.empty and cols_to_show:
            st.markdown("**Clusters d'interaction par paire**")
            st.dataframe(
                sub_all_pdb[cols_to_show].reset_index(drop=True),
                hide_index=True,
                width="stretch",
            )

    # ── Séquences avec résidus d'interface colorés ───────────────────────

    all_data_path = "data/filtered/filtered_all_data.csv"
    res3_path_seq = "data/filtered/details/3.interface_residues.csv"
    seq_inter = st.session_state.get("sel_inter")
    seq_node = st.session_state.get("sel_node")

    if (df_int1_g is not None
            and os.path.exists(all_data_path) and os.path.exists(res3_path_seq)
            and (seq_inter or seq_node)):

        def seq_html(seq, iface, col):
            W = 60
            rows = []
            for s in range(0, len(seq), W):
                chunk = seq[s:s + W]
                cells = ""
                for k, aa in enumerate(chunk):
                    if (s + k + 1) in iface:
                        cells += (
                            f'<span style="color:{col};font-weight:bold'
                            f';text-shadow:0 0 4px {col}88">{aa}</span>'
                        )
                    else:
                        cells += f'<span style="color:#aaa">{aa}</span>'
                num = str(s + 1)
                pad = "&nbsp;" * (5 - len(num))
                rows.append(
                    f'<span style="color:#555;user-select:none">'
                    f'{pad}{num}&nbsp;</span>{cells}'
                )
            return (
                '<div style="font-family:\'Courier New\',monospace;'
                'font-size:13px;line-height:2;background:#161b22;'
                'padding:12px 16px;border-radius:8px;overflow-x:auto">'
                + "<br>".join(rows) + "</div>"
            )

        def render_pair(ca, cb, iid, df_all_g, df_r3,
                        col_a="#FFD700", col_b="#39FF14", primary=None):
            ca_l, cb_l = ca.lower(), cb.lower()
            rseq = df_all_g[
                (df_all_g["subunit_1"].str.lower() == ca_l) &
                (df_all_g["subunit_2"].str.lower() == cb_l)
            ]
            if rseq.empty:
                st.info("Séquences non disponibles pour cette interaction.")
                return
            s1 = str(rseq.iloc[0].get("s1_sequence", ""))
            s2 = str(rseq.iloc[0].get("s2_sequence", ""))
            t1 = str(rseq.iloc[0].get("subunit_1_title", ca))
            t2 = str(rseq.iloc[0].get("subunit_2_title", cb))
            if iid is not None:
                sub_r3 = df_r3[df_r3["interaction_id"] == iid]
                r3_chain = sub_r3["chain"].str.lower()
                iface_a = set(
                    sub_r3[r3_chain == ca_l]["residue_number_sequence"]
                    .dropna().astype(int)
                )
                iface_b = set(
                    sub_r3[r3_chain == cb_l]["residue_number_sequence"]
                    .dropna().astype(int)
                )
            else:
                iface_a, iface_b = set(), set()
            # Si primary == cb, mettre cb à gauche
            swap = primary is not None and primary.lower() == cb.lower()
            if swap:
                pairs = [(cb, t2, s2, iface_b, col_a),
                         (ca, t1, s1, iface_a, col_b)]
            else:
                pairs = [(ca, t1, s1, iface_a, col_a),
                         (cb, t2, s2, iface_b, col_b)]
            sc1, sc2 = st.columns(2)
            for widget, (cid, title, seq, iface, col) in zip([sc1, sc2], pairs):
                with widget:
                    st.markdown(
                        f"**{title}** · `{cid}` · {len(seq)} aa · {len(iface)} résidus d'interface"
                    )
                    if seq and seq != "nan":
                        st.markdown(seq_html(seq, iface, col),
                                    unsafe_allow_html=True)

        df_all_g = read_csv(all_data_path)
        df_r3 = read_csv(res3_path_seq)

        # ── Arête sélectionnée : une paire ───────────────────────────────
        if seq_inter:
            row_i = df_int1_g[df_int1_g["interaction_id"] == seq_inter]
            if not row_i.empty:
                ca = str(row_i.iloc[0]["chain_A_id"])
                cb = str(row_i.iloc[0]["chain_B_id"])
                st.markdown("#### Séquences — résidus d'interface")
                render_pair(ca, cb, seq_inter, df_all_g, df_r3)

        # ── Nœud sélectionné : toutes ses interactions ───────────────────
        elif seq_node:
            # Source : filtered_all_data pour tout avoir (même sans données d'interface)
            node_rows = df_all_g[
                (df_all_g["subunit_1"].str.lower() == seq_node) |
                (df_all_g["subunit_2"].str.lower() == seq_node)
            ]
            if not node_rows.empty:
                n = len(node_rows)
                st.markdown(
                    f"#### Séquences — `{seq_node.split('_')[-1].upper()}` "
                    f"({n} interaction{'s' if n > 1 else ''})"
                )
                for _, row in node_rows.iterrows():
                    ca = str(row["subunit_1"])
                    cb = str(row["subunit_2"])
                    # Chercher l'interaction_id dans 1.interactions.csv (peut être absent)
                    match = df_int1_g[
                        ((df_int1_g["chain_A_id"].str.lower() == ca.lower()) & (df_int1_g["chain_B_id"].str.lower() == cb.lower())) |
                        ((df_int1_g["chain_A_id"].str.lower() == cb.lower()) & (
                            df_int1_g["chain_B_id"].str.lower() == ca.lower()))
                    ]
                    iid = int(match.iloc[0]["interaction_id"]
                              ) if not match.empty else None
                    render_pair(ca, cb, iid, df_all_g, df_r3,
                                col_a="#FFD700", col_b="#00E5FF",
                                primary=seq_node)

    # ── Protéines impliquées ──────────────────────────────────────────────
    proteins_path = "data/filtered/proteins_per_pdb.csv"
    if os.path.exists(proteins_path):
        st.markdown("**Protéines impliquées**")
        df_pp = read_csv(proteins_path)
        sub_pp = df_pp[df_pp["pdb_id"].str.upper() == selected_pdb]

        nb_actin = int(sub_pp["is_actin"].sum())
        nb_abp = int((~sub_pp["is_actin"]).sum())
        c1, c2, _ = st.columns([1, 1, 3])
        c1.metric("Actines", nb_actin)
        c2.metric("ABP", nb_abp)

        counts = (
            sub_pp.groupby(["protein"])
            .agg(
                Nombre=("chain", "count"),
                Chaînes=("chain", lambda x: ", ".join(
                    sorted(s.split("_")[-1].upper() for s in x)
                )),
            )
            .reset_index()
            .rename(columns={"protein": "Protéine"})
            .sort_values("Nombre", ascending=False)
            .reset_index(drop=True)
        )
        max_len = counts["Protéine"].str.len().max() if len(counts) else 0
        col_width = "large" if max_len > 40 else "medium" if max_len > 20 else "small"
        st.dataframe(counts, width="stretch", hide_index=True, column_config={
            "Protéine": st.column_config.TextColumn(width=col_width),
            "Chaînes": st.column_config.TextColumn(width="medium"),
            "Nombre": st.column_config.NumberColumn(width="small"),
        })
    else:
        st.info("Relancer le notebook `graphe_filter.ipynb` pour générer les données.")


# ---------------------------------------------------------------------------
# Section clusters
# ---------------------------------------------------------------------------

@st.cache_data(show_spinner="Chargement réseau binding sites…")
def _build_tripartite_graph_html(all_data_path: str) -> str:
    df_all = pd.read_csv(all_data_path)
    col_s1, col_c70, col_s2 = (
        "s1_binding_site_cluster_data_70",
        "cluster_data_70",
        "s2_binding_site_cluster_data_70",
    )
    df = df_all[[col_s1, col_c70, col_s2,
                 "subunit_2_title", "s2_actine",
                 "s2_sequence_cluster_70"]].copy()
    df = df.dropna(subset=[col_s1, col_c70, col_s2])
    df[col_s1] = df[col_s1].astype(str)
    df[col_c70] = df[col_c70].astype(str)
    df[col_s2] = df[col_s2].astype(str)
    df["s2_actine"] = df["s2_actine"].fillna(False).astype(bool)

    # Nœuds :
    #   S1 (rouge) = binding sites actine (S1 + S2 des interactions homo)
    #   C70 (violet) = clusters d'interface
    #   S2 (vert) = binding sites partenaire (S2 des interactions hétéro seulement)
    s1_counts: dict = {}
    c70_counts: dict = {}
    s2_counts: dict = {}
    s2_proteins: dict = {}  # s2_cluster → set of protein names (hétéro)
    s2_seq_cluster: dict = {}  # s2_binding_site_cluster → s2_sequence_cluster_70
    s1_c70: dict = {}        # (s1_cluster, c70) → weight
    c70_s2_hetero: dict = {}  # (c70, s2_cluster) → weight  (hétéro uniquement)
    # (c70, s2_actin_cluster) → weight (homo : S2 = actine)
    c70_s2_homo: dict = {}
    # c70 → list of (s1, s2, is_homo, w) pour le tooltip
    c70_pairs: dict = {}

    triplet_w: dict = {}  # (s1, c70, s2, is_homo) → weight

    for _, row in df.iterrows():
        s1, c70, s2 = row[col_s1], row[col_c70], row[col_s2]
        is_homo = bool(row["s2_actine"])
        s1_counts[s1] = s1_counts.get(s1, 0) + 1
        c70_counts[c70] = c70_counts.get(c70, 0) + 1
        s1_c70[(s1, c70)] = s1_c70.get((s1, c70), 0) + 1
        triplet_key = (s1, c70, s2, is_homo)
        triplet_w[triplet_key] = triplet_w.get(triplet_key, 0) + 1
        pair_key = (min(s1, s2), max(s1, s2),
                    is_homo) if is_homo else (s1, s2, is_homo)
        c70_pairs.setdefault(c70, {}).setdefault(pair_key, 0)
        c70_pairs[c70][pair_key] += 1
        if is_homo:
            # S2 est aussi de l'actine → même pool de nœuds rouges
            s1_counts[s2] = s1_counts.get(s2, 0) + 1
            c70_s2_homo[(c70, s2)] = c70_s2_homo.get((c70, s2), 0) + 1
        else:
            s2_counts[s2] = s2_counts.get(s2, 0) + 1
            prot = str(row["subunit_2_title"]) if pd.notna(
                row["subunit_2_title"]) else ""
            if prot:
                s2_proteins.setdefault(s2, set()).add(prot)
            seq_cl = row["s2_sequence_cluster_70"]
            if pd.notna(seq_cl):
                s2_seq_cluster[s2] = str(int(seq_cl))
            c70_s2_hetero[(c70, s2)] = c70_s2_hetero.get((c70, s2), 0) + 1

    all_counts = (list(s1_counts.values()) + list(c70_counts.values())
                  + list(s2_counts.values()))
    min_c, max_c = min(all_counts), max(all_counts)

    def node_size(cnt):
        if max_c == min_c:
            return 10
        return 8 + 18 * (cnt - min_c) / (max_c - min_c)

    net = Network(height="900px", width="100%", bgcolor="#ffffff")
    net.set_options("""{
      "layout": { "randomSeed": 42 },
      "physics": {
        "barnesHut": {
          "gravitationalConstant": -4000,
          "centralGravity": 1,
          "springLength": 10,
          "springConstant": 0.06,
          "damping": 0.6,
          "avoidOverlap": 1
        },
        "stabilization": { "enabled": true, "iterations": 2000, "fit": true },
        "minVelocity": 0.1
      },
      "edges": { "smooth": false }
    }""")

    # S1 / S2-actine → nœuds rouges (tous les binding sites actine)
    for nd, cnt in s1_counts.items():
        net.add_node(f"s1_{nd}", label=nd, color="#e05252",
                     size=node_size(cnt),
                     title=f"Actin binding site: {nd}\n{cnt} interactions",
                     font={"size": 11, "color": "#000000", "background": "white", "strokeWidth": 0})

    # C70 nodes (violet) — tooltip liste toutes les paires (S1→S2) du cluster
    for nd, cnt in c70_counts.items():
        pairs = c70_pairs.get(nd, {})
        # Trier par poids décroissant
        sorted_pairs = sorted(pairs.items(), key=lambda x: -x[1])
        lines = []
        for (s1, s2, is_homo), w in sorted_pairs[:20]:
            tag = "🔴↔🔴" if is_homo else "🔴→🟢"
            lines.append(f"  {tag}  {s1} ↔ {s2}  ({w}×)")
        tooltip = f"Cluster C70 : {nd}\n{cnt} interactions · {len(pairs)} paires uniques"
        if lines:
            tooltip += "\n\nInteractions :\n" + "\n".join(lines)
            if len(sorted_pairs) > 20:
                tooltip += f"\n  … +{len(sorted_pairs)-20} autres paires"
        net.add_node(f"c70_{nd}", label=nd, color="#7B52E0", shape="square",
                     size=node_size(cnt),
                     title=tooltip,
                     font={"size": 11, "color": "#000000", "background": "white", "strokeWidth": 0})

    # Table de recherche : node_id → termes recherchables (cluster + protéines)
    search_map: dict = {}

    # S2 binding site nodes (vert — partenaires hétéro uniquement)
    for nd, cnt in s2_counts.items():
        prots = sorted(s2_proteins.get(nd, set()))
        # Label = nom de protéine tronqué (ou cluster si plusieurs protéines)
        if len(prots) == 1:
            short = prots[0][:28] + ("…" if len(prots[0]) > 28 else "")
            display_label = short
        else:
            display_label = nd
        prot_str = "\n".join(f"  • {p}" for p in prots[:10])
        if len(prots) > 10:
            prot_str += f"\n  … +{len(prots)-10} autres"
        tooltip = f"S2 cluster: {nd}\n{cnt} interactions"
        if prot_str:
            tooltip += f"\n{prot_str}"
        search_map[f"s2_{nd}"] = " ".join([nd] + prots).lower()
        net.add_node(f"s2_{nd}", label=display_label, color="#52b788",
                     size=node_size(cnt),
                     title=tooltip,
                     font={"size": 6, "color": "#444444", "strokeWidth": 0})

    # Arêtes solides S1 ↔ C70 (toutes interactions)
    for (s1, c70), w in s1_c70.items():
        net.add_edge(f"s1_{s1}", f"c70_{c70}",
                     color="#e05252", width=1.5,
                     title=f"{w} interactions")

    # Arêtes solides C70 ↔ S2-partenaire (hétéro)
    for (c70, s2), w in c70_s2_hetero.items():
        net.add_edge(f"c70_{c70}", f"s2_{s2}",
                     color="#52b788", width=1.5,
                     title=f"{w} interactions")

    # Arêtes solides C70 ↔ S2-actine (homo — nœud rouge)
    for (c70, s2actin), w in c70_s2_homo.items():
        net.add_edge(f"c70_{c70}", f"s1_{s2actin}",
                     color="#e05252", width=1.5,
                     title=f"{w} interactions actine-actine")

    # C70 avec plus de 3 paires uniques → connexion directe S1↔S2 en pointillé
    busy_c70 = {c70 for c70, pairs in c70_pairs.items() if len(pairs) >= 2}
    direct_edges_added: set = set()

    # Pointillés par triplet (S1, C70, S2) — deux segments courbes qui passent par C70
    # Hétéro : gris · Homo : orange
    _dyn = {"enabled": True, "type": "dynamic"}
    for (s1, c70, s2, is_homo), w in triplet_w.items():
        if is_homo:
            dash_color = "#FF8C00"
            tip = f"Actine↔Actine · C70: {c70} · {w} interactions"
            s2_node = f"s1_{s2}"
        else:
            dash_color = "#888888"
            tip = f"S1: {s1} → C70: {c70} → S2: {s2} · {w} interactions"
            s2_node = f"s2_{s2}"
        net.add_edge(f"s1_{s1}", f"c70_{c70}",
                     dashes=True, color=dash_color, width=0.8,
                     title=tip, smooth=_dyn)
        net.add_edge(f"c70_{c70}", s2_node,
                     dashes=True, color=dash_color, width=0.8,
                     title=tip, smooth=_dyn)

        # Arête directe S1↔S2 pour les C70 très connectés (>3 paires)
        if c70 in busy_c70:
            s1_node = f"s1_{s1}"
            direct_key = (min(s1_node, s2_node), max(s1_node, s2_node))
            if direct_key not in direct_edges_added:
                direct_edges_added.add(direct_key)
                direct_tip = f"Via C70 {c70} · {w} interactions"
                if is_homo:
                    direct_color = "#FF6600"
                else:
                    direct_color = "#444444"
                net.add_edge(s1_node, s2_node,
                             dashes=True, color=direct_color, width=1.5,
                             title=direct_tip)

    # ── Arêtes pointillées vertes : nœuds S2 du même cluster de séquence ────
    from collections import defaultdict as _dd2

    # s2_sequence_cluster_70 → ensemble de nœuds s2_*
    _seqcl_to_s2: dict = _dd2(set)
    for _s2nd, _seqcl in s2_seq_cluster.items():
        _seqcl_to_s2[_seqcl].add(f"s2_{_s2nd}")

    _s2_done: set = set()
    for _seqcl, _s2nodes in _seqcl_to_s2.items():
        if len(_s2nodes) < 2:
            continue
        _s2nodes = sorted(_s2nodes)
        for _i in range(len(_s2nodes)):
            for _j in range(_i + 1, len(_s2nodes)):
                _u2, _v2 = _s2nodes[_i], _s2nodes[_j]
                _k2 = (min(_u2, _v2), max(_u2, _v2))
                if _k2 not in _s2_done:
                    _s2_done.add(_k2)
                    net.add_edge(_u2, _v2, color="#00aa44", width=2,
                                 dashes=True, title=f"Même cluster de séquence S2 : {_seqcl}")

    html = net.generate_html()

    import json as _json
    search_map_js = _json.dumps(search_map, ensure_ascii=False)

    inject = (
        "function _freeze(){\n"
        "  network.setOptions({physics:{enabled:false}});\n"
        "  network.fit({animation:{duration:600,easingFunction:'easeInOutQuad'}});\n"
        "}\n"
        "network.once('stabilizationIterationsDone', _freeze);\n"
        "network.once('stabilized', _freeze);\n"
        "(function(){\n"
        f"  var SEARCH_MAP = {search_map_js};\n"
        "  var cont = document.getElementById('mynetwork').parentElement;\n"
        "  cont.style.position = 'relative';\n"
        "  var box = document.createElement('input');\n"
        "  box.type = 'text';\n"
        "  box.placeholder = 'Rechercher n\\u0153ud ou prot\\u00e9ine...';\n"
        "  box.style.cssText = 'position:absolute;top:10px;right:10px;z-index:9999;'\n"
        "    + 'padding:6px 12px;border:2px solid #aaa;border-radius:20px;'\n"
        "    + 'font-size:13px;width:240px;outline:none;background:white;';\n"
        "  cont.appendChild(box);\n"
        "  var origColors = {};\n"
        "  network.body.data.nodes.get().forEach(function(n){ origColors[n.id] = n.color; });\n"
        "  function _matches(n, q){\n"
        "    if(!q) return false;\n"
        "    var extra = SEARCH_MAP[String(n.id)] || '';\n"
        "    var haystack = (extra + ' ' + String(n.label || '')).toLowerCase();\n"
        "    return haystack.indexOf(q) !== -1;\n"
        "  }\n"
        "  box.addEventListener('input', function(){\n"
        "    var q = this.value.trim().toLowerCase();\n"
        "    var updates = network.body.data.nodes.get().map(function(n){\n"
        "      return {id: n.id, color: _matches(n,q) ? '#FFD700' : origColors[n.id]};\n"
        "    });\n"
        "    network.body.data.nodes.update(updates);\n"
        "    if(q){\n"
        "      var matched = network.body.data.nodes.get().filter(function(n){\n"
        "        return _matches(n, q);\n"
        "      });\n"
        "      if(matched.length === 1)\n"
        "        network.focus(matched[0].id, {scale:1.8, animation:{duration:400}});\n"
        "    }\n"
        "  });\n"
        "})();\n"
    )

    for pat in [
        "network = new vis.Network(container, data, options);",
        "var network = new vis.Network(container, data, options);",
    ]:
        if pat in html:
            html = html.replace(pat, pat + "\n" + inject)
            break
    return html


@st.cache_data(show_spinner="Génération du graphe global…")
def _build_global_graph_html(all_data_path: str, summary_path: str) -> str:
    df_all = pd.read_csv(all_data_path)
    df_summary = pd.read_csv(summary_path)
    df_summary["Expect value"] = df_summary["Expect value"].astype(float)
    real_actin = set(
        df_summary[df_summary["Expect value"] == 0]["Result protein"].unique())

    col_s1, col_s2 = "s1_binding_site_cluster_data_70", "s2_binding_site_cluster_data_70"
    G = nx.Graph()
    edge_c70s: dict = {}  # (u,v) → {c70: count}
    for _, row in df_all.iterrows():
        patch = row[col_s1]
        s2_title = row.get("subunit_2_title")
        if pd.isna(patch) or pd.isna(s2_title):
            continue
        if not G.has_node(patch):
            G.add_node(patch, side="cluster", count=0)
        G.nodes[patch]["count"] += 1
        if s2_title in real_actin:
            s2_cluster = row[col_s2]
            if pd.isna(s2_cluster):
                continue
            if not G.has_node(s2_cluster):
                G.add_node(s2_cluster, side="cluster", count=0)
            G.nodes[s2_cluster]["count"] += 1
            target = s2_cluster
        else:
            if not G.has_node(s2_title):
                G.add_node(s2_title, side="protein", count=0)
            G.nodes[s2_title]["count"] += 1
            target = s2_title
        if G.has_edge(patch, target):
            G[patch][target]["weight"] += 1
        else:
            G.add_edge(patch, target, weight=1)
        c70 = row.get("cluster_data_70")
        if pd.notna(c70):
            d = edge_c70s.setdefault((patch, target), {})
            d[str(c70)] = d.get(str(c70), 0) + 1

    counts_all = [d["count"] for _, d in G.nodes(data=True)]
    min_c, max_c = min(counts_all), max(counts_all)

    net = Network(height="900px", width="100%", bgcolor="#ffffff")
    net.set_options("""{
      "layout": { "randomSeed": 42 },
      "physics": {
        "barnesHut": {
          "gravitationalConstant": -8000,
          "centralGravity": 0.3,
          "springLength": 120,
          "springConstant": 0.04,
          "damping": 0.6,
          "avoidOverlap": 1
        },
        "stabilization": { "enabled": true, "iterations": 1000, "fit": true },
        "minVelocity": 0.1
      },
      "edges": { "smooth": false }
    }""")

    for n, d in G.nodes(data=True):
        is_cluster = d["side"] == "cluster"
        size = 20 + 110 * (d["count"] - min_c) / \
            (max_c - min_c) if max_c > min_c else 50
        color = "#e05252" if is_cluster else "#52e07a"
        label = n
        font = {"size": 0} if is_cluster else {
            "size": 12, "background": "white", "strokeWidth": 0}
        net.add_node(n, label=label, color=color, size=size,
                     title=f"{n} - {d['count']} interactions", font=font)
    for u, v in G.edges():
        c70_counts = edge_c70s.get((u, v)) or edge_c70s.get((v, u)) or {}
        lines = "\n".join(
            f"- {c70} ({n})" for c70, n in sorted(c70_counts.items(), key=lambda x: -x[1])
        ) if c70_counts else "—"
        title = f"{lines}\n({G[u][v]['weight']} interactions)"
        net.add_edge(u, v, color="#000000FF", width=2, title=title)

    # ── Arêtes pointillées vertes : protéines de la même famille ─────────────
    import re as _re
    from collections import defaultdict as _dd

    _GENERIC = {
        "protein", "actin", "like", "type", "associated", "binding",
        "domain", "subunit", "complex", "chain", "alpha", "beta", "gamma",
        "and", "the", "isoform", "heavy", "muscle", "cardiac", "skeletal",
        "delta", "epsilon", "smooth", "enteric", "cytoplasmic",
        "cell", "nuclear", "membrane", "cytoskeletal",
    }

    def _protein_keys(name):
        name = _re.sub(r"^Isoform\s+\S+\s+of\s+", "",
                       name, flags=_re.IGNORECASE)
        name = name.split(",")[0]
        return {w for w in _re.findall(r"[A-Za-z]{3,}", name.lower()) if w not in _GENERIC}

    _protein_nodes = [n for n, d in G.nodes(
        data=True) if d["side"] == "protein"]
    _word_to_p: dict = _dd(set)
    for _pn in _protein_nodes:
        for _w in _protein_keys(_pn):
            _word_to_p[_w].add(_pn)

    _family_done: set = set()
    for _word, _proteins in _word_to_p.items():
        if len(_proteins) < 2:
            continue
        _proteins = sorted(_proteins)
        for _i in range(len(_proteins)):
            for _j in range(_i + 1, len(_proteins)):
                _u, _v = _proteins[_i], _proteins[_j]
                _key = (min(_u, _v), max(_u, _v))
                if _key not in _family_done:
                    _family_done.add(_key)
                    net.add_edge(_u, _v, color="#00aa44", width=2, dashes=True,
                                 title=f"Même famille : {_word}")

    html = net.generate_html()

    inject = (
        "function _freeze(){ network.setOptions({physics:{enabled:false}}); network.fit({offset:{x:-120,y:0},animation:{duration:600,easingFunction:'easeInOutQuad'}}); }\n"
        "network.once('stabilizationIterationsDone', _freeze);\n"
        "network.once('stabilized', _freeze);\n"
        "(function(){\n"
        "  var cont = document.getElementById('mynetwork').parentElement;\n"
        "  cont.style.position = 'relative';\n"
        "  var box = document.createElement('input');\n"
        "  box.type = 'text';\n"
        "  box.placeholder = 'Rechercher un n\\u0153ud...';\n"
        "  box.style.cssText = 'position:absolute;top:10px;right:10px;z-index:9999;'\n"
        "    + 'padding:6px 12px;border:2px solid #aaa;border-radius:20px;'\n"
        "    + 'font-size:13px;width:210px;outline:none;background:white;';\n"
        "  cont.appendChild(box);\n"
        "  var origColors = {};\n"
        "  network.body.data.nodes.get().forEach(function(n){ origColors[n.id] = n.color; });\n"
        "  box.addEventListener('input', function(){\n"
        "    var q = this.value.trim().toLowerCase();\n"
        "    var updates = network.body.data.nodes.get().map(function(n){\n"
        "      var hit = q && String(n.id).toLowerCase().indexOf(q) !== -1;\n"
        "      return {id: n.id, color: hit ? '#FFD700' : origColors[n.id]};\n"
        "    });\n"
        "    network.body.data.nodes.update(updates);\n"
        "    if(q){\n"
        "      var matched = network.body.data.nodes.get().filter(function(n){\n"
        "        return String(n.id).toLowerCase().indexOf(q) !== -1;\n"
        "      });\n"
        "      if(matched.length === 1)\n"
        "        network.focus(matched[0].id, {scale:1.8, animation:{duration:400}});\n"
        "    }\n"
        "  });\n"
        "})();\n"
    )

    for pat in [
        "network = new vis.Network(container, data, options);",
        "var network = new vis.Network(container, data, options);",
    ]:
        if pat in html:
            html = html.replace(pat, pat + "\n" + inject)
            break
    return html


st.divider()
st.header("Clusters d'intéractions")

PATCHES_S1_CSV = "data/filtered/patches_infos_s1_binding_site.csv"
PATCHES_C70_CSV = "data/filtered/patches_infos_cluster_data_70.csv"
GLOBAL_GRAPH_IMG = "visualisations/patch_graphs/global.png"
GLOBAL_HEATMAP_IMG = "visualisations/patch_graphs/heatmap_binary.png"
S1_GRAPHS_DIR = "visualisations/patch_graphs/s1_binding_site"
C70_GRAPHS_DIR = "visualisations/patch_graphs/cluster_data_70"
S1_HEATMAPS_F4_DIR = "visualisations/patch_heatmaps_s1_contacts"
C70_HEATMAPS_F4_DIR = "visualisations/patch_heatmaps_c70_contacts"
C70_CONTACTS_DIR = "visualisations/actin_c70_contacts"
C70_CONTACTS_SURFACE_AREA_DIR = "visualisations/actin_c70_contacts_surface_area"
C70_HEATMAP_SURFACE_AREA = "visualisations/actin_c70_heatmap_surface_area.png"
S1_HEATMAP_RAW = "visualisations/actin_s1_homo_used_heatmap.png"
S1_HEATMAP_EQUITABLE = "visualisations/actin_s1_all_equitable_heatmap.png"
S1_HEATMAP_ABSOLUTE = "visualisations/actin_s1_heatmap_absolute.png"
S1_INTERFACE_FREQ_CSV = "data/filtered/actin_s1_canon_area_by_cluster.csv"
S1_INTERFACE_CLUSTER_DIR = "visualisations/actin_s1_clusters"
S1_INTERFACE_CLUSTER_BY_C70_DIR = "visualisations/actin_s1_clusters_by_c70"
BFACTOR_CLUSTER_DIR = "data/filtered/details/structures_files/bfactor_cluster"


@st.cache_data
def _global_bfac_max(bfactor_dir: str) -> float:
    global_max = 1.0
    if not os.path.exists(bfactor_dir):
        return global_max
    for fname in os.listdir(bfactor_dir):
        if not fname.endswith(".pdb"):
            continue
        try:
            with open(os.path.join(bfactor_dir, fname)) as f:
                vals = [
                    float(line[60:66])
                    for line in f
                    if line.startswith("ATOM") and len(line) >= 66
                ]
            if vals:
                global_max = max(global_max, max(vals))
        except Exception:
            continue
    return global_max


if not os.path.exists(PATCHES_C70_CSV):
    st.info("Lancer le pipeline pour générer l'analyse des clusters.")
else:
    _all_data_path = "data/filtered/filtered_all_data.csv"
    _summary_path = "data/filtered/filtered_summary.csv"
    tab_s1, tab_c70 = st.tabs(
        ["Binding Site Cluster Data 70", "Interaction Cluster Data 70"])

    # --- S1 Binding Site ---
    with tab_s1:
        if os.path.exists(_all_data_path) and os.path.exists(_summary_path):
            st.components.v1.html(_build_global_graph_html(
                _all_data_path, _summary_path), height=780)
            st.divider()
        if os.path.exists(PATCHES_S1_CSV):

            df_s1 = read_csv(PATCHES_S1_CSV)
            st.caption(f"{len(df_s1)} patchs")

            df_s1_display = df_s1.drop(
                columns=["ids_interactions"], errors="ignore").copy()

            # Colonnes C70 : clusters d'interaction contenant les interactions du patch S1
            if os.path.exists(PATCHES_C70_CSV):
                import re as _re

                def parse_ids(s):
                    s = str(s)
                    s = _re.sub(r'np\.int64\((\d+)\)', r'\1', s)
                    return set(int(x) for x in _re.findall(r'\d+', s))

                df_c70_tmp = read_csv(PATCHES_C70_CSV)
                # Mapping interaction_id → C70 patch
                id_to_c70 = {}
                for _, r in df_c70_tmp.iterrows():
                    for iid in parse_ids(r['ids_interactions']):
                        id_to_c70[iid] = str(r['patch'])

                c70_clusters_col, n_c70_col = [], []
                for _, row in df_s1.iterrows():
                    iids = parse_ids(row['ids_interactions'])
                    c70s = sorted(set(id_to_c70[i]
                                      for i in iids if i in id_to_c70))
                    c70_clusters_col.append(", ".join(c70s))
                    n_c70_col.append(len(c70s))

                df_s1_display['n_c70_clusters'] = n_c70_col
                df_s1_display['c70_clusters'] = c70_clusters_col

            df_s1_display.index = range(1, len(df_s1_display) + 1)
            st.dataframe(df_s1_display, width="stretch")

            # --- Heatmap global homo + hétéro ---
            st.subheader(
                "Heatmap — positions canoniques actin à l'interface S1")
            if not os.path.exists(S1_INTERFACE_FREQ_CSV):
                st.info(
                    "Lancer le notebook `notebooks/interface_analysis_s1.ipynb` pour générer les données.")
            else:
                heatmap_mode = st.selectbox(
                    "Normalisation",
                    ["Relative (max cluster = 1)", "Absolue (max = 100 %)"],
                    key="heatmap_norm_mode",
                )
                if heatmap_mode == "Relative (max cluster = 1)":
                    st.caption(
                        "Chaque ligne normalisée par son propre max — compare les motifs d'interface entre clusters.")
                    if os.path.exists(S1_HEATMAP_EQUITABLE):
                        st.image(S1_HEATMAP_EQUITABLE, width='content')
                    else:
                        st.info(
                            "Relancer le notebook `interface_analysis_s1.ipynb`.")
                else:
                    st.caption(
                        "Valeurs absolues en % buried ASA — compare l'intensité d'enfouissement entre clusters.")
                    if os.path.exists(S1_HEATMAP_ABSOLUTE):
                        st.image(S1_HEATMAP_ABSOLUTE, width='content')
                    else:
                        st.info(
                            "Relancer le notebook `interface_analysis_s1.ipynb` (cellule 6 — heatmap absolu).")

            st.divider()

            # --- Cluster sélectionné ---
            all_s1 = df_s1.sort_values("n_interactions", ascending=False)[
                "patch"].astype(str).tolist()
            sel_s1 = st.selectbox("Patch S1 binding site", all_s1, key="sel_s1",
                                  format_func=lambda p: f"{p} — {int(df_s1[df_s1['patch'].astype(str) == p]['n_interactions'].values[0])} interactions")

            row_s1 = df_s1[df_s1["patch"].astype(str) == sel_s1].iloc[0]

            # ── Réseau bipartite interactif + Surface 3D ─────────────────
            col_net_s1, col_3d_s1 = st.columns([3, 3])

            with col_net_s1:
                st.markdown(
                    "**Réseau interactif — résidus actine ↔ partenaires**")
                _bip_ok = all(os.path.exists(f) for f in _BIPARTITE_FILES)
                if _bip_ok:
                    _html_bip, _n_r, _n_p, _n_t = _build_bipartite_html(
                        sel_s1, _BIP_CACHE_VERSION, *_bip_mtimes())
                    if _html_bip:
                        st.caption(
                            f"{_n_r} résidus actine · {_n_p} protéines partenaires"
                            f" · n={_n_t} interactions")
                        st.components.v1.html(
                            _html_bip, height=540, scrolling=False)
                    else:
                        st.info("Réseau non disponible pour ce cluster.")
                else:
                    st.info("Lancer le pipeline pour générer les données réseau.")

            with col_3d_s1:
                st.markdown("**Interface 3D — couple représentatif**")
                if _bip_ok:
                    _s1_3d = _build_s1_3d_html(sel_s1, _BIP_CACHE_VERSION, *_bip_mtimes())
                    if _s1_3d:
                        _html_3d_s1, _s1_max, _s2_max, _s1_is_homo = _s1_3d
                        st.components.v1.html(
                            _html_3d_s1, height=470, scrolling=False)
                        if _s1_is_homo:
                            st.caption(
                                f"S1 (actine) : beige → rouge · % ASA buried (0 → 100 %) · max observé : {_s1_max:.1f} %"
                            )
                        else:
                            st.caption(
                                f"S1 (actine) : beige → rouge · % ASA buried (0 → 100 %) · max : {_s1_max:.1f} % "
                                f"· S2 (ABP) : vert · % ASA, max {_s2_max:.1f} %"
                            )
                    else:
                        st.info("Vue 3D non disponible pour ce patch.")
                else:
                    st.info("Lancer le pipeline pour générer les données réseau.")

            # Heatmap 1 ligne pour ce cluster
            s1_cluster_img = os.path.join(
                S1_INTERFACE_CLUSTER_DIR, f"{sel_s1}.png")
            if os.path.exists(s1_cluster_img):
                st.image(s1_cluster_img, width='content')
            elif os.path.exists(S1_INTERFACE_FREQ_CSV):
                st.info("Heatmap non généré pour ce cluster (seuil non atteint).")

            # Heatmap séparé par sous-cluster C70
            s1_by_c70_img = os.path.join(
                S1_INTERFACE_CLUSTER_BY_C70_DIR, f"{sel_s1}.png")
            if os.path.exists(s1_by_c70_img):
                st.caption(
                    "Profil par sous-cluster d'interaction C70 (chaque ligne = 1 cluster d'interface)")
                st.image(s1_by_c70_img, width='content')

    # --- Cluster Data 70 ---
    with tab_c70:
        if os.path.exists(_all_data_path):
            st.caption(
                "Rouge : binding sites actine (S1 + S2 homo) · Violet : clusters C70 · "
                "Vert : binding sites partenaire (hétéro) · "
                "Pointillé gris : interaction actine↔partenaire · Pointillé orange : interaction actine↔actine"
            )
            st.components.v1.html(_build_tripartite_graph_html(
                _all_data_path), height=900)
            st.divider()
        if os.path.exists(PATCHES_C70_CSV):
            df_c70 = read_csv(PATCHES_C70_CSV)
            st.caption(f"{len(df_c70)} patchs")

            df_c70_display = df_c70.drop(
                columns=["ids_interactions"], errors="ignore").copy()

            df_c70_display.index = range(1, len(df_c70_display) + 1)
            st.dataframe(df_c70_display, width="stretch")

            st.divider()
            all_c70 = df_c70.sort_values("n_interactions", ascending=False)[
                "patch"].astype(str).tolist()
            sel_c70 = st.selectbox("Patch cluster_data_70", all_c70, key="sel_c70",
                                   format_func=lambda p: f"{p} — {int(df_c70[df_c70['patch'].astype(str) == p]['n_interactions'].values[0])} interactions")

            row_c70 = df_c70[df_c70["patch"].astype(str) == sel_c70].iloc[0]

            # Clusters S2 depuis df_all — source de vérité pour tout ce qui suit
            df_all_c70_sel = pd.DataFrame()
            if os.path.exists("data/filtered/filtered_all_data.csv"):
                df_all_c70_sel = read_csv(
                    "data/filtered/filtered_all_data.csv")
                df_all_c70_sel = df_all_c70_sel[
                    df_all_c70_sel["cluster_data_70"].astype(
                        str) == str(sel_c70)
                ]
            s2_counts = (
                df_all_c70_sel["s2_sequence_cluster_70"].dropna().astype(int)
                .value_counts().sort_values(ascending=False)
                if not df_all_c70_sel.empty else pd.Series(dtype=int)
            )

            st.markdown(
                "**Réseau interactif — résidus actine ↔ résidus ABP**")
            _bip_ok = all(os.path.exists(f) for f in _BIPARTITE_FILES)
            if _bip_ok:
                _use_bipartite = st.toggle(
                    "Vue bipartite  (actine S1 à gauche | partenaire S2 à droite)",
                    value=False, key=f"bip_toggle_{sel_c70}")
                _html_c70, _n_s1, _n_s2, _n_tot, _html_3d_c70 = _build_bipartite_c70_html(
                    sel_c70, _use_bipartite, _BIP_CACHE_VERSION, *_bip_mtimes())
                if _html_c70:
                    st.caption(
                        f"{_n_s1} résidus actine (S1) · {_n_s2} résidus partenaire (S2)"
                        f" · n={_n_tot} interactions")
                    _net_height = max(700, max(_n_s1, _n_s2) * 55 + 20) if _use_bipartite else 620
                    st.components.v1.html(
                        _html_c70, height=_net_height, scrolling=False)
                    if _html_3d_c70:
                        st.markdown("**Interface 3D — couple représentatif**")
                        st.components.v1.html(
                            _html_3d_c70, height=450, scrolling=False)
                else:
                    st.info("Réseau non disponible pour ce cluster.")
            else:
                st.info("Lancer le pipeline pour générer les données réseau.")

st.divider()
st.header("ABP")

# Vue globale des protéines non-actines
proteins_path = "data/filtered/proteins_per_pdb.csv"
_all_data_path = "data/filtered/filtered_all_data.csv"
if os.path.exists(proteins_path):
    st.subheader("Protéines non-actines (ABP) — global")
    df_pp_all = read_csv(proteins_path)
    df_abp = df_pp_all[~df_pp_all["is_actin"]]

    def _fmt_ids(series):
        vals = sorted(set(str(v).strip()
                      for v in series.dropna() if str(v) not in ("nan", "")))
        return ", ".join(vals)

    if os.path.exists(_all_data_path):
        df_all_g = read_csv(_all_data_path)
        hetero = df_all_g[df_all_g["s1_actine"] & ~df_all_g["s2_actine"]][
            ["subunit_2", "s1_binding_site_cluster_data_70", "cluster_data_70"]
        ].copy()
        hetero["subunit_2"] = hetero["subunit_2"].str.lower()
        df_abp = df_abp.copy()
        df_abp["chain_low"] = df_abp["chain"].str.lower()
        merged = df_abp.merge(hetero, left_on="chain_low",
                              right_on="subunit_2", how="inner")
        abp_global = (
            merged.groupby("protein")
            .agg(
                Nb_noeuds=("chain", "nunique"),
                PDB=("pdb_id", lambda x: ", ".join(sorted(x.unique()))),
                **{"Binding site S1": ("s1_binding_site_cluster_data_70", _fmt_ids)},
                **{"Cluster C70": ("cluster_data_70", _fmt_ids)},
            )
            .reset_index()
            .rename(columns={"protein": "Protéine", "Nb_noeuds": "Nb nœuds"})
            .sort_values("Nb nœuds", ascending=False)
            .reset_index(drop=True)
        )
    else:
        abp_global = (
            df_abp.groupby("protein")
            .agg(
                Nb_noeuds=("chain", "nunique"),
                PDB=("pdb_id", lambda x: ", ".join(sorted(x.unique()))),
            )
            .reset_index()
            .rename(columns={"protein": "Protéine", "Nb_noeuds": "Nb nœuds"})
            .sort_values("Nb nœuds", ascending=False)
            .reset_index(drop=True)
        )

    max_len = abp_global["Protéine"].str.len().max() if len(abp_global) else 0
    col_width = "large" if max_len > 40 else "medium" if max_len > 20 else "small"
    st.dataframe(abp_global, width="stretch", hide_index=True, column_config={
        "Protéine": st.column_config.TextColumn(width=col_width),
        "Nb nœuds": st.column_config.NumberColumn(width="small"),
        "PDB": st.column_config.TextColumn(width="large"),
        "Binding site S1": st.column_config.TextColumn(width="medium"),
        "Cluster C70": st.column_config.TextColumn(width="medium"),
    })

# ── Heatmap ABP × résidus actin ──────────────────────────────────────────────
st.subheader("Heatmap — résidus actin contactés par les ABP")
st.caption(
    "Chaque cellule = % ASA buried moyen du résidu actin (position canonique) "
    "sur toutes les interactions de l'ABP — les interactions sans contact comptent 0%."
)

_ABP_HM_FILES = [
    "data/filtered/details/3.interface_residues.csv",
    "data/filtered/details/1.interactions.csv",
    "data/filtered/filtered_all_data.csv",
    "data/filtered/proteins_per_pdb.csv",
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
        df3["buried_ASA_percent"].astype(str).str.replace("%", "", regex=False),
        errors="coerce",
    )
    df3["residue_number_canon_mafft"] = pd.to_numeric(
        df3["residue_number_canon_mafft"], errors="coerce"
    )
    df3 = df3[df3["residue_number_canon_mafft"].notna()
              & df3["buried_ASA_percent"].notna()].copy()

    _actin_ch = set(df_pp_h[df_pp_h["is_actin"]]["chain"].str.lower())
    homo_iids = set(
        df_int_h[df_int_h["chain_B_id"].str.lower().isin(_actin_ch)]["interaction_id"]
    )
    het_int = (
        df_int_h[~df_int_h["interaction_id"].isin(homo_iids)]
        .merge(
            df_all_h[["subunit_1", "subunit_2", "subunit_2_title", "s2_actine"]],
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

    het_ids = set(het_int["interaction_id"])
    df3_h = df3[df3["interaction_id"].isin(het_ids)].copy()
    df3_h["_s1c"] = df3_h["interaction_id"].map(_s1ch)
    df3_h = df3_h[df3_h["chain"].str.lower() == df3_h["_s1c"]].copy()
    df3_h["abp"] = df3_h["interaction_id"].map(_abpn)
    df3_h["canon"] = df3_h["residue_number_canon_mafft"].astype(int)

    abp_freq = (
        df3_h.groupby("abp")["interaction_id"].nunique()
        .sort_values(ascending=False)
    )
    # Moyenne en comptant 0 pour les interactions sans contact
    agg_sum = (
        df3_h.groupby(["abp", "canon"])["buried_ASA_percent"]
        .sum().reset_index(name="asa_sum")
    )
    agg_sum["buried_ASA_percent"] = (
        agg_sum["asa_sum"] / agg_sum["abp"].map(abp_freq)
    )
    agg = agg_sum[["abp", "canon", "buried_ASA_percent"]]
    n_abp = max(abp_freq.shape[0], 1)
    res_freq = agg.groupby("canon")["abp"].nunique() / n_abp

    pivot = agg.pivot(index="abp", columns="canon", values="buried_ASA_percent")
    pivot = pivot.loc[abp_freq.index.intersection(pivot.index)]
    return pivot, abp_freq, res_freq


if all(os.path.exists(f) for f in _ABP_HM_FILES):
    _hm_mtimes = tuple(
        os.path.getmtime(f) if os.path.exists(f) else 0.0
        for f in _ABP_HM_FILES
    )
    _hm_res = _build_abp_heatmap_data(*_hm_mtimes)
    if _hm_res is not None:
        _pivot_full_hm, _abp_freq_hm, _ = _hm_res

        _pivot_hm = _pivot_full_hm.copy()

        import seaborn as sns
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import pdist

        _n_rows = len(_pivot_hm)
        _n_cols = len(_pivot_hm.columns)
        _fig_w = max(14, _n_cols * 0.35 + 4)
        _fig_h = max(6, _n_rows * 0.45 + 2)

        # Clustering hiérarchique sur les lignes (ABP) uniquement
        _mat = _pivot_hm.fillna(0).values.astype(float)
        _row_order = leaves_list(linkage(pdist(_mat, metric="cosine"), method="average"))
        _pivot_clustered = _pivot_hm.iloc[_row_order].fillna(0)

        # Labels Y avec n interactions (après réordonnancement)
        _ylabels = [
            f"{r}  ({_abp_freq_hm.get(r, 0)})"
            for r in _pivot_clustered.index
        ]

        _fig_hm, _ax_hm = plt.subplots(figsize=(_fig_w, _fig_h))
        sns.heatmap(
            _pivot_clustered,
            ax=_ax_hm,
            cmap="YlOrRd",
            vmin=0, vmax=100,
            xticklabels=True,
            yticklabels=_ylabels,
            linewidths=0,
            cbar_kws={"label": "% ASA buried moy (0 si absent)", "shrink": 0.6},
        )
        _ax_hm.set_xlabel("Position canonique (MAFFT)", fontsize=11, labelpad=6)
        _ax_hm.set_ylabel("")
        _ax_hm.set_xticklabels(_ax_hm.get_xticklabels(), rotation=90, fontsize=8)
        _ax_hm.set_yticklabels(_ax_hm.get_yticklabels(), rotation=0, fontsize=22)
        _ax_hm.set_title(
            f"% ASA buried moyen — résidus actin × ABP  ({_n_rows} ABP, {_n_cols} résidus)",
            fontsize=11, pad=8,
        )
        plt.tight_layout()
        st.pyplot(_fig_hm)
        plt.close(_fig_hm)
