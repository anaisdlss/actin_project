"""Extracted from streamlit.py — s1_heatmaps view/build helpers (keeps streamlit.py light)."""
import os
import numpy as np
import pandas as pd
import streamlit as st


_ACTIN_LEN = 375


_S1_GLOBAL_FILES = [
    "data/filtered/actin_s1_canon_area_by_cluster.csv",
    "data/filtered/filtered_all_data.csv",
    "data/filtered/details/1.interactions.csv",
    "data/filtered/details/3.interface_residues.csv",
]


@st.cache_data(show_spinner="Computing S1 heatmap (fair C70)…")
def _build_s1_global_heatmap(_mtimes):
    """Réplique regenerate_s1_global_heatmap : profil équitable-C70 par cluster S1,
    séparé HOMO (actin/actin) et HÉTÉRO (actin/ABP). Renvoie
    (positions, homo_labels, homo_mat, hetero_labels, hetero_mat) en % ASA absolu."""
    from collections import defaultdict
    if not all(os.path.exists(f) for f in _S1_GLOBAL_FILES):
        return None
    area = pd.read_csv(_S1_GLOBAL_FILES[0], index_col="patch")
    # actin canonical = 375 résidus : on ignore les positions > 375 (artefacts MSA)
    positions = [int(p) for p in area.columns if int(p) <= _ACTIN_LEN]
    pos_index = {p: i for i, p in enumerate(positions)}
    df = pd.read_csv(_S1_GLOBAL_FILES[1], low_memory=False)
    di = pd.read_csv(_S1_GLOBAL_FILES[2])[
        ["interaction_id", "chain_A_id", "chain_B_id"]]
    res = pd.read_csv(_S1_GLOBAL_FILES[3])
    res["canon"] = pd.to_numeric(
        res["residue_number_canon_mafft"], errors="coerce")
    res["basa"] = pd.to_numeric(
        res["buried_ASA_percent"].astype(
            str).str.replace("%", "", regex=False),
        errors="coerce")
    m = df.merge(di, left_on=["subunit_1", "subunit_2"],
                 right_on=["chain_A_id", "chain_B_id"], how="left")
    m = m[m["s1_binding_site_cluster_data_70"].notna()
          & m["interaction_id"].notna()]

    def _profile(sub_rows):
        by_c70 = defaultdict(dict)   # c70 -> {iid: actin_chain}
        for _, r in sub_rows.iterrows():
            if pd.isna(r.cluster_data_70):
                continue
            by_c70[str(r.cluster_data_70)][int(r.interaction_id)] = r.subunit_1
        if not by_c70:
            return None
        c70_profiles = []
        for _c70, iid_ch in by_c70.items():
            iids = list(iid_ch)
            chains = set(iid_ch.values())
            s = res[(res.interaction_id.isin(iids)) & (res.chain.isin(chains))]
            if s.empty:
                continue
            prof = (s.groupby(["interaction_id", "canon"])["basa"].max()
                    .groupby(level="canon").sum() / len(iids))
            v = np.zeros(len(positions))
            for canon, val in prof.items():
                if canon in pos_index:
                    v[pos_index[int(canon)]] = val
            c70_profiles.append(v)
        if not c70_profiles:
            return None
        return np.mean(c70_profiles, axis=0)   # équitable entre C70

    homo_prof, hetero_prof = {}, {}
    for patch, g in m.groupby("s1_binding_site_cluster_data_70"):
        gh = g[g["s2_actine"].fillna(False).astype(bool)]
        ge = g[~g["s2_actine"].fillna(False).astype(bool)]
        if len(gh):
            p = _profile(gh)
            if p is not None:
                homo_prof[str(patch)] = p
        if len(ge):
            p = _profile(ge)
            if p is not None:
                hetero_prof[str(patch)] = p

    def _numkey(patch):
        try:
            return int(patch.split("_")[1])
        except (IndexError, ValueError):
            return 10 ** 9

    homo = sorted(homo_prof, key=_numkey)
    hetero = sorted(hetero_prof, key=_numkey)
    homo_mat = (np.array([homo_prof[p] for p in homo])
                if homo else np.empty((0, len(positions))))
    hetero_mat = (np.array([hetero_prof[p] for p in hetero])
                  if hetero else np.empty((0, len(positions))))
    return positions, homo, homo_mat, hetero, hetero_mat


def _render_s1_global_plotly(data, relative, valid_clusters=None):
    """HOMO + HÉTÉRO empilés + bande agrégée, positions à leur place (trous inclus).

    Un clic sur une cellule sélectionne le cluster (label Y) dans `sel_s1` et
    recharge la page pour afficher son détail plus bas."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    positions, homo, hm, hetero, em = data
    full = list(range(min(positions), max(positions) + 1))
    idx = {p: i for i, p in enumerate(full)}

    def _expand(mat):
        out = np.full((mat.shape[0], len(full)), np.nan)
        for j, p in enumerate(positions):
            out[:, idx[p]] = mat[:, j]
        return out

    def _norm(mat):
        if not relative or mat.shape[0] == 0:
            return mat
        mx = np.nanmax(mat, axis=1, keepdims=True) if mat.size else mat
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.divide(mat, mx, out=np.full_like(mat, np.nan), where=mx > 0)

    hE, eE = _expand(hm), _expand(em)
    # bande agrégée : nb de clusters (HOMO+HÉTÉRO) touchant chaque position
    cnt = ((np.nan_to_num(hE) > 0).sum(0) + (np.nan_to_num(eE) > 0).sum(0))
    hZ, eZ = _norm(hE), _norm(eE)

    nh, ne = len(homo), len(hetero)
    zmax = 1.0 if relative else float(
        max(np.nanmax(hE) if hE.size else 0, np.nanmax(eE) if eE.size else 0, 1))
    # titre de l'échelle de couleur ; l'infobulle montre toujours le %ASA absolu
    cbar_title = "relatif (0–1)" if relative else "% ASA"

    # hauteurs de rangées proportionnelles au nombre de clusters + bande
    track = max(6, (nh + ne) * 0.08)
    tot = nh + ne + track
    fig = make_subplots(
        rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.015,
        row_heights=[nh / tot, ne / tot, track / tot])

    # %ASA absolu affiché au survol (via customdata), même en mode relatif
    _ht = ("Cluster : %{y}<br>Position : %{x}"
           "<br>%ASA : %{customdata:.2f}<extra></extra>")
    if nh:
        fig.add_trace(go.Heatmap(
            z=hZ, x=full, y=homo, customdata=hE,
            colorscale="YlOrRd", zmin=0, zmax=zmax,
            showscale=False, hoverongaps=False, hovertemplate=_ht), row=1, col=1)
    if ne:
        fig.add_trace(go.Heatmap(
            z=eZ, x=full, y=hetero, customdata=eE,
            colorscale="YlOrRd", zmin=0, zmax=zmax,
            colorbar=dict(title=cbar_title, thickness=12, len=0.75, y=0.6),
            hoverongaps=False, hovertemplate=_ht), row=2, col=1)
    fig.add_trace(go.Bar(
        x=full, y=cnt, marker=dict(
            color=cnt, colorscale="Blues", showscale=False),
        hovertemplate="Position : %{x}<br>clusters touchant : %{y}<extra></extra>"),
        row=3, col=1)

    fig.update_layout(
        height=max(540, (nh + ne) * 11 + 170), bargap=0,
        margin=dict(l=120, r=4, t=28, b=44),
        plot_bgcolor="white",   # trous (NaN) en blanc
        title=dict(text="Global S1 heatmap — HOMO (actin/actin) & "
                        "HETERO (actin/ABP), fair C70", font=dict(size=13)))
    fig.update_yaxes(autorange="reversed", tickfont=dict(size=9),
                     automargin=True, row=1, col=1)
    fig.update_yaxes(autorange="reversed", tickfont=dict(size=8),
                     automargin=True, row=2, col=1)
    fig.update_yaxes(title_text="n clusters",
                     title_font=dict(size=10), row=3, col=1)
    fig.update_xaxes(dtick=25, row=3, col=1,
                     title_text="Position canonical de l'actin (MAFFT)",
                     title_font=dict(size=11))
    # annotations de section
    fig.add_annotation(text=f"HOMO — actin / actin ({nh})", xref="paper",
                       yref="paper", x=0, y=1.0, showarrow=False,
                       font=dict(size=11), xanchor="left")
    _ev = st.plotly_chart(fig, use_container_width=True, key="s1g_hm",
                          on_select="rerun", selection_mode="points")

    # Clic sur une cellule -> sélectionne ce cluster dans sel_s1 (détail plus bas)
    _labels = set(homo) | set(hetero)
    if valid_clusters is not None:
        _labels &= set(valid_clusters)
    try:
        _pts = _ev["selection"]["points"]
    except (KeyError, TypeError):
        _pts = []
    for _p in reversed(_pts):
        _y = _p.get("y")
        if _y in _labels:
            if st.session_state.get("sel_s1") != _y:
                st.session_state["sel_s1"] = _y
                st.session_state["_scroll_to_s1"] = True
                st.toast(f"Cluster {_y} selected — opening detail…")
                st.rerun()
            break


@st.cache_data(show_spinner=False)
def _s1_sources(_mtimes):
    """Charge/fusionne les sources S1 une fois : (positions, m, res)."""
    if not all(os.path.exists(f) for f in _S1_GLOBAL_FILES):
        return None
    area = pd.read_csv(_S1_GLOBAL_FILES[0], index_col="patch")
    positions = [int(p) for p in area.columns if int(p) <= _ACTIN_LEN]
    df = pd.read_csv(_S1_GLOBAL_FILES[1], low_memory=False)
    di = pd.read_csv(_S1_GLOBAL_FILES[2])[
        ["interaction_id", "chain_A_id", "chain_B_id"]]
    res = pd.read_csv(_S1_GLOBAL_FILES[3])
    res["canon"] = pd.to_numeric(
        res["residue_number_canon_mafft"], errors="coerce")
    res["basa"] = pd.to_numeric(
        res["buried_ASA_percent"].astype(
            str).str.replace("%", "", regex=False),
        errors="coerce")
    m = df.merge(di, left_on=["subunit_1", "subunit_2"],
                 right_on=["chain_A_id", "chain_B_id"], how="left")
    m = m[m["s1_binding_site_cluster_data_70"].notna()
          & m["interaction_id"].notna()]
    keep = ["interaction_id", "subunit_1", "subunit_2", "subunit_2_title",
            "s2_actine", "cluster_data_70", "s1_binding_site_cluster_data_70"]
    return positions, m[keep], res[["interaction_id", "chain", "canon", "basa"]]


@st.cache_data(show_spinner="Profil du cluster…")
def _build_s1_patch_detail(patch, _mtimes):
    """Profil équitable-C70 du patch + décomposition par sous-cluster C70.
    Recalculé depuis les données -> toujours complet (tous les C70) et à jour."""
    import re as _re_ab
    src = _s1_sources(_mtimes)
    if src is None:
        return None
    positions, m, res = src
    pos_index = {p: i for i, p in enumerate(positions)}
    sub = m[m["s1_binding_site_cluster_data_70"].astype(str) == str(patch)]
    if sub.empty:
        return None

    def _prof(rows):
        iids = set(int(i) for i in rows["interaction_id"])
        chains = set(rows["subunit_1"])
        s = res[(res["interaction_id"].isin(iids))
                & (res["chain"].isin(chains))]
        if s.empty:
            return None
        pr = (s.groupby(["interaction_id", "canon"])["basa"].max()
              .groupby(level="canon").sum() / len(iids))
        v = np.zeros(len(positions))
        for canon, val in pr.items():
            if canon in pos_index:
                v[pos_index[int(canon)]] = val
        return v

    c70_rows = []
    for c70, g in sub.groupby("cluster_data_70"):
        if pd.isna(c70):
            continue
        v = _prof(g)
        if v is None:
            continue
        het = g[~g["s2_actine"].fillna(False).astype(bool)]
        if len(het) and het["subunit_2_title"].notna().any():
            # TOUS les noms d'ABP distincts de ce sous-cluster C70 (pas un seul)
            _names = sorted({_re_ab.sub(r"\s*\(.*?\)", "", str(t)).strip()
                             for t in het["subunit_2_title"].dropna()})
            abp = " ; ".join(_names)
        else:
            abp = "actin (homo)"
        c70_rows.append((str(c70), int(g["interaction_id"].nunique()), abp, v))
    if not c70_rows:
        return None
    patch_profile = np.mean([r[3] for r in c70_rows], axis=0)
    c70_rows.sort(key=lambda r: -r[1])
    return positions, patch_profile, c70_rows


def _render_s1_patch_plotly(detail, patch):
    """Profil 1-ligne du patch + heatmap décomposée par C70, positions à leur place."""
    import plotly.graph_objects as go
    positions, patch_profile, c70_rows = detail
    full = list(range(min(positions), max(positions) + 1))
    idx = {p: i for i, p in enumerate(full)}

    def _exp(v1d):
        out = np.full(len(full), np.nan)
        for j, p in enumerate(positions):
            out[idx[p]] = v1d[j]
        return out

    _z = _exp(patch_profile).reshape(1, -1)
    # échelle ABSOLUE : 0 -> 100 % ASA (comparable entre patchs, pas relatif au max)
    fig1 = go.Figure(go.Heatmap(
        z=_z, x=full, y=[str(patch)], colorscale="YlOrRd", zmin=0, zmax=100,
        colorbar=dict(title="% ASA", thickness=12), hoverongaps=False,
        hovertemplate="Position : %{x}<br>%ASA : %{z:.2f}<extra></extra>"))
    fig1.update_layout(
        height=150, margin=dict(l=4, r=4, t=30, b=36), plot_bgcolor="white",
        title=dict(text=f"Patch {patch} — S1 interface profile (fair-C70 %ASA)",
                   font=dict(size=12)))
    fig1.update_yaxes(showticklabels=False)
    fig1.update_xaxes(dtick=25, title_text="Actin canonical position (MAFFT)",
                      title_font=dict(size=10))
    # Clic sur une position → met à jour le sélecteur « Position canonical » ci-dessus.
    _ev1 = st.plotly_chart(fig1, use_container_width=True, key=f"s1prof_{patch}",
                           on_select="rerun", selection_mode="points")
    try:
        _pts = (_ev1 or {}).get("selection", {}).get("points", [])
    except Exception:
        _pts = []
    if _pts and _pts[0].get("x") is not None:
        _cx = int(round(float(_pts[0]["x"])))
        # ne relance que si la sélection change réellement (évite la boucle)
        if st.session_state.get(f"s1posdet_{patch}") != _cx:
            st.session_state[f"_s1_click_{patch}"] = _cx
            st.rerun()

    labels = [f"C70={r[0]} (n={r[1]}) — {r[2]}" for r in c70_rows]
    zmat = np.array([_exp(r[3]) for r in c70_rows])
    fig2 = go.Figure(go.Heatmap(
        z=zmat, x=full, y=labels, colorscale="YlOrRd", zmin=0, zmax=100,
        colorbar=dict(title="% ASA", thickness=12), hoverongaps=False,
        hovertemplate="%{y}<br>Position: %{x}<br>%ASA: %{z:.2f}<extra></extra>"))
    fig2.update_layout(
        height=max(170, len(labels) * 34 + 96), margin=dict(l=4, r=4, t=30, b=36),
        plot_bgcolor="white",
        title=dict(text=f"Patch {patch} — breakdown by C70 sub-cluster "
                   f"({len(labels)} C70)", font=dict(size=12)))
    fig2.update_yaxes(autorange="reversed", tickfont=dict(size=9))
    fig2.update_xaxes(dtick=25, title_text="Actin canonical position (MAFFT)",
                      title_font=dict(size=10))
    st.plotly_chart(fig2, use_container_width=True)


_S1_TAXID_NAMES = {
    4577: "Zea mays (maize)", 4930: "Saccharomyces sp.",
    4932: "Saccharomyces cerevisiae (levure)", 5833: "Plasmodium falciparum",
    7955: "Danio rerio (zebrafish)", 9030: "Gallus sp.",
    9031: "Gallus gallus (poulet)", 9606: "Homo sapiens",
    9823: "Sus scrofa (porc)", 9913: "Bos taurus (bovin)",
    9986: "Oryctolagus cuniculus (lapin)", 10090: "Mus musculus (souris)",
    10116: "Rattus norvegicus (rat)", 36329: "Plasmodium falciparum 3D7",
    137071: "Plasmodium falciparum HB3",
}


@st.cache_data(show_spinner="Per-position detail…")
def _s1_position_detail(patch, _mtimes):
    """Pour un cluster S1 : par position canonical, aa d'actin par organisme et
    aa d'ABP en contact (+ %ASA). Renvoie (positions, res_actin, con_abp)."""
    if not all(os.path.exists(f) for f in _S1_GLOBAL_FILES):
        return None
    alld = pd.read_csv(_S1_GLOBAL_FILES[1], low_memory=False)
    di = pd.read_csv(_S1_GLOBAL_FILES[2])[
        ["interaction_id", "chain_A_id", "chain_B_id"]]
    sub = alld[alld["s1_binding_site_cluster_data_70"].astype(
        str) == str(patch)].copy()
    if sub.empty:
        return None
    sub["abp"] = (sub["subunit_2_title"].fillna("Unknown").astype(str)
                  .str.replace(r"\s*\(.*?\)", "", regex=True).str.strip().str[:50])
    m = sub.merge(di, left_on=["subunit_1", "subunit_2"],
                  right_on=["chain_A_id", "chain_B_id"], how="inner")
    if m.empty:
        return None
    iid2actinch = dict(zip(m["interaction_id"], m["chain_A_id"]))
    iid2abp = dict(zip(m["interaction_id"], m["abp"]))
    iid2tax = dict(zip(m["interaction_id"], m["s1_taxonomy_id"]))
    iid2pdb = dict(zip(m["interaction_id"], m["pdb_id"]))
    iid2s2act = dict(zip(m["interaction_id"],
                         m["s2_actine"].fillna(False).astype(bool)))
    iids = set(m["interaction_id"])

    # aa d'actin par (position, organisme) — côté actin (chaîne A)
    res = pd.read_csv(_S1_GLOBAL_FILES[3])
    res = res[res["interaction_id"].isin(iids)].copy()
    res["canon"] = pd.to_numeric(
        res["residue_number_canon_mafft"], errors="coerce")
    res = res.dropna(subset=["canon"])
    res["canon"] = res["canon"].astype(int)
    res["_actinch"] = res["interaction_id"].map(iid2actinch)
    res = res[res["chain"].str.lower() == res["_actinch"].str.lower()].copy()
    res["taxid"] = res["interaction_id"].map(iid2tax)
    res["pdb_id"] = res["interaction_id"].map(iid2pdb)
    res_actin = res[["canon", "residue_name", "taxid", "pdb_id"]].copy()

    # aa d'ABP en contact par (position, ABP) — contacts hétéro uniquement
    con = pd.read_csv("data/filtered/details/4.inter-residue_contacts.csv")
    con = con[con["interaction_id"].isin(iids)].copy()
    con = con[con["interaction_id"].map(iid2s2act) == False]
    con["canon"] = pd.to_numeric(con["residue_A_canon_mafft"], errors="coerce")
    con = con.dropna(subset=["canon"])
    con["canon"] = con["canon"].astype(int)
    con["abp"] = con["interaction_id"].map(iid2abp)
    con["asa"] = pd.to_numeric(con["asa_pct_A"], errors="coerce")
    con_abp = con[["canon", "abp", "residue_B_name", "residue_B_sequence",
                   "asa"]].copy()

    positions = sorted(set(res_actin["canon"]) | set(con_abp["canon"]))
    return positions, res_actin, con_abp


def _render_s1_position_detail(detail, patch):
    """Sélecteur de position + tableau (organisme -> aa actin) et (ABP -> aa ABP)."""
    positions, res_actin, con_abp = detail
    if not positions:
        return
    st.markdown("**Per-position detail — which aa by organism and ABP**")
    # Un clic sur le heatmap de profil (plus bas) a pu déposer une position ici :
    # on la reporte dans le sélecteur avant de l'instancier.
    _ckey = f"_s1_click_{patch}"
    if _ckey in st.session_state:
        _cp = st.session_state.pop(_ckey)
        if _cp in positions:
            st.session_state[f"s1posdet_{patch}"] = _cp
    sel_pos = st.selectbox(
        "Position canonical", positions, key=f"s1posdet_{patch}")

    c_org, c_abp = st.columns(2)
    with c_org:
        st.caption("Actin — aa by organism")
        _r = res_actin[res_actin["canon"] == sel_pos]
        if _r.empty:
            st.caption("—")
        else:
            _org = (_r.groupby("taxid")
                    .agg(aa=("residue_name",
                             lambda s: " / ".join(sorted(set(s.dropna())))),
                         nb_pdb=("pdb_id", "nunique"))
                    .reset_index())
            _org["Organism"] = _org["taxid"].map(
                lambda t: _S1_TAXID_NAMES.get(int(t), f"TaxID {int(t)}")
                if pd.notna(t) else "?")
            st.dataframe(
                _org.rename(columns={"aa": "actin aa", "nb_pdb": "# PDBs"})[
                    ["Organism", "actin aa", "# PDBs"]],
                hide_index=True, use_container_width=True)
    with c_abp:
        st.caption("ABP — contacting aa (+ %ASA on the actin side)")
        _c = con_abp[con_abp["canon"] == sel_pos]
        if _c.empty:
            st.caption("No atomic ABP contact at this position.")
        else:
            _c = _c.copy()
            _c["aa ABP"] = (_c["residue_B_name"].astype(str)
                            + _c["residue_B_sequence"].astype(str))
            _t = (_c.groupby(["abp", "aa ABP"])
                  .agg(asa=("asa", "max")).reset_index()
                  .sort_values("asa", ascending=False))
            _t["asa"] = _t["asa"].round(1)
            st.dataframe(
                _t.rename(columns={"abp": "ABP", "asa": "% ASA actin"})[
                    ["ABP", "aa ABP", "% ASA actin"]],
                hide_index=True, use_container_width=True)
