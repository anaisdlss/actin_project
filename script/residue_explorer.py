"""Explorateur interactif actin : résidu / cluster / ABP / paire d'ABP / séquence.

Cinq vues sur la même table « passeport résidu » (residue_passport.build_passport) :

  * Résidu       : heatmap Plotly cliquable (ABP × position) + fiche complète du résidu
  * Cluster      : résumé d'un cluster C70 + résidus utilisés et leurs infos
  * ABP          : résumé des clusters d'un ABP + empreinte sur l'actin
  * Paire d'ABP  : recouvrement des empreintes de deux ABP sur l'actin
  * Séquence     : la séquence query de l'utilisateur alignée sur l'actin canonical
"""

import os
import numpy as np
import pandas as pd
import streamlit as st

from residue_passport import pp_mtimes

_TAXID_NAMES = {
    4577: "Zea mays (maize)", 4930: "Saccharomyces sp.",
    4932: "Saccharomyces cerevisiae (levure)", 5833: "Plasmodium falciparum",
    7955: "Danio rerio (zebrafish)", 9030: "Gallus sp.",
    9031: "Gallus gallus (poulet)", 9606: "Homo sapiens",
    9823: "Sus scrofa (porc)", 9913: "Bos taurus (bovin)",
    9986: "Oryctolagus cuniculus (lapin)", 10090: "Mus musculus (souris)",
    10116: "Rattus norvegicus (rat)", 36329: "Plasmodium falciparum 3D7",
    137071: "Plasmodium falciparum HB3",
}


@st.cache_data(show_spinner=False)
def _actin_aa_by_organism(_mtimes):
    """Par position canonical : aa d'actin observé, avec organisme (taxid) et PDB.
    Global (tous clusters). Renvoie un DataFrame (canon, residue_name, taxid, pdb_id)."""
    f_all = "data/filtered/filtered_all_data.csv"
    f_int = "data/filtered/details/1.interactions.csv"
    f_res = "data/filtered/details/3.interface_residues.csv"
    if not all(os.path.exists(f) for f in (f_all, f_int, f_res)):
        return None
    alld = pd.read_csv(f_all, low_memory=False)
    di = pd.read_csv(f_int)[["interaction_id", "chain_A_id", "chain_B_id"]]
    sub = alld[alld["s1_actine"].fillna(False).astype(bool)]
    m = sub.merge(di, left_on=["subunit_1", "subunit_2"],
                  right_on=["chain_A_id", "chain_B_id"], how="inner")
    if m.empty:
        return None
    iid2ch = dict(zip(m["interaction_id"], m["chain_A_id"]))
    iid2tax = dict(zip(m["interaction_id"], m["s1_taxonomy_id"]))
    iid2pdb = dict(zip(m["interaction_id"], m["pdb_id"]))
    res = pd.read_csv(f_res)
    if "residue_number_canon_mafft" not in res.columns:
        res["residue_number_canon_mafft"] = pd.NA  # étape 3 pas encore faite
    res = res[res["interaction_id"].isin(set(m["interaction_id"]))].copy()
    res["canon"] = pd.to_numeric(res["residue_number_canon_mafft"], errors="coerce")
    res = res.dropna(subset=["canon"])
    res["canon"] = res["canon"].astype(int)
    res["_ch"] = res["interaction_id"].map(iid2ch)
    res = res[res["chain"].str.lower() == res["_ch"].str.lower()].copy()
    res["taxid"] = res["interaction_id"].map(iid2tax)
    res["pdb_id"] = res["interaction_id"].map(iid2pdb)
    return res[["canon", "residue_name", "taxid", "pdb_id"]]


# ── Utilitaires fiche résidu ────────────────────────────────────────────────

def _pos_row(pp, canon):
    pos = pp["pos"]
    r = pos[pos["canon"] == canon]
    return None if r.empty else r.iloc[0]


def _fmt(v, nd=2):
    try:
        if pd.isna(v):
            return "—"
        return f"{float(v):.{nd}f}"
    except (TypeError, ValueError):
        return "—" if v is None else str(v)


def render_residue_fiche(pp, canon):
    """Fiche complète d'un résidu d'actin (position canonical)."""
    row = _pos_row(pp, canon)
    aa = row["actin_aa"] if row is not None and pd.notna(row.get("actin_aa")) else "?"
    st.markdown(f"### Residue **{aa}{canon}**  (canonical position {canon})")

    # ── Classe ProteoCast (en gros) + sensibilité vs moyenne ────────────────
    if row is not None:
        _cls = row.get("residue_class")
        st.markdown(f"## ProteoCast class: {_cls}")

        # Conservation (positif, intuitif : PLUS HAUT = PLUS CONSERVÉ = plus
        # sensible aux mutations). Référence = moyenne des résidus de SURFACE
        # (RSA >= 0.20), pas le cœur enfoui (comparaison surface/surface).
        _SURF = 0.20
        _p = pp["pos"]
        _rsa = pd.to_numeric(_p["rsa"], errors="coerce")
        _cons = pd.to_numeric(_p["conservation"], errors="coerce")
        avg_surf = float(_cons[_rsa >= _SURF].mean())
        this_c = pd.to_numeric(pd.Series([row.get("conservation")]),
                               errors="coerce").iloc[0]
        c1, c2 = st.columns(2)
        c1.metric(
            "Residue conservation", _fmt(this_c),
            delta=(f"{float(this_c) - avg_surf:+.2f} vs surface"
                   if pd.notna(this_c) and pd.notna(avg_surf) else None))
        c2.metric("Mean actin surface", _fmt(avg_surf))
        if pd.notna(this_c) and pd.notna(avg_surf):
            _more = float(this_c) > avg_surf
            st.caption(
                ("→ **more** conserved than the surface average" if _more
                 else "→ **less** conserved than the surface average")
                + " (higher = more conserved = residue more sensitive to mutations).")
        st.caption(
            f"at interface: **{'yes' if bool(row.get('at_interface')) else 'no'}** · "
            f"ABPs that contact it: **{int(row.get('n_abp', 0))}**")
    else:
        st.caption("Position absent from the ProteoCast conservation table.")

    # ── aa d'actin à cette position (% + organismes) ───────────────────────
    aa_g = _actin_aa_by_organism(pp_mtimes())
    if aa_g is not None:
        _a = aa_g[aa_g["canon"] == canon]
        if not _a.empty:
            _tot = _a["pdb_id"].nunique()
            _agg = (_a.groupby("residue_name")
                    .agg(nb=("pdb_id", "nunique"),
                         orgs=("taxid", lambda s: ", ".join(sorted({
                             _TAXID_NAMES.get(int(t), f"TaxID {int(t)}")
                             for t in s.dropna()}))))
                    .reset_index().sort_values("nb", ascending=False))
            _agg["pct"] = (_agg["nb"] / max(_tot, 1) * 100).round(0).astype(int)
            st.markdown("**actin aa at this position (by organism)**")
            st.dataframe(
                _agg.rename(columns={
                    "residue_name": "actin aa", "pct": "% of structures",
                    "orgs": "Organisms"})[
                    ["actin aa", "% of structures", "Organisms"]],
                hide_index=True, use_container_width=True)

    # ── ABP qui utilisent ce résidu (+ aa d'ABP + % ASA côté ACTINE) ─────────
    ra = pp["res_abp"]
    sub = ra[ra["canon"] == canon].sort_values("asa_max", ascending=False).copy()
    # aa d'ABP en contact (depuis la table contacts), agrégé par ABP
    con = pp["contacts"]
    if con is not None and not sub.empty:
        cc = con[con["canon"] == canon].copy()
        if not cc.empty:
            cc["_aa"] = (cc["abp_aa"].astype(str) + cc["abp_resnum"].astype(str))
            _aa_map = cc.groupby("abp")["_aa"].apply(
                lambda s: ", ".join(sorted(set(s))))
            sub["ABP aa"] = sub["abp"].map(_aa_map).fillna("—")
        else:
            sub["ABP aa"] = "—"
    st.markdown("**ABPs that contact this residue**  "
                "_(% ASA = buried actin residue; ABP aa = contacting ABP residue)_")
    if sub.empty:
        st.caption("No ABP contacts this residue in the dataset.")
    else:
        _cols = ["ABP", "% ASA max", "# interactions", "Clusters C70", "Sites S1"]
        if "ABP aa" in sub.columns:
            _cols = ["ABP", "ABP aa", "% ASA max", "# interactions",
                     "Clusters C70", "Sites S1"]
        st.dataframe(
            sub.rename(columns={
                "abp": "ABP", "asa_max": "% ASA max",
                "n_int": "# interactions", "c70s": "Clusters C70", "sites": "Sites S1",
            })[_cols],
            hide_index=True, use_container_width=True,
        )



# ── Heatmap Plotly cliquable ────────────────────────────────────────────────

def _build_pivot(pp, value="asa_max"):
    ra = pp["res_abp"]
    order = (ra.groupby("abp")["canon"].nunique()
             .sort_values(ascending=False).index.tolist())
    canons = sorted(ra["canon"].unique())
    piv = ra.pivot_table(index="abp", columns="canon",
                         values=value, aggfunc="max")
    piv = piv.reindex(index=order, columns=canons)
    return piv


def render_residue_tab(pp):
    import plotly.graph_objects as go
    st.markdown(
        "Click a **cell** of the heatmap (or choose a position below) "
        "to show all the info of the corresponding actin residue.")

    piv = _build_pivot(pp, "asa_max")
    xcanon = list(piv.columns)
    xlabels = [str(c) for c in xcanon]
    z = piv.values.astype(float)

    fig = go.Figure(go.Heatmap(
        z=z, x=xlabels, y=list(piv.index),
        colorscale="YlOrRd", zmin=0, zmax=100,
        colorbar=dict(title="% ASA max", thickness=12),
        hovertemplate="ABP : %{y}<br>Position : %{x}<br>% ASA max : %{z:.0f}<extra></extra>",
    ))
    fig.update_layout(
        height=max(420, len(piv.index) * 15 + 120),
        margin=dict(l=4, r=4, t=10, b=40),
        xaxis=dict(title="Position canonical (MAFFT)", tickangle=90,
                   tickfont=dict(size=8)),
        yaxis=dict(tickfont=dict(size=9), autorange="reversed"),
    )

    ev = st.plotly_chart(fig, use_container_width=True, key="explo_hm",
                         on_select="rerun", selection_mode="points")

    # position sélectionnée : clic heatmap prioritaire, sinon selectbox
    clicked = None
    try:
        pts = ev["selection"]["points"]
        if pts:
            clicked = int(pts[-1]["x"])
    except (KeyError, TypeError, ValueError):
        clicked = None

    canon_opts = xcanon
    if clicked is not None and clicked in canon_opts:
        st.session_state["explo_res_sel"] = clicked
    default_idx = (canon_opts.index(st.session_state["explo_res_sel"])
                   if st.session_state.get("explo_res_sel") in canon_opts else 0)
    sel = st.selectbox(
        "Position canonical", canon_opts, index=default_idx,
        key="explo_res_selbox",
        format_func=lambda c: f"{c}"
        + (f"  ({_pos_row(pp, c)['actin_aa']})"
           if _pos_row(pp, c) is not None
           and pd.notna(_pos_row(pp, c).get('actin_aa')) else ""))
    st.session_state["explo_res_sel"] = sel

    st.divider()
    render_residue_fiche(pp, sel)
    _render_actin_3d(pp, sel)


# ── Structure 3D actin avec résidu surligné ────────────────────────────────

_ACTIN_PDB = "actin_human_pdb.pdb"
# Human actin (ACTB) UniProt accession — used as a fallback when the optional
# local landmark PDB above is absent, so the 3D views also work on a fresh clone.
_ACTIN_UNIPROT = "P60709"


@st.cache_data(show_spinner=False)
@st.cache_data(show_spinner="Fetching the actin reference structure…")
def _actin_pdb_text():
    """Reference actin structure for the 3D views.

    Prefers the local landmark file ``actin_human_pdb.pdb`` (project root); if it
    is absent (e.g. a fresh clone), falls back to the AlphaFold model of human
    actin (P60709) fetched from the EBI API and cached."""
    import os
    if os.path.exists(_ACTIN_PDB):
        with open(_ACTIN_PDB) as fh:
            return fh.read()
    import proteocast_view
    return proteocast_view._fetch_af_pdb(_ACTIN_UNIPROT)


def _render_actin_3d(pp, canon):
    """Vue 3D de l'actin, résidu sélectionné surligné (sphères magenta)."""
    txt = _actin_pdb_text()
    if txt is None:
        return
    row = _pos_row(pp, canon)
    if row is None or pd.isna(row.get("Residue")):
        return
    resnum = int(row["Residue"])
    with st.expander("3D location of the residue on actin", expanded=False):
        try:
            import py3Dmol
        except ImportError:
            st.caption("py3Dmol indisponible.")
            return
        # chaîne majoritaire du PDB
        v = py3Dmol.view(width="100%", height=420)
        v.addModel(txt, "pdb")
        v.setStyle({}, {"cartoon": {"color": "lightgrey"}})
        v.addStyle({"resi": resnum},
                   {"stick": {"colorscheme": "magentaCarbon"}})
        v.addStyle({"resi": resnum},
                   {"sphere": {"color": "magenta", "opacity": 0.6}})
        v.zoomTo({"resi": resnum})
        v.setBackgroundColor("white")
        st.components.v1.html(v._make_html(), height=430, scrolling=False)
        st.caption(f"Residue {row.get('actin_aa','')}{canon} "
                   f"(structure no. {resnum}) highlighted in magenta.")


# ── Vue cluster (#2) ────────────────────────────────────────────────────────

def _residues_summary(pp, df_long, title_extra=""):
    """Tableau agrégé des résidus d'actin d'un sous-ensemble res_long,
    enrichi de la conservation (pos)."""
    if df_long.empty:
        st.info("No actin residue in this selection.")
        return
    agg = (df_long.groupby("canon")
           .agg(actin_aa=("actin_aa", "first"),
                asa_max=("asa_pct", "max"),
                n_abp=("abp", "nunique"),
                abps=("abp", lambda s: " ; ".join(sorted(set(s)))))
           .reset_index())
    pos = pp["pos"][["canon", "conservation", "mean_vs",
                     "frac_impactful", "rsa", "residue_class"]]
    agg = agg.merge(pos, on="canon", how="left").sort_values(
        "asa_max", ascending=False)
    agg["asa_max"] = agg["asa_max"].round(1)
    agg["residu"] = agg["actin_aa"].fillna("?") + agg["canon"].astype(str)
    st.dataframe(
        agg.rename(columns={
            "residu": "Residue", "asa_max": "% ASA max",
            "n_abp": "# ABPs", "abps": "Other ABPs on this residue",
            "conservation": "Conservation", "mean_vs": "Mut. sens.",
            "frac_impactful": "Frac. impact.", "rsa": "RSA",
            "residue_class": "Class",
        })[["Residue", "% ASA max", "Conservation", "Mut. sens.",
            "Frac. impact.", "RSA", "Class", "# ABPs",
            "Other ABPs on this residue"]],
        hide_index=True, use_container_width=True,
    )


def render_cluster_tab(pp):
    rl = pp["res_long"]
    c70s = sorted(rl["c70"].dropna().unique(),
                  key=lambda x: -rl[rl["c70"] == x]["interaction_id"].nunique())
    if not c70s:
        st.info("No C70 cluster available.")
        return
    sel = st.selectbox("Cluster C70 (actin–ABP)", c70s, key="explo_c70")
    sub = rl[rl["c70"] == sel]
    abps = sorted(sub["abp"].unique())
    sites = sorted(sub["s1_site"].unique())
    c1, c2, c3 = st.columns(3)
    c1.metric("Interactions", sub["interaction_id"].nunique())
    c2.metric("Actin residues", sub["canon"].nunique())
    c3.metric("ABPs involved", len(abps))
    st.caption(f"ABP : {', '.join(abps)}  ·  Site(s) S1 : {', '.join(sites)}")
    st.markdown("**Actin residues used by this cluster**")
    _residues_summary(pp, sub)


# ── Vue ABP (#3) ────────────────────────────────────────────────────────────

def render_abp_tab(pp):
    rl = pp["res_long"]
    sel = st.selectbox("ABP", pp["abp_list"], key="explo_abp")
    sub = rl[rl["abp"] == sel]
    if sub.empty:
        st.info("No data for this ABP.")
        return
    c1, c2, c3 = st.columns(3)
    c1.metric("Clusters C70", sub["c70"].nunique())
    c2.metric("Actin residues contacted", sub["canon"].nunique())
    c3.metric("Interactions", sub["interaction_id"].nunique())

    # Résumé par cluster
    st.markdown("**Summary of this ABP's clusters**")
    csum = (sub.groupby("c70")
            .agg(n_int=("interaction_id", "nunique"),
                 n_res=("canon", "nunique"),
                 sites=("s1_site", lambda s: " ; ".join(sorted(set(s)))))
            .reset_index().sort_values("n_int", ascending=False))
    st.dataframe(
        csum.rename(columns={
            "c70": "Cluster C70", "n_int": "Nb interactions",
            "n_res": "# actin residues", "sites": "Site(s) S1"}),
        hide_index=True, use_container_width=True,
    )

    # Détail : tous les résidus de l'empreinte, avec quels autres ABP les partagent
    st.markdown("**Footprint on actin — residues and other ABPs sharing them**")
    # « autres ABP » calculé sur res_long complet (pas seulement cet ABP)
    others = (rl.groupby("canon")["abp"]
              .apply(lambda s: " ; ".join(sorted(set(s) - {sel}))))
    disp = sub.copy()
    _residues_summary_with_others(pp, disp, others)


def _residues_summary_with_others(pp, df_long, others_map):
    agg = (df_long.groupby("canon")
           .agg(actin_aa=("actin_aa", "first"),
                asa_max=("asa_pct", "max"))
           .reset_index())
    pos = pp["pos"][["canon", "conservation", "mean_vs",
                     "frac_impactful", "rsa", "residue_class"]]
    agg = agg.merge(pos, on="canon", how="left").sort_values(
        "asa_max", ascending=False)
    agg["asa_max"] = agg["asa_max"].round(1)
    agg["residu"] = agg["actin_aa"].fillna("?") + agg["canon"].astype(str)
    agg["autres"] = agg["canon"].map(others_map).fillna("")
    st.dataframe(
        agg.rename(columns={
            "residu": "Residue", "asa_max": "% ASA max",
            "conservation": "Conservation", "mean_vs": "Mut. sens.",
            "frac_impactful": "Frac. impact.", "rsa": "RSA",
            "residue_class": "Class", "autres": "Other ABPs sharing",
        })[["Residue", "% ASA max", "Conservation", "Mut. sens.",
            "Frac. impact.", "RSA", "Class", "Other ABPs sharing"]],
        hide_index=True, use_container_width=True,
    )


# ── Vue paire d'ABP (#4) ────────────────────────────────────────────────────

def _all_pairs_overlap(pp):
    """Recouvrement (résidus partagés + Jaccard) de TOUTES les paires d'ABP."""
    from itertools import combinations
    rl = pp["res_long"]
    sets = {abp: set(g["canon"]) for abp, g in rl.groupby("abp")}
    rows = []
    for a, b in combinations(sorted(sets), 2):
        sa, sb = sets[a], sets[b]
        inter = len(sa & sb)
        uni = len(sa | sb)
        rows.append({
            "ABP A": a, "ABP B": b,
            "Residues A": len(sa), "Residues B": len(sb),
            "Shared": inter,
            "Jaccard": round(inter / max(uni, 1), 2),
        })
    return pd.DataFrame(rows).sort_values(
        ["Jaccard", "Shared"], ascending=False).reset_index(drop=True)


def render_pair_tab(pp):
    import plotly.graph_objects as go
    st.markdown(
        "Compares **two ABPs**: which actin positions each one contacts, and "
        "which they **share**. Many shared positions = both ABPs "
        "aim at the same spot on actin → they **compete for the site** "
        "(potential competition).")

    # ── Tableau récapitulatif : recouvrement de toutes les paires (2 à 2) ────
    with st.expander("Overview — overlap of all ABP pairs (pairwise)"):
        _ov = _all_pairs_overlap(pp)
        _ov_pos = _ov[_ov["Shared"] > 0]
        st.caption(
            f"{len(_ov)} pairs in total · {len(_ov_pos)} with at least 1 shared "
            "residue. Sorted by descending overlap (Jaccard): at the top = the "
            "pairs competing most for the same actin site. Sortable columns.")
        st.dataframe(_ov, hide_index=True, use_container_width=True, height=320)

    rl = pp["res_long"]
    abps = pp["abp_list"]
    c1, c2 = st.columns(2)
    a1 = c1.selectbox("ABP A", abps, index=0, key="explo_pair_a")
    a2 = c2.selectbox("ABP B", abps,
                      index=min(1, len(abps) - 1), key="explo_pair_b")
    if a1 == a2:
        st.info("Choose two different ABPs.")
        return

    sub1 = rl[rl["abp"] == a1]
    sub2 = rl[rl["abp"] == a2]

    # ── Détail : restreindre chaque ABP à UN de ses clusters côté actin ────
    _ALL = "All clusters"

    def _sites(sub):
        return [_ALL] + sorted(
            sub["s1_site"].dropna().unique(),
            key=lambda s: -sub[sub["s1_site"] == s]["canon"].nunique())
    d1, d2 = st.columns(2)
    cl1 = d1.selectbox(f"Cluster of {a1[:22]} (actin side)", _sites(sub1),
                       key="explo_pair_cl_a")
    cl2 = d2.selectbox(f"Cluster of {a2[:22]} (actin side)", _sites(sub2),
                       key="explo_pair_cl_b")
    if cl1 != _ALL:
        sub1 = sub1[sub1["s1_site"] == cl1]
    if cl2 != _ALL:
        sub2 = sub2[sub2["s1_site"] == cl2]
    lab1 = a1[:20] + (f" · {cl1}" if cl1 != _ALL else "")
    lab2 = a2[:20] + (f" · {cl2}" if cl2 != _ALL else "")

    s1 = set(sub1["canon"])
    s2 = set(sub2["canon"])
    inter = s1 & s2
    only1 = s1 - s2
    only2 = s2 - s1
    union = s1 | s2
    jac = len(inter) / max(len(union), 1)

    m1, m2, m3, m4 = st.columns(4)
    m1.metric(f"Residues {lab1}", len(s1))
    m2.metric(f"Residues {lab2}", len(s2))
    m3.metric("Shared residues", len(inter),
              help="Positions contacted by BOTH (current selections).")
    m4.metric("Jaccard (overlap)", f"{jac:.2f}",
              help="shared / union. 0 = no overlap, 1 = identical footprints.")

    st.caption(
        f"**How to read the band**: each stroke = an actin position. "
        f"**shared** line (purple) = touched by BOTH · **{lab1}** (blue) · "
        f"**{lab2}** (orange). The denser the purple line, the more they overlap.")

    allc = sorted(union)

    def _cat(c):
        return ("shared" if c in inter else lab1 if c in only1 else lab2)
    colors = {"shared": "#6a3d9a", lab1: "#1f77b4", lab2: "#ff7f0e"}
    fig = go.Figure()
    for cat in ["shared", lab1, lab2]:
        cs = [c for c in allc if _cat(c) == cat]
        fig.add_trace(go.Scatter(
            x=cs, y=[cat] * len(cs), mode="markers",
            marker=dict(size=16, color=colors[cat], symbol="line-ns",
                        line=dict(color=colors[cat], width=2)),
            name=cat, hovertemplate="position %{x}<extra>" + cat + "</extra>"))
    fig.update_layout(
        height=240, margin=dict(l=4, r=4, t=10, b=30),
        xaxis=dict(title="Actin canonical position", dtick=25),
        yaxis=dict(autorange="reversed"), showlegend=False)
    st.plotly_chart(fig, use_container_width=True)

    # ── Actin 3D : lieux de contact, mêmes couleurs (partagé / A / B) ──────
    _render_pair_3d(pp, inter, only1, only2, lab1, lab2, colors)

    if inter:
        st.markdown("**Shared residues (detail)**")
        _residues_summary(pp, pd.concat([sub1, sub2])[
            pd.concat([sub1, sub2])["canon"].isin(inter)])
    else:
        st.info("These two selections share no actin residue.")


def _render_pair_3d(pp, inter, only1, only2, a1, a2, colors):
    """Surface d'actin : positions partagées / propres à A / propres à B, en couleur."""
    txt = _actin_pdb_text()
    if txt is None:
        return
    try:
        import py3Dmol
    except ImportError:
        return
    canon2resi = {int(c): int(r) for c, r in
                  zip(pp["pos"]["canon"], pp["pos"]["Residue"])
                  if pd.notna(c) and pd.notna(r)}

    def _resis(cset):
        return [canon2resi[c] for c in cset if c in canon2resi]

    v = py3Dmol.view(width="100%", height=460)
    v.addModel(txt, "pdb")
    v.setStyle({}, {"cartoon": {"color": "white"}})
    v.addSurface(py3Dmol.SES, {"opacity": 1.0, "color": "#eeeeee"})  # base grise
    for cset, col in [(only1, colors[a1]), (only2, colors[a2]),
                      (inter, colors["shared"])]:
        rl = _resis(cset)
        if rl:
            v.addSurface(py3Dmol.SES, {"opacity": 1.0, "color": col},
                         {"resi": rl})
    v.zoomTo()
    v.setBackgroundColor("white")
    st.components.v1.html(v._make_html(), height=470, scrolling=False)
    st.markdown(
        f"<span style='color:{colors['shared']};font-size:18px'>■</span> shared "
        f"&nbsp; <span style='color:{colors[a1]};font-size:18px'>■</span> {a1} only "
        f"&nbsp; <span style='color:{colors[a2]};font-size:18px'>■</span> {a2} only "
        "&nbsp; (grey = rest of the actin)", unsafe_allow_html=True)


# ── Vue séquence query (#1) ──────────────────────────────────────────────────

_ACTIN_REF_FASTA = "data/P60709_ref.fasta"
_CONS_CSV = "data/proteocast/conservation_vs_asa_per_position.csv"


@st.cache_data(show_spinner="Fetching the UniProt sequence…")
def _fetch_uniprot(acc):
    """Récupère la séquence FASTA d'un accession UniProt. Renvoie (header, seq) ou None."""
    import urllib.request
    acc = acc.strip().split()[0]
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        with urllib.request.urlopen(url, timeout=20) as r:
            txt = r.read().decode()
    except Exception:
        return None
    lines = txt.splitlines()
    if not lines or not lines[0].startswith(">"):
        return None
    header = lines[0][1:].strip()
    seq = "".join(l.strip() for l in lines[1:] if not l.startswith(">"))
    return (header, seq) if seq else None


def _actin_bfactor_pdb(pp, col):
    """PDB de l'actin avec b-factor = valeur de la colonne `col` du passeport
    (ex. 'conservation'). Renvoie (pdb_text, vmin, vmax)."""
    txt = _actin_pdb_text()
    if txt is None:
        return None, 0.0, 1.0
    vmap = {}
    for _, r in pp["pos"].iterrows():
        if pd.notna(r.get("Residue")) and pd.notna(r.get(col)):
            vmap[int(r["Residue"])] = float(r[col])
    if not vmap:
        return None, 0.0, 1.0
    vmin, vmax = min(vmap.values()), max(vmap.values())
    out = []
    for ln in txt.splitlines():
        if ln.startswith(("ATOM", "HETATM")) and len(ln) >= 66:
            try:
                resi = int(ln[22:26])
            except ValueError:
                out.append(ln)
                continue
            out.append(ln[:60] + f"{vmap.get(resi, vmin):6.2f}" + ln[66:])
        else:
            out.append(ln)
    return "\n".join(out), vmin, vmax


@st.cache_data(show_spinner=False)
def _canonical_actin():
    """Séquence de référence P60709 (ACTB humaine) + mapping position P60709 -> canon."""
    if not (os.path.exists(_ACTIN_REF_FASTA) and os.path.exists(_CONS_CSV)):
        return None
    seq = "".join(l.strip() for l in open(_ACTIN_REF_FASTA)
                  if not l.startswith(">"))
    cons = pd.read_csv(_CONS_CSV)
    res2canon = {int(r): int(c) for r, c in
                 zip(cons["Residue"], cons["canon"]) if pd.notna(r) and pd.notna(c)}
    return seq, res2canon


@st.cache_data(show_spinner="Alignement MAFFT…")
def _mafft_align(query, ref_seq):
    """Aligne (query, ref) avec MAFFT. Renvoie (aln_ref, aln_query) ou None."""
    import subprocess
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as f:
        f.write(f">ref\n{ref_seq}\n>query\n{query}\n")
        path = f.name
    try:
        out = subprocess.run(["mafft", "--auto", "--quiet", path],
                             capture_output=True, text=True, timeout=120)
    except (subprocess.TimeoutExpired, FileNotFoundError):
        os.unlink(path)
        return None
    os.unlink(path)
    if out.returncode != 0:
        return None
    seqs, name = {}, None
    for line in out.stdout.splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    aln = {k: "".join(v).upper() for k, v in seqs.items()}
    if "ref" not in aln or "query" not in aln:
        return None
    return aln["ref"], aln["query"]


def render_sequence_tab(pp):
    st.markdown(
        "Paste an **actin sequence** (variant, isoform, ortholog). "
        "It is **aligned by MAFFT** to the canonical human actin (ACTB, P60709), "
        "and for each position **that varies** we show the ABPs known to use it.")

    ref = _canonical_actin()
    if ref is None:
        st.info("Canonical reference unavailable (data/P60709_ref.fasta).")
        return
    ref_seq, res2canon = ref

    with st.form("explo_seq_form"):
        acc = st.text_input(
            "UniProt accession no. (e.g. P68133 = skeletal α-actin) — "
            "fetches the sequence automatically",
            key="explo_seq_acc")
        seq_in = st.text_area(
            "… or paste an actin sequence (one letter per aa; digits ignored)",
            height=120, key="explo_seq_in",
            placeholder="MDDDIAALVVDNGSGMCKAGF...")
        st.form_submit_button("Analyse sequence")
    acc = acc.strip()

    import re
    if acc:
        _fetched = _fetch_uniprot(acc)
        if _fetched is None:
            st.error(f"Could not fetch “{acc}” from UniProt "
                     "(invalid accession or no network).")
            return
        _header, q = _fetched
        q = re.sub(r"[^A-Za-z]", "", q).upper()
        st.caption(f"Sequence fetched — **{acc}** · {_header[:90]} "
                   f"({len(q)} residues).")
    elif seq_in.strip():
        q = re.sub(r"[^A-Za-z]", "", seq_in).upper()
    else:
        st.caption(f"Reference: human ACTB P60709 ({len(ref_seq)} residues). "
                   "Enter a UniProt accession **or** paste a sequence, then click "
                   "**Analyse sequence**. MAFFT handles the alignment (any length).")
        return

    if len(q) < 30:
        st.warning("Sequence too short (< 30 residues).")
        return

    aln = _mafft_align(q, ref_seq)
    if aln is None:
        st.error("MAFFT alignment failed.")
        return
    aln_ref, aln_q = aln

    # Parcours de l'alignement : substitutions, insertions, délétions
    subs, indels = [], 0
    rpos = 0                      # position P60709 (1-based)
    identical = 0
    for cr, cq in zip(aln_ref, aln_q):
        ref_res = cr != "-"
        q_res = cq != "-"
        if ref_res:
            rpos += 1
        if ref_res and q_res:
            if cr == cq:
                identical += 1
            else:
                subs.append((rpos, cr, cq))
        elif ref_res != q_res:
            indels += 1

    _pid = 100 * identical / max(sum(1 for c in aln_ref if c != "-"), 1)
    c1, c2, c3 = st.columns(3)
    c1.metric("Identity vs ACTB", f"{_pid:.1f}%")
    c2.metric("Substitutions", len(subs))
    c3.metric("Indels (colonnes)", indels)

    if not subs:
        st.success("No substitution — sequence identical to canonical actin.")
        return

    # Pour chaque substitution : canon -> ABP concernés
    ra = pp["res_abp"]
    rows = []
    for rpos, wt, mut in subs:
        canon = res2canon.get(rpos)
        sub = ra[ra["canon"] == canon] if canon is not None else ra.iloc[0:0]
        abps = " ; ".join(sub.sort_values("asa_max", ascending=False)["abp"])
        pr = _pos_row(pp, canon) if canon is not None else None
        rows.append({
            "Variation (P60709)": f"{wt}{rpos}{mut}",
            "canon": canon,
            "Conservation": None if pr is None else pr.get("conservation"),
            "Mut. sens.": None if pr is None else pr.get("mean_vs"),
            "at interface": "yes" if (pr is not None and bool(pr.get("at_interface"))) else "no",
            "# ABPs": len(sub),
            "ABPs involved": abps,
        })
    df = pd.DataFrame(rows).sort_values("# ABPs", ascending=False)
    _n_iface = int((df["# ABPs"] > 0).sum())
    st.markdown(f"**{len(subs)} substitutions**, of which **{_n_iface}** at a position "
                "known at an ABP interface:")
    st.dataframe(df, hide_index=True, use_container_width=True)
    st.caption("The variations at the top (most ABPs) are the most likely "
               "to alter the actin–ABP interactions. “canon” = MAFFT "
               "canonical position used everywhere else in the app.")

    # ── Visuels : carte positionnelle + actin 3D avec variations surlignées ──
    import plotly.graph_objects as go
    _dd = df.dropna(subset=["canon"]).copy()
    _allp = sorted(pp["res_abp"]["canon"].unique())
    figt = go.Figure()
    figt.add_trace(go.Scatter(
        x=_allp, y=[0] * len(_allp), mode="markers", showlegend=False,
        marker=dict(size=6, color="#e2e2e2", symbol="square"),
        hoverinfo="skip"))
    if not _dd.empty:
        _cons = pd.to_numeric(_dd["Conservation"], errors="coerce")
        figt.add_trace(go.Scatter(
            x=_dd["canon"], y=[0] * len(_dd), mode="markers+text",
            marker=dict(size=15, color=_cons, colorscale="Blues",
                        line=dict(color="#333", width=1),
                        colorbar=dict(title="mut. sensitivity", thickness=12)),
            text=_dd["Variation (P60709)"], textposition="top center",
            textfont=dict(size=9), customdata=_cons,
            hovertemplate="%{text}<br>sensitivity: %{customdata:.2f}<extra></extra>",
            showlegend=False))
    figt.update_layout(
        height=190, margin=dict(l=4, r=4, t=34, b=34),
        xaxis=dict(title="Actin canonical position", dtick=25),
        yaxis=dict(visible=False),
        title=dict(text="Where your variations fall (grey = interface positions; "
                        "colour = mutational sensitivity, dark = more sensitive)",
                   font=dict(size=11)))
    st.plotly_chart(figt, use_container_width=True)

    # actin 3D : PAS de gradient — actin neutre, seules les positions
    # substituées sont mises en avant (grosses sphères vertes bien visibles).
    txt = _actin_pdb_text()
    if txt is not None:
        try:
            import py3Dmol
            v = py3Dmol.view(width="100%", height=500)
            v.addModel(txt, "pdb")
            v.setStyle({}, {"cartoon": {"color": "#c9ced6"}})
            v.addSurface(py3Dmol.SES, {"opacity": 0.85, "color": "#e6e8ec"})
            for _rp, _wt, _mut in subs:      # variations = sphères + surface vertes
                v.addStyle({"resi": _rp},
                           {"sphere": {"color": "#16d100", "radius": 2.2}})
                v.addSurface(py3Dmol.SES, {"opacity": 1.0, "color": "#16d100"},
                             {"resi": _rp})
            v.zoomTo()
            v.setBackgroundColor("white")
            st.components.v1.html(v._make_html(), height=510, scrolling=False)
            st.caption(
                "Neutral actin (grey). **Green** = your substituted positions, "
                "highlighted. (Their sensitivity and involved ABPs are in the table "
                "and the track above.)")
        except ImportError:
            pass


# ── Point d'entrée ──────────────────────────────────────────────────────────

def render_explorer(pp):
    # Deux entrées ici. Par résidu -> section « Empreinte ABP sur l'actin » ;
    # par cluster -> section « Clusters d'interactions » ; par ABP -> « Détail
    # par ABP ». Elles existent déjà, on ne les duplique pas.
    st.caption(
        "Two entry points: paste a query **actin sequence** (positions that "
        "vary × involved ABPs), or compare an **ABP pair** (overlap).")
    mode = st.radio(
        "Exploration mode",
        ["Query sequence", "ABP pair"],
        horizontal=True, key="explo_mode")
    st.divider()
    if mode == "ABP pair":
        render_pair_tab(pp)
    elif mode == "Query sequence":
        render_sequence_tab(pp)


# ── Vue « grand actin » : empreinte ABP colorée + clic résidu ────────────────

def _actin_nabp_pdb(pp):
    """PDB de l'actin, b-factor = nombre d'ABP contactant chaque résidu."""
    txt = _actin_pdb_text()
    if txt is None:
        return None, 1.0
    pos = pp["pos"]
    nmap = {}
    for _, r in pos.iterrows():
        if pd.notna(r.get("Residue")):
            nmap[int(r["Residue"])] = float(r.get("n_abp") or 0)
    mx = max(nmap.values(), default=1.0) or 1.0
    out = []
    for ln in txt.splitlines():
        if ln.startswith(("ATOM", "HETATM")) and len(ln) >= 66:
            try:
                resi = int(ln[22:26])
            except ValueError:
                out.append(ln)
                continue
            out.append(ln[:60] + f"{nmap.get(resi, 0.0):6.2f}" + ln[66:])
        else:
            out.append(ln)
    return "\n".join(out), mx


def _render_actin_overview_3d(pp, sel):
    """Surface de l'actin colorée par nb d'ABP + résidu sélectionné surligné en vert."""
    pdb_nabp, mx = _actin_nabp_pdb(pp)
    if pdb_nabp is None:
        return
    try:
        import py3Dmol
    except ImportError:
        st.info("py3Dmol indisponible.")
        return
    row = _pos_row(pp, sel)
    resi = (int(row["Residue"]) if row is not None
            and pd.notna(row.get("Residue")) else None)
    v = py3Dmol.view(width="100%", height=460)
    v.addModel(pdb_nabp, "pdb")
    v.setStyle({}, {"cartoon": {"color": "white"}})
    _grad = ["#FFFFCC", "#FEB24C", "#FC4E2A", "#B10026"]
    # surface pleine (aucune transparence)
    v.addSurface(py3Dmol.SES, {"opacity": 1.0, "colorscheme": {
        "prop": "b", "gradient": "linear", "colors": _grad,
        "min": 0, "max": mx}})
    if resi is not None:
        # résidu sélectionné : patch de surface en JAUNE FLUO, opaque
        v.addSurface(py3Dmol.SES, {"opacity": 1.0, "color": "#eaff00"},
                     {"resi": resi})
    v.zoomTo()   # vue d'ensemble stable — on ne re-zoome pas sur le résidu
    v.setBackgroundColor("white")
    st.components.v1.html(v._make_html(), height=470, scrolling=False)
    _aa = row["actin_aa"] if row is not None and pd.notna(row.get("actin_aa")) else ""
    st.caption(f"Surface: white = 0 ABP → red = {int(mx)} ABPs. "
               f"Selected residue **{_aa}{sel}** = fluo yellow.")


def render_actin_overview(pp):
    """Grand actin coloré par nb d'ABP + barre cliquable -> résidu surligné + fiche."""
    import plotly.graph_objects as go
    st.caption(
        "Each actin residue is coloured by the **number of different ABPs** "
        "that contact it. Click a bar (or choose a position): the residue "
        "is highlighted on the structure and all its info shows below.")

    pos = pp["pos"].copy()
    pos = pos[pd.to_numeric(pos["canon"], errors="coerce").notna()]
    pos = pos[pos["canon"] <= 375].sort_values("canon")
    x = pos["canon"].astype(int).tolist()
    y = pos["n_abp"].fillna(0).astype(int).tolist()
    aa = pos["actin_aa"].fillna("?").tolist()

    col3d, colbar = st.columns([1, 2])

    # Barre cliquable (rendue avant le 3D pour lire le clic, mais affichée à droite)
    with colbar:
        fig = go.Figure(go.Bar(
            x=x, y=y, customdata=aa,
            marker=dict(color=y, colorscale="YlOrRd",
                        colorbar=dict(title="n ABP", thickness=12)),
            hovertemplate=("Residue %{customdata}%{x}"
                           "<br>ABPs in contact: %{y}<extra></extra>")))
        fig.update_layout(
            height=460, margin=dict(l=6, r=6, t=10, b=44), bargap=0.1,
            xaxis=dict(title="Actin canonical position (MAFFT)", dtick=25),
            yaxis=dict(title="number of ABPs in contact"))
        ev = st.plotly_chart(fig, use_container_width=True, key="actin_ov",
                             on_select="rerun", selection_mode="points")

    # Clic -> on pilote le selectbox PAR SA CLÉ (sinon Streamlit ignore l'update)
    clicked = None
    try:
        pts = ev["selection"]["points"]
        if pts:
            clicked = int(pts[-1]["x"])
    except (KeyError, TypeError, ValueError):
        clicked = None
    if clicked is not None and clicked in x:
        st.session_state["actin_ov_selbox"] = clicked

    sel = st.session_state.get("actin_ov_selbox")
    if sel not in x:
        sel = x[0]

    # 3D à gauche, avec le résidu sélectionné surligné
    with col3d:
        _render_actin_overview_3d(pp, sel)

    # Sélecteur (repli / choix manuel) — piloté par la même clé que le clic
    st.selectbox(
        "Position to detail", x, key="actin_ov_selbox",
        format_func=lambda c: f"{c}"
        + (f"  ({_pos_row(pp, c)['actin_aa']})"
           if _pos_row(pp, c) is not None
           and pd.notna(_pos_row(pp, c).get('actin_aa')) else ""))
    sel = st.session_state.get("actin_ov_selbox", sel)
    render_residue_fiche(pp, sel)
