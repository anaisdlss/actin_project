"""Extracted from streamlit.py — proteocast_abp view/build helpers (keeps streamlit.py light)."""
import os
import csv
import re as _re_sd
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path as _Path
import proteocast_view


_PROTEOCAST_ABP_DIR = _Path("data/proteocast/abp")


_AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")


def _find_proteocast_csv(slug):
    """Cherche le ProteoCast d'un ABP : data/proteocast/abp/<slug>.csv
    ou data/proteocast/abp/<slug>/4.query_ProteoCast.csv."""
    for c in (_PROTEOCAST_ABP_DIR / f"{slug}.csv",
              _PROTEOCAST_ABP_DIR / slug / "4.query_ProteoCast.csv"):
        if c.exists():
            return c
    return None


@st.cache_data(show_spinner=False)
def _load_proteocast(path, mtime):
    pc = pd.read_csv(path)
    pc["aa"] = pc["Mutation"].astype(str).str[-1]
    piv = pc.pivot_table(index="aa", columns="Residue", values="Variant_score",
                         aggfunc="first").reindex(_AA_ORDER)
    return piv


def _render_proteocast_mutland(csv_path, iface_asa=None, title="", domains=None,
                               surface_only=False, rsa=None, focus=None):
    """Paysage mutationnel ProteoCast INTERACTIF (Plotly) : heatmap 20 subs ×
    positions (vert), piste %ASA d'interface (vert) et piste domaines (couleur
    par domaine). Tout partage l'axe X (zoom/survol comme les autres heatmaps).

    surface_only : ne garde que les résidus EXPOSÉS (RSA >= 0.2) et recalcule le
    gradient sur eux (c'est là qu'est l'interface ABP-actin). `rsa` : dict
    {position: RSA} calculé sur la vraie structure ABP (compute_abp_rsa).
    focus : (lo, hi) pour recadrer l'axe X sur la zone de liaison à l'actin
    (gros ABP multi-domaines : évite d'afficher 2500 résidus hors-sujet)."""
    import numpy as np
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    piv = _load_proteocast(str(csv_path), _Path(csv_path).stat().st_mtime)
    positions = [int(p) for p in piv.columns]
    aa = list(piv.index)
    _dom = list(domains or [])
    _ndom = len(_dom)

    # matrice des scores (copie modifiable)
    _z = piv.values.astype(float).copy()
    _z_bur = None                                # couche « cœur enfoui » en gris
    if surface_only:
        _rsa = dict(rsa or {})
        _exposed = np.array([_rsa.get(p, 1.0) >= 0.2 for p in positions])
        if _exposed.any():
            _z_bur = _z.copy()
            _z_bur[:, _exposed] = np.nan         # gris : seulement l'enfoui
            _z[:, ~_exposed] = np.nan            # vert : seulement l'exposé
    # gradient RELATIF : bornes = min/max des scores affichés. Si zoom sur la zone
    # de liaison, on calcule les bornes sur cette zone (sinon un pic ailleurs dans
    # les 2500 résidus hors-sujet écraserait le contraste de l'interface).
    _zg = _z
    if focus:
        _fm = np.array([focus[0] <= p <= focus[1] for p in positions])
        if _fm.any():
            _zg = _z[:, _fm]
    _zmin = float(np.nanmin(_zg)) if np.isfinite(np.nanmin(_zg)) else -8
    _zmax = float(np.nanmax(_zg)) if np.isfinite(np.nanmax(_zg)) else 1

    # colorscale ProteoCast vert, ancres RÉGULIÈRES (clair/foncé ~50/50) :
    # noir(-8, très délétère) -> vert -> blanc(+1, toléré)
    _green = [
        [0.000, "#000000"], [0.125, "#08260f"], [0.250, "#00441b"],
        [0.375, "#116b32"], [0.500, "#238b45"], [0.625, "#41ab5d"],
        [0.750, "#74c476"], [0.875, "#c7e9c0"], [1.000, "#ffffff"]]
    _asa_green = [[0, "#ffffff"], [0.25, "#c7e9c0"], [0.5, "#74c476"],
                  [0.75, "#238b45"], [1.0, "#00441b"]]

    _rows = 3 if _ndom else 2
    _heights = ([0.72, 0.08, 0.20] if _ndom else [0.9, 0.1])
    fig = make_subplots(rows=_rows, cols=1, shared_xaxes=True,
                        vertical_spacing=0.03, row_heights=_heights)

    # 1a) couche « cœur enfoui » en GRIS (blanc toléré -> noir délétère), affichée
    #     sous la couche verte pour garder le contexte sans la confondre avec la
    #     surface exposée. Même échelle que le vert pour rester comparable.
    if _z_bur is not None:
        _grey = [[0.0, "#000000"], [0.5, "#9e9e9e"], [1.0, "#ffffff"]]
        fig.add_trace(go.Heatmap(
            z=_z_bur, x=positions, y=aa, colorscale=_grey,
            zmin=_zmin, zmax=_zmax, showscale=False, hoverongaps=False,
            hovertemplate="aa %{y} · position %{x} (buried core)<br>"
                          "Variant_score: %{z:.2f}<extra></extra>",
        ), row=1, col=1)

    # 1) heatmap des scores (surface exposée, vert)
    fig.add_trace(go.Heatmap(
        z=_z, x=positions, y=aa, colorscale=_green, zmin=_zmin, zmax=_zmax,
        colorbar=dict(title="Variant_score", thickness=12, len=0.55, y=0.78),
        hoverongaps=False,
        hovertemplate="aa %{y} · position %{x}<br>Variant_score : %{z:.2f}<extra></extra>",
    ), row=1, col=1)

    # 2) piste ABP en contact à l'actin (%ASA, vert)
    #    En mode surface_only : on masque aussi les contacts ENFOUIS (RSA < 0.2)
    #    pour que la piste coïncide avec la heatmap (sinon barres de contact sous
    #    des colonnes masquées -> "l'ASA est là où c'est enfoui").
    _asa = {int(k): float(v) for k, v in (iface_asa or {}).items()}
    if surface_only:
        _rsa2 = dict(rsa or {})
        _asa_row = [(_asa.get(p, 0.0) if _rsa2.get(p, 1.0) >= 0.2 else np.nan)
                    for p in positions]
    else:
        _asa_row = [_asa.get(p, 0.0) for p in positions]
    fig.add_trace(go.Heatmap(
        z=[_asa_row], x=positions, y=["ABP contact actin"], colorscale=_asa_green,
        zmin=0, zmax=100, showscale=False, hoverongaps=False,
        hovertemplate="position %{x}<br>%ASA (ABP au contact) : %{z:.0f}<extra></extra>",
    ), row=2, col=1)

    # 3) piste domaines : un segment coloré par domaine (hover = nom + bornes)
    if _ndom:
        _pal = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5"]
        for _i, _d in enumerate(_dom):
            _c = _pal[_i % len(_pal)]
            _lab = f"{_d['name'][:30]} ({_d['db']})"
            for _j, (_s, _e) in enumerate(_d["spans"]):
                fig.add_trace(go.Scatter(
                    x=[_s, _e], y=[_lab, _lab], mode="lines",
                    line=dict(color=_c, width=12), showlegend=False,
                    hovertemplate=f"{_d['name']} ({_d['db']})<br>"
                                  f"{_s}–{_e}<extra></extra>"),
                    row=3, col=1)

    fig.update_layout(
        height=(460 + (_ndom * 26 if _ndom else 60)),
        margin=dict(l=4, r=4, t=34, b=40),
        title=dict(text=title, font=dict(size=13)),
        plot_bgcolor="white", showlegend=False)
    fig.update_yaxes(autorange="reversed", tickfont=dict(size=8), row=1, col=1)
    fig.update_yaxes(tickfont=dict(size=8), row=2, col=1)
    if _ndom:
        fig.update_yaxes(tickfont=dict(size=8), row=3, col=1,
                         autorange="reversed")
        fig.update_xaxes(title_text="position (ABP residue)", row=3, col=1)
    else:
        fig.update_xaxes(title_text="position (ABP residue)", row=2, col=1)
    if focus:                                    # recadrage sur la zone de liaison
        fig.update_xaxes(range=[focus[0] - 0.5, focus[1] + 0.5])
    st.plotly_chart(fig, use_container_width=True)


def _abp_actin_focus(sel_abp, slug):
    """(lo, hi, qlen) de la zone de liaison à l'actin si l'ABP est gros ET que
    l'interface y est localisée, sinon None. Sert à zoomer paysage + 3D sur la
    région pertinente (ex. Filamin-A : ~30–146 sur 2647 aa). Marge 5 % (≥ 10)."""
    iface = proteocast_view.abp_interface_asa_on_query(sel_abp, slug, 0.0)
    if not iface:
        return None
    _qs = proteocast_view._query_seq(slug)
    qlen = len(_qs) if _qs else max(iface)
    lo, hi = min(iface), max(iface)
    span = hi - lo + 1
    if qlen and span < 0.6 * qlen and (qlen - span) > 50:
        mar = max(10, int(0.05 * span))
        return (max(1, lo - mar), min(qlen, hi + mar), qlen)
    return None


def _crop_pdb_residues(pdb_text, lo, hi):
    """Ne garde que les atomes des résidus [lo, hi] (numérotation = query/UniProt
    pour le modèle AlphaFold complet). Sert à restreindre la 3D à l'ABD."""
    out = []
    for ln in pdb_text.splitlines():
        if ln.startswith(("ATOM", "HETATM", "TER")) and len(ln) >= 26:
            try:
                r = int(ln[22:26])
            except ValueError:
                out.append(ln)
                continue
            if lo <= r <= hi:
                out.append(ln)
        else:
            out.append(ln)
    return "\n".join(out)


def _render_abp_proteocast(sel_abp, abp_subunits):
    """Panneau ProteoCast pour l'ABP sélectionné dans 'Détail par ABP'."""
    slug = _re_sd.sub(r"[^a-zA-Z0-9]+", "_", str(sel_abp)).strip("_")[:60]
    csv = _find_proteocast_csv(slug)
    # empreinte de l'ABP DÉJÀ recalée sur la numérotation query ProteoCast :
    # {position query : %ASA MAX}. L'alignement par chaîne corrige les décalages
    # de numérotation PDB vs UniProt soumise (ex. Adducin).
    iface_asa = proteocast_view.abp_interface_asa_on_query(sel_abp, slug, 0.0)
    if csv is not None:
        # domaines Pfam/InterPro de l'ABP (via UniProt du manifest) — toujours affichés
        _doms = None
        _pcs = proteocast_view.load_status(0.0)
        _uni = None
        if _pcs is not None:
            _mm = _pcs[_pcs["abp_title"] == sel_abp]
            _uni = _mm.iloc[0]["uniprot"] if len(_mm) else None
        if _uni:
            _doms = proteocast_view.fetch_domains(_uni)
        # Gros ABP multi-domaines : seule une région contacte l'actin. On zoome
        # dessus par défaut (contacts ± marge) plutôt que d'afficher 2500 résidus
        # hors-sujet (ex. Filamin-A : interface ~30–146 sur 2647 aa).
        _focus = None
        _fz = _abp_actin_focus(sel_abp, slug)
        if _fz:
            _flo, _fhi, _qlen = _fz
            _whole = st.toggle(
                f"Whole protein ({_qlen} aa) — otherwise zoom on the "
                f"actin-binding region (~{_flo}–{_fhi})",
                value=False, key=f"pc_whole_{slug}",
                help="This ABP is multi-domain but only one region contacts "
                     "actin. We zoom on it by default; enable to see "
                     "the whole protein.")
            if not _whole:
                _focus = (_flo, _fhi)
        _surface = st.toggle(
            "Exposed surface only (RSA ≥ 0.2) — recompute the gradient on these positions",
            value=False, key=f"pc_surface_{slug}",
            help="Keeps only the EXPOSED residues of the ABP (RSA computed on the "
                 "real Arp/ABP structure from our data) and re-spreads the colours "
                 "over them: that is where the ABP-actin interface is.")
        _rsa_struct = (proteocast_view.abp_rsa_on_query(sel_abp, slug, 0.0)
                       if _surface else None)
        _render_proteocast_mutland(csv, iface_asa, f"ProteoCast — {sel_abp}",
                                   domains=_doms, surface_only=_surface,
                                   rsa=_rsa_struct, focus=_focus)
        st.caption(
            "**Green** heatmap: conservation (black/dark = deleterious/conserved). "
            "**ABP in contact with actin** track (green): positions of **the ABP** that "
            "touch actin, max buried %ASA. **Domains** track (Pfam/InterPro). "
            "In *exposed surface* mode, the **buried core** (RSA < 0.2) stays visible "
            "in **greyscale**; only the exposed surface is green. "
            "Interactive: zoom + hover."
        )
    else:
        st.info(
            f"ProteoCast not computed yet for **{sel_abp}**. "
            f"Drop the result into `data/proteocast/abp/{slug}.csv` "
            "(format `4.query_ProteoCast.csv`) and it will show up here."
        )


@st.cache_data(show_spinner=False)
def _load_actin_conservation():
    p = _Path("data/proteocast/conservation_vs_asa_per_position.csv")
    return pd.read_csv(p) if p.exists() else None


@st.cache_data(show_spinner=False)
def _abp_actin_footprint(sel_abp, mtime):
    """Positions canonical de l'actin contactées par cet ABP (toutes interactions)."""
    df = pd.read_csv("data/filtered/filtered_all_data.csv", low_memory=False)
    di = pd.read_csv("data/filtered/details/1.interactions.csv")[
        ["interaction_id", "chain_A_id", "chain_B_id"]]
    res = pd.read_csv("data/filtered/details/3.interface_residues.csv")
    res["canon"] = pd.to_numeric(
        res["residue_number_canon_mafft"], errors="coerce")
    m = df.merge(di, left_on=["subunit_1", "subunit_2"],
                 right_on=["chain_A_id", "chain_B_id"], how="left")
    sub = m[(m.s1_actine) & (~m.s2_actine) & (m.subunit_2_title == sel_abp)]
    canon = set()
    for _, r in sub.iterrows():
        rr = res[(res.interaction_id == r.interaction_id)
                 & (res.chain == r.subunit_1)]
        canon |= set(rr.canon.dropna().astype(int))
    return sorted(canon)


def _render_abp_actin_conservation(sel_abp):
    """Côté actin : la conservation (ProteoCast actin) des résidus que cet ABP touche."""
    import matplotlib.pyplot as plt
    cons = _load_actin_conservation()
    if cons is None:
        st.info(
            "Actin conservation unavailable (data/proteocast/conservation_vs_asa_per_position.csv).")
        return
    _cp = _Path("data/filtered/details/3.interface_residues.csv")
    fp = _abp_actin_footprint(
        sel_abp, _cp.stat().st_mtime if _cp.exists() else 0)
    if not fp:
        st.info("No actin footprint for this ABP.")
        return
    cons = cons.dropna(subset=["conservation"]).copy()
    cons["is_fp"] = cons["canon"].isin(fp)
    # « reste » = SURFACE de l'actin uniquement (résidus exposés, RSA ≥ 0.20),
    # hors empreinte. On EXCLUT le cœur enfoui (très conservé) : sinon la
    # comparaison est biaisée (surface vs cœur au lieu de surface vs surface).
    _SURF_RSA = 0.20
    cons["is_surface"] = pd.to_numeric(
        cons.get("rsa"), errors="coerce") >= _SURF_RSA
    fpv = cons[cons.is_fp]["conservation"]
    other = cons[cons.is_surface & ~cons.is_fp]["conservation"]

    c1, c2, c3 = st.columns(3)
    c1.metric("Actin residues contacted", len(fpv))
    c2.metric("Mean conservation (footprint)", f"{fpv.mean():.2f}")
    c3.metric("Mean conservation (rest of surface)", f"{other.mean():.2f}")
    st.caption(
        f"“rest of surface” = {len(other)} **exposed** actin residues (RSA ≥ {_SURF_RSA:g}) "
        "outside the footprint — the **buried core** (highly conserved) is excluded to compare "
        "surface vs surface.")
    try:
        from scipy.stats import mannwhitneyu
        _p = mannwhitneyu(fpv, other, alternative="two-sided").pvalue
        verdict = ("more conserved" if fpv.mean() >
                   other.mean() else "less conserved")
        st.caption(f"The actin residues contacted by {sel_abp} are **{verdict}** than the rest "
                   f"of the surface (Mann-Whitney p = {_p:.1e}).")
    except Exception:
        pass

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(13, 3.4),
                                   gridspec_kw={"width_ratios": [3, 1]})
    axa.plot(cons["canon"], cons["conservation"],
             color="lightgrey", lw=0.7, zorder=1)
    axa.scatter(cons.loc[cons.is_fp, "canon"], cons.loc[cons.is_fp, "conservation"],
                color="#e63946", s=14, zorder=3, label="footprint of this ABP")
    axa.set_xlabel("actin canonical position")
    axa.set_ylabel("conservation")
    axa.legend(fontsize=8)
    axa.set_title(
        f"Actin conservation — residues contacted by {sel_abp}", fontsize=10)
    axb.boxplot([other.dropna(), fpv.dropna()],
                tick_labels=["rest\n(surface)", "footprint"],
                patch_artist=True, boxprops=dict(facecolor="#cfe8e0"))
    axb.set_ylabel("conservation")
    axb.set_title("Footprint vs surface", fontsize=10)
    fig.tight_layout()
    st.pyplot(fig)
    plt.close(fig)
