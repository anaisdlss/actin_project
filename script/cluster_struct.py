"""Extracted from streamlit.py — cluster_struct view/build helpers (keeps streamlit.py light)."""
import re as _re_sd
import numpy as np
import pandas as pd
import streamlit as st
from pathlib import Path as _Path
import folddisco_view


_ABP_SD = _Path("data/exports/abp_site_domain")


def _sd_slug(s):
    return _re_sd.sub(r"[^A-Za-z0-9]+", "_", str(s)).strip("_")[:50]


def _struct_mtime():
    ts = []
    for n in ["actin_footprint_overlap.csv", "interface_tm_sweep.csv",
              "abp_master.csv", "familles.csv", "whole_pairs_all.tsv",
              "interface_secondary_structure.csv", "interface_chemistry.csv",
              "interface_domain_by_abp.csv"]:
        p = _ABP_SD / n
        if p.exists():
            ts.append(p.stat().st_mtime)
    return tuple(ts)


@st.cache_data(show_spinner=False)
def _load_struct_tables(mtime):
    def rd(n, **kw):
        p = _ABP_SD / n
        return pd.read_csv(p, **kw) if p.exists() else None
    wp_p = _ABP_SD / "whole_pairs_all.tsv"
    return {
        "fp":     rd("actin_footprint_overlap.csv"),
        "sweep":  rd("interface_tm_sweep.csv"),
        "master": rd("abp_master.csv"),
        "fam":    rd("familles.csv"),
        "ss":     rd("interface_secondary_structure.csv"),
        "chem":   rd("interface_chemistry.csv"),
        "idom":   rd("interface_domain_by_abp.csv"),
        "wp":     (pd.read_csv(wp_p, sep="\t", header=None,
                               names=["q", "tt", "qtm", "ttm", "rmsd", "fident", "alnlen"])
                   if wp_p.exists() else None),
    }


def _render_cluster_struct(cid):
    """Panneau récap structural des ABP d'un site : TM interface / TM entière / %id /
    recouvrement d'empreinte sur l'actin / familles / domaines Pfam."""
    cid = str(cid)
    T = _load_struct_tables(_struct_mtime())
    fp = T["fp"]
    if fp is None:
        st.info("Structural analyses unavailable "
                "(run the scripts in `script/abp_site_domain/`).")
        return
    sub = fp[fp["site"].astype(str) == cid]
    if sub.empty:
        _avail = sorted(fp["site"].astype(str).unique(),
                        key=lambda x: [int(p) if p.isdigit() else p for p in x.split("_")])
        st.info(
            f"No structural comparison for **{cid}**: this site has only **one "
            "ABP** (or it is a homo actin-actin site). The comparison only exists "
            "for sites with **≥ 2 partner ABPs**.")
        st.caption(
            f"**{len(_avail)} clusters have a structural comparison** — pick "
            f"one of them in “S1 binding-site patch”: {', '.join(_avail)}")
        return

    master, fam, sweep, wp = T["master"], T["fam"], T["sweep"], T["wp"]
    ss = T["ss"]
    ssinfo = ss.set_index("abp_title") if ss is not None else None
    chem = T["chem"]
    cheminfo = chem.set_index("abp_title") if chem is not None else None

    fam_of = {}
    if fam is not None:
        for _, r in fam.iterrows():
            for a in str(r["membres"]).split(" · "):
                fam_of[a.strip()] = r["famille"]

    # Domaine Pfam RÉELLEMENT utilisé au contact de l'actin (par ABP)
    idom = T.get("idom")
    idom_of = (dict(zip(idom["abp_title"], idom["domaine_interface"]))
               if idom is not None else {})

    abps = sorted(set(sub["a"]) | set(sub["b"]))
    minfo = master.set_index("abp_title") if master is not None else None
    n_fam = len({fam_of.get(a, a) for a in abps})

    # Verdict « famille de domaine » basé sur le DOMAINE À L'INTERFACE
    _idoms = [idom_of.get(a, "—") for a in abps]
    _distinct = {d for d in _idoms if d and d != "—"}
    if len(_distinct) == 1 and all(d and d != "—" for d in _idoms):
        _verdict = (f"All ABPs bind actin through the **same domain** "
                    f"“{list(_distinct)[0]}” → **same domain family**.")
    elif len(_distinct) >= 2:
        _verdict = ("These ABPs bind actin through **different domains** "
                    f"({', '.join(sorted(_distinct))}) → **true convergence** "
                    "(distinct domain families at the same site).")
    elif len(_distinct) == 1:
        _verdict = (f"Common interface domain “{list(_distinct)[0]}” for the resolved "
                    "ABPs (some “—” = region outside any Pfam domain).")
    else:
        _verdict = ("Interface domain unresolved (“—”): binding regions "
                    "outside any annotated Pfam domain.")
    st.caption(
        f"**{len(abps)} ABPs** · **{n_fam} families** at this site. {_verdict}")

    # Table par ABP : protéine, famille, domaine À L'INTERFACE, SS, domaines Pfam
    rows = []
    for a in abps:
        dom = uni = ""
        if minfo is not None and a in minfo.index:
            dom = minfo.loc[a, "pfam_domains"]
            uni = minfo.loc[a, "uniprot"]
        ss_txt = "—"
        if ssinfo is not None and a in ssinfo.index:
            sr = ssinfo.loc[a]
            ss_txt = (f"{sr['SS_dominante']} "
                      f"(H{sr['pct_helice']}/β{sr['pct_brin']}/L{sr['pct_boucle']})")
        rows.append({"ABP (protein)": a, "Family": fam_of.get(a, "—"),
                     "Interface domain": idom_of.get(a, "—"),
                     "Interface (secondary struct.)": ss_txt,
                     "Pfam domains": dom, "UniProt": uni})
    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)

    # Lookups
    def _sw(a, b, level):
        if sweep is None:
            return None
        m = sweep[(sweep.level == level) &
                  (((sweep.a == a) & (sweep.b == b)) | ((sweep.a == b) & (sweep.b == a)))]
        if not len(m):
            return None
        v = m.tm.iloc[0]
        return None if pd.isna(v) else round(float(v), 2)

    stem = {}
    if minfo is not None:
        for a in abps:
            if a in minfo.index:
                stem[a] = f"{_sd_slug(a)}__{minfo.loc[a, 'pdb']}_{minfo.loc[a, 'abp_chain']}"
    wpidx = {}
    if wp is not None:
        for _, r in wp.iterrows():
            wpidx[(r.q, r.tt)] = (r.qtm, r.fident)

    def _whole(a, b):
        c = [wpidx.get((stem.get(a), stem.get(b))),
             wpidx.get((stem.get(b), stem.get(a)))]
        c = [x for x in c if x]
        return max(c, key=lambda x: x[0]) if c else (None, None)

    # FoldDisco (motif d'interface discontinu) par paire — fusionné dans le tableau
    _fd_map = folddisco_view.folddisco_pair_map(abps)

    prows = []
    for _, r in sub.iterrows():
        a, b = r["a"], r["b"]
        tmw, pid = _whole(a, b)
        tm_ent = _sw(a, b, "entière")
        if tm_ent is None and tmw is not None:
            tm_ent = round(float(tmw), 2)
        _da, _db = idom_of.get(a, "—"), idom_of.get(b, "—")
        _same_dom = ("Oui" if (_da != "—" and _da == _db)
                     else ("Non" if (_da != "—" and _db != "—") else "—"))
        # empreinte : % + nombre brut de résidus d'actin partagés (partages / plus
        # petit des deux jeux), ex. « 82% (14/17) »
        _emp = f"{100*r['recouvrement_min']:.0f}%"
        if {"partages", "n_a", "n_b"} <= set(r.index):
            _emp += f" ({int(r['partages'])}/{int(min(r['n_a'], r['n_b']))})"
        _row = {
            "ABP A": a, "ABP B": b,
            "Same domain family": _same_dom,
            "TM whole structure": tm_ent,
            "%id sequence": (f"{100*pid:.0f}%" if pid is not None else "—"),
            "Actin footprint (% shared)": _emp,
        }
        _fd = (_fd_map or {}).get(frozenset((a, b)), {})
        _sc = _fd.get("FoldDisco score")
        _sn = _fd.get("Norm. score")
        _rm = _fd.get("RMSD (Å)")
        _row["Matched residues"] = _fd.get("Matched residues", "—")
        _row["% shared motif"] = _fd.get("% shared motif", "—")
        _row["FoldDisco score"] = "—" if _sc is None else f"{_sc:.1f}"
        _row["Norm. score (0-1)"] = "—" if _sn is None else f"{_sn:.2f}"
        _row["RMSD (Å)"] = "—" if _rm is None else f"{_rm:.1f}"
        _row["Interpretation"] = _fd.get("Interpretation", "—")
        prows.append(_row)
    pair = pd.DataFrame(prows).sort_values(
        ["Same domain family", "TM whole structure"], ascending=[False, False])
    st.markdown("**Pairwise comparison of the site's ABPs**")
    # deux tableaux (sinon 11 colonnes -> les colonnes FoldDisco passent hors écran)
    _cols_struct = ["ABP A", "ABP B", "Same domain family",
                    "TM whole structure", "%id sequence",
                    "Actin footprint (% shared)"]
    _cols_fd = ["ABP A", "ABP B", "Matched residues", "% shared motif",
                "FoldDisco score", "Norm. score (0-1)", "RMSD (Å)",
                "Interpretation"]
    st.markdown("*Whole structure & footprint on actin*")
    st.dataframe(pair[_cols_struct], hide_index=True, use_container_width=True)
    st.markdown("*Interface region — 3D motif (FoldDisco)*")
    st.dataframe(pair[_cols_fd], hide_index=True, use_container_width=True)
    st.caption(
        "**Whole structure** — **Same domain family** = same Pfam domain at "
        "the contact (Yes/No; “—” = unresolved) · **TM whole structure** = same "
        "global fold (≥ 0.5) · **%id** sequence · **Actin footprint** = fraction "
        "of shared actin residues.  \n"
        "**Zone d'interface (FoldDisco)** — les deux ABP agrippent-ils l'actin avec le "
        "**same 3D residue motif**? No clear-cut verdict (no threshold separates "
        "proprement les cas) : on montre toutes les sorties de `folddisco query` (pour le "
        "**smaller of the two motifs** searched in the other — one direction, consistent). "
        "**Matched residues** = raw `n_match` (motif residues found / motif "
        "size) · **% shared motif** = same in % → **coverage** · **FoldDisco score** "
        "= raw score (rarity, no absolute scale) · **Norm. score (0-1)** = score ÷ "
        "self-match → **geometric quality** (close to 1 = very similar motifs) · "
        "**RMSD (Å)** = superposition of matched residues: **low = true 3D overlap**, "
        "**high (> ~5 Å) = loose/coincidental** (⚠ trivially low if few residues, e.g. 3/12 at "
        "0.2 Å). Truly same motif = **high % + high Norm. score + low RMSD on enough "
        "residues**. **Interpretation** = graded reading (a guide, not a verdict) combining "
        "these 4 signals: *same motif (strong)* / *largely shared but loose* / *partial* "
        "/ *different* / *few residues — inconclusive* (< 5 matched residues).  \n"
        "Reading convergence: **high footprint + different fold (low TM) + rarely "
        "shared motif** → same actin site, distinct binding solutions "
        "(**convergence positionnelle**)."
    )

    # Pooled PDB discovery for this site (motifs of all its ABPs)
    folddisco_view.render_discovery_cluster(cid)

    # ── Interface chemistry ───────────────────────────────────────────────────
    if cheminfo is not None:
        crows = []
        for a in abps:
            if a not in cheminfo.index:
                continue
            cr = cheminfo.loc[a]
            crows.append({
                "ABP": a,
                "ABP charge": f"{int(cr['charge_nette']):+d}",
                "% hydrophobic": f"{int(cr['pct_hydrophobe'])}%",
                "% aromatic": f"{int(cr['pct_aromatique'])}%",
                "Actin patch charge": ("—" if pd.isna(cr["charge_actine_patch"])
                                        else f"{int(cr['charge_actine_patch']):+d}"),
                "Area (Å²)": ("—" if pd.isna(cr["interface_area"]) else str(int(cr["interface_area"]))),
                "H-bonds": ("—" if pd.isna(cr["num_hbonds"]) else str(int(cr["num_hbonds"]))),
                "Salt bridges": ("—" if pd.isna(cr["num_salt_bridges"]) else str(int(cr["num_salt_bridges"]))),
            })
        if crows:
            st.markdown("**Interface chemistry**")
            st.dataframe(pd.DataFrame(crows), hide_index=True,
                         use_container_width=True)
            st.caption(
                "ABP / actin-patch charge = (# K,R residues) − (# D,E residues). "
                "General trend: **charged + ABP facing a charged − actin** "
                "(electrostatic complementarity) + hydrophobic core → "
                "a shared chemical strategy even across different families."
            )
