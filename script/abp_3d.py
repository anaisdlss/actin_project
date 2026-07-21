"""Extracted from streamlit.py — abp_3d view/build helpers (keeps streamlit.py light)."""
import os
import numpy as np
import pandas as pd
import streamlit as st
from pathlib import Path as _Path


_BFAC_C70_DIR = _Path(
    "data/filtered/details/structures_files/bfactor_c70_interface")


@st.cache_data(show_spinner=False)
def _s1_abp_3d_options(patch, mtime):
    """Pour un cluster S1 : liste (label, pdb) des ABP (sous-clusters C70) disponibles en 3D."""
    df = pd.read_csv("data/filtered/filtered_all_data.csv", low_memory=False)
    sub = df[df["s1_binding_site_cluster_data_70"].astype(str) == str(patch)].dropna(
        subset=["cluster_data_70"])
    opts, seen = [], set()
    for _, r in sub.iterrows():
        c70 = str(r["cluster_data_70"])
        if c70 in seen:
            continue
        seen.add(c70)
        p = _BFAC_C70_DIR / f"{c70}.pdb"
        if p.exists():
            opts.append({"label": f"{r['subunit_2_title']}  ·  C70 {c70}",
                         "pdb": str(p)})
    return sorted(opts, key=lambda o: o["label"])


def _build_abp_actin_3d(pdb_path):
    """Vue 3D actin (chaîne majoritaire, gradient rouge = %ASA) + ABP (autre chaîne,
    gradient vert = intensité de contact). b-factors des deux côtés."""
    import py3Dmol
    from collections import Counter, defaultdict
    txt = _Path(pdb_path).read_text()
    ca = Counter()
    bmax = defaultdict(float)
    for ln in txt.splitlines():
        if ln.startswith("ATOM") and len(ln) > 65:
            ch = ln[21]
            if ln[12:16].strip() == "CA":
                ca[ch] += 1
            try:
                bmax[ch] = max(bmax[ch], float(ln[60:66]))
            except ValueError:
                pass
    if len(ca) < 2:
        return None
    # convention des PDB bfactor_c70_interface : chaîne A = actin, chaîne B = partenaire
    actin_ch = "A" if "A" in ca else ca.most_common(1)[0][0]
    abp_ch = max((c for c in ca if c != actin_ch), key=lambda c: ca[c])
    _red = ["#FFFFCC", "#FEE186", "#FDAA48", "#FC5A2D", "#D30F20", "#800026"]
    _grn = ["#FFFFFF", "#E5F5E0", "#A1D99B", "#41AB5D", "#238B45", "#00441B"]
    v = py3Dmol.view(width=560, height=470)
    v.addModel(txt, "pdb")
    v.setStyle({}, {})
    v.addSurface(py3Dmol.SES,
                 {"opacity": 1, "colorscheme": {"prop": "b", "gradient": "linear",
                  "colors": _red, "min": 0, "max": max(bmax[actin_ch], 1.0)}},
                 {"chain": actin_ch})
    v.addSurface(py3Dmol.SES,
                 {"opacity": 1, "colorscheme": {"prop": "b", "gradient": "linear",
                  "colors": _grn, "min": 0, "max": max(bmax[abp_ch], 1.0)}},
                 {"chain": abp_ch})
    v.setBackgroundColor("white")
    v.zoomTo()
    import re as _re_fog
    html = v._make_html()
    html2 = _re_fog.sub(r'(viewer_\w+)\.render\(\)',
                        r'\1.setFogParameters({density:0});\1.render()', html, count=1)
    return (html2 if html2 != html else html), bmax[actin_ch], bmax[abp_ch]


_ABP_MULTI_COLORS = ["#1f9e3a", "#4363d8", "#f58231", "#911eb4", "#e6194B",
                     "#42d4f4", "#f032e6", "#9A6324", "#808000", "#000075"]


@st.cache_data(show_spinner="Superimposing ABPs on actin…")
def _build_all_abp_pdb(pdb_paths, labels):
    """Superpose tous les complexes ABP sur une actin commune (biopython) et renvoie
    (pdb_combiné, actin_chain='A', {chain_id: label})."""
    from Bio.PDB import PDBParser, Superimposer, PDBIO
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model
    from io import StringIO
    parser = PDBParser(QUIET=True)

    def _actin_chain(model):
        # convention : chaîne A = actin ; sinon plus grande chaîne
        for ch in model:
            if ch.id == "A":
                return ch
        best, bn = None, -1
        for ch in model:
            n = sum(1 for r in ch if "CA" in r)
            if n > bn:
                bn, best = n, ch
        return best

    def _abp_chain(model, actin):
        best, bn = None, -1
        for ch in model:
            if ch.id == actin.id:
                continue
            n = sum(1 for r in ch if "CA" in r)
            if n > bn:
                bn, best = n, ch
        return best

    comb = Structure("c")
    cm = Model(0)
    comb.add(cm)
    ref = parser.get_structure("r", pdb_paths[0])
    rm = next(iter(ref))
    ract = _actin_chain(rm)
    ract.id = "A"
    cm.add(ract.copy())
    ref_ca = {r.id[1]: r["CA"] for r in ract if "CA" in r}
    ids = iter("BCDEFGHIJKLMNOPQRSTUVWXYZ0123456789")
    chain_map = {}
    for p, lbl in zip(pdb_paths, labels):
        try:
            m = next(iter(parser.get_structure("m", p)))
            ac = _actin_chain(m)
            pairs = [(r["CA"], ref_ca[r.id[1]]) for r in ac
                     if "CA" in r and r.id[1] in ref_ca]
            if len(pairs) < 3:
                continue
            sup = Superimposer()
            sup.set_atoms([b for _, b in pairs], [a for a, _ in pairs])
            sup.apply(list(m.get_atoms()))
            abp = _abp_chain(m, ac)
            if abp is None:
                continue
            nc = next(ids)
            ab2 = abp.copy()
            ab2.id = nc
            cm.add(ab2)
            chain_map[nc] = lbl
        except Exception:
            continue
    io = PDBIO()
    sio = StringIO()
    io.set_structure(comb)
    io.save(sio)
    return sio.getvalue(), "A", chain_map


def _build_all_abp_3d(pdb_content, actin_ch, chain_map):
    """Actin (%ASA rouge) + tous les ABP superposés, une couleur distincte par ABP."""
    import py3Dmol
    _red = ["#FFFFCC", "#FEE186", "#FDAA48", "#FC5A2D", "#D30F20", "#800026"]
    bmax = 0.0
    for ln in pdb_content.splitlines():
        if ln.startswith("ATOM") and len(ln) > 65 and ln[21] == actin_ch:
            try:
                bmax = max(bmax, float(ln[60:66]))
            except ValueError:
                pass
    v = py3Dmol.view(width=560, height=470)
    v.addModel(pdb_content, "pdb")
    v.setStyle({}, {})
    v.addSurface(py3Dmol.SES,
                 {"opacity": 1, "colorscheme": {"prop": "b", "gradient": "linear",
                  "colors": _red, "min": 0, "max": max(bmax, 1.0)}}, {"chain": actin_ch})
    for i, ch in enumerate(chain_map):
        v.addSurface(py3Dmol.SES,
                     {"opacity": 0.85,
                         "color": _ABP_MULTI_COLORS[i % len(_ABP_MULTI_COLORS)]},
                     {"chain": ch})
    v.setBackgroundColor("white")
    v.zoomTo()
    import re as _re_fog
    html = v._make_html()
    html2 = _re_fog.sub(r'(viewer_\w+)\.render\(\)',
                        r'\1.setFogParameters({density:0});\1.render()', html, count=1)
    return (html2 if html2 != html else html), bmax
