from display_helpers import viewer_html
"""Extracted from streamlit.py — abp_3d view/build helpers (keeps streamlit.py light)."""
import os
import numpy as np
import pandas as pd
import streamlit as st
from pathlib import Path as _Path


_BFAC_C70_DIR = _Path(
    "data/filtered/details/structures_files/bfactor_c70_interface")


def _s1_abp_3d_options(patch, mtime):
    sources = [_Path("data/filtered/filtered_all_data.csv")] + [
        _Path("data/filtered/details") / name for name in
        ("1.interactions.csv", "3.interface_residues.csv", "8.structures.csv")]
    for folder in ("assembly", "pairwise"):
        sources.extend(sorted((_Path("data/filtered/details/structures_files") / folder).glob("*.pdb")))
    signature = tuple((str(p), p.stat().st_mtime_ns, p.stat().st_size) for p in sources if p.exists())
    return _cached_s1_abp_options(patch, (mtime, signature))


@st.cache_data(show_spinner=False)
def _cached_s1_abp_options(patch, source_signature):
    """Use observed chain pairs even when pre-generated cluster structures are absent."""
    import hashlib
    import tempfile
    from structure_utils import pairwise_text
    data = pd.read_csv("data/filtered/filtered_all_data.csv", low_memory=False)
    pairs = pd.read_csv("data/filtered/details/1.interactions.csv")
    structures = pd.read_csv("data/filtered/details/8.structures.csv").set_index("interaction_id")
    residues = pd.read_csv("data/filtered/details/3.interface_residues.csv", dtype={"residue_number_structure": str})
    first = data.s1_binding_site_cluster_data_70.astype(str).eq(str(patch)) & data.s1_actine.eq(True)
    second = data.s2_binding_site_cluster_data_70.astype(str).eq(str(patch)) & data.s2_actine.eq(True)
    subset = data[first | second].copy().sort_values("area", ascending=False)
    options, seen = [], set()
    cache = _Path(tempfile.gettempdir()) / "actin-observed-pairs"
    cache.mkdir(exist_ok=True)
    for _, row in subset.iterrows():
        actin_side = 1 if str(row.s1_binding_site_cluster_data_70) == str(patch) and row.s1_actine else 2
        partner_side = 3 - actin_side
        if row[f"s{partner_side}_actine"]:
            continue
        actin, partner = row[f"subunit_{actin_side}"], row[f"subunit_{partner_side}"]
        label = str(row[f"subunit_{partner_side}_title"])
        group = (str(row.cluster_data_70), label)
        if group in seen:
            continue
        matches = pairs[((pairs.chain_A_id == actin) & (pairs.chain_B_id == partner)) |
                        ((pairs.chain_B_id == actin) & (pairs.chain_A_id == partner))]
        for _, pair in matches.iterrows():
            iid = pair.interaction_id
            if iid not in structures.index:
                continue
            assembly = structures.loc[iid, "biological_assembly_pdb_file"]
            if isinstance(assembly, pd.Series):
                assembly = assembly.iloc[0]
            text = pairwise_text(iid, pair.pdb_id, pair.chain_A_id.split("_", 1)[1],
                                 pair.chain_B_id.split("_", 1)[1],
                                 assembly_pdb=assembly if pd.notna(assembly) else None)
            if not text:
                continue
            contact = residues[residues.interaction_id == iid]
            values = {}
            for chain, output in [(actin, "A"), (partner, "B")]:
                sub = contact[contact.chain == chain]
                for pos, val in zip(sub.residue_number_structure, pd.to_numeric(sub.buried_ASA_percent, errors="coerce")):
                    if pd.notna(val):
                        values[(output, str(pos).strip())] = float(val)
            actin_original = "A" if pair.chain_A_id == actin else "B"
            output, present = [], set()
            for line in text.splitlines():
                if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
                    chain = "A" if line[21] == actin_original else "B"
                    position = line[22:27].strip()
                    value = values.get((chain, position), 0.0)
                    line = line[:21] + chain + line[22:60] + f"{value:6.2f}" + line[66:]
                    present.add(chain)
                    output.append(line)
            if present != {"A", "B"}:
                continue
            content = "\n".join(output + ["END", ""])
            target = cache / (hashlib.sha256(content.encode()).hexdigest() + ".pdb")
            if not target.exists():
                target.write_text(content)
            options.append({"label": f"{label} · {str(pair.pdb_id).upper()} · C70 {group[0]}",
                            "pdb": str(target), "interaction_id": int(iid)})
            seen.add(group)
            break
    return sorted(options, key=lambda item: item["label"].casefold())


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
    _grn = ["#F7FBFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#2171B5", "#08306B"]
    v = py3Dmol.view(width="100%", height=470)
    v.addModel(txt, "pdb")
    v.setStyle({}, {})
    v.addSurface(py3Dmol.SES,
                 {"opacity": 1, "colorscheme": {"prop": "b", "gradient": "linear",
                  "colors": _red, "min": 0, "max": 100.0}},
                 {"chain": actin_ch})
    v.addSurface(py3Dmol.SES,
                 {"opacity": 1, "colorscheme": {"prop": "b", "gradient": "linear",
                  "colors": _grn, "min": 0, "max": 100.0}},
                 {"chain": abp_ch})
    v.setBackgroundColor("white")
    v.zoomTo()
    import re as _re_fog
    html = viewer_html(v)
    html2 = _re_fog.sub(r'(viewer_\w+)\.render\(\)',
                        r'\1.setFogParameters({density:0});\1.render()', html, count=1)
    return (html2 if html2 != html else html), bmax[actin_ch], bmax[abp_ch]


_ABP_MULTI_COLORS = ["#0072B2", "#E69F00", "#CC79A7", "#56B4E9", "#332288",
                     "#AA4499", "#666666", "#DDCC77", "#1177AA", "#AA7744"]


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
    from scientific_analysis import isoform_map
    from Bio.SeqUtils import seq1
    def resolved(chain):
        return [res for res in chain if "CA" in res and seq1(res.resname, custom_map={"HIC": "H"}) != "X"]
    ref_res = resolved(ract)
    ref_seq = "".join(seq1(res.resname, custom_map={"HIC": "H"}) for res in ref_res)
    ids = iter("BCDEFGHIJKLMNOPQRSTUVWXYZ0123456789")
    chain_map = {}
    for p, lbl in zip(pdb_paths, labels):
        try:
            m = next(iter(parser.get_structure("m", p)))
            ac = _actin_chain(m)
            moving = resolved(ac)
            sequence = "".join(seq1(res.resname, custom_map={"HIC": "H"}) for res in moving)
            mapping = isoform_map(sequence, ref_seq)
            pairs = [(moving[i-1]["CA"], ref_res[j-1]["CA"]) for i,j in mapping.items()
                     if sequence[i-1] == ref_seq[j-1]]
            if len(pairs) < 100 or len(pairs) / min(len(moving), len(ref_res)) < .8:
                continue
            sup = Superimposer()
            sup.set_atoms([b for _, b in pairs], [a for a, _ in pairs])
            if sup.rms > 5.0:
                continue
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
    v = py3Dmol.view(width="100%", height=470)
    v.addModel(pdb_content, "pdb")
    v.setStyle({}, {})
    # Neutral reference actin: its first pair's ASA must not imply a pooled footprint.
    v.addSurface(py3Dmol.SES, {"opacity": 1, "color": "#D9D9D9"}, {"chain": actin_ch})
    for i, ch in enumerate(chain_map):
        v.addSurface(py3Dmol.SES,
                     {"opacity": 0.85,
                         "color": _ABP_MULTI_COLORS[i % len(_ABP_MULTI_COLORS)]},
                     {"chain": ch})
    v.setBackgroundColor("white")
    v.zoomTo()
    import re as _re_fog
    html = viewer_html(v)
    html2 = _re_fog.sub(r'(viewer_\w+)\.render\(\)',
                        r'\1.setFogParameters({density:0});\1.render()', html, count=1)
    return (html2 if html2 != html else html), bmax
