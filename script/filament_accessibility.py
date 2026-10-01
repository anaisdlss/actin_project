"""Separate, reproducible accessibility audit; never replaces legacy RSA data."""
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Data.PDBData import protein_letters_3to1, protein_letters_3to1_extended, residue_sasa_scales
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue


MAX_ASA = dict(residue_sasa_scales["Wilke"])
EXPECTED_ATOMS = {
    "ALA": "N CA C O CB", "ARG": "N CA C O CB CG CD NE CZ NH1 NH2",
    "ASN": "N CA C O CB CG OD1 ND2", "ASP": "N CA C O CB CG OD1 OD2",
    "CYS": "N CA C O CB SG", "GLN": "N CA C O CB CG CD OE1 NE2",
    "GLU": "N CA C O CB CG CD OE1 OE2", "GLY": "N CA C O",
    "HIS": "N CA C O CB CG ND1 CD2 CE1 NE2", "ILE": "N CA C O CB CG1 CG2 CD1",
    "LEU": "N CA C O CB CG CD1 CD2", "LYS": "N CA C O CB CG CD CE NZ",
    "MET": "N CA C O CB CG SD CE", "PHE": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ",
    "PRO": "N CA C O CB CG CD", "SER": "N CA C O CB OG",
    "THR": "N CA C O CB OG1 CG2", "TRP": "N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
    "TYR": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH", "VAL": "N CA C O CB CG1 CG2",
}


def exact_reference_map(observed, reference):
    """Map observed parent amino acids; reject mismatches/insertions/ambiguity.

    Deletions in coordinates remain missing; author residue numbers are never
    assumed to be reference positions. Modified residues are normalized to their
    Bio.PDB parent letter for mapping only, not for RSA normalization.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -.5
    alignments = aligner.align(reference, observed)
    if len(alignments) > 100:
        raise ValueError("Reference mapping has too many equally optimal alignments")
    maps = []
    for alignment in alignments:
        mapping = {}
        for target, query in zip(*alignment.aligned):
            for t, q in zip(range(*target), range(*query)):
                if reference[t] != observed[q]:
                    raise ValueError(f"Observed amino acid {observed[q]} does not match P60709 {reference[t]}{t + 1}")
                mapping[q] = t + 1
        maps.append(mapping)
    if not maps or any(mapping != maps[0] for mapping in maps) or len(maps[0]) != len(observed):
        raise ValueError("Observed residues do not have a complete unambiguous identity mapping to P60709")
    return maps[0]


def protein_heavy_model(source_model, chain_ids):
    """Heavy atoms of amino-acid residues, including modified parent residues.

    Selected alternate locations are those chosen by Bio.PDB (highest occupancy).
    Waters, ions, nucleotides and non-protein ligands are absent. Standard amino
    acids from unselected ligand/peptide chains are also absent.
    """
    model = Model(0)
    for chain_id in chain_ids:
        if chain_id not in source_model:
            raise ValueError(f"Selected chain {chain_id!r} is absent from the assembly")
        chain = Chain(chain_id)
        model.add(chain)
        for residue in source_model[chain_id]:
            if residue.resname not in protein_letters_3to1_extended or "CA" not in residue:
                continue
            clean = Residue(residue.id, residue.resname, residue.segid)
            # Attach the hierarchy before adding atoms: ShrakeRupley uses the
            # cached atom.full_id during residue aggregation, not get_full_id().
            # Otherwise equivalent residue IDs from different chains collide.
            chain.add(clean)
            for wrapped in residue:
                atom = wrapped.selected_child if wrapped.is_disordered() == 2 else wrapped
                if atom.element.upper() in {"H", "D"} or atom.occupancy == 0:
                    continue
                clean.add(atom.copy())
            if not len(clean):
                chain.detach_child(clean.id)
        if not len(chain):
            model.detach_child(chain.id)
    return model


def complete_standard_residue(residue):
    expected = EXPECTED_ATOMS.get(residue.resname)
    return expected is not None and set(expected.split()).issubset(set(residue.child_dict))


def residue_profile(reference, residues, mapping, areas, pdb_id, chain_id):
    """One row per reference position, NaN for missing/modified/incomplete RSA."""
    rows = [{"position": p, "aa": aa, "pdb_id": pdb_id, "chain": chain_id,
             "pdb_residue_number": None, "pdb_insertion_code": "", "pdb_residue_name": "",
             "coordinate_status": "missing_coordinates", "rsa_eligible": False,
             **{f"sasa_{context}_A2": np.nan for context in areas},
             **{f"rsa_{context}": np.nan for context in areas}} for p, aa in enumerate(reference, 1)]
    for i, residue in enumerate(residues):
        position = mapping[i]
        row = rows[position - 1]
        standard = residue.resname in protein_letters_3to1
        eligible = complete_standard_residue(residue)
        row.update(pdb_residue_number=residue.id[1], pdb_insertion_code=residue.id[2].strip(),
                   pdb_residue_name=residue.resname, rsa_eligible=eligible,
                   coordinate_status="complete_standard_residue" if eligible else "incomplete_atoms" if standard else "modified_residue")
        for context, values in areas.items():
            asa = values.get(residue.id, np.nan)
            row[f"sasa_{context}_A2"] = asa
            if eligible:
                row[f"rsa_{context}"] = asa / MAX_ASA[residue.resname]
    frame = pd.DataFrame(rows)
    frame["pdb_residue_number"] = frame.pdb_residue_number.astype("Int64")
    if {"isolated", "actin_fragment"}.issubset(areas):
        frame["delta_sasa_actin_A2"] = frame.sasa_isolated_A2 - frame.sasa_actin_fragment_A2
        frame["delta_rsa_actin"] = frame.rsa_isolated - frame.rsa_actin_fragment
    if {"actin_fragment", "with_abp"}.issubset(areas):
        frame["delta_sasa_abp_A2"] = frame.sasa_actin_fragment_A2 - frame.sasa_with_abp_A2
        frame["delta_rsa_abp"] = frame.rsa_actin_fragment - frame.rsa_with_abp
    return frame


def render_filament_accessibility(root=".", key_prefix="filament_rsa"):
    import streamlit as st
    import plotly.graph_objects as go
    from plot_interaction import position_hover

    base = Path(root) / "reports/scientific_audit"
    table = base / "filament_accessibility_7pdz_I.csv"
    with st.expander("Reference accessibility: isolated actin, actin fragment and capping proteins"):
        if not table.exists():
            st.info("The reference accessibility calculation is not installed.")
            return
        frame = pd.read_csv(table)
        st.caption("7PDZ chain I, same experimental coordinates in all three calculations. Isolated means this chain "
                   "removed from its neighbors, not a separately determined or relaxed G-actin. The finite fragment contains "
                   "six actins; the additional ABP context includes the two capping-protein chains. This capped-end "
                   "reference does not represent every position in an infinite filament.")
        st.caption("Shrake–Rupley, 1.4 Å probe, 960 sphere points per atom; heavy protein atoms only. "
                   "RSA = SASA / Tien et al. theoretical maximum for the residue (Biopython Wilke scale), without clipping. "
                   "Missing coordinates, incomplete standard residues and modified residues have no RSA. "
                   "H73 is modified (HIC): its atoms contribute to occlusion and its raw SASA is retained. "
                   "Nucleotides, ions, waters and phalloidin are excluded. This separate audit does not replace the old RSA table.")
        fig = go.Figure()
        for col, label, color in [("rsa_isolated", "Isolated chain, same conformation", "#777777"),
                                  ("rsa_actin_fragment", "Six-actin fragment only", "#E69F00"),
                                  ("rsa_with_abp", "Fragment + capping proteins", "#0072B2")]:
            fig.add_trace(go.Scatter(x=frame.position, y=frame[col], mode="lines", name=label,
                                    line=dict(color=color), connectgaps=False, customdata=frame.aa,
                                    hovertemplate="%{customdata}%{x}: %{y:.3f}<extra>%{fullData.name}</extra>"))
        fig.update_layout(xaxis_title="P60709 position", yaxis_title="Relative solvent accessibility", height=400)
        st.plotly_chart(position_hover(fig), use_container_width=True, key=f"{key_prefix}_profile")
        cutoff = st.slider("Surface RSA threshold for this reference comparison", 0.0, 1.0, .2, .01, key=f"{key_prefix}_cutoff")
        valid = frame[frame.rsa_eligible]
        st.dataframe(pd.DataFrame([
            {"Context": label, "Residues with RSA": int(valid[col].notna().sum()), "Residues at or above threshold": int(valid[col].ge(cutoff).sum())}
            for col, label in [("rsa_isolated", "Isolated"), ("rsa_actin_fragment", "Actin fragment"), ("rsa_with_abp", "Fragment + capping proteins")]
        ]), hide_index=True, width="stretch")
        st.dataframe(frame, hide_index=True, width="stretch")
        st.download_button("Download separate reference SASA/RSA calculation", table.read_bytes(), file_name=table.name, mime="text/csv", key=f"{key_prefix}_csv")
        for name, label in [("filament_accessibility_manifest.json", "Download accessibility methods and provenance"),
                            ("filament_accessibility_convergence.csv", "Download numerical convergence check")]:
            path = base / name
            if path.exists():
                st.download_button(label, path.read_bytes(), file_name=name, key=f"{key_prefix}_{name}")


def render_representative_geometry(root=".", key_prefix="geometry_audit"):
    import streamlit as st

    base = Path(root) / "reports/scientific_audit"
    summary = base / "representative_geometry_summary.csv"
    with st.expander("Representative actin-pair geometry audit"):
        st.caption("Illustrative comparison of 3J8A, 5YU8 and 6VAO. One actin is fitted by its mapped Cα atoms; "
                   "the neighbor is measured under the SAME transformation, without a second fit. At least 300 Cα "
                   "per subunit are required. Both reference-chain assignments are tried and the smaller neighbor RMSD is retained. "
                   "RMSD is in Å. These examples do not assign biological major/minor labels, test steric clashes, "
                   "or validate every cluster. 3J8A also contains tropomyosin.")
        if not summary.exists():
            st.info("The representative geometry audit is not installed.")
            return
        frame = pd.read_csv(summary)
        st.dataframe(frame, hide_index=True, width="stretch")
        for name in (summary.name, "representative_geometry.csv", "representative_geometry_manifest.json"):
            path = base / name
            if path.exists():
                st.download_button(f"Download {name}", path.read_bytes(), file_name=name, key=f"{key_prefix}_{name}")
