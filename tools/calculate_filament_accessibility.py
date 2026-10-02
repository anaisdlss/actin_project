#!/usr/bin/env python
"""Offline, separate SASA/RSA reference audit for local experimental 7PDZ chain I.

Keeps the legacy RSA table unchanged. Uses the same atom coordinates in three
contexts and exports observed positions plus explicit missing values to P60709.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time

import streamlit  # import package before script/ shadows its name
import Bio
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley, ATOMIC_RADII
from Bio.Data.PDBData import protein_letters_3to1_extended
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "script"))
from filament_accessibility import (MAX_ASA, exact_reference_map, protein_heavy_model,
                                    residue_profile)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def calculate_contexts(model, contexts, target, points):
    areas, stats = {}, {}
    for context, chains in contexts.items():
        clean = protein_heavy_model(model, chains)
        start = time.perf_counter()
        ShrakeRupley(probe_radius=1.4, n_points=points).compute(clean, level="R")
        elapsed = time.perf_counter() - start
        areas[context] = {residue.id: float(residue.sasa) for residue in clean[target]}
        stats[context] = {"chain_ids": chains, "atom_count": sum(1 for _ in clean.get_atoms()),
                          "elapsed_seconds": elapsed}
        print(f"{context}, {points} points: {stats[context]['atom_count']} atoms, {elapsed:.2f} s", flush=True)
    return areas, stats


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args(argv)
    root = args.root.resolve()
    pdb_id, target = "7pdz", "I"
    assembly = root / "data/filtered/details/structures_files/assembly/7pdz.pdb"
    reference_path = root / "data/P60709_ref.fasta"
    proteins_path = root / "data/filtered/proteins_per_pdb.csv"
    out = root / "reports/scientific_audit"
    out.mkdir(parents=True, exist_ok=True)
    reference = str(SeqIO.read(reference_path, "fasta").seq)
    if len(reference) != 375:
        raise ValueError("Expected the 375-residue P60709 reference")
    model = PDBParser(QUIET=True).get_structure(pdb_id, assembly)[0]
    proteins = pd.read_csv(proteins_path)
    proteins = proteins[proteins.pdb_id.str.lower().eq(pdb_id)].copy()
    proteins["author_chain"] = proteins.chain.str.split("_", n=1).str[1]
    actin_flags = proteins.is_actin.astype(str).str.lower()
    actin = sorted(set(proteins.loc[actin_flags.eq("true"), "author_chain"]))
    abp = sorted(set(proteins.loc[actin_flags.eq("false"), "author_chain"]))
    if target not in actin or not abp or len(actin) < 2:
        raise ValueError("Expected target actin, neighboring actins and separately classified capping proteins")
    contexts = {"isolated": [target], "actin_fragment": actin, "with_abp": sorted(set(actin + abp))}
    target_model = protein_heavy_model(model, [target])
    residues = list(target_model[target])
    observed = "".join(protein_letters_3to1_extended[r.resname] for r in residues)
    mapping = exact_reference_map(observed, reference)
    if len(mapping) < 300:
        raise ValueError("Insufficient sequence-mapped reference coverage")
    areas, stats = calculate_contexts(model, contexts, target, 960)
    frame = residue_profile(reference, residues, mapping, areas, pdb_id.upper(), target)
    for col in ("delta_sasa_actin_A2", "delta_sasa_abp_A2"):
        if (frame[col].dropna() < -1e-8).any():
            raise ValueError("Adding occluding protein atoms unexpectedly increased SASA")
    # A second discretization documents numerical sensitivity; it is not an
    # experimental uncertainty or an independent biological replicate.
    coarse, coarse_stats = calculate_contexts(model, contexts, target, 480)
    lower = residue_profile(reference, residues, mapping, coarse, pdb_id.upper(), target)
    convergence = []
    for context in contexts:
        a = frame[f"rsa_{context}"]
        b = lower[f"rsa_{context}"]
        delta = (a - b).abs().dropna()
        convergence.append({"context": context, "n_positions": len(delta), "points_fine": 960, "points_coarse": 480,
                            "mean_absolute_RSA_difference": float(delta.mean()), "max_absolute_RSA_difference": float(delta.max()),
                            "surface_threshold": .2, "surface_classification_disagreements": int(((a.ge(.2) != b.ge(.2)) & a.notna() & b.notna()).sum())})
    paths = [out / "filament_accessibility_7pdz_I.csv", out / "filament_accessibility_convergence.csv"]
    frame.to_csv(paths[0], index=False)
    pd.DataFrame(convergence).to_csv(paths[1], index=False)
    modified = frame[frame.coordinate_status.eq("modified_residue")]
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(), "biopython_version": Bio.__version__,
        "pdb_id": pdb_id.upper(), "model_index": 0, "target_author_chain": target,
        "reference": {"id": "UniProt P60709", "file": str(reference_path.relative_to(root)), "sha256": digest(reference_path)},
        "structure": {"file": str(assembly.relative_to(root)), "sha256": digest(assembly), "url": "https://www.rcsb.org/structure/7PDZ"},
        "chain_classification_source": {"file": str(proteins_path.relative_to(root)), "sha256": digest(proteins_path),
                                        "proteins": proteins[["chain", "protein", "is_actin"]].to_dict("records")},
        "method": "Bio.PDB.SASA.ShrakeRupley; same coordinates in all contexts; no structural relaxation, no added atoms",
        "probe_radius_A": 1.4, "sphere_points_per_atom": 960, "convergence_check_points": 480,
        "atomic_radii_A": {element: float(ATOMIC_RADII[element]) for element in sorted({a.element for a in target_model.get_atoms()})},
        "rsa_normalization": "SASA / Tien et al. (2013) theoretical maximum; Bio.Data.PDBData.residue_sasa_scales['Wilke']; not clipped to [0,1]",
        "maximum_ASA_A2": MAX_ASA,
        "mapping": "Global BLOSUM62, gap open -10, extend -0.5; every observed parent amino acid must match P60709 exactly and map identically across all optimal alignments. Author numbers are retained separately.",
        "observed_parent_residues": len(mapping), "rsa_eligible_residues": int(frame.rsa_eligible.sum()),
        "missing_coordinate_positions": frame.loc[frame.coordinate_status.eq("missing_coordinates"), "position"].tolist(),
        "modified_residues_without_RSA": modified[["position", "aa", "pdb_residue_name"]].to_dict("records"),
        "atom_policy": "Selected highest-occupancy alternate location from Bio.PDB; exclude occupancy=0, H and D. Keep amino-acid heavy atoms, including modified-parent residues, only in selected chains. Exclude waters, ions, nucleotides, non-protein ligands and phalloidin chains. RSA absent for modified or incomplete standard residues; raw observed-atom SASA retained.",
        "contexts_960": stats, "contexts_480": coarse_stats,
        "excluded_chain_ids": sorted(set(model.child_dict) - set(contexts["with_abp"])),
        "limitations": ["7PDZ is a finite capped-end experimental fragment, not a universal infinite-filament/interior reference.",
                        "Isolated chain is extracted from the filament in the same conformation, not relaxed or separately measured G-actin.",
                        "Missing atoms and excluded ligands can affect accessibility of neighboring residues; coordinates are not completed.",
                        "Discretization convergence is numerical sensitivity, not biological uncertainty.",
                        "Current app RSA uses these explicit contexts; the legacy table remains archived and is not used. This calculation does not validate biological major/minor interface classifications."],
        "primary_method_sources": ["https://biopython.org/docs/latest/api/Bio.PDB.SASA.html", "https://doi.org/10.1016/0022-2836(73)90011-9", "https://doi.org/10.1371/journal.pone.0080635"],
        "code_sha256": {name: digest(ROOT / name) for name in ("script/filament_accessibility.py", "tools/calculate_filament_accessibility.py")},
        "outputs_sha256": {str(p.relative_to(root)): digest(p) for p in paths},
        "reproduce": "python tools/calculate_filament_accessibility.py",
    }
    (out / "filament_accessibility_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(f"Wrote {len(frame)} reference positions, {len(mapping)} observed parent residues, {int(frame.rsa_eligible.sum())} RSA-eligible residues")
    print(pd.DataFrame(convergence).to_string(index=False))


if __name__ == "__main__":
    main()
