#!/usr/bin/env python
"""Offline SASA/RSA calculation for human beta-actin, 8DNH chain B.

The same deposited coordinates are used alone and in the four-actin fragment.
Acquire the identified sources with tools/fetch_rsa_reference.py before rebuilding.
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
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from fetch_rsa_reference import validate_entity
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
        clean = protein_heavy_model(model, chains, include_acetyl_caps=True)
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
    pdb_id, target = "8dnh", "B"
    assembly = root / "data/reference_structures/8dnh.cif"
    entity_path = root / "data/reference_structures/8dnh_entity.json"
    source_path = root / "data/reference_structures/8dnh_source.json"
    reference_path = root / "data/P60709_ref.fasta"
    out = root / "reports/scientific_audit"
    out.mkdir(parents=True, exist_ok=True)
    source = json.loads(source_path.read_text())
    for name, item in source['files'].items():
        if digest(assembly.parent / name) != item['sha256']:
            raise ValueError(f"Reference source changed since download: {name}")
    entity = json.loads(entity_path.read_text())
    validate_entity(entity)
    reference = str(SeqIO.read(reference_path, "fasta").seq)
    if len(reference) != 375:
        raise ValueError("Expected the 375-residue P60709 reference")
    model = MMCIFParser(QUIET=True).get_structure(pdb_id, assembly)[0]
    actin = sorted(entity['rcsb_polymer_entity_container_identifiers']['auth_asym_ids'])
    # Four deposited subunits, no ABP and no extrapolation by helical symmetry.
    contexts = {"isolated": [target], "actin_fragment": actin}
    target_model = protein_heavy_model(model, [target], include_acetyl_caps=True)
    residues = [r for r in target_model[target] if r.resname in protein_letters_3to1_extended and "CA" in r]
    observed = "".join(protein_letters_3to1_extended[r.resname] for r in residues)
    # The deposition's mutation flag describes N-terminal acetylation. Preserve
    # that annotation and require the expected modification; do not call it WT.
    cif = MMCIF2Dict(str(assembly))
    if set(cif.get('_struct_ref_seq_dif.details', [])) != {'acetylation'}:
        raise ValueError("Review changed sequence-difference annotations before recalculating")
    if not any(r.resname == 'ACE' for r in target_model[target]):
        raise ValueError("Expected the experimentally modelled N-terminal acetyl cap")
    mapping = exact_reference_map(observed, reference)
    if len(mapping) < 300:
        raise ValueError("Insufficient sequence-mapped reference coverage")
    if mapping[0] != 2 or residues[0].resname != 'ASP':
        raise ValueError("Expected acetylated N-terminal Asp2 in P60709 numbering")
    for chain in actin:
        other = [r for r in model[chain] if r.resname in protein_letters_3to1_extended and 'CA' in r]
        exact_reference_map(''.join(protein_letters_3to1_extended[r.resname] for r in other), reference)
    areas, stats = calculate_contexts(model, contexts, target, 960)
    frame = residue_profile(reference, residues, mapping, areas, pdb_id.upper(), target, modified_positions={2})
    for col in ("delta_sasa_actin_A2",):
        if (frame[col].dropna() < -1e-8).any():
            raise ValueError("Adding occluding protein atoms unexpectedly increased SASA")
    # A second discretization documents numerical sensitivity; it is not an
    # experimental uncertainty or an independent biological replicate.
    coarse, coarse_stats = calculate_contexts(model, contexts, target, 480)
    lower = residue_profile(reference, residues, mapping, coarse, pdb_id.upper(), target, modified_positions={2})
    convergence = []
    for context in contexts:
        a = frame[f"rsa_{context}"]
        b = lower[f"rsa_{context}"]
        delta = (a - b).abs().dropna()
        convergence.append({"context": context, "n_positions": len(delta), "points_fine": 960, "points_coarse": 480,
                            "mean_absolute_RSA_difference": float(delta.mean()), "max_absolute_RSA_difference": float(delta.max()),
                            "surface_threshold": .2, "surface_classification_disagreements": int(((a.ge(.2) != b.ge(.2)) & a.notna() & b.notna()).sum())})
    paths = [out / "filament_accessibility_8dnh_B.csv", out / "filament_accessibility_convergence.csv"]
    frame.to_csv(paths[0], index=False)
    pd.DataFrame(convergence).to_csv(paths[1], index=False)
    modified = frame[frame.coordinate_status.eq("modified_residue")]
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(), "biopython_version": Bio.__version__,
        "pdb_id": pdb_id.upper(), "model_index": 0, "target_author_chain": target,
        "reference": {"id": "UniProt P60709", "file": str(reference_path.relative_to(root)), "sha256": digest(reference_path)},
        "structure": {"file": str(assembly.relative_to(root)), "sha256": digest(assembly), "url": "https://www.rcsb.org/structure/8DNH"},
        "entity_metadata": {"file": str(entity_path.relative_to(root)), "sha256": digest(entity_path)},
        "download_receipt": {"file": str(source_path.relative_to(root)), "sha256": digest(source_path)},
        "organism": "Homo sapiens", "taxonomy_id": 9606, "protein": "ACTB / beta-actin", "uniprot": "P60709",
        "deposition_mutation_count": entity['entity_poly']['rcsb_mutation_count'],
        "sequence_difference_annotation": "N-terminal acetylation (ACE); every observed parent amino acid matches P60709",
        "target_chain_choice": "B is one of the two middle subunits in the deposited four-actin fragment; no missing neighbor is generated.",
        "method": "Bio.PDB.SASA.ShrakeRupley; same coordinates in all contexts; no structural relaxation, no added atoms",
        "probe_radius_A": 1.4, "sphere_points_per_atom": 960, "convergence_check_points": 480,
        "atomic_radii_A": {element: float(ATOMIC_RADII[element]) for element in sorted({a.element for a in target_model.get_atoms()})},
        "rsa_normalization": "SASA / Tien et al. (2013) theoretical maximum; Bio.Data.PDBData.residue_sasa_scales['Wilke']; not clipped to [0,1]",
        "maximum_ASA_A2": MAX_ASA,
        "mapping": "Global BLOSUM62, gap open -10, extend -0.5; every observed parent amino acid must match P60709 exactly and map identically across all optimal alignments. Author numbers are retained separately.",
        "observed_parent_residues": len(mapping), "rsa_eligible_residues": int(frame.rsa_eligible.sum()),
        "missing_coordinate_positions": frame.loc[frame.coordinate_status.eq("missing_coordinates"), "position"].tolist(),
        "modified_residues_without_RSA": modified[["position", "aa", "pdb_residue_name"]].to_dict("records"),
        "atom_policy": "Selected highest-occupancy alternate location from Bio.PDB; exclude occupancy=0, H and D. Keep amino-acid heavy atoms, HIC and covalently attached ACE caps in selected chains. Exclude waters, ions and nucleotides. RSA absent for incomplete residues, acetylated Asp2 and methylated His73; raw observed-residue SASA retained. ACE occludes solvent but is not mapped to a reference residue or assigned RSA.",
        "contexts_960": stats, "contexts_480": coarse_stats,
        "excluded_chain_ids": sorted(set(model.child_dict) - set(contexts["actin_fragment"])),
        "limitations": ["8DNH is human recombinant beta-actin expressed in yeast, with deposited post-translational modifications. Its four-subunit fragment is not an infinite-filament reference.",
                        "Isolated chain is extracted from the filament in the same conformation, not relaxed or separately measured G-actin.",
                        "Missing atoms and excluded ligands can affect accessibility of neighboring residues; coordinates are not completed.",
                        "Discretization convergence is numerical sensitivity, not biological uncertainty.",
                        "Current app RSA uses these explicit contexts; the legacy table remains archived and is not used. This calculation does not validate biological major/minor interface classifications."],
        "primary_method_sources": ["https://biopython.org/docs/latest/api/Bio.PDB.SASA.html", "https://doi.org/10.1016/0022-2836(73)90011-9", "https://doi.org/10.1371/journal.pone.0080635"],
        "code_sha256": {name: digest(ROOT / name) for name in ("script/filament_accessibility.py", "tools/calculate_filament_accessibility.py", "tools/fetch_rsa_reference.py")},
        "outputs_sha256": {str(p.relative_to(root)): digest(p) for p in paths},
        "reproduce": "python tools/calculate_filament_accessibility.py",
    }
    (out / "filament_accessibility_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(f"Wrote {len(frame)} reference positions, {len(mapping)} observed parent residues, {int(frame.rsa_eligible.sum())} RSA-eligible residues")
    print(pd.DataFrame(convergence).to_string(index=False))


if __name__ == "__main__":
    main()
