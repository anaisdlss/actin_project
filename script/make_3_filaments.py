"""
make_3_filaments.py
Génère 3 filaments F-actine 13 sous-unités correspondant aux 3 géométries
hélicales longitudinales distinctes trouvées dans les données :

  filament_c70_0_7797_50  — cluster Cofilin-1   (dist=42.2Å, angle=164°)
  filament_c70_0_7797_10  — cluster Cofilin-2   (dist=40.7Å, angle=164°)
  filament_c70_0_7797_11  — cluster PHOSPHATASE (dist=44.4Å, angle=158°)

Usage :
    python -m script.make_3_filaments
"""

from pathlib import Path
import numpy as np
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
C70_DIR  = PROJECT_ROOT / "data/filtered/details/structures_files/bfactor_c70_interface"
OUT_DIR  = PROJECT_ROOT / "data/filtered/details/structures_files/filament/3_geometries"
OUT_DIR.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(PROJECT_ROOT))
from script.filament_by_abp import (
    build_filament_from_pairwise,
    _helix_angle_and_longitudinal,
    _get_ca,
    _load_pairwise_atoms,
    _kabsch,
    CHAIN_LABELS,
    N_SUBUNITS,
)

CLUSTERS = [
    ("0_7797_50", "Cofilin-1",    "cofilin"),
    ("0_7797_10", "Cofilin-2 / WD-repeat", "cofilin2_wd"),
    ("0_7797_11", "PHOSPHATASE",  "phosphatase"),
]


def measure_helix(path: Path) -> tuple[float, float]:
    a_lines, b_lines = _load_pairwise_atoms(path)
    ca_A, rA = _get_ca(a_lines)
    ca_B, rB = _get_ca(b_lines)
    common = sorted(set(rA) & set(rB))
    cA, cB = ca_A.mean(0), ca_B.mean(0)
    dist = float(np.linalg.norm(cB - cA))
    rA_idx = {r: i for i, r in enumerate(rA)}
    rB_idx = {r: i for i, r in enumerate(rB)}
    P = ca_A[[rA_idx[r] for r in common]]
    Q = ca_B[[rB_idx[r] for r in common]]
    R, t = _kabsch(P, Q)
    angle = float(np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))))
    return dist, angle


pdb_paths = []
for cluster_id, label, short in CLUSTERS:
    src = C70_DIR / f"{cluster_id}.pdb"
    out = OUT_DIR / f"filament_{short}.pdb"

    if not src.exists():
        print(f"[MANQUANT] {src.name} — ignoré")
        continue

    dist, angle = measure_helix(src)
    print(f"{cluster_id}  ({label})")
    print(f"  dist centres = {dist:.1f} Å   angle = {angle:.1f}°")

    if not out.exists():
        print(f"  → construction {out.name} …")
        try:
            build_filament_from_pairwise(src, None, out, n_subunits=14)
            print(f"  ✓ {out.name}")
        except Exception as e:
            print(f"  [ERREUR] {e}")
            continue
    else:
        print(f"  {out.name} déjà existant")

    pdb_paths.append((out, label, short))
    print()

# ── Script PyMOL ─────────────────────────────────────────────────────────────
NORMAL_PDB = PROJECT_ROOT / "data/filtered/details/structures_files/filament/filament_global_base.pdb"

pml_path = OUT_DIR / "3_filaments.pml"
lines = [
    "# 3 filaments F-actine — géométries hélicales distinctes",
    "# filament_normal      : cluster 0_7797_0   (bare actin,  166.9°, rise=27.8Å)",
    "# filament_cofilin    : cluster 0_7797_50  (Cofilin-1,   dist=42.2Å, 164°)",
    "# filament_cofilin2_wd: cluster 0_7797_10  (Cofilin-2,   dist=40.7Å, 164°)",
    "# filament_phosphatase: cluster 0_7797_11  (PHOSPHATASE, dist=44.4Å, 158°)",
    "",
    f"load {NORMAL_PDB.as_posix()}, filament_normal",
    "hide everything, filament_normal",
    "show surface, filament_normal",
    "color white, filament_normal",
    "",
]
colors = ["blue", "green", "red"]
for (pdb, label, short), color in zip(pdb_paths, colors):
    obj = f"filament_{short}"
    lines += [
        f"load {pdb.as_posix()}, {obj}",
        f"hide everything, {obj}",
        f"show surface, {obj}",
        f"color {color}, {obj}",
        "",
    ]

lines += [
    "set surface_quality, 1",
    "bg_color white",
    "zoom all",
]
pml_path.write_text("\n".join(lines))
print(f"Script PyMOL : {pml_path}")
print("\nPour visualiser :")
print(f"  pymol {pml_path}")
