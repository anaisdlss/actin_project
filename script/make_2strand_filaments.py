"""
make_2strand_filaments.py
Construit des filaments F-actine à 2 brins pour chaque combinaison unique
de clusters longitudinal + latéral observée dans les données ABP.

Combinaisons (PHOSPHATASE ignorée) :
  cofilin1   : long=0_7797_50  lat=0_7797_12
  cofilin2_wd: long=0_7797_10  lat=0_7797_12
  coronin    : long=0_7797_50  lat=0_7797_1

Usage : python -m script.make_2strand_filaments
"""

from pathlib import Path
import numpy as np
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
C70_DIR = PROJECT_ROOT / "data/filtered/details/structures_files/bfactor_c70_interface"
OUT_DIR = PROJECT_ROOT / "data/filtered/details/structures_files/filament/2strand"
OUT_DIR.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(PROJECT_ROOT))
from script.filament_by_abp import (
    _load_pairwise_atoms, _get_ca, _kabsch,
    _atom_xyz, _set_atom_xyz, _apply_helical,
    CHAIN_LABELS,
)

N_SUBUNITS_TOTAL = 14  # 7 par brin = 14 sous-unités
CHAIN_S1 = "ACEGIKM"  # indices pairs  (brin 1)
CHAIN_S2 = "BDFHJLN"  # indices impairs (brin 2)

COMBINATIONS = [
    ("cofilin1",    "Cofilin-1",          "0_7797_50"),
    ("cofilin2_wd", "Cofilin-2 / WD-rep", "0_7797_10"),
    ("coronin",     "Coronin-1B",         "0_7797_50"),
]


def build_two_strand(long_path: Path, out_path: Path,
                     n_total: int = N_SUBUNITS_TOTAL) -> None:
    """
    Construit un filament 2-brins en propageant le transform hélicoïdal n_total fois.
    Brin 1 = sous-unités paires (0,2,4,...), brin 2 = sous-unités impaires (1,3,5,...).
    C'est la vraie géométrie double-brin de F-actine : les deux brins sont
    les sous-unités alternées du même hélice.
    """
    a_long, b_long = _load_pairwise_atoms(long_path)
    ca_A, rA = _get_ca(a_long)
    ca_B, rB = _get_ca(b_long)
    common = sorted(set(rA) & set(rB))
    rA_idx = {r: i for i, r in enumerate(rA)}
    rB_idx = {r: i for i, r in enumerate(rB)}
    P = ca_A[[rA_idx[r] for r in common]]
    Q = ca_B[[rB_idx[r] for r in common]]
    R_long, t_long = _kabsch(P, Q)
    xyz_0 = np.array([_atom_xyz(l) for l in a_long])

    # Construire 14 sous-unités dans l'hélice
    subunits = []
    xyz_k = xyz_0.copy()
    for k in range(n_total):
        subunits.append(xyz_k.copy())
        xyz_k = _apply_helical(xyz_k, R_long, t_long, 1)

    # Assigner : brin 1 = 0,2,4,6,8,10,12  brin 2 = 1,3,5,7,9,11,13
    n_per_strand = n_total // 2
    out_lines = []
    serial = 1
    for brin, start in [(CHAIN_S1, 0), (CHAIN_S2, 1)]:
        for j in range(n_per_strand):
            chain_label = brin[j]
            xyz_k = subunits[start + j * 2]
            for i, line in enumerate(a_long):
                new_line = (line[:6]
                            + f"{serial:5d}"
                            + line[11:21]
                            + chain_label
                            + line[22:60]
                            + f"{0.0:6.2f}"
                            + line[66:])
                new_line = _set_atom_xyz(new_line, xyz_k[i])
                out_lines.append(new_line)
                serial += 1
            out_lines.append(f"TER   {serial:5d}      {chain_label}\n")
            serial += 1

    out_lines.append("END\n")
    out_path.write_text("".join(out_lines))


# ── main ─────────────────────────────────────────────────────────────────────
pdb_paths = []
for short, label, long_id in COMBINATIONS:
    long_f = C70_DIR / f"{long_id}.pdb"
    out    = OUT_DIR / f"filament_2s_{short}.pdb"

    if not long_f.exists():
        print(f"[MANQUANT] {long_id} — ignoré")
        continue

    print(f"{label}  (long={long_id})")
    try:
        build_two_strand(long_f, out)
        print(f"  ✓ {out.name}")
        pdb_paths.append((out, label, short))
    except Exception as e:
        print(f"  [ERREUR] {e}")
    print()

# ── Script PyMOL ─────────────────────────────────────────────────────────────
pml = OUT_DIR / "2strand_filaments.pml"
colors = {"cofilin1": "blue", "cofilin2_wd": "green", "coronin": "red"}
lines = [
    "# Filaments F-actine 2 brins — géométries distinctes par combinaison de clusters",
    "# Brin 1 = chaînes A-G  |  Brin 2 = chaînes H-N",
    "",
]
for out, label, short in pdb_paths:
    obj = f"fil_{short}"
    c1, c2 = colors.get(short, "cyan"), colors.get(short, "cyan")
    lines += [
        f"load {out.as_posix()}, {obj}",
        f"hide everything, {obj}",
        f"show surface, {obj}",
        f"spectrum chain, rainbow, {obj}",
        "",
    ]
lines += ["set surface_quality, 1", "bg_color white", "zoom all"]
pml.write_text("\n".join(lines))
print(f"PyMOL : {pml}")
print(f"\npymol {pml}")
