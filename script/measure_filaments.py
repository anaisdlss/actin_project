"""
measure_filaments.py
Charge pymol_session/filaments.pse et mesure la longueur axiale (Z) de chaque objet.

Usage :
    python -m script.measure_filaments
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SESSION_PATH = PROJECT_ROOT / "pymol_session" / "filaments.pse"
OUT_PATH     = PROJECT_ROOT / "pymol_session" / "filament_lengths.txt"

try:
    import pymol
    from pymol import cmd
except ImportError:
    sys.exit("PyMOL non trouvé. Active l'environnement pixi : pixi run python -m script.measure_filaments")

import numpy as np

pymol.finish_launching(["pymol", "-cq"])   # mode headless, sans fenêtre
cmd.load(str(SESSION_PATH))

results = []
for obj in sorted(cmd.get_object_list("all")):
    coords = cmd.get_coords("%s and name CA" % obj, 1)
    if coords is None or len(coords) == 0:
        continue
    z_vals = coords[:, 2]
    length_A = z_vals.max() - z_vals.min()
    results.append((obj, length_A))

if not results:
    print("Aucun objet avec des CA trouvé dans la session.")
    sys.exit(1)

header = "%-35s  %10s  %8s" % ("Objet", "Longueur (Å)", "Longueur (nm)")
separator = "-" * len(header)
print(header)
print(separator)
for obj, la in results:
    print("%-35s  %10.1f  %8.2f" % (obj, la, la / 10))

with open(OUT_PATH, "w") as f:
    f.write(header + "\n" + separator + "\n")
    for obj, la in results:
        f.write("%-35s  %10.1f  %8.2f\n" % (obj, la, la / 10))

print("\nRésultats sauvegardés : %s" % OUT_PATH)
cmd.quit()
