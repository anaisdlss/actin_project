"""
Génère un script PyMOL pour visualiser les partenaires C70 avec surface complète.

Logique :
  - Actine de référence = chaîne A du PDB de 8iah_L (l'actine la plus fréquente dans all data)
  - Pour chaque cluster C70 :
      * Chargement du PDB bfactor (chain A = actin S1, chain B = partenaire S2)
      * Alignement chain A → actine de référence
      * Extraction chain B seule dans un objet séparé → surface entière non masquée
      * Coloration B-factor : vert (ABP, hétéro) ou rose (actine, homo)

Usage :
  cd /path/to/actin_project
  python script/bfactor_c70_pymol_full_surface.py

Sortie :
  data/filtered/details/structures_files/bfactor_c70_interface/view_full_surface.pml
"""

import re
import pandas as pd
from pathlib import Path

ALL_DATA = "data/filtered/filtered_all_data.csv"
C70_CSV = "data/filtered/patches_infos_cluster_data_70.csv"
BFAC_DIR = Path("data/filtered/details/structures_files/bfactor_c70_interface")
PAIRWISE_DIR = Path("data/filtered/details/structures_files/pairwise")
OUT_PML = BFAC_DIR / "view_full_surface.pml"

df_all = pd.read_csv(ALL_DATA)
df_c70 = pd.read_csv(C70_CSV)

# ── Actine de référence : 8iah_L est la plus fréquente en s1 ──────────────────
#    On extrait chain A de son PDB pairwise (interaction 182)
ref_pdb_path = PAIRWISE_DIR / "182_8iah_L_0.pdb"
if not ref_pdb_path.exists():
    # Fallback : chercher n'importe quel PDB contenant 8iah_L
    candidates = sorted(PAIRWISE_DIR.glob("*_8iah_L_*.pdb"))
    if not candidates:
        raise FileNotFoundError("PDB de référence 8iah_L introuvable")
    ref_pdb_path = candidates[0]

ref_pdb_abs = ref_pdb_path.resolve()
print(f"Actine de référence : {ref_pdb_path}")

# ── Homo vs hétéro par patch ──────────────────────────────────────────────────
patch_type = dict(zip(df_c70["patch"].astype(str), df_c70["interaction_type"]))

# ── Construction du script PyMOL ──────────────────────────────────────────────
lines = [
    "# PyMOL — clusters C70, surface complète des partenaires",
    "# Chain B de chaque cluster extraite seule → surface non masquée",
    "# Vert = ABP (hétéro) | Rose = actine (homo) | Intensité = B-factor (% ASA interface)",
    "",
    "# ── Actine de référence ──────────────────────────────────────────────────",
    f"load {ref_pdb_abs}, _ref_complex",
    "create base_actin, _ref_complex and chain A",
    "delete _ref_complex",
    "hide everything, base_actin",
    "show surface, base_actin",
    "color grey70, base_actin",
    "set transparency, 0.55, base_actin",
    "",
    "# ── Partenaires par cluster ──────────────────────────────────────────────",
]

n_loaded = 0
for _, row in df_c70.sort_values("n_interactions", ascending=False).iterrows():
    patch = str(row["patch"])
    itype = str(row["interaction_type"])
    pdb_path = (BFAC_DIR / f"{patch}.pdb").resolve()
    if not pdb_path.exists():
        print(f"  SKIP {patch} (pas de PDB)")
        continue

    # Nom d'objet PyMOL valide (pas de tirets, commence par lettre)
    obj_tmp = "tmp_" + re.sub(r"[^A-Za-z0-9_]", "_", patch)
    obj_partner = "p_" + re.sub(r"[^A-Za-z0-9_]", "_", patch)

    # Max B-factor réel de chain B (même logique que le viewer Streamlit)
    bmax = 1.0
    for line in open(pdb_path):
        if line.startswith("ATOM") and len(line) > 66 and line[21] == "B":
            try:
                v = float(line[60:66].strip())
                if v > bmax:
                    bmax = v
            except ValueError:
                pass
    bmax = round(bmax, 2)

    # Couleur et surface selon homo/hétéro
    # Homo : on charge directement sous obj_partner (pas de create) pour conserver
    # la transformation d'alignement ; on ne montre que chain B et b > 0 (résidus
    # d'interface). Sans create, PyMOL utilise les coordonnées post-alignement.
    # Hétéro : create isole chain B (ABP) puis surface complète.
    if itype == "homo":
        lines += [
            f"# {patch} ({itype}, {row['n_interactions']} interactions)",
            f"load {pdb_path}, {obj_partner}",
            f"align {obj_partner} and chain A, base_actin",
            f"hide everything, {obj_partner}",
            f"show surface, {obj_partner} and chain B and b > 0.001",
            f"spectrum b, white_hotpink, {obj_partner} and chain B, minimum=0, maximum={bmax}",
            "",
        ]
    else:
        lines += [
            f"# {patch} ({itype}, {row['n_interactions']} interactions)",
            f"load {pdb_path}, {obj_tmp}",
            f"align {obj_tmp} and chain A, base_actin",
            f"create {obj_partner}, {obj_tmp} and chain B",
            f"delete {obj_tmp}",
            f"hide everything, {obj_partner}",
            f"show surface, {obj_partner}",
            f"spectrum b, white_green,   {obj_partner}, minimum=0, maximum={bmax}",
            "",
        ]
    n_loaded += 1

lines += [
    "# ── Rendu global ─────────────────────────────────────────────────────────",
    "set surface_quality, 1",
    "bg_color white",
    "set ray_opaque_background, off",
    "zoom base_actin",
]

OUT_PML.write_text("\n".join(lines))
print(f"\n{n_loaded} clusters chargés")
print(f"Script PyMOL : {OUT_PML}")
