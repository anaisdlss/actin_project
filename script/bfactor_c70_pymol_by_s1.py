"""
Génère un script PyMOL par cluster S1 (site de liaison actine C70).

Chaque script montre :
  - L'actine de référence (8iah_L, chain A semi-transparente)
  - Tous les partenaires S2 dont le C70 patch interagit avec ce cluster S1
    (surface complète chain B, colorée par B-factor)

Usage :
  cd /path/to/actin_project
  python -m script.bfactor_c70_pymol_by_s1

Sortie :
  data/filtered/details/structures_files/bfactor_c70_interface/by_s1_cluster/{s1}.pml
"""

import re
import pandas as pd
from pathlib import Path

C70_CSV   = "data/filtered/patches_infos_cluster_data_70.csv"
ROLES_CSV = "data/filtered/cluster_s1_s2_roles.csv"
BFAC_DIR  = Path("data/filtered/details/structures_files/bfactor_c70_interface")
PAIRWISE_DIR = Path("data/filtered/details/structures_files/pairwise")
OUT_DIR   = BFAC_DIR / "by_s1_cluster"
OUT_DIR.mkdir(parents=True, exist_ok=True)

df_c70   = pd.read_csv(C70_CSV)
df_roles = pd.read_csv(ROLES_CSV)

# Index interaction_type par patch
patch_type = dict(zip(df_c70["patch"].astype(str), df_c70["interaction_type"]))
patch_n    = dict(zip(df_c70["patch"].astype(str), df_c70["n_interactions"]))

# Actine de référence : 8iah_L
ref_pdb_path = PAIRWISE_DIR / "182_8iah_L_0.pdb"
if not ref_pdb_path.exists():
    candidates = sorted(PAIRWISE_DIR.glob("*_8iah_L_*.pdb"))
    if not candidates:
        raise FileNotFoundError("PDB de référence 8iah_L introuvable")
    ref_pdb_path = candidates[0]
ref_pdb_abs = ref_pdb_path.resolve()

# Construire le mapping actin_site → liste de patches C70
# On inclut :
#   - s1_binding_sites pour tous les patches
#   - s2_binding_sites pour les patches homo (S2 = actine aussi)
actin_to_patches: dict[str, list[str]] = {}

def _add(site: str, patch: str):
    site = site.strip()
    if not site:
        return
    actin_to_patches.setdefault(site, [])
    if patch not in actin_to_patches[site]:
        actin_to_patches[site].append(patch)

for _, row in df_roles.iterrows():
    patch = str(row["patch"])
    itype = patch_type.get(patch, "hetero")

    for s1 in str(row["s1_binding_sites"]).split(";"):
        _add(s1, patch)

    # Pour les homos, S2 est aussi actine → même groupement
    if itype == "homo":
        for s2 in str(row["s2_binding_sites"]).split(";"):
            _add(s2, patch)

print(f"Clusters actine uniques : {len(actin_to_patches)}")

n_written = 0
for s1, patches in sorted(actin_to_patches.items()):
    # Trier les patches par n_interactions décroissant
    patches_sorted = sorted(patches, key=lambda p: patch_n.get(p, 0), reverse=True)

    lines = [
        f"# PyMOL — cluster actine : {s1}",
        f"# {len(patches_sorted)} partenaires C70",
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
        "# ── Partenaires ──────────────────────────────────────────────────────────",
    ]

    n_loaded = 0
    for patch in patches_sorted:
        pdb_path = (BFAC_DIR / f"{patch}.pdb").resolve()
        if not pdb_path.exists():
            continue

        itype = patch_type.get(patch, "hetero")
        obj_tmp     = "tmp_" + re.sub(r"[^A-Za-z0-9_]", "_", patch)
        obj_partner = "p_"   + re.sub(r"[^A-Za-z0-9_]", "_", patch)

        # Max B-factor de chain B pour l'autoscale
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

        if itype == "homo":
            color_cmd = f"spectrum b, white_hotpink, {obj_partner}, minimum=0, maximum={bmax}"
        else:
            color_cmd = f"spectrum b, white_green,   {obj_partner}, minimum=0, maximum={bmax}"

        lines += [
            f"# {patch} ({itype}, {patch_n.get(patch, '?')} interactions)",
            f"load {pdb_path}, {obj_tmp}",
            f"align {obj_tmp} and chain A, base_actin",
            f"create {obj_partner}, {obj_tmp} and chain B",
            f"delete {obj_tmp}",
            f"hide everything, {obj_partner}",
            f"show surface, {obj_partner}",
            color_cmd,
            "",
        ]
        n_loaded += 1

    if n_loaded == 0:
        continue

    lines += [
        "# ── Rendu global ─────────────────────────────────────────────────────────",
        "set surface_quality, 1",
        "bg_color white",
        "set ray_opaque_background, off",
        "zoom base_actin",
    ]

    # Nom de fichier : remplacer caractères non-valides
    safe_s1 = re.sub(r"[^A-Za-z0-9_\-]", "_", s1)
    out_pml = OUT_DIR / f"{safe_s1}.pml"
    out_pml.write_text("\n".join(lines))
    n_written += 1

print(f"{n_written} scripts PyMOL générés dans {OUT_DIR}/")
