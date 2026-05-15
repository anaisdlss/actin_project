"""
Génère un script PyMOL par cluster de site de liaison actine (Binding Site Cluster Data 70).

Chaque script correspond à un nœud rouge du réseau "Binding Site Cluster Data 70"
et montre :
  - L'actine de référence (8iah_L, surface grise semi-transparente)
  - Tous les partenaires C70 associés à ce site actine (surface chain B colorée par B-factor)

Logique de groupement (identique au réseau Streamlit) :
  - Pour chaque interaction de df_all : s1_binding_site_cluster_data_70 → cluster_data_70
  - Pour les interactions homo (s2_actine=True) : s2_binding_site_cluster_data_70 aussi actine

Usage :
  cd /path/to/actin_project
  python -m script.bfactor_c70_pymol_by_s1

Sortie :
  data/filtered/details/structures_files/bfactor_c70_interface/by_s1_cluster/{actin_site}.pml
"""

import re
import pandas as pd
from pathlib import Path
from collections import defaultdict

C70_CSV  = "data/filtered/patches_infos_cluster_data_70.csv"
ALL_CSV  = "data/filtered/filtered_all_data.csv"
BFAC_DIR = Path("data/filtered/details/structures_files/bfactor_c70_interface")
PAIRWISE_DIR = Path("data/filtered/details/structures_files/pairwise")
OUT_DIR  = BFAC_DIR / "by_s1_cluster"
OUT_DIR.mkdir(parents=True, exist_ok=True)

df_c70 = pd.read_csv(C70_CSV)
df_all = pd.read_csv(ALL_CSV)
df_all["s2_actine"] = df_all["s2_actine"].fillna(False)

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

# ── Groupement actin_site → set de C70 patches (même logique que le réseau) ───
actin_to_patches: dict[str, set] = defaultdict(set)

for _, row in df_all.iterrows():
    s1_bs = row.get("s1_binding_site_cluster_data_70")
    c70   = row.get("cluster_data_70")
    if pd.isna(s1_bs) or pd.isna(c70):
        continue
    s1_bs = str(s1_bs)
    c70   = str(c70)
    actin_to_patches[s1_bs].add(c70)

    # Pour les homos, S2 est aussi actine → même groupement
    if row["s2_actine"]:
        s2_bs = row.get("s2_binding_site_cluster_data_70")
        if pd.notna(s2_bs):
            actin_to_patches[str(s2_bs)].add(c70)

print(f"Clusters actine uniques : {len(actin_to_patches)}")

# ── Génération des scripts PyMOL ───────────────────────────────────────────────
n_written = 0
for actin_site, patches_set in sorted(actin_to_patches.items()):
    patches_sorted = sorted(patches_set, key=lambda p: patch_n.get(p, 0), reverse=True)

    lines = [
        f"# PyMOL — cluster actine : {actin_site}",
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

    safe_name = re.sub(r"[^A-Za-z0-9_\-]", "_", actin_site)
    (OUT_DIR / f"{safe_name}.pml").write_text("\n".join(lines))
    n_written += 1

print(f"{n_written} scripts PyMOL générés dans {OUT_DIR}/")
