"""
Génère un script PyMOL par cluster de site de liaison actine (Binding Site Cluster Data 70).

Logique identique au réseau Streamlit "Binding Site Cluster Data 70" :
  - Un script par nœud rouge (S1 binding site cluster, ou S2 actin pour homo)
  - À l'intérieur : un objet PyMOL par S2 binding site connecté (arête dans le réseau)
    → représenté par le C70 patch avec le plus d'interactions pour cette paire (S1, S2)

Chaque script montre :
  - L'actine de référence (8iah_L, surface grise semi-transparente)
  - Un objet par S2 binding site (surface chain B, colorée par B-factor)

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

# ── Groupement : actin_site → { s2_binding_site → meilleur C70 patch } ────────
# Identique au réseau : une arête par (S1 binding site, S2 binding site)
# Si plusieurs C70 patches pour la même paire → on garde celui avec le plus
# d'interactions (représentant le plus statistiquement solide).
#
# Structure : actin_to_s2 [actin_site][s2_site] = best_c70_patch
actin_to_s2: dict[str, dict[str, str]] = defaultdict(dict)

df_valid = df_all[df_all["cluster_data_70"].notna()].copy()
df_valid["s1_bs"] = df_valid["s1_binding_site_cluster_data_70"].astype(str)
df_valid["s2_bs"] = df_valid["s2_binding_site_cluster_data_70"].astype(str)
df_valid["c70"]   = df_valid["cluster_data_70"].astype(str)

for _, row in df_valid.iterrows():
    s1_bs = row["s1_bs"]
    s2_bs = row["s2_bs"]
    c70   = row["c70"]
    if s1_bs == "nan" or s2_bs == "nan" or c70 == "nan":
        continue

    # Mettre à jour si ce C70 patch a plus d'interactions que le précédent
    current = actin_to_s2[s1_bs].get(s2_bs)
    if current is None or patch_n.get(c70, 0) > patch_n.get(current, 0):
        actin_to_s2[s1_bs][s2_bs] = c70

    # Pour les homos, S2 est aussi actine → même logique depuis le côté S2
    if row["s2_actine"]:
        current2 = actin_to_s2[s2_bs].get(s1_bs)
        if current2 is None or patch_n.get(c70, 0) > patch_n.get(current2, 0):
            actin_to_s2[s2_bs][s1_bs] = c70

print(f"Clusters actine uniques : {len(actin_to_s2)}")

# ── Génération des scripts PyMOL ───────────────────────────────────────────────
n_written = 0
for actin_site, s2_to_c70 in sorted(actin_to_s2.items()):
    # Trier les S2 binding sites par n_interactions du représentant (décroissant)
    s2_sorted = sorted(s2_to_c70.items(), key=lambda x: patch_n.get(x[1], 0), reverse=True)

    lines = [
        f"# PyMOL — cluster actine : {actin_site}",
        f"# {len(s2_sorted)} S2 binding sites (un objet par arête du réseau Binding Site Cluster Data 70)",
        "# Vert = ABP (hétéro) | Rose = actine (homo) | Intensité = B-factor (% ASA interface)",
        "# Représentant par S2 binding site = C70 patch avec le plus d'interactions",
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
        "# ── S2 binding sites ─────────────────────────────────────────────────────",
    ]

    n_loaded = 0
    for s2_site, c70_patch in s2_sorted:
        pdb_path = (BFAC_DIR / f"{c70_patch}.pdb").resolve()
        if not pdb_path.exists():
            continue

        itype = patch_type.get(c70_patch, "hetero")
        # Nom d'objet PyMOL basé sur le S2 binding site (pas le C70 patch)
        obj_tmp     = "tmp_" + re.sub(r"[^A-Za-z0-9_]", "_", s2_site)
        obj_partner = "s2_"  + re.sub(r"[^A-Za-z0-9_]", "_", s2_site)

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
            f"# S2 site : {s2_site} ({itype}) — représentant C70 : {c70_patch} ({patch_n.get(c70_patch, '?')} interactions)",
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
