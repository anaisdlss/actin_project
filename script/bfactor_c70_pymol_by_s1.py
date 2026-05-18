"""
Génère un script PyMOL par cluster de site de liaison actine (Binding Site Cluster Data 70).

Chaque script montre :
  - L'actine S1 de référence colorée par B-factor du cluster (cyan = résidus d'interface)
  - Un objet par S2 binding site connecté, nommé {s2_cluster}_{nom_proteine}
    Vert = ABP (hétéro) | Rose = actine (homo) | Intensité = % ASA enfouie

Usage :
  python -m script.bfactor_c70_pymol_by_s1

Sortie :
  data/filtered/details/structures_files/bfactor_c70_interface/by_s1_cluster/{actin_site}.pml
"""

import re
import pandas as pd
from pathlib import Path
from collections import defaultdict

C70_CSV      = "data/filtered/patches_infos_cluster_data_70.csv"
ALL_CSV      = "data/filtered/filtered_all_data.csv"
BFAC_DIR     = Path("data/filtered/details/structures_files/bfactor_c70_interface")
BFAC_S1_DIR  = Path("data/filtered/details/structures_files/bfactor_cluster")
PAIRWISE_DIR = Path("data/filtered/details/structures_files/pairwise")
OUT_DIR      = BFAC_DIR / "by_s1_cluster"
OUT_DIR.mkdir(parents=True, exist_ok=True)

df_c70 = pd.read_csv(C70_CSV)
df_all = pd.read_csv(ALL_CSV)
df_all["s2_actine"] = df_all["s2_actine"].fillna(False)

patch_type = dict(zip(df_c70["patch"].astype(str), df_c70["interaction_type"]))
patch_n    = dict(zip(df_c70["patch"].astype(str), df_c70["n_interactions"]))

# ── Actine de référence pour l'alignement (8iah_L) ────────────────────────────
ref_pdb_path = PAIRWISE_DIR / "182_8iah_L_0.pdb"
if not ref_pdb_path.exists():
    candidates = sorted(PAIRWISE_DIR.glob("*_8iah_L_*.pdb"))
    if not candidates:
        raise FileNotFoundError("PDB de référence 8iah_L introuvable")
    ref_pdb_path = candidates[0]
ref_pdb_abs = ref_pdb_path.resolve()

# ── Mapping s2_binding_site → nom de protéine nettoyé ─────────────────────────
def _clean_name(title: str) -> str:
    if not title or str(title).lower() in ("nan", ""):
        return ""
    part = str(title).split(",")[0].strip()
    part = re.sub(r"[^A-Za-z0-9]", "_", part)
    part = re.sub(r"_+", "_", part).strip("_").lower()
    return part[:35]

df_valid = df_all[df_all["cluster_data_70"].notna()].copy()
df_valid["s1_bs"] = df_valid["s1_binding_site_cluster_data_70"].astype(str)
df_valid["s2_bs"] = df_valid["s2_binding_site_cluster_data_70"].astype(str)
df_valid["c70"]   = df_valid["cluster_data_70"].astype(str)

# Titre le plus fréquent par s2_binding_site
_s2bs_name: dict[str, str] = {}
for s2_bs, grp in df_valid.groupby("s2_bs"):
    titles = grp["subunit_2_title"].dropna()
    if len(titles):
        _s2bs_name[s2_bs] = _clean_name(titles.value_counts().index[0])

# ── Coloration physicochimique sur résidus d'interface (b > seuil) ───────────
_BTHRESH = 10  # % ASA buried minimum

def _restype_surf_cmds(obj: str, bthresh: float = _BTHRESH) -> list[str]:
    """Surface blanche + couleur physicochimique sur les résidus d'interface."""
    sel = f"{obj} and b > {bthresh}"
    return [
        f"color white, {obj}",
        f"color grey80,    {sel} and (resn ALA+GLY+ILE+LEU+MET+VAL)",
        f"color lightpink, {sel} and (resn PHE+TRP+TYR)",
        f"color palecyan,  {sel} and (resn HIS+ASN+GLN+SER+THR)",
        f"color blue,      {sel} and (resn LYS+ARG)",
        f"color red,       {sel} and (resn ASP+GLU)",
        f"color paleyellow,{sel} and resn CYS",
        f"color palegreen, {sel} and resn PRO",
    ]

# ── Lecture du bmax d'un PDB sur une chaîne donnée ────────────────────────────
def _bmax_chain(pdb_path: Path, ch: str) -> float:
    bmax = 0.0
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and len(line) > 66 and line[21] == ch:
                try:
                    v = float(line[60:66].strip())
                    if v > bmax:
                        bmax = v
                except ValueError:
                    pass
    return round(bmax, 2)

# ── Groupement actin_site → { s2_site → (c70_patch, is_swapped) } ─────────────
actin_to_s2: dict[str, dict[str, tuple[str, bool]]] = defaultdict(dict)

for _, row in df_valid.iterrows():
    s1_bs = row["s1_bs"]
    s2_bs = row["s2_bs"]
    c70   = row["c70"]
    if s1_bs == "nan" or s2_bs == "nan" or c70 == "nan":
        continue

    current = actin_to_s2[s1_bs].get(s2_bs)
    if current is None or patch_n.get(c70, 0) > patch_n.get(current[0], 0):
        actin_to_s2[s1_bs][s2_bs] = (c70, False)

    if row["s2_actine"]:
        current2 = actin_to_s2[s2_bs].get(s1_bs)
        if current2 is None or patch_n.get(c70, 0) > patch_n.get(current2[0], 0):
            actin_to_s2[s2_bs][s1_bs] = (c70, True)

print(f"Clusters actine uniques : {len(actin_to_s2)}")

# ── Génération des scripts PyMOL ───────────────────────────────────────────────
n_written = 0
for actin_site, s2_to_c70 in sorted(actin_to_s2.items()):
    s2_sorted = sorted(s2_to_c70.items(), key=lambda x: patch_n.get(x[1][0], 0), reverse=True)

    # ── Actine S1 : B-factor du cluster si disponible, sinon gris ─────────────
    # Le PDB bfactor_cluster contient déjà chain A avec les B-factors calculés.
    # On le charge directement comme base_actin — les S2 s'alignent dessus.
    s1_bfac_pdb = (BFAC_S1_DIR / f"{actin_site}.pdb").resolve()
    if s1_bfac_pdb.exists():
        ref_section = [
            "# ── Actine S1 — couleurs physicochimiques sur résidus d'interface ───────",
            f"load {s1_bfac_pdb}, base_actin",
            "hide everything, base_actin",
            "show surface, base_actin",
        ] + _restype_surf_cmds("base_actin")
    else:
        ref_section = [
            "# ── Actine S1 — référence grise (pas de B-factor pour ce cluster) ────────",
            f"load {ref_pdb_abs}, _ref_complex",
            "create base_actin, _ref_complex and chain A",
            "delete _ref_complex",
            "hide everything, base_actin",
            "show surface, base_actin",
            "color grey70, base_actin",
        ]

    lines = [
        f"# PyMOL — cluster actine S1 : {actin_site}",
        f"# {len(s2_sorted)} partenaires (un objet par arête du réseau Binding Site Cluster Data 70)",
        "# Interface (b > 10%) : gris=hydrophobe · rose=aromatique · cyan=polaire",
        "#                       bleu=positif · rouge=négatif · jaune=Cys · vert=Pro",
        "# Vert (partenaire) = ABP hétéro | Rose (partenaire) = actine homo",
        "",
    ] + ref_section + ["", "# ── Partenaires S2 ───────────────────────────────────────────────────────"]

    n_loaded = 0
    for s2_site, (c70_patch, is_swapped) in s2_sorted:
        pdb_path = (BFAC_DIR / f"{c70_patch}.pdb").resolve()
        if not pdb_path.exists():
            continue

        itype = patch_type.get(c70_patch, "hetero")

        # Nom d'objet : {s2_site}_{nom_proteine}
        prot_name = _s2bs_name.get(s2_site, "")
        safe_s2   = re.sub(r"[^A-Za-z0-9_]", "_", s2_site)
        obj_base  = safe_s2 + ("_" + prot_name if prot_name else "")
        obj_tmp     = "tmp_" + safe_s2
        obj_partner = obj_base

        align_ch   = "B" if is_swapped else "A"
        partner_ch = "A" if is_swapped else "B"

        bmax = _bmax_chain(pdb_path, partner_ch)

        swap_note = " [vue swappée]" if is_swapped else ""
        label = f"# {s2_site} ({itype}{swap_note}) — C70 : {c70_patch} ({patch_n.get(c70_patch, '?')} interactions)"

        if itype == "homo":
            color_cmd = (f"spectrum b, white_hotpink, {obj_partner}, minimum=0, maximum={bmax}"
                         if bmax > 1.0 else f"color hotpink, {obj_partner}")
        else:
            color_cmd = f"spectrum b, white_green, {obj_partner}, minimum=0, maximum={bmax}"

        lines += [
            label,
            f"load {pdb_path}, {obj_tmp}",
            f"align {obj_tmp} and chain {align_ch}, base_actin",
            f"create {obj_partner}, {obj_tmp} and chain {partner_ch}",
            f"delete {obj_tmp}",
            f"hide everything, {obj_partner}",
            f"show surface, {obj_partner}",
        ] + _restype_surf_cmds(obj_partner) + [""]
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
