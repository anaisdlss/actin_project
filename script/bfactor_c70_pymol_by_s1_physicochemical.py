"""
Génère un script PyMOL par cluster S1 — couleurs physicochimiques + intensité B-factor.

Chaque résidu d'interface (b > 10) est coloré selon :
  - la teinte  → type physicochimique (bleu=+, rouge=-, gris=hydrophobe…)
  - l'intensité → valeur du B-factor (blanc=faible, couleur saturée=élevé)

Usage :
  python -m script.bfactor_c70_pymol_by_s1_physicochemical

Sortie :
  data/filtered/details/structures_files/bfactor_c70_interface/by_s1_physicochemical/{actin_site}.pml
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
OUT_DIR      = BFAC_DIR / "by_s1_physicochemical"
OUT_DIR.mkdir(parents=True, exist_ok=True)

df_c70 = pd.read_csv(C70_CSV)
df_all = pd.read_csv(ALL_CSV)
df_all["s2_actine"] = df_all["s2_actine"].fillna(False)

patch_type = dict(zip(df_c70["patch"].astype(str), df_c70["interaction_type"]))
patch_n    = dict(zip(df_c70["patch"].astype(str), df_c70["n_interactions"]))

ref_pdb_path = PAIRWISE_DIR / "182_8iah_L_0.pdb"
if not ref_pdb_path.exists():
    candidates = sorted(PAIRWISE_DIR.glob("*_8iah_L_*.pdb"))
    if not candidates:
        raise FileNotFoundError("PDB de référence 8iah_L introuvable")
    ref_pdb_path = candidates[0]
ref_pdb_abs = ref_pdb_path.resolve()


def _clean_name(title: str) -> str:
    if not title or str(title).lower() in ("nan", ""):
        return ""
    part = str(title).split(",")[0].strip()
    part = re.sub(r"[^A-Za-z0-9]", "_", part)
    part = re.sub(r"_+", "_", part).strip("_").lower()
    return part[:35]


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


_BTHRESH = 10

def _restype_surf_cmds(obj: str, bthresh: float = _BTHRESH,
                       base_color: str = "orange", bmax: float = 0.0) -> list[str]:
    """Physicochimique + intensité B-factor.

    Résidus non-interface (b <= bthresh) : base_color.
    Résidus interface (b > bthresh)      : teinte physicochimique, intensité ∝ B-factor.
    """
    sel = f"{obj} and b > {bthresh}"
    cmds = [f"color {base_color}, {obj}"]
    if bmax > 1.0:
        m = f"minimum=0, maximum={bmax}"
        cmds += [
            f"spectrum b, white_grey60,   {sel} and (resn ALA+GLY+ILE+LEU+MET+VAL), {m}",
            f"spectrum b, white_hotpink,  {sel} and (resn PHE+TRP+TYR), {m}",
            f"spectrum b, white_cyan,     {sel} and (resn HIS+ASN+GLN+SER+THR), {m}",
            f"spectrum b, white_blue,     {sel} and (resn LYS+ARG), {m}",
            f"spectrum b, white_red,      {sel} and (resn ASP+GLU), {m}",
            f"spectrum b, white_yellow,   {sel} and resn CYS, {m}",
            f"spectrum b, white_limegreen,{sel} and resn PRO, {m}",
        ]
    else:
        cmds += [
            f"color grey60,    {sel} and (resn ALA+GLY+ILE+LEU+MET+VAL)",
            f"color hotpink,   {sel} and (resn PHE+TRP+TYR)",
            f"color cyan,      {sel} and (resn HIS+ASN+GLN+SER+THR)",
            f"color blue,      {sel} and (resn LYS+ARG)",
            f"color red,       {sel} and (resn ASP+GLU)",
            f"color yellow,    {sel} and resn CYS",
            f"color limegreen, {sel} and resn PRO",
        ]
    return cmds


df_valid = df_all[df_all["cluster_data_70"].notna()].copy()
df_valid["s1_bs"] = df_valid["s1_binding_site_cluster_data_70"].astype(str)
df_valid["s2_bs"] = df_valid["s2_binding_site_cluster_data_70"].astype(str)
df_valid["c70"]   = df_valid["cluster_data_70"].astype(str)

_s2bs_name: dict[str, str] = {}
for s2_bs, grp in df_valid.groupby("s2_bs"):
    titles = grp["subunit_2_title"].dropna()
    if len(titles):
        _s2bs_name[s2_bs] = _clean_name(titles.value_counts().index[0])

actin_to_s2: dict[str, dict[str, tuple[str, bool]]] = defaultdict(dict)
for _, row in df_valid.iterrows():
    s1_bs, s2_bs, c70 = row["s1_bs"], row["s2_bs"], row["c70"]
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

n_written = 0
for actin_site, s2_to_c70 in sorted(actin_to_s2.items()):
    s2_sorted = sorted(s2_to_c70.items(), key=lambda x: patch_n.get(x[1][0], 0), reverse=True)

    s1_bfac_pdb = (BFAC_S1_DIR / f"{actin_site}.pdb").resolve()
    if s1_bfac_pdb.exists():
        bmax_s1 = _bmax_chain(s1_bfac_pdb, "A")
        ref_section = [
            "# ── Actine S1 — physicochimique + intensité B-factor, semi-transparente ─",
            f"load {s1_bfac_pdb}, base_actin",
            "hide everything, base_actin",
            "show surface, base_actin",
        ] + _restype_surf_cmds("base_actin", base_color="orange", bmax=bmax_s1)
    else:
        ref_section = [
            "# ── Actine S1 — orange (pas de B-factor) ────────────────────────────",
            f"load {ref_pdb_abs}, _ref_complex",
            "create base_actin, _ref_complex and chain A",
            "delete _ref_complex",
            "hide everything, base_actin",
            "show surface, base_actin",
            "color orange, base_actin",
        ]

    lines = [
        f"# PyMOL — physicochimique + B-factor — cluster actine S1 : {actin_site}",
        f"# {len(s2_sorted)} partenaires",
        "# Résidus interface (b > 10) : teinte=type physicochimique, intensité=B-factor",
        "# gris=hydrophobe · rose=aromatique · cyan=polaire · bleu=+ · rouge=- · jaune=Cys",
        "# Base : orange=actine · vert pâle=ABP",
        "",
    ] + ref_section + ["", "# ── Partenaires S2 ─────────────────────────────────────────────────────"]

    n_loaded = 0
    for s2_site, (c70_patch, is_swapped) in s2_sorted:
        pdb_path = (BFAC_DIR / f"{c70_patch}.pdb").resolve()
        if not pdb_path.exists():
            continue

        itype       = patch_type.get(c70_patch, "hetero")
        prot_name   = _s2bs_name.get(s2_site, "")
        safe_s2     = re.sub(r"[^A-Za-z0-9_]", "_", s2_site)
        obj_tmp     = "tmp_" + safe_s2
        obj_partner = safe_s2 + ("_" + prot_name if prot_name else "")
        align_ch    = "B" if is_swapped else "A"
        partner_ch  = "A" if is_swapped else "B"

        bmax        = _bmax_chain(pdb_path, partner_ch)
        base_col    = "orange" if itype == "homo" else "palegreen"

        swap_note = " [vue swappée]" if is_swapped else ""
        label = (f"# {s2_site} ({itype}{swap_note}) — C70 : {c70_patch} "
                 f"({patch_n.get(c70_patch, '?')} interactions)")

        lines += [
            label,
            f"load {pdb_path}, {obj_tmp}",
            f"align {obj_tmp} and chain {align_ch}, base_actin",
            f"create {obj_partner}, {obj_tmp} and chain {partner_ch}",
            f"delete {obj_tmp}",
            f"hide everything, {obj_partner}",
            f"show surface, {obj_partner}",
        ] + _restype_surf_cmds(obj_partner, base_color=base_col, bmax=bmax) + [""]
        n_loaded += 1

    if n_loaded == 0:
        continue

    lines += [
        "# ── Rendu global ───────────────────────────────────────────────────────",
        "set surface_quality, 1",
        "bg_color white",
        "set ray_opaque_background, off",
        "zoom base_actin",
    ]

    safe_name = re.sub(r"[^A-Za-z0-9_\-]", "_", actin_site)
    (OUT_DIR / f"{safe_name}.pml").write_text("\n".join(lines))
    n_written += 1

print(f"{n_written} scripts PyMOL générés dans {OUT_DIR}/")
