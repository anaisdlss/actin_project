# ──────────────────────────────────────────────────────────────────────────────
# Génère un PDB par cluster C70 avec b-factors = moyenne équitable % ASA interface.
# Chain A = S1 (actine), Chain B = S2 (ABP ou actine pour les homos).
# B-factor = 0 pour les résidus hors interface du représentatif.
#
# Usage :
#   cd /path/to/actin_project
#   python script/bfactor_c70_interface.py
#
# Sorties :
#   data/filtered/details/structures_files/bfactor_c70_interface/<patch>.pdb
#   data/filtered/details/structures_files/bfactor_c70_interface/load_all.pml
# ──────────────────────────────────────────────────────────────────────────────

import re
import os
import numpy as np
import pandas as pd
from pathlib import Path

TABLE1   = "data/filtered/details/1.interactions.csv"
TABLE3   = "data/filtered/details/3.interface_residues.csv"
TABLE4   = "data/filtered/details/4.inter-residue_contacts.csv"
TABLE8   = "data/filtered/details/8.structures.csv"
ALL_DATA = "data/filtered/filtered_all_data.csv"
C70_CSV  = "data/filtered/patches_infos_cluster_data_70.csv"
PP_CSV   = "data/filtered/proteins_per_pdb.csv"
OUT_DIR  = Path("data/filtered/details/structures_files/bfactor_c70_interface")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Chargement ─────────────────────────────────────────────────────────────────
df_int  = pd.read_csv(TABLE1)
df4_raw = pd.read_csv(TABLE4)
df4_raw["contact_area"] = pd.to_numeric(df4_raw["contact_area"], errors="coerce").fillna(0)
df4_raw["residue_B_canon_mafft"] = pd.to_numeric(df4_raw["residue_B_canon_mafft"], errors="coerce")
df4_raw["residue_A_canon_mafft"] = pd.to_numeric(df4_raw["residue_A_canon_mafft"], errors="coerce")
df3_raw  = pd.read_csv(TABLE3)
df3_raw["residue_number_canon_mafft"] = pd.to_numeric(
    df3_raw["residue_number_canon_mafft"], errors="coerce")
df8      = pd.read_csv(TABLE8)
df_all   = pd.read_csv(ALL_DATA)
df_c70   = pd.read_csv(C70_CSV)
df_pp    = pd.read_csv(PP_CSV)

def parse_ids(s):
    s = re.sub(r"np\.int64\((\d+)\)", r"\1", str(s))
    return {int(x) for x in re.findall(r"\d+", s)}

patch_to_iids = {str(r["patch"]): parse_ids(r["ids_interactions"])
                 for _, r in df_c70.iterrows()}
all_iids_global = set().union(*patch_to_iids.values())

# ── Détermination des swaps S1↔S2 ─────────────────────────────────────────────
# Ancres chargées depuis data/filtered/binding_site_roles.csv (éditez ce fichier
# pour ajouter de nouveaux patches similaires aux 4 sites canoniques).
# Propagation locale puis fallback vote majoritaire global.
_bsr_path = Path("data/filtered/binding_site_roles.csv")
if _bsr_path.exists():
    _bsr = pd.read_csv(str(_bsr_path))
    _FIXED_S1 = set(_bsr[_bsr["role"] == "S1"]["patch"].astype(str))
    _FIXED_S2 = set(_bsr[_bsr["role"] == "S2"]["patch"].astype(str))
else:
    _FIXED_S1 = {"6685_1", "6685_3"}
    _FIXED_S2 = {"6685_2", "6685_4"}

df_all["s2_actine"] = df_all["s2_actine"].fillna(False)
_homo_df = df_all[df_all["s2_actine"] == True]
_gs1 = _homo_df["s1_binding_site_cluster_data_70"].astype(str).value_counts()
_gs2 = _homo_df["s2_binding_site_cluster_data_70"].astype(str).value_counts()
_is_canonical_s1 = {
    p: int(_gs1.get(p, 0)) >= int(_gs2.get(p, 0))
    for p in set(_gs1.index) | set(_gs2.index)
}

_cmap_raw = (
    df_int.merge(
        df_all[["subunit_1", "subunit_2",
                "s1_binding_site_cluster_data_70",
                "s2_binding_site_cluster_data_70"]],
        left_on=["chain_A_id", "chain_B_id"],
        right_on=["subunit_1", "subunit_2"], how="left"
    )
    .drop_duplicates("interaction_id")
    .set_index("interaction_id")
)
_s1p_map = _cmap_raw["s1_binding_site_cluster_data_70"].astype(str)
_s2p_map = _cmap_raw["s2_binding_site_cluster_data_70"].astype(str)

_patch_role: dict = {p: "S1" for p in _FIXED_S1}
_patch_role.update({p: "S2" for p in _FIXED_S2})

_changed = True
while _changed:
    _changed = False
    for _iid in all_iids_global:
        if _iid not in _cmap_raw.index:
            continue
        _cs1, _cs2 = _s1p_map[_iid], _s2p_map[_iid]
        if _cs1 in _patch_role and _cs2 not in _patch_role:
            _patch_role[_cs2] = "S2" if _patch_role[_cs1] == "S1" else "S1"
            _changed = True
        elif _cs2 in _patch_role and _cs1 not in _patch_role:
            _patch_role[_cs1] = "S1" if _patch_role[_cs2] == "S2" else "S2"
            _changed = True

for _iid in all_iids_global:
    if _iid not in _cmap_raw.index:
        continue
    for _pp in (_s1p_map[_iid], _s2p_map[_iid]):
        if _pp not in _patch_role:
            _patch_role[_pp] = "S1" if _is_canonical_s1.get(_pp, True) else "S2"

# Interactions homo uniquement (actin-actin) — le swap sémantique S1/S2 ne
# s'applique qu'aux homos ; pour les hétéro, chain_A = actin par convention PPI3D.
_actin_ch = set(df_pp[df_pp["is_actin"]]["chain"].str.lower())
_homo_iids = set(
    df_int[df_int["chain_B_id"].str.lower().isin(_actin_ch)]["interaction_id"]
)

swap_iids: set = set()
for _iid in all_iids_global:
    if _iid not in _cmap_raw.index or _iid not in _homo_iids:
        continue
    if (_patch_role.get(_s1p_map[_iid]) == "S2" and
            _patch_role.get(_s2p_map[_iid]) == "S1"):
        swap_iids.add(_iid)

print(f"Interactions swappées (S1/S2 sémantique, homo seulement) : {len(swap_iids)} / {len(_homo_iids & all_iids_global)}")

# ── Correction swap PPI3D table4 (Jaccard overlap, homo seulement) ─────────────
# Table4 inverse parfois residue_A et residue_B par rapport à la convention
# chain_A_id / chain_B_id de table1. On détecte via l'overlap avec table3.
_int_ch4 = df_int.set_index("interaction_id")
_t4_swapped_ppi3d: set = set()
for _iid in _homo_iids & all_iids_global:
    if _iid not in _int_ch4.index:
        continue
    _cA = str(_int_ch4.at[_iid, "chain_A_id"]).lower()
    _cB = str(_int_ch4.at[_iid, "chain_B_id"]).lower()
    _sub3 = df3_raw[df3_raw["interaction_id"] == _iid]
    _t3A = set(_sub3[_sub3["chain"].str.lower() == _cA]["residue_number_canon_mafft"].dropna().astype(int))
    _t3B = set(_sub3[_sub3["chain"].str.lower() == _cB]["residue_number_canon_mafft"].dropna().astype(int))
    if not _t3A and not _t3B:
        continue
    _t4A = set(df4_raw[df4_raw["interaction_id"] == _iid]["residue_A_canon_mafft"].dropna().astype(int))
    if not _t4A:
        continue
    _ovA = len(_t4A & _t3A) / max(len(_t4A | _t3A), 1)
    _ovB = len(_t4A & _t3B) / max(len(_t4A | _t3B), 1)
    if _ovB > _ovA:
        _t4_swapped_ppi3d.add(_iid)

print(f"Interactions avec swap PPI3D table4 (Jaccard) : {len(_t4_swapped_ppi3d)}")

# ── Appliquer les swaps dans df4 ───────────────────────────────────────────────
# Ordre : 1) swap PPI3D table4 → residue_A = chain_A
#          2) swap sémantique S1/S2 → residue_A = S1
df4 = df4_raw[df4_raw["interaction_id"].isin(all_iids_global)].copy()

if _t4_swapped_ppi3d:
    _sw_ppi3d = df4["interaction_id"].isin(_t4_swapped_ppi3d)
    df4.loc[_sw_ppi3d, ["residue_A_canon_mafft", "residue_B_canon_mafft"]] = \
        df4.loc[_sw_ppi3d, ["residue_B_canon_mafft", "residue_A_canon_mafft"]].values
    df4.loc[_sw_ppi3d, ["asa_pct_A", "asa_pct_B"]] = \
        df4.loc[_sw_ppi3d, ["asa_pct_B", "asa_pct_A"]].values
    df4.loc[_sw_ppi3d, ["residue_A_name", "residue_B_name"]] = \
        df4.loc[_sw_ppi3d, ["residue_B_name", "residue_A_name"]].values
    df4.loc[_sw_ppi3d, ["residue_A_structure", "residue_B_structure"]] = \
        df4.loc[_sw_ppi3d, ["residue_B_structure", "residue_A_structure"]].values

_sw = df4["interaction_id"].isin(swap_iids)
df4.loc[_sw, ["asa_pct_A", "asa_pct_B"]] = \
    df4.loc[_sw, ["asa_pct_B", "asa_pct_A"]].values
df4.loc[_sw, ["residue_A_canon_mafft", "residue_B_canon_mafft"]] = \
    df4.loc[_sw, ["residue_B_canon_mafft", "residue_A_canon_mafft"]].values
df4.loc[_sw, ["residue_A_name", "residue_B_name"]] = \
    df4.loc[_sw, ["residue_B_name", "residue_A_name"]].values

# Attacher le patch C70 à chaque contact
iid_to_patch = {iid: p for p, iids in patch_to_iids.items() for iid in iids}
df4["patch"] = df4["interaction_id"].map(iid_to_patch)

# ── Moyenne équitable ASA par position canonique ────────────────────────────────
def equitable_mean_asa(df_patch, side: str) -> dict:
    """Moyenne équitable par interaction (max par résidu si contacts multiples)."""
    if side == "S1":
        canon_col, asa_col = "residue_A_canon_mafft", "asa_pct_A"
    else:
        canon_col, asa_col = "residue_B_canon_mafft", "asa_pct_B"
    sub = df_patch[[canon_col, asa_col, "interaction_id"]].dropna(
        subset=[canon_col, asa_col])
    if sub.empty:
        return {}
    n_iids = float(sub["interaction_id"].nunique())
    profile = (
        sub.rename(columns={canon_col: "canon", asa_col: "asa"})
        .groupby(["interaction_id", "canon"])["asa"].max()
        .groupby(level="canon").sum()
        / n_iids
    )
    return profile.to_dict()

# ── Index structures (pairwise PDB par interaction) ─────────────────────────────
df8_valid = df8[
    df8["pairwise_pdb_file"].notna() &
    (df8["pairwise_pdb_file"].astype(str).str.lower() != "nan")
].copy()

ok, skip = 0, 0
pymol_lines: list[str] = ["# PyMOL : clusters C70 interface — b-factors = % ASA interface équitable"]
ref_obj: str | None = None

for patch, iids in patch_to_iids.items():
    df_patch = df4[df4["patch"] == patch]
    if df_patch.empty:
        skip += 1
        continue

    asa_s1 = equitable_mean_asa(df_patch, "S1")
    asa_s2 = equitable_mean_asa(df_patch, "S2")
    if not asa_s1 and not asa_s2:
        skip += 1
        continue

    # Représentatif : premier interaction_id avec un pairwise PDB valide
    rep_iid, pdb_file, is_sw_rep = None, None, False
    for iid in sorted(iids):
        rows = df8_valid[df8_valid["interaction_id"] == iid]
        for _, r in rows.iterrows():
            pp = str(r["pairwise_pdb_file"])
            if os.path.exists(pp):
                rep_iid, pdb_file = iid, pp
                # Lookup direct binding site (même logique que le réseau)
                _int_row = df_int[df_int["interaction_id"] == iid]
                if not _int_row.empty:
                    _ch_A = _int_row.iloc[0]["chain_A_id"]
                    _ch_B = _int_row.iloc[0]["chain_B_id"]
                    _row = df_all[
                        (df_all["subunit_1"] == _ch_A) &
                        (df_all["subunit_2"] == _ch_B)
                    ]
                    if not _row.empty:
                        _s1 = str(_row.iloc[0]["s1_binding_site_cluster_data_70"])
                        _s2 = str(_row.iloc[0]["s2_binding_site_cluster_data_70"])
                        if _patch_role.get(_s1) == "S1":
                            is_sw_rep = False   # subunit_1 = S1 → chain A physique = S1
                        elif _patch_role.get(_s2) == "S1":
                            is_sw_rep = True    # subunit_2 = S1 → chain B physique = S1 → swap
                        else:
                            is_sw_rep = iid in swap_iids  # fallback propagation
                    else:
                        is_sw_rep = iid in swap_iids
                else:
                    is_sw_rep = iid in swap_iids
                break
        if pdb_file:
            break

    if pdb_file is None:
        skip += 1
        continue

    # Mapping canon → résidu structure dans le PDB représentatif (données BRUTES)
    df_rep = df4_raw[df4_raw["interaction_id"] == rep_iid]
    if is_sw_rep:
        # Chaîne physique B = S1
        s1_c2s = (df_rep[["residue_B_canon_mafft", "residue_B_structure"]]
                  .dropna().drop_duplicates("residue_B_canon_mafft")
                  .set_index("residue_B_canon_mafft")["residue_B_structure"]
                  .astype(int).to_dict())
        s2_c2s = (df_rep[["residue_A_canon_mafft", "residue_A_structure"]]
                  .dropna().drop_duplicates("residue_A_canon_mafft")
                  .set_index("residue_A_canon_mafft")["residue_A_structure"]
                  .astype(int).to_dict())
        ch_s1, ch_s2 = "B", "A"
    else:
        # Chaîne physique A = S1
        s1_c2s = (df_rep[["residue_A_canon_mafft", "residue_A_structure"]]
                  .dropna().drop_duplicates("residue_A_canon_mafft")
                  .set_index("residue_A_canon_mafft")["residue_A_structure"]
                  .astype(int).to_dict())
        s2_c2s = (df_rep[["residue_B_canon_mafft", "residue_B_structure"]]
                  .dropna().drop_duplicates("residue_B_canon_mafft")
                  .set_index("residue_B_canon_mafft")["residue_B_structure"]
                  .astype(int).to_dict())
        ch_s1, ch_s2 = "A", "B"

    bfac_s1 = {s1_c2s[c]: v for c, v in asa_s1.items() if c in s1_c2s}
    bfac_s2 = {s2_c2s[c]: v for c, v in asa_s2.items() if c in s2_c2s}

    # Réécrire le PDB : chaîne A = S1, chaîne B = S2, b-factor = ASA équitable
    out_lines: list[str] = []
    for line in open(pdb_file):
        if line.startswith(("ATOM", "HETATM")) and len(line) > 26:
            ch_phys = line[21]
            try:
                resnum = int(line[22:26].strip())
            except ValueError:
                out_lines.append(line)
                continue
            if ch_phys == ch_s1:
                bfac = bfac_s1.get(resnum, 0.0)
                line = line[:21] + "A" + line[22:]   # renomme → A (S1)
            elif ch_phys == ch_s2:
                bfac = bfac_s2.get(resnum, 0.0)
                line = line[:21] + "B" + line[22:]   # renomme → B (S2)
            else:
                out_lines.append(line)
                continue
            line = f"{line[:60]}{bfac:6.2f}{line[66:]}"
        out_lines.append(line)

    out_path = OUT_DIR / f"{patch}.pdb"
    with open(out_path, "w") as f:
        f.write("".join(out_lines))
    ok += 1

    # Ligne PyMOL pour ce cluster
    obj = patch.replace("-", "_")
    if ref_obj is None:
        ref_obj = obj
    pymol_lines.append(f"\nload {out_path}, {obj}")
    if ref_obj != obj:
        pymol_lines.append(f"align {obj} and chain A, {ref_obj} and chain A")
    pymol_lines.append(f"show surface, {obj}")
    pymol_lines.append(f"spectrum b, white_red,   {obj} and chain A, minimum=0, maximum=100")
    pymol_lines.append(f"spectrum b, white_blue,  {obj} and chain B, minimum=0, maximum=100")

# ── Script PyMOL global ────────────────────────────────────────────────────────
pymol_lines += [
    "",
    "set surface_quality, 1",
    "bg_color white",
    "set ray_opaque_background, off",
]
pymol_path = OUT_DIR / "load_all.pml"
with open(pymol_path, "w") as f:
    f.write("\n".join(pymol_lines))

print(f"Terminé — {ok} clusters sauvegardés, {skip} ignorés")
print(f"Répertoire : {OUT_DIR}/")
print(f"Script PyMOL : {pymol_path}")

# ── Export CSV rôles S1/S2 par cluster ────────────────────────────────────────
_rows_roles = []
for _patch, _iids in patch_to_iids.items():
    _s1_sites, _s2_sites = set(), set()
    for _iid in _iids:
        if _iid not in _cmap_raw.index:
            continue
        _cs1 = _s1p_map[_iid]
        _cs2 = _s2p_map[_iid]
        if _patch_role.get(_cs1) == "S1":
            _s1_sites.add(_cs1)
            _s2_sites.add(_cs2)
        elif _patch_role.get(_cs2) == "S1":
            _s1_sites.add(_cs2)
            _s2_sites.add(_cs1)
    _rows_roles.append({
        "patch": _patch,
        "s1_binding_sites": ";".join(sorted(_s1_sites)),
        "s2_binding_sites": ";".join(sorted(_s2_sites)),
    })

_roles_csv = Path("data/filtered/cluster_s1_s2_roles.csv")
pd.DataFrame(_rows_roles).to_csv(_roles_csv, index=False)
print(f"CSV rôles S1/S2 : {_roles_csv}")
