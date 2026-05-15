# ------------------------------------------------------------
# B-factor = moyenne équitable % buried ASA par position canonique MAFFT
# Même méthode que le heatmap equitable : moyenne-de-moyennes par sous-cluster C70
# S1 et S2 actin confondus (même patch = même cluster)
# Un PDB par cluster (chaîne actine de référence uniquement)
# Sortie : data/filtered/details/structures_files/bfactor_cluster/
# ------------------------------------------------------------

import io
import re
import numpy as np
import pandas as pd
import os
from pathlib import Path
from Bio import PDB as _PDB

ALL_DATA = "data/filtered/filtered_all_data.csv"
TABLE1 = "data/filtered/details/1.interactions.csv"
TABLE3 = "data/filtered/details/3.interface_residues.csv"
TABLE8 = "data/filtered/details/8.structures.csv"
REF_CSV = "data/filtered/s1_cluster_reference.csv"
C70_CSV = "data/filtered/patches_infos_cluster_data_70.csv"
OUT_DIR = Path("data/filtered/details/structures_files/bfactor_cluster")
OUT_DIR.mkdir(parents=True, exist_ok=True)

df_all = pd.read_csv(ALL_DATA)
df_int = pd.read_csv(TABLE1)
df3 = pd.read_csv(TABLE3)
df8 = pd.read_csv(TABLE8)
df_ref = pd.read_csv(REF_CSV)
df_c70 = pd.read_csv(C70_CSV)

df3["buried_ASA_percent"] = pd.to_numeric(
    df3["buried_ASA_percent"].str.replace("%", "", regex=False),
    errors="coerce").fillna(0)

# Mapping interaction_id → C70 patch


def parse_ids(s):
    s = re.sub(r'np\.int64\((\d+)\)', r'\1', str(s))
    return set(int(x) for x in re.findall(r'\d+', s))


id_to_c70 = {}
for _, row in df_c70.iterrows():
    for iid in parse_ids(row["ids_interactions"]):
        id_to_c70[iid] = row["patch"]

# Mapping interaction_id → patch S1 et S2
merge_cols = df_int.merge(
    df_all[["subunit_1", "subunit_2",
            "s1_binding_site_cluster_data_70",
            "s2_binding_site_cluster_data_70"]],
    left_on=["chain_A_id", "chain_B_id"],
    right_on=["subunit_1", "subunit_2"],
    how="inner"
)[["interaction_id",
   "s1_binding_site_cluster_data_70",
   "s2_binding_site_cluster_data_70"]].drop_duplicates("interaction_id").set_index("interaction_id")

iid_to_s1_patch = merge_cols["s1_binding_site_cluster_data_70"]
iid_to_s2_patch = merge_cols["s2_binding_site_cluster_data_70"]

# --- Résidus S1 (chain_A) ---
id_to_s1chain = df_int.set_index("interaction_id")["chain_A_id"].str.lower()
df3_tmp = df3.copy()
df3_tmp["s1_chain"] = df3_tmp["interaction_id"].map(id_to_s1chain)
df3_s1 = df3_tmp[
    df3_tmp["residue_number_canon_mafft"].notna() &
    (df3_tmp["chain"].str.lower() == df3_tmp["s1_chain"])
].copy()
df3_s1["patch"] = df3_s1["interaction_id"].map(iid_to_s1_patch)

# --- Résidus S2 (chain_B) pour interactions homo uniquement ---
homo_iids = set(
    df_all[df_all["s2_actine"] == True].merge(
        df_int[["interaction_id", "chain_A_id", "chain_B_id"]],
        left_on=["subunit_1", "subunit_2"],
        right_on=["chain_A_id", "chain_B_id"],
        how="left"
    )["interaction_id"].dropna().astype(int)
)
id_to_s2chain = df_int.set_index("interaction_id")["chain_B_id"].str.lower()
df3_s2 = df3_tmp[
    df3_tmp["interaction_id"].isin(homo_iids) &
    df3_tmp["residue_number_canon_mafft"].notna() &
    (df3_tmp["chain"].str.lower() ==
     df3_tmp["interaction_id"].map(id_to_s2chain))
].copy()
df3_s2["patch"] = df3_s2["interaction_id"].map(iid_to_s2_patch)

# --- Combinaison S1 + S2 ---
merged = pd.concat([df3_s1, df3_s2], ignore_index=True)
merged = merged.dropna(subset=["patch"])
merged["c70"] = merged["interaction_id"].map(id_to_c70)

all_positions = sorted(merged["residue_number_canon_mafft"].dropna().unique())

# Moyenne équitable par cluster : moyenne-de-moyennes par C70
# max par (iid, canon) d'abord pour éviter le double-comptage si plusieurs
# résidus physiques mappent au même canon dans la même interaction


def equitable_mean(sub_cluster):
    """Pour un cluster, retourne Series(canon_pos → moyenne équitable)."""
    c70_profiles = []
    for _, c70_sub in sub_cluster.groupby("c70"):
        n_iids = c70_sub["interaction_id"].nunique()
        profile = (
            c70_sub.groupby(["interaction_id", "residue_number_canon_mafft"])[
                "buried_ASA_percent"]
            .max()
            .groupby(level="residue_number_canon_mafft").sum()
            .reindex(all_positions, fill_value=0)
            .fillna(0)
            / n_iids
        )
        c70_profiles.append(profile.values)
    if not c70_profiles:
        return pd.Series(0.0, index=all_positions)
    return pd.Series(np.mean(c70_profiles, axis=0), index=all_positions)


mean_per_cluster = {}
for patch, sub in merged.groupby("patch"):
    mean_per_cluster[patch] = equitable_mean(sub)

print(f"{len(mean_per_cluster)} clusters avec B-factor calculé")


def rewrite_bfactor(pdb_path, chain_letter, bfactor_map, pdb_content=None):
    """Garde uniquement chain_letter, remplace B-factor par bfactor_map.
    bfactor_map : {residue_number_structure: valeur}
    pdb_content : contenu PDB en mémoire (str) si pdb_path est un sentinel CIF converti.
    """
    lines_out = []
    lines = pdb_content.splitlines(
        keepends=True) if pdb_content else open(pdb_path)
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            if len(line) <= 26:
                continue
            if line[21] != chain_letter:
                continue
            try:
                resnum = int(line[22:26].strip())
            except ValueError:
                lines_out.append(line)
                continue
            bfac = bfactor_map.get(resnum, 0.0)
            line = f"{line[:60]}{bfac:6.2f}{line[66:]}"
        lines_out.append(line)
    if not pdb_content:
        lines.close()
    return "".join(lines_out)


ok, skip = 0, 0
for _, ref_row in df_ref.iterrows():
    patch = ref_row["patch"]
    pdb_id = ref_row["pdb_id"]
    chain_full = str(ref_row["chain"])
    chain_letter = chain_full.split("_")[-1]

    # Fichier PDB via table 8 (assembly) — fallback CIF → PDB via BioPython
    sub8 = df8[df8["pdb_id"] == pdb_id]
    pdb_path = None
    cif_path = None
    for _, r8 in sub8.iterrows():
        p = str(r8.get("biological_assembly_pdb_file", "")).strip()
        if p and p.lower() != "nan" and os.path.exists(p):
            pdb_path = p
            break
        c = str(r8.get("biological_assembly_cif_file", "")).strip()
        if c and c.lower() != "nan" and os.path.exists(c) and cif_path is None:
            cif_path = c

    if pdb_path is None and cif_path is not None:
        try:
            parser = _PDB.MMCIFParser(QUIET=True)
            struct = parser.get_structure(pdb_id, cif_path)
            buf = io.StringIO()
            io_obj = _PDB.PDBIO()
            io_obj.set_structure(struct)
            io_obj.save(buf)
            pdb_path = f"__cif_converted_{pdb_id}__"
            _cif_pdb_content = buf.getvalue()
        except Exception:
            pdb_path = None

    if pdb_path is None:
        skip += 1
        continue

    # Mapping canon_mafft → residue_structure pour cette chaîne de référence
    ref_chain_df = df3[df3["chain"] == chain_full][
        ["residue_number_canon_mafft", "residue_number_structure"]
    ].dropna().drop_duplicates()
    canon_to_struct = dict(zip(
        ref_chain_df["residue_number_canon_mafft"].astype(int),
        ref_chain_df["residue_number_structure"].astype(int)
    ))

    if patch not in mean_per_cluster:
        skip += 1
        continue

    bfactor_map = {}
    for canon, mean_asa in mean_per_cluster[patch].items():
        struct_res = canon_to_struct.get(int(canon))
        if struct_res is not None:
            bfactor_map[struct_res] = float(mean_asa)

    cif_content = _cif_pdb_content if pdb_path.startswith(
        "__cif_converted_") else None
    new_pdb = rewrite_bfactor(pdb_path, chain_letter,
                              bfactor_map, pdb_content=cif_content)

    out_path = OUT_DIR / f"{patch}.pdb"
    with open(out_path, "w") as f:
        f.write(new_pdb)
    ok += 1

print(f"Terminé — {ok} clusters sauvegardés, {skip} ignorés → {OUT_DIR}/")
