#!/usr/bin/env python
# coding: utf-8

import os as _os
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # cwd = racine projet (robuste, peu importe d'où on lance)


# # Contruction de graphe et filtrage

# ## Import librairie et data

# In[1]:


import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import math
from textwrap import fill

PDB_PATH     = "data/raw/pdb_entry_results.csv"
SUM_PATH     = "data/raw/ppi3d_actin_summary.csv"
ALL_DATA_RAW = "data/raw/all_data.csv"
ALL_DATA     = "data/filtered/filtered_all_data.csv"

df = pd.read_csv(PDB_PATH, sep=";")
df_summary = pd.read_csv(SUM_PATH, sep=";")


# ## 1er filtre : garder uniquement les interactinos proteines-proteines

# In[2]:


print("Nombre de lignes dans df pdb :", len(df))
print("Nombre de lignes dans df_summary :", len(df_summary))

print(df["Interface type"].unique())
print(df_summary["Interface type"].unique())

# garder seulement protéine-protéine
df = df[df["Interface type"].isin(["hetero", "homo"])]
df_summary = df_summary[df_summary["Interface type"].isin(["hetero", "homo"])]
print(df["Interface type"].unique())
print(df_summary["Interface type"].unique())

print("Nombre de lignes dans df avec juste interface protéine-protéine:",
      len(df))
print("Nombre de lignes dans df_summary avec juste interface", 
      " protéine-protéine:", len(df_summary))


# ## Création du graphe + 2eme filtrage (>= 5 actines connectées)

# In[3]:


graphs = {}
node_titles_dict = {}
node_labels_dict = {}
actin_nodes_dict = {}

valid_pdb = []
invalid_pdb = []

df_summary["Expect value"] = pd.to_numeric(df_summary["Expect value"], errors="coerce")
df["pdb_id"] = df["pdb_id"].astype(str).str.upper()

def normalize_title(value):
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().lower().split())

def display_label(node_id):
    suffix = str(node_id).split("_")[-1]
    return suffix[-1] if suffix else str(node_id)

target_proteins = set(
    df_summary.loc[df_summary["Expect value"] == 0.0, "Result protein"]
    .dropna()
    .map(normalize_title)
)

for pdb in df["pdb_id"].unique():
    sub = df[df["pdb_id"] == pdb]

    G = nx.Graph()
    node_titles = {}
    node_labels = {}

    for _, row in sub.iterrows():
        a = row["Interactor 1"]
        b = row["Interactor 2"]

        title_a = row["Interactor 1 title"]
        title_b = row["Interactor 2 title"]

        node_titles[a] = title_a
        node_titles[b] = title_b

        node_labels[a] = fill(title_a, width=20)
        node_labels[b] = fill(title_b, width=20)

        G.add_edge(a, b)

    graphs[pdb] = G
    node_titles_dict[pdb] = node_titles
    node_labels_dict[pdb] = node_labels

    actin_nodes = [
        node for node, title in node_titles.items()
        if normalize_title(title) in target_proteins
    ]
    actin_nodes_dict[pdb] = actin_nodes

    G_actin = G.subgraph(actin_nodes)

    keep = False
    if len(actin_nodes) > 0:
        keep = any(len(comp) >= 5 for comp in nx.connected_components(G_actin))

    if keep:
        valid_pdb.append(pdb)
    else:
        invalid_pdb.append(pdb)

print("PDB retenus :", len(valid_pdb))
print("PDB rejetés :", len(invalid_pdb))


# ## Visualisation des graphes retenues

# In[4]:


pdb_ids = valid_pdb

n_cols = 3
n_rows = math.ceil(len(pdb_ids) / n_cols)

fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 5 * n_rows))
axes = axes.flatten()

fig.suptitle(
    'PDB structures avec >= 5 actines',
    fontsize=16
)

for i, pdb in enumerate(pdb_ids):

    G = graphs[pdb]
    labels = node_labels_dict[pdb]
    actin_nodes = set(actin_nodes_dict[pdb])

    colors = [
        "red" if node in actin_nodes else "lightblue"
        for node in G.nodes()
    ]

    pos = nx.spring_layout(G, seed=42)

    nx.draw(
        G,
        pos,
        labels=labels,
        with_labels=True,
        node_color=colors,
        node_size=600,
        font_size=7,
        ax=axes[i]
    )

    axes[i].set_title(pdb)

for j in range(i + 1, len(axes)):
    axes[j].axis("off")

plt.tight_layout(rect=[0, 0, 1, 0.98])

plt.show()


# ## Créer les df filtrés

# In[5]:


# actines réelles (Expect value == 0)
real_actin = df_summary.loc[df_summary["Expect value"] == 0.0, "Result protein"].dropna().unique()

# filtrer df pdb
df_filtered = df[df["pdb_id"].isin(valid_pdb)]
print(f"df pdb: {len(df)} -> {len(df_filtered)}")

# filtrer all_data : valid_pdb + S1 = actine réelle + domain/fragment + excl. cluster 55649
valid_pdb_lower = [p.lower() for p in valid_pdb]
# --- Exclusion permanente de structures (PHACTR1 / 4b1z retiré de TOUTES les analyses) ---
EXCLUDE_PDBS = {"4b1z"}
valid_pdb = [p for p in valid_pdb if p.lower() not in EXCLUDE_PDBS]
valid_pdb_lower = [p for p in valid_pdb_lower if p not in EXCLUDE_PDBS]
df_filtered = df_filtered[~df_filtered["pdb_id"].astype(str).str.lower().isin(EXCLUDE_PDBS)]
df_all = pd.read_csv(ALL_DATA_RAW, low_memory=False)
df_all = df_all[df_all["pdb_id"].isin(valid_pdb_lower)]
df_all_pdb = df_all.copy()  # sauvegarde avant filtre S1/S2
df_all = df_all[df_all["subunit_1_title"].isin(real_actin)]

mask_all = (
    df_all["subunit_1_title"].str.contains("domain|fragment", case=False, na=False) |
    df_all["subunit_2_title"].str.contains("domain|fragment", case=False, na=False)
)
df_all = df_all[~mask_all]
df_all = df_all[df_all['s1_sequence_cluster_70'] != 55649]

# Récupérer les interactions inversées (ABP en S1, actine en S2) et swapper S1↔S2
df_reversed = df_all_pdb[
    ~df_all_pdb["subunit_1_title"].isin(real_actin) &
    df_all_pdb["subunit_2_title"].isin(real_actin)
].copy()
mask_rev = (
    df_reversed["subunit_1_title"].str.contains("domain|fragment", case=False, na=False) |
    df_reversed["subunit_2_title"].str.contains("domain|fragment", case=False, na=False)
)
df_reversed = df_reversed[~mask_rev]
df_reversed = df_reversed[df_reversed["s2_sequence_cluster_70"] != 55649]

swap_pairs = [
    ("subunit_1",                           "subunit_2"),
    ("subunit_1_title",                     "subunit_2_title"),
    ("s1_protein_type",                     "s2_protein_type"),
    ("scop_family_1",                       "scop_family_2"),
    ("s1_taxonomy_id",                      "s2_taxonomy_id"),
    ("s1_number_of_residues",               "s2_number_of_residues"),
    ("s1_sequence",                         "s2_sequence"),
    ("s1_sequence_cluster_95",              "s2_sequence_cluster_95"),
    ("s1_sequence_cluster_70",              "s2_sequence_cluster_70"),
    ("s1_sequence_cluster_40",              "s2_sequence_cluster_40"),
    ("s1_number_of_visible_residues",       "s2_number_of_visible_residues"),
    ("s1_binding_site_cluster_data_95",     "s2_binding_site_cluster_data_95"),
    ("s1_binding_site_cluster_data_70",     "s2_binding_site_cluster_data_70"),
    ("s1_binding_site_cluster_data_40",     "s2_binding_site_cluster_data_40"),
    ("s1_binding_site_cluster_data_40_area","s2_binding_site_cluster_data_40_area"),
]
for c1, c2 in swap_pairs:
    if c1 in df_reversed.columns and c2 in df_reversed.columns:
        df_reversed[c1], df_reversed[c2] = df_reversed[c2].copy(), df_reversed[c1].copy()

df_all = pd.concat([df_all, df_reversed], ignore_index=True)
df_all = df_all.drop_duplicates(subset=["pdb_id", "subunit_1", "subunit_2"])
print(f"  dont {len(df_reversed)} interactions inversées récupérées (actine replacée en S1)")

# normaliser les séquences en majuscules
df_all["s1_sequence"] = df_all["s1_sequence"].str.upper()
df_all["s2_sequence"] = df_all["s2_sequence"].str.upper()

df_all["s1_actine"] = df_all["subunit_1_title"].isin(real_actin)
df_all["s2_actine"] = df_all["subunit_2_title"].isin(real_actin)

print(f"all_data: {sum(1 for _ in open(ALL_DATA_RAW)) - 1} -> {len(df_all)} (filtre PDB + S1 actine + domain/fragment + excl. cluster 55649 + inversées récupérées)")

# filtrer summary : valid_pdb + Result protein = actine réelle + domain/fragment + excl. cluster 55649 (via df_all)
df_summary_filtered = df_summary[df_summary["PDB ID"].isin(valid_pdb_lower)]
df_summary_filtered = df_summary_filtered[df_summary_filtered["Result protein"].isin(real_actin)]

mask_summary = (
    df_summary_filtered["Result protein"].str.contains("domain|fragment", case=False, na=False) |
    df_summary_filtered["Interacts with"].str.contains("domain|fragment", case=False, na=False)
)
df_summary_filtered = df_summary_filtered[~mask_summary]
df_summary_filtered = df_summary_filtered[df_summary_filtered["PDB ID"].isin(df_all["pdb_id"])]

# pairs_all inclut maintenant les inversées swappées → mask_dir les garde dans le summary
pairs_all = set(zip(df_all["subunit_1"], df_all["subunit_2"]))
mask_dir = df_summary_filtered.apply(
    lambda r: (r["Protein ID"], r["Interactor ID"]) in pairs_all, axis=1
)
df_summary_filtered = df_summary_filtered[mask_dir]

print(f"df summary: {len(df_summary)} -> {len(df_summary_filtered)} (filtre PDB + real actin + domain/fragment + excl. cluster 55649 + sens cohérent df_all)")


# In[6]:


OUTPUT_SUMMARY   = "data/filtered/filtered_summary.csv"
OUTPUT_PDB_ENTRY = "data/filtered/filtered_pdb_entry.csv"

# Normalisation : PDB ID en minuscules, lettre de chaîne préservée
# Ex : "6VEC_A" → "6vec_A"  |  "6VEC_a" → "6vec_a"  (chaînes distinctes)
def norm_chain_id(s):
    if pd.isna(s): return s
    parts = str(s).split("_", 1)
    return parts[0].lower() + "_" + parts[1] if len(parts) == 2 else str(s).lower()

df_filtered = df_filtered.copy()
df_filtered["pdb_id"] = df_filtered["pdb_id"].str.lower()
for col in ["Interactor 1", "Interactor 2"]:
    if col in df_filtered.columns:
        df_filtered[col] = df_filtered[col].map(norm_chain_id)

df_summary_filtered = df_summary_filtered.copy()
df_summary_filtered["PDB ID"]        = df_summary_filtered["PDB ID"].str.lower()
df_summary_filtered["Protein ID"]    = df_summary_filtered["Protein ID"].map(norm_chain_id)
df_summary_filtered["Interactor ID"] = df_summary_filtered["Interactor ID"].map(norm_chain_id)

df_all["subunit_1"] = df_all["subunit_1"].map(norm_chain_id)
df_all["subunit_2"] = df_all["subunit_2"].map(norm_chain_id)

df_filtered.to_csv(OUTPUT_PDB_ENTRY, index=False)

df_summary_filtered["interaction_id"] = df_summary_filtered["Link to details"].str.split(" ").str[0].astype(int)
df_summary_filtered.to_csv(OUTPUT_SUMMARY, index=False)

df_all.to_csv(ALL_DATA, index=False)

print(f"Sauvegardé : {OUTPUT_PDB_ENTRY}")
print(f"Sauvegardé : {OUTPUT_SUMMARY}  ({len(df_summary_filtered)} lignes)")
print(f"Sauvegardé : {ALL_DATA}  ({len(df_all)} lignes)")


# In[7]:


# Sauvegarder les graphes individuels + CSV protéines par PDB
import os

output_dir_pdb = "data/visualisations/pdb_graphs"
os.makedirs(output_dir_pdb, exist_ok=True)

for pdb in valid_pdb:
    G = graphs[pdb]
    labels = node_labels_dict[pdb]
    actin_nodes = set(actin_nodes_dict[pdb])

    colors = ["red" if node in actin_nodes else "lightblue" for node in G.nodes()]

    fig, ax = plt.subplots(figsize=(6, 5))
    pos = nx.spring_layout(G, seed=42)
    nx.draw(
        G, pos,
        labels=labels,
        with_labels=True,
        node_color=colors,
        node_size=600,
        font_size=7,
        font_color="black",
        ax=ax,
    )
    ax.set_title(pdb, fontsize=12)
    plt.tight_layout()
    plt.savefig(f"{output_dir_pdb}/{pdb}.png", dpi=120, bbox_inches="tight")
    plt.close(fig)

print(f"Graphes individuels sauvegardés pour {len(valid_pdb)} PDB → {output_dir_pdb}/")

# CSV : une ligne par chaîne unique, avec flag is_actin
# norm_chain_id est défini dans la cellule c58344f7
rows = []
for pdb in valid_pdb:
    actin_nodes = set(actin_nodes_dict[pdb])
    for chain, title in node_titles_dict[pdb].items():
        rows.append({
            "pdb_id": pdb.lower(),
            "chain": norm_chain_id(chain),
            "protein": title,
            "is_actin": chain in actin_nodes,
        })

df_proteins_pdb = pd.DataFrame(rows)
df_proteins_pdb.to_csv("data/filtered/proteins_per_pdb.csv", index=False)
print(f"proteins_per_pdb.csv sauvegardé ({len(df_proteins_pdb)} lignes)")


# In[8]:


# Normalisation des IDs de chaînes dans les fichiers de détails
# (données brutes téléchargées → norm_chain_id appliqué ici, pas dans get_interaction_details.py)

import pandas as pd
from pathlib import Path

DETAILS = Path("data/filtered/details")

def norm_col(df, col):
    if col in df.columns:
        df[col] = df[col].apply(
            lambda s: (str(s).split("_")[0].lower() + "_" + str(s).split("_", 1)[1])
            if pd.notna(s) and "_" in str(s) else (str(s).lower() if pd.notna(s) else s)
        )
    return df

files_cols = {
    "1.interactions.csv":          ["chain_A_id", "chain_B_id", "pdb_id"],
    "2.proteins.csv":              ["chain_id"],
    "3.interface_residues.csv":    ["chain"],
    "4.inter-residue_contacts.csv":["chain_A_id", "chain_B_id"],
}

for fname, cols in files_cols.items():
    path = DETAILS / fname
    if not path.exists():
        print(f"Non trouvé : {path}")
        continue
    df = pd.read_csv(path, dtype=str)
    for col in cols:
        df = norm_col(df, col)
    df.to_csv(path, index=False)
    print(f"Normalisé : {fname}")


# ## Ajout % ASA buried dans la table 4 (inter-residue contacts)

# In[9]:


from pathlib import Path

TABLE3 = "data/filtered/details/3.interface_residues.csv"
TABLE4 = "data/filtered/details/4.inter-residue_contacts.csv"

_tables_ready = Path(TABLE3).exists() and Path(TABLE4).exists()

if not _tables_ready:
    print("Fichiers non disponibles")
    print("   -> Sera complete a l'etape 5.7 du pipeline (apres get_interaction_details)")
    df4 = None
else:
    df3 = pd.read_csv(TABLE3)
    df4 = pd.read_csv(TABLE4)

    df3["buried_ASA_percent"] = pd.to_numeric(
        df3["buried_ASA_percent"].astype(str).str.replace("%", "", regex=False),
        errors="coerce"
    )

    df3_dedup = (
        df3.groupby(["interaction_id", "chain", "residue_number_sequence"],
                    sort=False)["buried_ASA_percent"]
        .max()
        .reset_index()
    )

    df4 = df4.drop(columns=["asa_pct_A", "asa_pct_B"], errors="ignore")

    df4 = df4.merge(
        df3_dedup.rename(columns={
            "chain": "chain_A_id",
            "residue_number_sequence": "residue_A_sequence",
            "buried_ASA_percent": "asa_pct_A"
        }),
        on=["interaction_id", "chain_A_id", "residue_A_sequence"],
        how="left"
    )

    df4 = df4.merge(
        df3_dedup.rename(columns={
            "chain": "chain_B_id",
            "residue_number_sequence": "residue_B_sequence",
            "buried_ASA_percent": "asa_pct_B"
        }),
        on=["interaction_id", "chain_B_id", "residue_B_sequence"],
        how="left"
    )

    ok_A = df4["asa_pct_A"].notna().sum()
    ok_B = df4["asa_pct_B"].notna().sum()
    n = len(df4)
    print(f"Lignes table 4 : {n}")
    print(f"asa_pct_A renseigne : {ok_A} ({ok_A/n*100:.1f} %)")
    print(f"asa_pct_B renseigne : {ok_B} ({ok_B/n*100:.1f} %)")
    print()
    print(df4[["interaction_id", "residue_A_name", "residue_A_sequence", "asa_pct_A",
               "residue_B_name", "residue_B_sequence", "asa_pct_B", "contact_area"]].head(6))


# In[10]:


# Sauvegarder table 4 enrichie
if df4 is None:
    print("⚠️  Enrichissement ASA ignoré (fichiers non disponibles)")
else:
    df4.to_csv(TABLE4, index=False)
    print(f"Sauvegardé : {TABLE4}  ({len(df4)} lignes, {len(df4.columns)} colonnes)")
    print(f"Nouvelles colonnes : asa_pct_A, asa_pct_B")

    # Retirer les colonnes ASA de filtered_all_data.csv si elles y ont été ajoutées
    df_all = pd.read_csv(ALL_DATA)
    asa_cols = ["mean_asa_pct_S1", "median_asa_pct_S1", "n_residues_S1",
                "mean_asa_pct_S2", "median_asa_pct_S2", "n_residues_S2"]
    cols_to_drop = [c for c in asa_cols if c in df_all.columns]
    if cols_to_drop:
        df_all = df_all.drop(columns=cols_to_drop)
        df_all.to_csv(ALL_DATA, index=False)
        print(f"Retiré de filtered_all_data.csv : {cols_to_drop}")

