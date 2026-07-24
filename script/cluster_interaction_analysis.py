#!/usr/bin/env python
# coding: utf-8

import os as _os
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # cwd = racine projet (robuste, peu importe d'où on lance)


# # Analyse des clusters d'interaction actine
# 
# **Objectif** : à partir de `filtered_all_data.csv`, grouper les interactions par cluster de séquence, calculer des stats de fréquence des résidus d'interface, et représenter chaque cluster sous forme de graphe NetworkX.
# 
# **Rappel** : chaque ligne de `filtered_all_data.csv` est un **représentant d'interaction** (pas une moyenne). Les résidus d'interface viennent de `details/` et les contacts sont toujours **S1 (actine) ↔ S2 (partenaire)** — les deux côtés sont donc nécessaires pour les arêtes du graphe.

# ## 0. Imports et chemins

# In[1]:


import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import textwrap
from matplotlib.lines import Line2D
import os
import re

DETAILS_PATH    = "data/filtered/details"

ALL_DATA        = "data/filtered/filtered_all_data.csv"
SUMMARY         = "data/filtered/filtered_summary.csv"
INTERACTIONS_1  = os.path.join(DETAILS_PATH, "1.interactions.csv")
INTER_RES_4     = os.path.join(DETAILS_PATH, "4.inter-residue_contacts.csv")


# In[2]:


# ═══════════════════════════════════════════════════════════════════════════
# Fonctions utilitaires — partagées par toutes les cellules d'analyse
# ═══════════════════════════════════════════════════════════════════════════

ACTIN_LABEL = "Actine"   # nœud unifié pour toutes les isoformes d'actine

# ── Graphes par patch ──────────────────────────────────────────────────────

def build_patch_graphs(DF, col_patch, real_actin):
    """
    Construit un graphe NetworkX par patch.
    Toutes les isoformes d'actine sont fusionnées en un nœud unique "Actine".
    Nœuds : side='S1'/'S2', is_actin=bool.
    Arêtes pondérées par nombre d'interactions.
    """
    real_actin_set = set(real_actin)

    def _label(title):
        return ACTIN_LABEL if title in real_actin_set else title

    graphs = {}
    for patch in DF[col_patch].dropna().unique():
        patch_data = DF[DF[col_patch] == patch]
        G, n_homo, n_total = nx.Graph(), 0, 0
        for _, row in patch_data.iterrows():
            s1_raw, s2_raw = row['subunit_1_title'], row['subunit_2_title']
            if pd.isna(s1_raw) or pd.isna(s2_raw):
                continue
            s1_is_actin = s1_raw in real_actin_set
            s2_is_actin = s2_raw in real_actin_set
            s1 = _label(s1_raw)
            s2 = _label(s2_raw)
            if not G.has_node(s1):
                G.add_node(s1, side='S1', is_actin=s1_is_actin)
            if not G.has_node(s2):
                G.add_node(s2, side='S2', is_actin=s2_is_actin)
            if G.has_edge(s1, s2):
                G[s1][s2]['weight'] += 1
            else:
                G.add_edge(s1, s2, weight=1)
            n_total += 1
            if s1_is_actin and s2_is_actin:
                n_homo += 1
        if n_homo == 0:
            itype = 'hetero'
        elif n_homo == n_total:
            itype = 'homo'
        else:
            itype = 'mixed'
        G.graph['interaction_type'] = itype
        graphs[patch] = G
    return graphs


def build_patch_info(DF, col_patch, patch_graphs, df_res_contacts=None):
    """
    Construit df_patches_infos depuis patch_graphs.
    df_res_contacts : si fourni, calcule avg_contacts.
    Si DF contient s2_sequence_cluster_70, ajoute s2_seq_clusters et n_s2_seq_clusters.
    """
    has_s2_seq = 's2_sequence_cluster_70' in DF.columns
    rows = []
    for patch, G in patch_graphs.items():
        sub = DF[DF[col_patch] == patch]
        ids   = list(sub['interaction_id'].dropna().astype(int).unique())
        s2_np = set(n for n, d in G.nodes(data=True) if d['side'] == 'S2' and not d['is_actin'])
        row   = {
            'patch':            patch,
            'interaction_type': G.graph['interaction_type'],
            'n_noeuds':         G.number_of_nodes(),
            'n_arretes':        G.number_of_edges(),
            'n_interactions':   len(sub),
            's2_partners':      list(s2_np),
            'n_s2_partners':    len(s2_np),
            'ids_interactions': ids,
        }
        if df_res_contacts is not None and ids:
            c = df_res_contacts[df_res_contacts['interaction_id'].isin(ids)]
            row['avg_contacts'] = round(len(c) / len(ids), 2)
        if has_s2_seq:
            s2_cls = sorted(sub['s2_sequence_cluster_70'].dropna().astype(int).astype(str).unique())
            row['s2_seq_clusters']   = ', '.join(s2_cls)
            row['n_s2_seq_clusters'] = len(s2_cls)
        rows.append(row)

    df = pd.DataFrame(rows).sort_values('n_interactions', ascending=False)
    s1_map = {p: sorted(set(n for n, d in G.nodes(data=True) if d['side'] == 'S1'))
              for p, G in patch_graphs.items()}
    s2_map = {p: sorted(set(n for n, d in G.nodes(data=True) if d['side'] == 'S2' and not d['is_actin']))
              for p, G in patch_graphs.items()}
    df['s2_partners']   = df['patch'].map(lambda p: s2_map.get(p, []))
    df['n_s2_partners'] = df['s2_partners'].apply(len)
    return df


def _node_colors(G, real_actin):
    """S1 → rouge | S2 non actine → vert."""
    return [
        'red' if d['side'] == 'S1' else 'lightgreen' for _, d in G.nodes(data=True)
    ]


def _patch_legend():
    """Légende standard patches."""
    return [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='red',        markersize=10, label='S1 actine'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='lightgreen', markersize=10, label='S2 non actine'),
    ]


def plot_top_patches(DF, col_patch, patch_graphs, real_actin, n=10):
    """Affiche les n patches les plus peuplés en grille verticale."""
    top = [(p, patch_graphs[p]) for p in DF[col_patch].value_counts().head(n).index
           if p in patch_graphs]
    fig, axes = plt.subplots(nrows=len(top), ncols=1, figsize=(10, 5 * len(top)))
    if len(top) == 1:
        axes = [axes]
    for ax, (patch, G) in zip(axes, top):
        pos    = nx.spring_layout(G, seed=42)
        labels = {n: textwrap.fill(n, width=20) for n in G.nodes()}
        widths = [min(G[u][v]['weight'], 10) for u, v in G.edges()]
        nx.draw(G, pos, node_color=_node_colors(G, real_actin), labels=labels,
                font_size=8, node_size=600, ax=ax, width=widths)
        n_inter = G.graph.get('n_interactions', G.number_of_edges())
        ax.set_title(f"Patch {patch} — {G.graph['interaction_type']} — {n_inter} interactions", fontsize=12)
    fig.legend(handles=_patch_legend(), loc='upper right', fontsize=10)
    plt.tight_layout()
    plt.show()


def save_patch_graphs(DF, col_patch, patch_graphs, real_actin, output_dir):
    """Sauvegarde un PNG par patch dans output_dir."""
    os.makedirs(output_dir, exist_ok=True)
    for patch in DF[col_patch].value_counts().index:
        if patch not in patch_graphs:
            continue
        G      = patch_graphs[patch]
        pos    = nx.spring_layout(G, seed=42)
        labels = {n: textwrap.fill(n, width=20) for n in G.nodes()}
        widths = [min(G[u][v]['weight'], 10) for u, v in G.edges()]
        fig, ax = plt.subplots(figsize=(10, 6))
        nx.draw(G, pos, node_color=_node_colors(G, real_actin), labels=labels,
                font_size=8, node_size=600, ax=ax, width=widths)
        n_inter = G.graph.get('n_interactions', G.number_of_edges())
        ax.set_title(f"Patch {patch} — {G.graph['interaction_type']} — {n_inter} interactions", fontsize=12)
        fig.legend(handles=_patch_legend(), loc='upper right', fontsize=10)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{patch}.png", dpi=120, bbox_inches='tight')
        plt.close(fig)
    print(f"Sauvegardé : {len(patch_graphs)} graphes → {output_dir}/")


print("Fonctions utilitaires définies.")


# ## 1. Chargement données

# In[3]:


df_all     = pd.read_csv(ALL_DATA)
df_summary = pd.read_csv(SUMMARY)  # contient déjà interaction_id
df_interactions = pd.read_csv(INTERACTIONS_1)
df_res = pd.read_csv(INTER_RES_4)

print(f"df_all    : {len(df_all)} lignes")
print(f"df_summary: {len(df_summary)} lignes")
print(f"df_res    : {df_res.shape}")


# ## Tag des vrais actines

# In[4]:


# actine reelle = expect value 0
df_summary["Expect value"] = df_summary["Expect value"].astype(float)
real_actin = df_summary[df_summary["Expect value"] == 0]["Result protein"].unique()
for actine in real_actin:
    print(actine)



# ## Quel seuil de cluster de similarité d'interface regroupe tout les actines ?

# In[5]:


actine_df = df_all[df_all["s1_actine"] == True]

DF = actine_df

seuils = [95, 70, 40]

for seuil in seuils:
    cluster_col = f"s1_sequence_cluster_{seuil}"

    print(f"\n--- Cluster level {seuil}% ---")

    n_clusters = DF[cluster_col].nunique()

    print(f"Nombre de clusters contenant des actines : {n_clusters}")


# In[6]:


df_merged = df_all.merge(
    df_summary,
    left_on=["subunit_1", "subunit_2"],
    right_on=["Protein ID", "Interactor ID"],
    how="left"
)
DF = df_merged
print(f"df_merged : {len(DF)} interactions")


# ## Un graphe par patch

# In[7]:


DF           = df_merged
col_patch    = 's1_binding_site_cluster_data_70'
col_patch_s2 = 's2_binding_site_cluster_data_70'
df_homo      = DF[DF['subunit_2_title'].isin(real_actin)]

# Étape 1 : graphes S1 de base (toutes isoformes actine → nœud "Actine")
patch_graphs_s1 = build_patch_graphs(DF, col_patch, real_actin)

# Étape 2 : enrichir chaque graphe avec les interactions homo où ce patch est S2
s2_only_patches = {}
for patch2, sub_s2 in df_homo.groupby(col_patch_s2):
    G = patch_graphs_s1.get(patch2)
    is_new = G is None
    if is_new:
        G = nx.Graph()
        G.graph['interaction_type'] = 'homo'

    for _, row in sub_s2.iterrows():
        s1t_raw, s2t_raw = row['subunit_1_title'], row['subunit_2_title']
        if pd.isna(s1t_raw) or pd.isna(s2t_raw):
            continue
        s1t = ACTIN_LABEL
        s2t = ACTIN_LABEL
        if not G.has_node(s1t):
            G.add_node(s1t, side='S1', is_actin=True)
        if not G.has_node(s2t):
            G.add_node(s2t, side='S2', is_actin=True)
        if G.has_edge(s1t, s2t):
            G[s1t][s2t]['weight'] += 1
        else:
            G.add_edge(s1t, s2t, weight=1)

    if not is_new:
        if G.graph.get('interaction_type') == 'hetero':
            G.graph['interaction_type'] = 'mixed'
        patch_graphs_s1[patch2] = G
    else:
        s2_only_patches[patch2] = G

patch_graphs = {**patch_graphs_s1, **s2_only_patches}

# Étape 3 : infos de base
df_infos_s1 = build_patch_info(DF, col_patch, patch_graphs_s1)
df_infos_s2 = build_patch_info(df_homo, col_patch_s2, s2_only_patches)
df_patches_infos = pd.concat([df_infos_s1, df_infos_s2], ignore_index=True)

# Étape 4 : n_interactions = nb total d'apparitions dans all_data (S1 + S2)
s1_counts    = DF[col_patch].value_counts()
s2_counts    = DF[col_patch_s2].value_counts()
total_counts = s1_counts.add(s2_counts, fill_value=0).astype(int)
df_patches_infos['n_interactions'] = (
    df_patches_infos['patch'].astype(str).map(total_counts).fillna(0).astype(int)
)

# Étape 5 : stocker le compte dans chaque graphe (utilisé dans les titres des PNG)
for patch, G in patch_graphs.items():
    G.graph['n_interactions'] = int(total_counts.get(patch, 0))

# Étape 6 : PDB impliqués dans les interactions homo pour chaque patch (S1 + S2)
pdbs_s1_homo = df_homo.groupby(col_patch)['pdb_id'].apply(lambda x: set(x.dropna()))
pdbs_s2_homo = df_homo.groupby(col_patch_s2)['pdb_id'].apply(lambda x: set(x.dropna()))
all_homo_patches = set(pdbs_s1_homo.index) | set(pdbs_s2_homo.index)
homo_pdb_map = {
    p: sorted((set(pdbs_s1_homo.get(p) or [])) | (set(pdbs_s2_homo.get(p) or [])))
    for p in all_homo_patches
}
df_patches_infos['homo_pdb_ids'] = df_patches_infos['patch'].astype(str).map(homo_pdb_map)
df_patches_infos['n_homo_pdbs']  = df_patches_infos['homo_pdb_ids'].apply(
    lambda x: len(x) if isinstance(x, list) else 0
)

print(f"Patchs S1 base     : {len(patch_graphs_s1)}")
print(f"Patchs S2 nouveaux : {len(s2_only_patches)}")
print(f"Total              : {len(patch_graphs)}")
print(df_patches_infos['interaction_type'].value_counts())
print()
for p in ['6685_109', '6685_274', '6685_2']:
    row = df_patches_infos[df_patches_infos['patch'] == p]
    if not row.empty:
        print(f"{p} → n={row['n_interactions'].values[0]}  type={row['interaction_type'].values[0]}  n_homo_pdbs={row['n_homo_pdbs'].values[0]}")

plot_top_patches(DF, col_patch, patch_graphs_s1, real_actin)


# In[8]:


save_patch_graphs(DF, col_patch, patch_graphs, real_actin, 'data/visualisations/patch_graphs/s1_binding_site')
df_patches_infos.to_csv('data/filtered/patches_infos_s1_binding_site.csv', index=False)
print(f"Sauvegardé : patches_infos_s1_binding_site.csv ({len(df_patches_infos)} patchs)")


# ## Graphe global interactif (pyvis)

# In[9]:


import math
from pyvis.network import Network
from IPython.display import IFrame
from collections import defaultdict

DF     = df_merged
col_s1 = "s1_binding_site_cluster_data_70"
col_s2 = "s2_binding_site_cluster_data_70"

G = nx.Graph()
edge_c70s = {}

for _, row in DF.iterrows():
    patch    = row[col_s1]
    s2_title = row['subunit_2_title']
    if pd.isna(patch) or pd.isna(s2_title):
        continue
    if not G.has_node(patch):
        G.add_node(patch, side="cluster", count=0)
    G.nodes[patch]['count'] += 1
    if s2_title in real_actin:
        s2_cluster = row[col_s2]
        if pd.isna(s2_cluster):
            continue
        if not G.has_node(s2_cluster):
            G.add_node(s2_cluster, side="cluster", count=0)
        G.nodes[s2_cluster]['count'] += 1
        target = s2_cluster
    else:
        if not G.has_node(s2_title):
            G.add_node(s2_title, side="protein", count=0)
        G.nodes[s2_title]['count'] += 1
        target = s2_title
    if G.has_edge(patch, target):
        G[patch][target]['weight'] += 1
    else:
        G.add_edge(patch, target, weight=1)
    c70 = row.get('cluster_data_70')
    if pd.notna(c70):
        key = (patch, target)
        edge_c70s.setdefault(key, set()).add(str(c70))

# ── Calcul des familles de protéines (avant le layout) ────────────────────
GENERIC_WORDS = {
    'protein', 'actin', 'like', 'type', 'associated', 'binding',
    'domain', 'subunit', 'complex', 'chain', 'alpha', 'beta', 'gamma',
    'and', 'the', 'isoform', 'heavy', 'muscle', 'cardiac', 'skeletal',
    'delta', 'epsilon', 'smooth', 'enteric', 'cytoplasmic',
    'cell', 'nuclear', 'membrane', 'cytoskeletal',
}

def protein_key_words(name):
    name = re.sub(r'^Isoform\s+\S+\s+of\s+', '', name, flags=re.IGNORECASE)
    name = name.split(',')[0]
    words = re.findall(r'[A-Za-z]{3,}', name.lower())
    return {w for w in words if w not in GENERIC_WORDS}

protein_nodes_list = [n for n, d in G.nodes(data=True) if d['side'] == 'protein']
word_to_proteins   = defaultdict(set)
for pname in protein_nodes_list:
    for word in protein_key_words(pname):
        word_to_proteins[word].add(pname)

# Paires de protéines de la même famille + mot-clé associé.
# ETOILE (chaque membre relié à UN représentant) plutôt que graphe complet :
# meme regroupement au layout, mais beaucoup moins de traits verts -> centre lisible.
family_pairs = {}  # (u, v) → word
for word, proteins in word_to_proteins.items():
    if len(proteins) < 2:
        continue
    proteins = sorted(proteins)
    hub = proteins[0]
    for v in proteins[1:]:
        key = (min(hub, v), max(hub, v))
        if key not in family_pairs:
            family_pairs[key] = word

# ── Layout ORGANIQUE : spring = graine, physique barnesHut = clusters separes ──
import json as _json
G_layout = G.copy()
for (u, v) in family_pairs:
    if not G_layout.has_edge(u, v):
        G_layout.add_edge(u, v, weight=0.4)
pos = nx.spring_layout(G_layout, seed=42, k=2.5, scale=2000)

# taille pilotee par le nombre d'INTERACTIONS contenant ce cluster (count)
cluster_counts = [d['count'] for _, d in G.nodes(data=True) if d['side'] == 'cluster']
protein_counts = [d['count'] for _, d in G.nodes(data=True) if d['side'] == 'protein']
min_cc, max_cc = (min(cluster_counts), max(cluster_counts)) if cluster_counts else (0, 1)
min_pc, max_pc = (min(protein_counts), max(protein_counts)) if protein_counts else (0, 1)
def _node_size(d):
    if d['side'] == 'cluster':
        return 24 + 66 * (d['count'] - min_cc) / (max_cc - min_cc) if max_cc > min_cc else 32
    return 20 + 42 * (d['count'] - min_pc) / (max_pc - min_pc) if max_pc > min_pc else 26
node_sizes = {n: _node_size(d) for n, d in G.nodes(data=True)}
print(f"Nœuds : {G.number_of_nodes()}  |  Arêtes : {G.number_of_edges()}")

net = Network(height="900px", width="100%", bgcolor="#ffffff")
_phys_opts = {
    "physics": {
        "enabled": True, "solver": "barnesHut",
        "barnesHut": {"gravitationalConstant": -38000, "centralGravity": 0.30,
                      "springLength": 130, "springConstant": 0.10,
                      "damping": 0.09, "avoidOverlap": 1},
        "stabilization": {"enabled": True, "iterations": 400, "fit": True},
        "minVelocity": 0.75,
    },
    "edges": {"smooth": False},
}
net.set_options(_json.dumps(_phys_opts))

for n, d in G.nodes(data=True):
    is_cluster = d['side'] == 'cluster'
    color = '#e05252' if is_cluster else '#52e07a'
    label = ' ' if is_cluster else n
    font  = {'size': 0} if is_cluster else {'size': 12, 'background': 'white', 'strokeWidth': 0}
    x, y  = pos[n]
    net.add_node(n, label=label, color=color, size=node_sizes[n],
                 title=f"{n}<br>{d['count']} interactions", font=font,
                 x=float(x), y=float(y), physics=True)

for u, v in G.edges():
    c70s = edge_c70s.get((u, v)) or edge_c70s.get((v, u)) or set()
    c70_str = "<br>".join(sorted(c70s)) if c70s else "—"
    title = f"cluster_data_70 :<br>{c70_str}<br>({G[u][v]['weight']} interactions)"
    net.add_edge(u, v, color='#00000088', width=3.5, title=title)

for (u, v), word in family_pairs.items():
    net.add_edge(u, v, color='#00aa44', width=4.0, dashes=True,
                 title=f"Même famille : {word}")
print(f"Arêtes de même famille : {len(family_pairs)}")

# ── Anti-chevauchement DANS LE NAVIGATEUR (rayon injecte explicitement) ──────
_sizes_js = _json.dumps({str(n): float(node_sizes[n]) for n in G.nodes()})
_deoverlap_js = """
  var NODE_SIZES = __SIZES__;
  network.once('stabilizationIterationsDone', function(){
    network.setOptions({physics:{enabled:false}});
    setTimeout(function(){
      var ids = network.body.nodeIndices, B = network.body.nodes;
      function R(id){ return (NODE_SIZES[id] || 10); }
      // on continue TANT QU'il reste au moins un chevauchement (pas un seuil de mouvement,
      // qui s'arretait trop tot dans les amas denses).
      for (var it=0; it<4000; it++){
        var overlaps = 0;
        for (var a=0; a<ids.length; a++){
          var na = B[ids[a]];
          for (var b=a+1; b<ids.length; b++){
            var nb = B[ids[b]];
            var dx = nb.x-na.x, dy = nb.y-na.y;
            var dist = Math.sqrt(dx*dx + dy*dy) || 0.01;
            var mind = R(ids[a]) + R(ids[b]) + 12;
            if (dist < mind){
              overlaps++;
              var pp = (mind-dist)/2 + 0.5, ux = dx/dist, uy = dy/dist;  // +0.5 : depasse le contact
              na.x -= ux*pp; na.y -= uy*pp;
              nb.x += ux*pp; nb.y += uy*pp;
            }
          }
        }
        if (overlaps === 0) break;
      }
      var upd = [];
      ids.forEach(function(id){ upd.push({id:id, x:B[id].x, y:B[id].y}); });
      network.body.data.nodes.update(upd);
      network.fit();
    }, 50);
  });
""".replace("__SIZES__", _sizes_js)

_html = net.generate_html()
for _pat in ["network = new vis.Network(container, data, options);",
             "var network = new vis.Network(container, data, options);"]:
    if _pat in _html:
        _html = _html.replace(_pat, _pat + _deoverlap_js)
        break

nx.write_graphml(G, "data/visualisations/patch_graphs/global.graphml")
out = "data/visualisations/patch_graphs/global_pyvis.html"
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w", encoding="utf-8") as _f:
    _f.write(_html)
print(f"Sauvegardé : {out}")
IFrame(out, width="100%", height=900)


# ## Graphe global interactif — SUPERCLUSTERS
# 
# Même réseau que ci-dessus mais les nœuds rouges sont les **superclusters** (`s1_supercluster` / `s2_supercluster`, regroupement plus grossier des binding-site clusters_70) au lieu des clusters_70. Verts = ABP. Sortie : `global_supercluster_pyvis.html`.

# In[10]:


import math
from pyvis.network import Network
from IPython.display import IFrame
from collections import defaultdict

DF     = df_merged
col_s1 = "s1_supercluster"
col_s2 = "s2_supercluster"

G = nx.Graph()
edge_c70s = {}

for _, row in DF.iterrows():
    patch    = row[col_s1]
    s2_title = row['subunit_2_title']
    if pd.isna(patch) or pd.isna(s2_title):
        continue
    if not G.has_node(patch):
        G.add_node(patch, side="cluster", count=0)
    G.nodes[patch]['count'] += 1
    if s2_title in real_actin:
        s2_cluster = row[col_s2]
        if pd.isna(s2_cluster):
            continue
        if not G.has_node(s2_cluster):
            G.add_node(s2_cluster, side="cluster", count=0)
        G.nodes[s2_cluster]['count'] += 1
        target = s2_cluster
    else:
        if not G.has_node(s2_title):
            G.add_node(s2_title, side="protein", count=0)
        G.nodes[s2_title]['count'] += 1
        target = s2_title
    if G.has_edge(patch, target):
        G[patch][target]['weight'] += 1
    else:
        G.add_edge(patch, target, weight=1)
    c70 = row.get('cluster_data_70')
    if pd.notna(c70):
        key = (patch, target)
        edge_c70s.setdefault(key, set()).add(str(c70))

# ── Calcul des familles de protéines (avant le layout) ────────────────────
GENERIC_WORDS = {
    'protein', 'actin', 'like', 'type', 'associated', 'binding',
    'domain', 'subunit', 'complex', 'chain', 'alpha', 'beta', 'gamma',
    'and', 'the', 'isoform', 'heavy', 'muscle', 'cardiac', 'skeletal',
    'delta', 'epsilon', 'smooth', 'enteric', 'cytoplasmic',
    'cell', 'nuclear', 'membrane', 'cytoskeletal',
}

def protein_key_words(name):
    name = re.sub(r'^Isoform\s+\S+\s+of\s+', '', name, flags=re.IGNORECASE)
    name = name.split(',')[0]
    words = re.findall(r'[A-Za-z]{3,}', name.lower())
    return {w for w in words if w not in GENERIC_WORDS}

protein_nodes_list = [n for n, d in G.nodes(data=True) if d['side'] == 'protein']
word_to_proteins   = defaultdict(set)
for pname in protein_nodes_list:
    for word in protein_key_words(pname):
        word_to_proteins[word].add(pname)

# Paires de protéines de la même famille + mot-clé associé.
# ETOILE (chaque membre relié à UN représentant) plutôt que graphe complet :
# meme regroupement au layout, mais beaucoup moins de traits verts -> centre lisible.
family_pairs = {}  # (u, v) → word
for word, proteins in word_to_proteins.items():
    if len(proteins) < 2:
        continue
    proteins = sorted(proteins)
    hub = proteins[0]
    for v in proteins[1:]:
        key = (min(hub, v), max(hub, v))
        if key not in family_pairs:
            family_pairs[key] = word

# ── Layout ORGANIQUE : spring = graine, physique barnesHut = clusters separes ──
import json as _json
G_layout = G.copy()
for (u, v) in family_pairs:
    if not G_layout.has_edge(u, v):
        G_layout.add_edge(u, v, weight=0.4)
pos = nx.spring_layout(G_layout, seed=42, k=2.5, scale=2000)

# taille pilotee par le nombre d'INTERACTIONS contenant ce cluster (count)
cluster_counts = [d['count'] for _, d in G.nodes(data=True) if d['side'] == 'cluster']
protein_counts = [d['count'] for _, d in G.nodes(data=True) if d['side'] == 'protein']
min_cc, max_cc = (min(cluster_counts), max(cluster_counts)) if cluster_counts else (0, 1)
min_pc, max_pc = (min(protein_counts), max(protein_counts)) if protein_counts else (0, 1)
def _node_size(d):
    if d['side'] == 'cluster':
        return 24 + 66 * (d['count'] - min_cc) / (max_cc - min_cc) if max_cc > min_cc else 32
    return 20 + 42 * (d['count'] - min_pc) / (max_pc - min_pc) if max_pc > min_pc else 26
node_sizes = {n: _node_size(d) for n, d in G.nodes(data=True)}
print(f"Nœuds : {G.number_of_nodes()}  |  Arêtes : {G.number_of_edges()}")

net = Network(height="900px", width="100%", bgcolor="#ffffff")
_phys_opts = {
    "physics": {
        "enabled": True, "solver": "barnesHut",
        "barnesHut": {"gravitationalConstant": -38000, "centralGravity": 0.30,
                      "springLength": 130, "springConstant": 0.10,
                      "damping": 0.09, "avoidOverlap": 1},
        "stabilization": {"enabled": True, "iterations": 400, "fit": True},
        "minVelocity": 0.75,
    },
    "edges": {"smooth": False},
}
net.set_options(_json.dumps(_phys_opts))

for n, d in G.nodes(data=True):
    is_cluster = d['side'] == 'cluster'
    color = '#e05252' if is_cluster else '#52e07a'
    label = ' ' if is_cluster else n
    font  = {'size': 0} if is_cluster else {'size': 12, 'background': 'white', 'strokeWidth': 0}
    x, y  = pos[n]
    net.add_node(n, label=label, color=color, size=node_sizes[n],
                 title=f"{n}<br>{d['count']} interactions", font=font,
                 x=float(x), y=float(y), physics=True)

for u, v in G.edges():
    c70s = edge_c70s.get((u, v)) or edge_c70s.get((v, u)) or set()
    c70_str = "<br>".join(sorted(c70s)) if c70s else "—"
    title = f"cluster_data_70 :<br>{c70_str}<br>({G[u][v]['weight']} interactions)"
    net.add_edge(u, v, color='#00000088', width=3.5, title=title)

for (u, v), word in family_pairs.items():
    net.add_edge(u, v, color='#00aa44', width=4.0, dashes=True,
                 title=f"Même famille : {word}")
print(f"Arêtes de même famille : {len(family_pairs)}")

# ── Anti-chevauchement DANS LE NAVIGATEUR (rayon injecte explicitement) ──────
_sizes_js = _json.dumps({str(n): float(node_sizes[n]) for n in G.nodes()})
_deoverlap_js = """
  var NODE_SIZES = __SIZES__;
  network.once('stabilizationIterationsDone', function(){
    network.setOptions({physics:{enabled:false}});
    setTimeout(function(){
      var ids = network.body.nodeIndices, B = network.body.nodes;
      function R(id){ return (NODE_SIZES[id] || 10); }
      // on continue TANT QU'il reste au moins un chevauchement (pas un seuil de mouvement,
      // qui s'arretait trop tot dans les amas denses).
      for (var it=0; it<4000; it++){
        var overlaps = 0;
        for (var a=0; a<ids.length; a++){
          var na = B[ids[a]];
          for (var b=a+1; b<ids.length; b++){
            var nb = B[ids[b]];
            var dx = nb.x-na.x, dy = nb.y-na.y;
            var dist = Math.sqrt(dx*dx + dy*dy) || 0.01;
            var mind = R(ids[a]) + R(ids[b]) + 12;
            if (dist < mind){
              overlaps++;
              var pp = (mind-dist)/2 + 0.5, ux = dx/dist, uy = dy/dist;  // +0.5 : depasse le contact
              na.x -= ux*pp; na.y -= uy*pp;
              nb.x += ux*pp; nb.y += uy*pp;
            }
          }
        }
        if (overlaps === 0) break;
      }
      var upd = [];
      ids.forEach(function(id){ upd.push({id:id, x:B[id].x, y:B[id].y}); });
      network.body.data.nodes.update(upd);
      network.fit();
    }, 50);
  });
""".replace("__SIZES__", _sizes_js)

_html = net.generate_html()
for _pat in ["network = new vis.Network(container, data, options);",
             "var network = new vis.Network(container, data, options);"]:
    if _pat in _html:
        _html = _html.replace(_pat, _pat + _deoverlap_js)
        break

nx.write_graphml(G, "data/visualisations/patch_graphs/global_supercluster.graphml")
out = "data/visualisations/patch_graphs/global_supercluster_pyvis.html"
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w", encoding="utf-8") as _f:
    _f.write(_html)
print(f"Sauvegardé : {out}")
IFrame(out, width="100%", height=900)


# # Cluster_data_70 : patch par interaction

# ## Vue générale

# In[11]:


DF        = df_merged
col_patch = 'cluster_data_70'

patch_graphs     = build_patch_graphs(DF, col_patch, real_actin)
df_patches_infos = build_patch_info(DF, col_patch, patch_graphs, df_res_contacts=df_res)

# Binding site clusters (S1 + S2) inclus dans chaque patch C70
bs_cols = [c for c in ['s1_binding_site_cluster_data_70', 's2_binding_site_cluster_data_70'] if c in DF.columns]
if bs_cols:
    iid_to_bs = DF[['interaction_id'] + bs_cols].drop_duplicates('interaction_id').set_index('interaction_id')
    def _binding_sites(ids):
        sites = set()
        for col in bs_cols:
            sites |= set(iid_to_bs.loc[iid_to_bs.index.isin(ids), col].dropna().unique())
        return sorted(sites)
    df_patches_infos['binding_site_clusters'] = df_patches_infos['ids_interactions'].apply(_binding_sites)
    df_patches_infos['n_binding_site_clusters'] = df_patches_infos['binding_site_clusters'].apply(len)

print(f"Patchs : {len(patch_graphs)}")
print(df_patches_infos.head(10))
print(df_patches_infos['interaction_type'].value_counts())

plot_top_patches(DF, col_patch, patch_graphs, real_actin)


# In[12]:


save_patch_graphs(DF, col_patch, patch_graphs, real_actin, 'data/visualisations/patch_graphs/cluster_data_70')
df_patches_infos.to_csv('data/filtered/patches_infos_cluster_data_70.csv', index=False)
print(f"Sauvegardé : patches_infos_cluster_data_70.csv ({len(df_patches_infos)} patchs)")


# ## Aligner les sequences

# ## Export pour pipeline MAFFT / PyMOL
# 
# Exporte depuis `df_merged` :
# - **`data/exports/interactions_for_mafft.csv`** — table de correspondance complète (interaction_id, pdb_id, chaînes, séquences, clusters)
# - **`data/exports/fasta_s1/`** — FASTA S1 (actine) par groupe `s1_sequence_cluster_70` (1–2 fichiers)
# - **`data/exports/fasta_s2_patches/s1patch/`** — FASTA S2 par patch S1 × sous-groupe `s2_sequence_cluster_70`
# - **`data/exports/fasta_s2_patches/c70patch/`** — FASTA S2 par patch C70 × sous-groupe `s2_sequence_cluster_70`
# 
# Ces fichiers sont l'input de `mafft_alignment.ipynb`.

# In[13]:


import re
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ── Répertoires d'export ───────────────────────────────────────────────────
EXPORT_DIR = Path('data/exports')
FASTA_S1_DIR  = EXPORT_DIR / 'fasta_s1'
FASTA_S2_DIR  = EXPORT_DIR / 'fasta_s2_patches'

for d in [
    FASTA_S1_DIR,
    FASTA_S2_DIR / 's1patch',
    FASTA_S2_DIR / 'c70patch',
]:
    d.mkdir(parents=True, exist_ok=True)

# ── Table de correspondance complète ──────────────────────────────────────
cols_export = [
    'interaction_id', 'pdb_id', 'subunit_1', 'subunit_2',
    's1_sequence', 's2_sequence',
    's1_sequence_cluster_70', 's2_sequence_cluster_70',
    's1_binding_site_cluster_data_70', 'cluster_data_70'
]
df_export = df_merged[[c for c in cols_export if c in df_merged.columns]].copy()
df_export.to_csv(EXPORT_DIR / 'interactions_for_mafft.csv', index=False)
print(f"Table exportée : {len(df_export)} interactions → {EXPORT_DIR}/interactions_for_mafft.csv")

# ── Helper FASTA ───────────────────────────────────────────────────────────
def write_fasta_group(group_df, seq_col, id_cols, fasta_path):
    """Écrit un FASTA en dédupliquant par séquence."""
    records, seen = [], set()
    for _, row in group_df.iterrows():
        seq = str(row.get(seq_col, '')).strip()
        if seq and seq != 'nan' and seq not in seen:
            header = '_'.join(str(row[c]) for c in id_cols if c in row.index)
            header += f"_iid{int(row['interaction_id'])}"
            records.append(SeqRecord(Seq(seq), id=header, description=''))
            seen.add(seq)
    if len(records) >= 2:
        SeqIO.write(records, fasta_path, 'fasta')
        return len(records)
    return 0

# ── FASTA S1 par groupe s1_sequence_cluster_70 ────────────────────────────
n_s1 = 0
for cluster_id, group in df_export.groupby('s1_sequence_cluster_70', dropna=True):
    path = FASTA_S1_DIR / f"s1_cluster_{cluster_id}.fasta"
    n = write_fasta_group(group, 's1_sequence', ['pdb_id', 'subunit_1'], path)
    if n:
        n_s1 += 1
print(f"S1 FASTA : {n_s1} groupe(s) → {FASTA_S1_DIR}/")

# ── FASTA S2 par patch S1 × sous-groupe s2_sequence_cluster_70 ────────────
n_s2_s1 = 0
for patch in df_export['s1_binding_site_cluster_data_70'].dropna().unique():
    patch_grp = df_export[df_export['s1_binding_site_cluster_data_70'] == patch]
    for s2_cl, sub in patch_grp.groupby('s2_sequence_cluster_70', dropna=True):
        safe = re.sub(r'[^\w\-]', '_', str(patch))
        path = FASTA_S2_DIR / 's1patch' / f"patch_{safe}_s2cl_{s2_cl}.fasta"
        n = write_fasta_group(sub, 's2_sequence', ['pdb_id', 'subunit_2'], path)
        if n:
            n_s2_s1 += 1
print(f"S2 FASTA (patches S1)  : {n_s2_s1} fichiers → {FASTA_S2_DIR}/s1patch/")

# ── FASTA S2 par patch C70 × sous-groupe s2_sequence_cluster_70 ───────────
n_s2_c70 = 0
for patch in df_export['cluster_data_70'].dropna().unique():
    patch_grp = df_export[df_export['cluster_data_70'] == patch]
    for s2_cl, sub in patch_grp.groupby('s2_sequence_cluster_70', dropna=True):
        safe = re.sub(r'[^\w\-]', '_', str(patch))
        path = FASTA_S2_DIR / 'c70patch' / f"patch_{safe}_s2cl_{s2_cl}.fasta"
        n = write_fasta_group(sub, 's2_sequence', ['pdb_id', 'subunit_2'], path)
        if n:
            n_s2_c70 += 1
print(f"S2 FASTA (patches C70) : {n_s2_c70} fichiers → {FASTA_S2_DIR}/c70patch/")

