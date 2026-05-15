# Projet interaction actine-ABP

Récupération et analyse des interactions actine–protéines de liaison à l'actine (ABP) depuis la base de données [PPI3D](https://bioinformatics.lt/ppi3d), à partir de l'identifiant UniProt **P60709** (Actin, cytoplasmic 1).

---

## 1. Installation

Assurez-vous que [pixi](https://pixi.sh) est installé.

```bash
git clone https://github.com/anaisdlss/actin_project.git
cd actin_project
pixi install
pixi shell
```

---

## 2. Pipeline complet — générer toutes les données

Le pipeline se déroule en **deux phases** à exécuter dans l'ordre.

### Phase 1 — Téléchargement et filtrage des données PPI3D

```bash
python -m script.data_extract.pipeline_data
```

> Sur macOS, pour éviter la mise en veille :
> ```bash
> caffeinate -i python -m script.data_extract.pipeline_data
> ```

Cette commande enchaîne automatiquement les étapes suivantes :

| Étape | Description |
|-------|-------------|
| 1/11 | Récupération du summary des interactions (BLAST PPI3D) |
| 2/11 | Récupération des entrées PDB associées |
| 3/11 | Récupération de la table globale des clusters PPI3D |
| 4/11 | Filtrage des structures (≥ 5 actines) via notebook |
| 5/11 | Téléchargement des détails d'interface (résidus de contact) |
| 6/11 | Alignements MAFFT des séquences |
| 7/11 | Calcul des patches S1 binding site |
| 8/11 | Calcul des B-factors C70 interface |
| 9/11 | Calcul des B-factors cluster |
| 10/11 | Génération du heatmap S1 |
| 11/11 | Calcul des B-factors cluster (finalisation) |

Les étapes déjà à jour sont automatiquement ignorées.

### Phase 2 — Analyse des interfaces (notebooks)

```bash
jupyter nbconvert --to notebook --execute --inplace notebooks/interface_analysis_s1.ipynb
jupyter nbconvert --to notebook --execute --inplace notebooks/interface_analysis_c70.ipynb
```

Ces notebooks génèrent les visualisations dans `visualisations/`.

---

## 3. Lancer l'interface Streamlit

```bash
streamlit run script/streamlit.py
```

L'application permet de visualiser les données filtrées et d'explorer les clusters d'interactions actine–ABP.

---

## Structure des données générées

```
data/
├── raw/
│   ├── ppi3d_actin_summary.csv
│   ├── pdb_entry_results.csv
│   ├── all_data.csv
│   └── metadata.json
└── filtered/
    ├── filtered_pdb_entry.csv
    ├── filtered_summary.csv
    ├── filtered_all_data.csv
    ├── patches_infos_s1_binding_site.csv
    ├── patches_infos_cluster_data_70.csv
    ├── actin_s1_canon_area_by_cluster.csv
    └── details/
        ├── 1.interactions.csv
        ├── 2.proteins.csv
        ├── 3.interface_residues.csv
        ├── 4.inter-residue_contacts.csv
        ├── 5.ligands.csv
        ├── 6.meta_alignement.csv
        ├── 7.alignment_sequences.csv
        └── 8.structures.csv

visualisations/
├── actin_s1_all_equitable_heatmap.png
├── actin_c70_heatmap_surface_area.png
├── actin_c70_contacts/
└── actin_c70_contacts_surface_area/
```

> `metadata.json` stocke les identifiants de jobs PPI3D pour la reproductibilité. Relancer le pipeline si la base PPI3D est mise à jour.
