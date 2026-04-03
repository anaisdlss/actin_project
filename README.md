# Projet interaction actine-ABP

Récupération et analyse des interactions actine–protéines de liaison à l'actine (ABP) depuis la base de données [PPI3D](https://bioinformatics.lt/ppi3d), à partir de l'identifiant UniProt **P60709** (Actin, cytoplasmic 1).

## 1. Installation

Assurez-vous que [pixi](https://pixi.sh) et Python ≥ 3.11 sont installés.

```bash
git clone https://github.com/anaisdlss/actin_project.git
cd actin_project
pixi install
pixi shell
```

## 2. Générer les données

### Option A — Pipeline complet (recommandé)

Le pipeline enchaîne les 6 étapes automatiquement :

```bash
python script/data_extract/pipeline_data.py
```

> Sur macOS, pour éviter la mise en veille pendant l'exécution :
> ```bash
> caffeinate -i python script/data_extract/pipeline_data.py
> ```

Les étapes sont :

| Étape | Description |
|-------|-------------|
| 1/6 | Récupération du summary des interactions (BLAST PPI3D) |
| 2/6 | Récupération des entrées PDB associées |
| 3/6 | Récupération de la table globale des clusters PPI3D |
| 4/6 | Filtrage des structures (≥ 5 actines) via notebook |
| 5/6 | Téléchargement des détails d'interface (résidus de contact) |
| 6/6 | Filtrage de la table all_data selon les PDB retenus |

Les étapes déjà à jour sont automatiquement ignorées (vérification par date de modification).

### Option B — Étapes individuelles

```bash
python -m script.data_extract.get_summary_results
python -m script.data_extract.get_pdb_entries
python -m script.data_extract.get_cluster_table
python -m script.data_extract.get_interaction_details
```

### Structure des données générées

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
    └── details/
        ├── 1.interactions.csv
        ├── 2.proteins.csv
        ├── 3.interface_residues.csv
        ├── 4.inter-residue_contacts.csv
        ├── 5.ligands.csv
        ├── 6.meta_alignement.csv
        ├── 7.alignment_sequences.csv
        └── 8.structures.csv
```

> Le fichier `metadata.json` stocke les identifiants de jobs PPI3D pour assurer la reproductibilité. Quand la base PPI3D est mise à jour, relancer le pipeline pour obtenir les données les plus récentes.

## 3. Interface Streamlit

L'application Streamlit permet de lancer le pipeline, visualiser les données filtrées et explorer les clusters d'interactions depuis un navigateur.

```bash
streamlit run script/streamlit.py
```

## 4. Notebooks d'analyse

Les notebooks sont dans le dossier `notebooks/` :

- `graphe_filter.ipynb` — filtrage des structures PDB (étape 4 du pipeline)
- `cluster_interaction_analysis.ipynb` — analyse des clusters d'interactions actine
