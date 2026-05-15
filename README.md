# Projet interaction actine-ABP

Récupération et analyse des interactions actine–protéines de liaison à l'actine (ABP) depuis la base de données [PPI3D](https://bioinformatics.lt/ppi3d), à partir de l'identifiant UniProt **P60709** (Actin, cytoplasmic 1).

---

## Étape 1 — Installer les outils nécessaires

### 1.1 Installer pixi

Pixi est le gestionnaire d'environnement utilisé par ce projet. Il installe automatiquement Python et toutes les dépendances.

**Sur macOS :**
Ouvre le terminal (cherche "Terminal" dans Spotlight avec `Cmd + Espace`), puis colle cette commande et appuie sur Entrée :

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

Ferme le terminal et rouvre-le pour que pixi soit pris en compte.

---

## Étape 2 — Récupérer le projet

Choisis un dossier où tu veux mettre le projet (par exemple ton Bureau ou tes Documents), puis dans le terminal :

```bash
# Exemple : aller dans le dossier Documents
cd ~/Documents

# Cloner le projet depuis GitHub
git clone https://github.com/anaisdlss/actin_project.git

# Entrer dans le dossier du projet
cd actin_project
```

> Si `git` n'est pas installé, macOS te proposera de l'installer automatiquement au premier usage — clique OK.

---

## Étape 3 — Installer l'environnement

Toujours dans le terminal, dans le dossier `actin_project` :

```bash
pixi install
```

Cette commande télécharge et installe Python, Jupyter, et toutes les librairies nécessaires. Cela peut prendre quelques minutes la première fois.

---

## Étape 4 — Lancer le pipeline complet

Le pipeline télécharge toutes les données depuis PPI3D et génère l'ensemble des analyses automatiquement.

```bash
pixi run python -m script.data_extract.pipeline_data
```

> **Important :** ne mets pas ton ordinateur en veille pendant l'exécution. Sur Mac, tu peux éviter la mise en veille automatique avec :
> ```bash
> caffeinate -i pixi run python -m script.data_extract.pipeline_data
> ```

Le pipeline enchaîne les étapes suivantes (peut durer 30–60 minutes selon la connexion) :

| Étape | Description |
|-------|-------------|
| 1/11 | Téléchargement du summary des interactions (PPI3D) |
| 2/11 | Téléchargement des entrées PDB |
| 3/11 | Téléchargement de la table des clusters |
| 4/11 | Filtrage des structures (≥ 4 actines connectées) |
| 5/11 | Téléchargement des détails d'interface |
| 6/11 | Alignements MAFFT des séquences |
| 7/11 | Analyse des clusters d'interaction C70 |
| 8/11 | Calcul des B-factors interface C70 |
| 9/11 | Analyse détaillée interface par cluster C70 |
| 10/11 | Heatmap S1 binding site |
| 11/11 | Calcul des B-factors S1 par cluster |

Les étapes déjà réalisées sont automatiquement ignorées si les données sont à jour.

---

## Étape 5 — Lancer l'interface web

Une fois le pipeline terminé, lance l'interface Streamlit :

```bash
pixi run streamlit run script/streamlit.py
```

Une page s'ouvre automatiquement dans ton navigateur. Si ce n'est pas le cas, ouvre manuellement l'adresse affichée dans le terminal (généralement `http://localhost:8501`).

---

## Structure des données générées

Après le pipeline, le dossier `data/` contiendra :

```
data/
├── raw/                        ← données brutes téléchargées depuis PPI3D
│   ├── ppi3d_actin_summary.csv
│   ├── pdb_entry_results.csv
│   ├── all_data.csv
│   └── metadata.json
└── filtered/                   ← données filtrées et analysées
    ├── filtered_pdb_entry.csv
    ├── filtered_summary.csv
    ├── filtered_all_data.csv
    ├── patches_infos_s1_binding_site.csv
    ├── patches_infos_cluster_data_70.csv
    └── details/
        ├── 1.interactions.csv
        ├── 2.proteins.csv
        ├── 3.interface_residues.csv
        ├── 4.inter-residue_contacts.csv
        ├── 5.ligands.csv
        ├── 6.meta_alignement.csv
        ├── 7.alignment_sequences.csv
        └── 8.structures.csv

visualisations/                 ← graphiques générés par le pipeline
```
