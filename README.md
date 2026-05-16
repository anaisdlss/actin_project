# Projet interaction actine–ABP

Analyse automatique des interactions entre l'actine et les protéines de liaison à l'actine (ABP) à partir de la base de données [PPI3D](https://bioinformatics.lt/ppi3d), en utilisant l'identifiant UniProt **P60709** (Actin, cytoplasmic 1).

---

## Ce que fait ce projet

1. **Télécharge** toutes les structures 3D d'interactions actine–ABP disponibles dans PPI3D
2. **Filtre** les structures pertinentes (assemblages avec au moins 4 actines connectées)
3. **Aligne** les séquences pour comparer les résidus d'interface entre structures
4. **Regroupe** les interactions en clusters selon leur similarité d'interface (C70)
5. **Calcule** des B-factors représentant le pourcentage de surface enfouie à l'interface
6. **Génère** des visualisations interactives (interface web Streamlit) et des sessions PyMOL

---

## Prérequis

Avoir un Mac avec accès à Internet. Rien d'autre n'est nécessaire — tout sera installé automatiquement.

---

## Étape 1 — Installer pixi

Pixi est l'outil qui installe Python et toutes les dépendances du projet automatiquement.

Ouvre le **Terminal** (appuie sur `Cmd + Espace`, tape "Terminal", appuie sur Entrée), puis colle cette ligne et appuie sur Entrée :

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

**Ferme le terminal et rouvre-le** avant de continuer.

---

## Étape 2 — Récupérer le projet

Dans le terminal, choisis où mettre le projet (par exemple ton Bureau) :

```bash
cd ~/Desktop
git clone https://github.com/anaisdlss/actin_project.git
cd actin_project
```

> Si macOS te demande d'installer les outils de développement, clique **Installer**.

---

## Étape 3 — Installer l'environnement

```bash
pixi install
```

Cette commande télécharge Python, Jupyter et toutes les librairies. Cela prend quelques minutes la première fois.

---

## Étape 4 — Lancer le pipeline complet

Le pipeline télécharge les données et génère toutes les analyses en une seule commande.

```bash
caffeinate -i pixi run python -m script.data_extract.pipeline_data
```

> `caffeinate` empêche l'ordinateur de se mettre en veille pendant le calcul.  
> Le pipeline peut durer **30 à 60 minutes** selon la connexion Internet.

### Étapes du pipeline

| Étape | Description |
|-------|-------------|
| 1/12 | Téléchargement du résumé des interactions (PPI3D BLAST) |
| 2/12 | Téléchargement des entrées PDB pour chaque structure |
| 3/12 | Téléchargement de la table complète des clusters PPI3D |
| 4/12 | Filtrage des structures (≥ 4 actines connectées) |
| 5/12 | Téléchargement des détails d'interface résidu par résidu |
| 6/12 | Alignement MAFFT des séquences par cluster |
| 7/12 | Analyse des clusters d'interaction C70 |
| 8/12 | Calcul des B-factors d'interface C70 par cluster |
| 9/12 | Génération des scripts PyMOL (vue globale + vue par site) |
| 10/12 | Analyse détaillée de l'interface par cluster C70 |
| 11/12 | Heatmap S1 binding site et références clusters |
| 12/12 | Calcul des B-factors S1 par cluster |

Les étapes déjà réalisées sont automatiquement ignorées si les données sont à jour.

---

## Étape 5 — Lancer l'interface web

```bash
pixi run streamlit run script/streamlit.py
```

Une page s'ouvre dans ton navigateur. Si ce n'est pas le cas, ouvre manuellement l'adresse affichée dans le terminal (généralement `http://localhost:8501`).

---

## Étape 6 — Visualiser les structures 3D dans PyMOL

Les scripts PyMOL sont générés automatiquement dans :

```
data/filtered/details/structures_files/bfactor_c70_interface/
```

### Vue globale — tous les clusters

Ouvre **`view_full_surface.pml`** dans PyMOL pour voir en une seule session :
- L'actine de référence (8iah_L, surface grise semi-transparente)
- Les 255 clusters C70, chacun représenté par la surface de son partenaire colorée par B-factor :
  - **Vert → blanc** : partenaires ABP (interactions hétéro)
  - **Rose → blanc** : actines partenaires (interactions homo)
  - L'intensité de couleur = pourcentage de surface enfouie à l'interface

### Vue par site de liaison actine

Le dossier **`by_s1_cluster/`** contient 172 scripts PyMOL, un par cluster de site de liaison actine (nœud rouge dans le réseau Streamlit).

Chaque fichier (ex : `6685_0.pml`) montre :
- L'actine de référence
- Un objet par partenaire connecté à ce site dans le réseau (une surface par arête)

**Comment ouvrir un script PyMOL :**
1. Ouvre PyMOL
2. Fais `File > Run Script…` et sélectionne le fichier `.pml`  
   ou dans la console PyMOL : `@chemin/vers/le/fichier.pml`

---

## Structure des fichiers générés

Après le pipeline, les dossiers principaux sont :

```
data/
├── raw/                            ← données brutes PPI3D
└── filtered/                       ← données filtrées et analysées
    ├── filtered_all_data.csv           données complètes filtrées
    ├── patches_infos_cluster_data_70.csv  clusters C70
    ├── patches_infos_s1_binding_site.csv  clusters binding site S1
    ├── actin_s1_canon_area_by_cluster.csv profils d'interface par cluster
    ├── s1_cluster_reference.csv        actines de référence par cluster
    └── details/
        ├── 1.interactions.csv
        ├── 3.interface_residues.csv
        ├── 4.inter-residue_contacts.csv
        └── structures_files/
            ├── pairwise/               PDB des paires d'interactions
            └── bfactor_c70_interface/  PDB avec B-factors + scripts PyMOL
                ├── view_full_surface.pml   vue globale (255 clusters)
                └── by_s1_cluster/          172 vues par site de liaison actine

visualisations/                     ← graphiques et heatmaps générés
```
