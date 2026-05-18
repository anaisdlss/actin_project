# Pipeline d'analyse des interactions actine–ABP

Analyse structurale automatisée des interactions entre l'actine et ses protéines de liaison (ABP) extraites de la base de données [PPI3D](https://bioinformatics.lt/ppi3d) (entrée UniProt **P60709**, Actin cytoplasmic 1).

Le pipeline récupère l'ensemble des co-structures 3D disponibles, filtre les assemblages contenant au moins quatre sous-unités d'actine connectées, regroupe les résidus d'interface à 70 % de similarité (C70), calcule la surface enfouie à l'interface (encodée en B-factor), et produit des visualisations interactives via une interface web Streamlit et des scripts de session PyMOL.

---

## Prérequis

| Outil | Rôle |
|-------|------|
| [pixi](https://pixi.sh) | Gestion de l'environnement et des dépendances (installé ci-dessous) |
| [PyMOL](https://pymol.org) | Visualisation des structures 3D |
| Accès Internet | Récupération des données PPI3D et RCSB PDB |

Python et toutes les librairies Python sont gérés automatiquement par pixi — aucune installation manuelle n'est requise.

---

## Installation

### 1 — Installer pixi

**macOS / Linux**

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

**Windows** (PowerShell)

```powershell
iwr -useb https://pixi.sh/install.ps1 | iex
```

Fermer et rouvrir le terminal avant de continuer.

### 2 — Cloner le dépôt

```bash
git clone https://github.com/anaisdlss/actin_project.git
cd actin_project
```

> Windows : utiliser Git Bash ou PowerShell avec [Git for Windows](https://git-scm.com/download/win).

### 3 — Installer l'environnement

```bash
pixi install
```

Quelques minutes sont nécessaires lors de la première exécution.

---

## Exécution du pipeline

**macOS**

```bash
caffeinate -i pixi run python -m script.data_extract.pipeline_data
```

**Linux / Windows**

```bash
pixi run python -m script.data_extract.pipeline_data
```

> `caffeinate` (macOS uniquement) empêche la mise en veille du système pendant l'exécution.  
> Durée estimée : **30 à 60 minutes** selon la vitesse de connexion.

Le pipeline est incrémental : les étapes dont les résultats sont déjà à jour sont automatiquement ignorées lors des exécutions suivantes.

| Étape | Description |
|-------|-------------|
| 1/13 | Téléchargement du résumé PPI3D |
| 2/13 | Téléchargement des métadonnées PDB |
| 3/13 | Téléchargement de la table des clusters |
| 4/13 | Filtrage des structures (≥ 4 sous-unités d'actine connectées) |
| 5/13 | Téléchargement des détails d'interface résidu par résidu |
| 6/13 | Alignement de séquences MAFFT par cluster |
| 7/13 | Analyse des clusters d'interaction C70 |
| 8/13 | Calcul des B-factors d'interface C70 |
| 9/13 | Génération du script PyMOL vue globale |
| 10/13 | Génération des scripts PyMOL par site de liaison |
| 11/13 | Analyse détaillée de l'interface C70 |
| 12/13 | Heatmap du site de liaison S1 |
| 13/13 | Calcul des B-factors par cluster S1 |

---

## Interface web

```bash
pixi run streamlit run script/streamlit.py
```

L'interface est accessible à l'adresse `http://localhost:8501`. Elle propose des réseaux bipartites interactifs, des statistiques par cluster, des alignements de séquences et des visualisations 3D des B-factors.

---

## Visualisations PyMOL

Les scripts sont générés dans :

```
data/filtered/details/structures_files/bfactor_c70_interface/
```

**Vue globale** — tous les représentants de clusters C70 superposés sur l'actine de référence (8iah chaîne L) :

```
view_full_surface.pml
```

**Vue par site de liaison** — un script par cluster de site de liaison S1 (nœuds rouges dans le réseau Streamlit) :

```
by_s1_cluster/<nom_du_cluster>.pml
```

Pour exécuter un script, utiliser `File > Run Script…` dans PyMOL, ou taper dans la console PyMOL :

```
@/chemin/vers/le/script.pml
```

**Code couleur :** dégradé vert = partenaire ABP (interaction hétérologue) · dégradé rose = partenaire actine (interaction homologue) · l'intensité encode le pourcentage de surface accessible enfouie à l'interface.
