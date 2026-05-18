# Projet interaction actine–ABP

Analyse automatisée des interactions entre l'actine et ses partenaires protéiques (ABP) à partir de la base de données [PPI3D](https://bioinformatics.lt/ppi3d), identifiant UniProt **P60709** (Actin, cytoplasmic 1).

Le pipeline télécharge les données, filtre les structures, aligne les séquences, regroupe les interfaces en clusters et génère des visualisations PyMOL et une interface web Streamlit.

---

## Prérequis

- Accès à Internet
- [PyMOL](https://pymol.org) installé (pour les visualisations 3D)

Tout le reste (Python, dépendances) s'installe automatiquement via **pixi**.

---

## Installation

### 1. Installer pixi

**macOS / Linux** — ouvre un terminal et colle :

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

**Windows** — ouvre **PowerShell** et colle :

```powershell
iwr -useb https://pixi.sh/install.ps1 | iex
```

**Ferme et rouvre le terminal** avant de continuer.

### 2. Récupérer le projet

```bash
cd ~/Desktop
git clone https://github.com/anaisdlss/actin_project.git
cd actin_project
```

> macOS : si on te propose d'installer les outils de développement, accepte.  
> Windows : utilise Git Bash ou PowerShell avec [Git for Windows](https://git-scm.com/download/win).

### 3. Installer l'environnement

```bash
pixi install
```

Quelques minutes la première fois.

---

## Générer les données

**macOS** :
```bash
caffeinate -i pixi run python -m script.data_extract.pipeline_data
```

**Linux / Windows** :
```bash
pixi run python -m script.data_extract.pipeline_data
```

> Sur macOS, `caffeinate` empêche la mise en veille pendant le calcul.  
> Durée : **30 à 60 minutes** selon la connexion.

Le pipeline reprend automatiquement là où il s'est arrêté si tu le relances.

| Étape | Description |
|-------|-------------|
| 1/13 | Téléchargement du résumé PPI3D |
| 2/13 | Téléchargement des entrées PDB |
| 3/13 | Téléchargement de la table des clusters |
| 4/13 | Filtrage des structures (≥ 4 actines connectées) |
| 5/13 | Téléchargement des détails d'interface |
| 6/13 | Alignement MAFFT des séquences |
| 7/13 | Analyse des clusters C70 |
| 8/13 | Calcul des B-factors d'interface C70 |
| 9/13 | Génération du script PyMOL vue globale |
| 10/13 | Génération des scripts PyMOL par site de liaison |
| 11/13 | Analyse détaillée de l'interface C70 |
| 12/13 | Heatmap S1 binding site |
| 13/13 | Calcul des B-factors S1 par cluster |

---

## Lancer l'interface web

```bash
pixi run streamlit run script/streamlit.py
```

Une page s'ouvre dans le navigateur (sinon, va sur `http://localhost:8501`).

---

## Visualiser dans PyMOL

Les scripts sont générés dans :

```
data/filtered/details/structures_files/bfactor_c70_interface/
```

**Vue globale** — tous les clusters C70 superposés sur l'actine de référence :

```
view_full_surface.pml
```

**Vue par site de liaison** — un fichier par cluster de site actine (nœud rouge dans le réseau Streamlit) :

```
by_s1_cluster/<nom_du_cluster>.pml
```

Pour ouvrir un script : dans PyMOL, `File > Run Script…` ou dans la console PyMOL :

```
@chemin/vers/le/fichier.pml
```

**Code couleur :** vert = ABP (hétéro) · rose = actine partenaire (homo) · intensité = % surface enfouie à l'interface
