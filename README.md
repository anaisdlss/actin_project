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

### 0 — Ouvrir un terminal

Toutes les commandes ci-dessous s'exécutent dans un terminal :

- **macOS** : `Cmd + Espace` → taper `Terminal` → Entrée
- **Linux** : `Ctrl + Alt + T`
- **Windows** : `Windows + R` → taper `powershell` → Entrée

### 1 — Installer pixi

Copier-coller la commande correspondante dans le terminal, puis appuyer sur Entrée.

**macOS / Linux**

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

**Windows** (PowerShell)

```powershell
iwr -useb https://pixi.sh/install.ps1 | iex
```

**Fermer et rouvrir le terminal** avant de continuer.

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

### 4 — Lancer l'interface web

```bash
pixi run streamlit run script/streamlit.py
```

Une page s'ouvre automatiquement dans le navigateur (sinon, aller sur `http://localhost:8501`).

---

## Générer les données

Une fois l'interface ouverte, se rendre dans la section **Téléchargement des données** et cliquer sur **Lancer le téléchargement**.

Le pipeline s'exécute entièrement depuis l'interface — aucune commande supplémentaire n'est nécessaire. La durée est de **30 à 60 minutes** selon la vitesse de connexion. Les étapes déjà réalisées sont automatiquement ignorées lors des exécutions suivantes.

---

## Visualisations PyMOL

Les scripts sont générés automatiquement par le pipeline dans :

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
