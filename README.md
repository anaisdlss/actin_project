# Projet interaction entre actin et actin binding protein

## Installation du projet et environnement

Assurez vous que pixi et python >3.11 (anaconda) soient installés.

Clonez le repertoire, allez au dossier du projet, et activez l'environnement pixi:

```bash
git clone https://github.com/anaisdlss/actin_project.git
cd actin_project
pixi install
pixi shell
```
## Data générées

Les données ont été recupéré à partir du script `script/ppi3d_query.py` et sont stockées dans `data`: `ppi3d_{NAME_PROT}_results.csv` est le resultat du BLAST de la protiene d'entrée (ici l'actine) et `ppi3d_actin_details_dataset.csv` est le resultats en détails des interactions qui en sont issues.
