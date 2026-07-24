# Déployer la version partagée (Streamlit Community Cloud)

Version **slim** de l'app (données réduites embarquées, ~600 Mo avec la 3D).
Le pipeline « Run / update » et le calcul ProteoCast sont **désactivés** sur cette
version (impossibles sur Cloud) — l'app détecte le marqueur `data/.slim_deploy`.

## 1. Générer le build slim

Avec les données complètes présentes en local :

```bash
pixi run python script/make_slim_deploy.py
```

→ crée le dossier **`deploy/`** = app autonome :
`deploy/script/`, `deploy/data/` (slim), `deploy/requirements.txt`,
`deploy/.slim_deploy` (dans data/). Ton dossier `data/` complet n'est **pas**
touché.

## 2. Publier sur GitHub (une seule fois)

Crée un dépôt GitHub vide (ex. `actin-abp-app`), puis :

```bash
cd deploy
git init && git add -A && git commit -m "slim build"
git branch -M main
git remote add origin https://github.com/<toi>/actin-abp-app.git
git push -u origin main
```

Le dépôt fait ~600 Mo (assembly 3D embarqué) — GitHub l'accepte (aucun fichier
> 100 Mo).

## 3. Déployer sur Streamlit Cloud

1. Va sur **share.streamlit.io** → **New app**.
2. Repository : `actin-abp-app` · Branch : `main` · Main file : `script/streamlit.py`.
3. **Deploy**. (Streamlit installe `requirements.txt` puis lance l'app.)

## 4. Mettre à jour

Oui, c'est **à toi** de mettre à jour, mais c'est simple : Streamlit Cloud
**redéploie tout seul à chaque push**.

- **Code seulement** : modifie sur le projet principal, relance
  `make_slim_deploy.py`, puis dans `deploy/` : `git add -A && git commit -m "maj" && git push`.
- **Données** : régénère en local (le pipeline ne tourne pas sur Cloud) →
  `make_slim_deploy.py` → commit + push dans `deploy/`.

L'app se recharge automatiquement en ~1 min.

## Notes

- Si l'app échoue à l'import sur Cloud (dépendance manquante), ajoute-la à
  `deploy/requirements.txt` et push.
- Pour alléger davantage (< 300 Mo) : on peut retirer l'assembly 3D et charger
  les structures à la demande depuis RCSB (à faire ensemble si besoin).
