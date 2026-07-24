#!/usr/bin/env python
"""Produit une copie SLIM des données pour un déploiement Streamlit Cloud.

Ne garde que ce que l'app AFFICHE : les CSV, les PNG de visualisation, les
assembly PDB/CIF (pour la 3D) et les `4.query_ProteoCast.csv`. On jette tout le
lourd non affiché (pairwise, filament, bfactor, pymol, MSA/structures ProteoCast,
sessions .pse, gros .a3m/.fasta).

Sortie : `deploy/data/` (arbre slim prêt à committer sur la branche `deploy`),
plus un marqueur `deploy/data/.slim_deploy` qui met l'app en mode déployé (cache
le pipeline « Run / update » et le calcul ProteoCast, impossibles sur Cloud).

Usage :  pixi run python script/make_slim_deploy.py
"""
import os
import shutil
from pathlib import Path

SRC = Path("data")
DST = Path("deploy/data")

# Dossiers entiers à EXCLURE (lourds, non affichés par l'app).
# On GARDE bfactor_cluster + bfactor_c70_interface : ils fournissent les PDB de
# la vue « Interface 3D — actin ↔ ABP » (S1). On jette filament/pymol/raw
# (aucun viewer de l'app ne les lit) et pairwise (reconstruit depuis l'assembly).
EXCLUDE_DIRS = {
    "data/filtered/details/structures_files/pairwise",
    "data/filtered/details/structures_files/filament",
    "data/filtered/details/structures_files/pymol",
    "data/raw",
}
# Extensions lourdes à exclure partout.
EXCLUDE_EXT = {".pse", ".a3m"}
# Taille max pour un .fasta (les MSA ProteoCast sont énormes).
FASTA_MAX = 3 * 1024 * 1024


def _keep(path: Path) -> bool:
    p = str(path).replace(os.sep, "/")
    for d in EXCLUDE_DIRS:
        if p.startswith(d + "/") or p == d:
            return False
    # ProteoCast : garder les CSV (conservation/substitution) + les PDB (3D) ;
    # jeter MSA (.a3m/.fasta), images (.jpg/.png), et le dossier GMM/.
    if "/proteocast/abp/" in p:
        if "/GMM/" in p:
            return False
        if path.is_file() and path.suffix.lower() not in (".csv", ".pdb"):
            return False
    if path.suffix.lower() in EXCLUDE_EXT:
        return False
    if path.suffix.lower() == ".fasta" and path.is_file() \
            and path.stat().st_size > FASTA_MAX:
        return False
    return True


def main():
    if DST.exists():
        shutil.rmtree(DST)
    DST.mkdir(parents=True)
    n_files = 0
    total = 0
    for root, dirs, files in os.walk(SRC):
        rootp = Path(root)
        # élaguer les dossiers exclus
        dirs[:] = [d for d in dirs if _keep(rootp / d)]
        for f in files:
            src = rootp / f
            if not _keep(src):
                continue
            rel = src.relative_to(SRC)
            dst = DST / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            n_files += 1
            total += src.stat().st_size
    # marqueur mode déployé
    (DST / ".slim_deploy").write_text("1\n")

    # ── Copier le CODE et écrire requirements.txt → deploy/ = app autonome ────
    deploy = DST.parent
    if (deploy / "script").exists():
        shutil.rmtree(deploy / "script")
    shutil.copytree("script", deploy / "script")
    (deploy / "requirements.txt").write_text(REQUIREMENTS)
    (deploy / ".gitignore").write_text("__pycache__/\n*.pyc\n")
    (deploy / "README.md").write_text(DEPLOY_README)

    print(f"deploy/data/ : {n_files} fichiers, {total/1e6:.0f} Mo")
    print("deploy/ = app autonome (code + data slim + requirements.txt).")
    print("Marqueur .slim_deploy créé (l'app masque pipeline + calcul ProteoCast).")


# requirements PIP pour Streamlit Cloud (affichage seulement, pas de pipeline :
# on n'inclut ni mafft/foldseek/pymol/tmtools ni jupyter).
REQUIREMENTS = """\
streamlit>=1.55,<2
pandas>=2,<3
numpy>=1.26,<3
plotly>=5.24,<6
networkx>=3.6,<4
pyvis>=0.3.2,<1
matplotlib>=3.10,<4
scipy>=1.15,<2
biopython>=1.87,<2
py3Dmol
seaborn>=0.13.2,<0.14
requests
openpyxl
lxml
beautifulsoup4
adjustText>=1.3,<2
streamlit-agraph>=0.0.45
"""

DEPLOY_README = """\
# Actin–ABP interaction analysis — shared (slim) build

Read-only Streamlit build (slim data pre-bundled). Deployed on Streamlit
Community Cloud. The data pipeline and ProteoCast computation are disabled here
(they run only on the full local project).

Run locally: `streamlit run script/streamlit.py`
"""


if __name__ == "__main__":
    main()
