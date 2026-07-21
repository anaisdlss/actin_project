#!/usr/bin/env python
"""Génère le manifest ProteoCast des ABP à partir de abp_master.csv (étape 9).

Un ABP = un titre de protéine. On récupère son accession UniProt (déjà résolue
dans abp_master), d'où ProteoCast dérive tout (proteocast.ijm.fr va chercher la
séquence + le modèle AlphaFold lui-même). Colonnes gardées compatibles avec
proteocast_view.load_status : abp_title, slug, uniprot, alphafold_id, longueur,
fasta_entree, sortie_attendue, fusion.
"""
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MASTER = ROOT / "data/exports/abp_site_domain/abp_master.csv"
OUT = ROOT / "data/proteocast/abp_inputs/manifest.csv"


def slug(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(s)).strip("_")[:50]


def main():
    if not MASTER.exists():
        raise FileNotFoundError(
            f"{MASTER} manquant (l'étape 9 doit tourner avant).")
    m = pd.read_csv(MASTER)

    rows = []
    for title, g in m.groupby("abp_title", sort=False):
        uni = g["uniprot"].dropna().astype(str)
        uni = uni.iloc[0].strip() if len(uni) else ""
        sl = slug(title)
        rows.append({
            "abp_title": title,
            "slug": sl,
            "uniprot": uni or None,
            "alphafold_id": f"AF-{uni}-F1" if uni else "",
            "longueur": "",
            "fasta_entree": f"data/proteocast/abp_inputs/{sl}.fasta",
            "sortie_attendue": f"data/proteocast/abp/{sl}.csv",
            "fusion": "," in str(title),
        })

    df = pd.DataFrame(rows).drop_duplicates(subset=["slug"])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    n_uni = int(df["uniprot"].notna().sum())
    print(f"Manifest ProteoCast : {len(df)} ABP ({n_uni} avec UniProt) -> {OUT}")


if __name__ == "__main__":
    main()
