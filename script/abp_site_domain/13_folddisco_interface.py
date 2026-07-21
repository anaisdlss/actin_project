#!/usr/bin/env python
"""
Comparaison des MOTIFS d'interface actine des ABP avec **FoldDisco**.

Pourquoi : le TM-align sur les seuls résidus de contact (±0) s'effondre pour des
fragments courts et **discontinus** (~20 résidus éparpillés) — sa normalisation
par la longueur (d0 ≈ 0,3 Å) le rend absurde, et il superpose mal un nuage de
points. FoldDisco (Steinegger lab) est **fait pour les motifs discontinus** :
il encode des features géométriques position-indépendantes entre résidus et
score par rareté, sans superposition de chaîne.

Méthode :
  1. Indexer les chaînes d'ABP (data/exports/abp_site_domain/abp_chains/*.pdb).
  2. Pour chaque ABP, interroger SES résidus de contact (motif d'interface) contre
     l'index → score FoldDisco / nb de résidus appariés / RMSD vs chaque autre ABP.
  3. Agréger en une table paire-à-paire.

Sortie :
  data/exports/abp_site_domain/folddisco_interface.csv
    (query_abp, target_abp, n_match, score, rmsd, self)
"""
import re
import subprocess
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/exports/abp_site_domain"
CHAINS = OUT / "abp_chains"
INDEX = OUT / "folddisco_index" / "idx"
INDEX.parent.mkdir(parents=True, exist_ok=True)


def slug(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(s)).strip("_")[:50]


def main():
    rep = pd.read_csv(OUT / "abp_representatives.csv")
    res = pd.read_csv(ROOT / "data/filtered/details/3.interface_residues.csv")
    # stem <-> titre (nom de fichier des chaînes = {slug}__{pdb}_{chain})
    rep["stem"] = rep.apply(
        lambda r: f"{slug(r.abp_title)}__{r.pdb}_{r.abp_chain}", axis=1)
    stem2title = dict(zip(rep["stem"], rep["abp_title"]))

    # 1) indexer le dossier des chaînes ABP
    print(f"[1/3] index FoldDisco de {CHAINS} …")
    subprocess.run(
        ["folddisco", "index", "-p", str(CHAINS), "-i", str(INDEX),
         "-t", "4", "--id", "relpath"],
        check=True, capture_output=True, text=True)

    # 2) pour chaque ABP : motif d'interface -> query
    rows = []
    for _, r in rep.iterrows():
        pdb_chain = f"{r.pdb}_{r.abp_chain}"
        sub = res[(res.interaction_id == int(r.interaction_id))
                  & (res.chain == pdb_chain)]
        contacts = sorted(pd.to_numeric(
            sub.residue_number_structure, errors="coerce").dropna().astype(int))
        if len(contacts) < 3:
            continue
        query = ",".join(f"{r.abp_chain}{c}" for c in contacts)
        qpdb = CHAINS / f"{r.stem}.pdb"
        if not qpdb.exists():
            continue
        print(f"[2/3] query {r.abp_title} ({len(contacts)} résidus)…")
        out = subprocess.run(
            ["folddisco", "query", "-p", str(qpdb), "-q", query,
             "-i", str(INDEX), "-t", "4"],
            capture_output=True, text=True)
        for line in out.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            tgt_path, n_match, score, rmsd = parts[0], parts[1], parts[2], parts[3]
            tstem = Path(tgt_path).stem
            ttitle = stem2title.get(tstem, tstem)
            try:
                rows.append({
                    "query_abp": r.abp_title,
                    "target_abp": ttitle,
                    "query_size": len(contacts),   # taille réelle du motif requête
                    "n_match": int(n_match),
                    "score": float(score),
                    "rmsd": float(rmsd),
                    "self": r.abp_title == ttitle,
                })
            except ValueError:
                continue

    df = pd.DataFrame(rows)
    # meilleur hit par (query, target) : score max
    df = (df.sort_values("score", ascending=False)
          .drop_duplicates(["query_abp", "target_abp"]))
    # normaliser le score par le self-score de la requête (self = 1.0)
    self_score = (df[df["self"]].set_index("query_abp")["score"].to_dict())
    df["score_norm"] = df.apply(
        lambda x: round(x["score"] / self_score.get(x["query_abp"], x["score"]), 3)
        if self_score.get(x["query_abp"]) else None, axis=1)
    out_csv = OUT / "folddisco_interface.csv"
    df.to_csv(out_csv, index=False)
    print(f"\n[3/3] écrit : {out_csv}  ({len(df)} paires query-target)")


if __name__ == "__main__":
    main()
