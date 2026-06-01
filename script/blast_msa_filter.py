#!/usr/bin/env python3
"""
Recalcule e-value / bit-score / %identité pour chaque séquence d'un MSA .a3m
en lançant blastp (query = P60709) contre les séquences du MSA, puis (option)
écrit un .a3m filtré selon un seuil.

Usage :
  # 1) Juste la table (e-value, bitscore, identité) + distribution :
  python script/blast_msa_filter.py

  # 2) Écrire un a3m filtré (garde identité >= 50%) :
  python script/blast_msa_filter.py --min-pident 50
  # ou par e-value :
  python script/blast_msa_filter.py --max-evalue 1e-5
  # ou par bit-score :
  python script/blast_msa_filter.py --min-bitscore 200

L'e-value dépend de la taille de la base. Par défaut elle est calculée contre
les séquences du MSA. Utilise --dbsize 150000000 pour imiter UniRef90.
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
A3M_DEFAULT = ROOT / "data" / "AF-P60709-F1-msa_v6.a3m"
QUERY_DEFAULT = ROOT / "data" / "P60709_ref.fasta"
OUTDIR = ROOT / "data" / "blast"


def parse_a3m(path):
    """Retourne (heads, seqs). seqs[0] = query."""
    heads, seqs = [], []
    h, buf = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if h is not None:
                    heads.append(h)
                    seqs.append("".join(buf))
                h, buf = line[1:], []
            else:
                buf.append(line)
    if h is not None:
        heads.append(h)
        seqs.append("".join(buf))
    return heads, seqs


def degap(seq):
    """Séquence acide aminé réelle : enlève gaps, met les insertions en majuscule."""
    return seq.replace("-", "").replace(".", "").upper()


def build_db(heads, seqs, outdir):
    """Écrit la base FASTA (IDs numériques) + mapping id->header. Retourne (db_fasta, id2head)."""
    outdir.mkdir(parents=True, exist_ok=True)
    db_fasta = outdir / "msa_db.fasta"
    id2head = {}
    with open(db_fasta, "w") as out:
        for i, (h, s) in enumerate(zip(heads, seqs)):
            sid = f"seq{i:05d}"
            id2head[sid] = h
            aa = degap(s)
            if aa:
                out.write(f">{sid}\n{aa}\n")
    return db_fasta, id2head


def run_blast(query, db_fasta, outdir, dbsize=None, evalue_cut=10.0, threads=4):
    db_prefix = outdir / "msa_db"
    subprocess.run(
        ["makeblastdb", "-in", str(db_fasta), "-dbtype", "prot",
         "-out", str(db_prefix)],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
    )
    out_tsv = outdir / "blast_results.tsv"
    cols = ["sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs"]
    cmd = [
        "blastp", "-query", str(query), "-db", str(db_prefix),
        "-outfmt", "6 " + " ".join(cols),
        "-evalue", str(evalue_cut),
        "-max_target_seqs", "100000",
        "-num_threads", str(threads),
        "-out", str(out_tsv),
    ]
    if dbsize:
        cmd += ["-dbsize", str(int(dbsize))]
    subprocess.run(cmd, check=True, stderr=subprocess.PIPE)
    df = pd.read_csv(out_tsv, sep="\t", names=cols)
    # garder le meilleur HSP par séquence (plus haut bitscore)
    df = df.sort_values("bitscore", ascending=False).drop_duplicates("sseqid")
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a3m", default=str(A3M_DEFAULT))
    ap.add_argument("--query", default=str(QUERY_DEFAULT))
    ap.add_argument("--outdir", default=str(OUTDIR))
    ap.add_argument("--dbsize", type=float, default=None,
                    help="Taille de base pour l'e-value (ex 150000000 ~ UniRef90)")
    ap.add_argument("--min-pident", type=float, default=None)
    ap.add_argument("--max-evalue", type=float, default=None)
    ap.add_argument("--min-bitscore", type=float, default=None)
    args = ap.parse_args()

    a3m = Path(args.a3m)
    outdir = Path(args.outdir)
    heads, seqs = parse_a3m(a3m)
    print(f"MSA : {len(seqs)-1} séquences (+ query)")

    # base = tous les hits (on saute la query seqs[0])
    db_fasta, id2head = build_db(heads[1:], seqs[1:], outdir)
    print("blastp en cours...")
    df = run_blast(Path(args.query), db_fasta, outdir, dbsize=args.dbsize)
    df["header"] = df["sseqid"].map(id2head)
    print(f"Hits significatifs (e-value <= 10) : {len(df)} / {len(seqs)-1}")

    # table complète triée par bit-score
    table = df[["header", "pident", "evalue", "bitscore", "qcovs", "length"]].copy()
    table = table.sort_values("bitscore", ascending=False)
    out_csv = outdir / "blast_table.csv"
    table.to_csv(out_csv, index=False)
    print(f"Table écrite : {out_csv}")

    # distribution e-value
    print("\nDistribution e-value (-log10) :")
    nl = -np.log10(df["evalue"].clip(lower=1e-300))
    for lo in [0, 5, 10, 20, 50, 100, 200]:
        hi = {0: 5, 5: 10, 10: 20, 20: 50, 50: 100, 100: 200, 200: 1e9}[lo]
        n = ((nl >= lo) & (nl < hi)).sum()
        print(f"  e-value 1e-{lo:<3d}..1e-{hi if hi<1e9 else '∞':<4} : {n:6d}")

    # où tombent les actin-related ?
    import re
    pat = re.compile(r"actin[- ]?related|actin[- ]?like|\bArp\d?\b|ARP\d", re.I)
    arp = df[df["header"].str.contains(pat)]
    if len(arp):
        print(f"\n'actin-related/like/Arp' : {len(arp)} séq")
        print(f"  bitscore  médiane={arp['bitscore'].median():.0f}  "
              f"[{arp['bitscore'].min():.0f} .. {arp['bitscore'].max():.0f}]")
        print(f"  pident    médiane={arp['pident'].median():.0f}%  "
              f"[{arp['pident'].min():.0f} .. {arp['pident'].max():.0f}]")

    # filtrage optionnel -> a3m filtré
    if any(v is not None for v in (args.min_pident, args.max_evalue, args.min_bitscore)):
        keep = pd.Series(True, index=df.index)
        crit = []
        if args.min_pident is not None:
            keep &= df["pident"] >= args.min_pident
            crit.append(f"pident>={args.min_pident}")
        if args.max_evalue is not None:
            keep &= df["evalue"] <= args.max_evalue
            crit.append(f"evalue<={args.max_evalue}")
        if args.min_bitscore is not None:
            keep &= df["bitscore"] >= args.min_bitscore
            crit.append(f"bitscore>={args.min_bitscore}")
        kept_ids = set(df[keep]["sseqid"])
        kept_heads = {id2head[s] for s in kept_ids}

        suffix = "_".join(crit).replace(">=", "ge").replace("<=", "le").replace(".", "p")
        out_a3m = a3m.with_name(a3m.stem + f"_filt_{suffix}.a3m")
        n_written = 0
        with open(out_a3m, "w") as out:
            # toujours garder la query en tête
            out.write(f">{heads[0]}\n{seqs[0]}\n")
            for h, s in zip(heads[1:], seqs[1:]):
                if h in kept_heads:
                    out.write(f">{h}\n{s}\n")
                    n_written += 1
        print(f"\nFiltre [{', '.join(crit)}] -> {n_written} séquences gardées (+query)")
        print(f"a3m filtré écrit : {out_a3m}")


if __name__ == "__main__":
    main()
