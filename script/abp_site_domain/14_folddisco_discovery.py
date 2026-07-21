#!/usr/bin/env python
"""
DÉCOUVERTE : autres protéines du PDB partageant le motif d'interface d'un ABP,
**décomposé par cluster de site d'actine**.

Complément de 13_folddisco_interface.py. Le script 13 compare tes ABP ENTRE EUX
(index local `abp_chains/`). Ici on cherche le motif d'interface dans TOUT le PDB,
via le serveur web FoldDisco (search.foldseek.com) — l'API par « ticket ».

DÉCOMPOSITION PAR CLUSTER : un ABP qui touche l'actine à plusieurs sites a
PLUSIEURS motifs (un par cluster). On ne prend donc pas un motif unique par ABP,
mais un motif par PAIRE (ABP, cluster). Source du découpage :
`abp_site_table.csv` (toutes les interactions + leur `actin_site_cluster`).
Pour chaque (ABP, cluster) on prend l'interaction de **meilleure résolution**,
ses résidus de contact ABP = le motif, sa structure = la requête.

Méthode, pour chaque (ABP, cluster) :
  1. motif = résidus de contact actine de l'ABP dans CE cluster.
  2. structure = chaîne ABP extraite de l'assembly {iid}_{pdb} (cache local).
  3. POST structure + motif -> /api/ticket/folddisco (base pdb_folddisco).
  4. polling /api/ticket/{id} jusqu'à COMPLETE.
  5. GET /api/result/folddisco/{id} -> hits (target PDB, score, rmsd…).

Sortie : data/exports/abp_site_domain/folddisco_discovery.csv
  (query_abp, query_cluster, query_size, is_source, pdb_hit, target_chain,
   nodecount, coverage, idfscore, rmsd)

Reprise-able : saute les paires (ABP, cluster) déjà présentes dans le CSV.
--force pour tout refaire.

Usage :
  pixi run python script/abp_site_domain/14_folddisco_discovery.py
  pixi run python script/abp_site_domain/14_folddisco_discovery.py --limit 3   # test
"""
import argparse
import os
import re
import time
from pathlib import Path

import pandas as pd
import requests
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/exports/abp_site_domain"
ASM = ROOT / "data/filtered/details/structures_files/assembly"
CHAINS = OUT / "abp_chains_disco"          # chaînes ABP extraites par (pdb, chaîne)
CHAINS.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUT / "folddisco_discovery.csv"

BASE = "https://search.foldseek.com"
SUBMIT = f"{BASE}/api/ticket/folddisco"
# Deux bases interrogées en UNE requête : PDB (structures réelles, hits = chaîne
# d'un complexe) + AlphaFold proteome (1 protéine prédite par hit, nom via UniProt).
DBS = ["pdb_folddisco", "afdb-proteome_folddisco"]

POLL_EVERY = 3               # s entre deux vérifications de statut
POLL_MAX = 120               # s max d'attente par requête
PAUSE_BETWEEN = 5            # s de politesse entre deux requêtes
RETRY_TRIES = 6              # tentatives sur 429/5xx (backoff exponentiel)


# ── extraction de la chaîne ABP depuis l'assembly (cache) ────────────────────
class _ChainSel(Select):
    def __init__(self, cid):
        self.cid = cid

    def accept_chain(self, chain):
        return chain.id == self.cid

    def accept_residue(self, residue):
        return residue.id[0] == " "        # acides aminés seulement


_io = PDBIO()


def chain_pdb(iid, pdb, chain):
    """Renvoie le chemin de la chaîne ABP (extraite/caché) ; None si échec."""
    dest = CHAINS / f"{pdb}_{chain}.pdb"
    if dest.exists():
        return dest
    base = f"{pdb}"
    f = ASM / (base + ".pdb")
    parser = PDBParser(QUIET=True)
    if not f.exists():
        f = ASM / (base + ".cif")
        parser = MMCIFParser(QUIET=True)
    if not f.exists():
        return None
    try:
        struct = parser.get_structure(base, str(f))
        model = next(iter(struct))
        if chain not in [c.id for c in model]:
            return None
        _io.set_structure(struct)
        _io.save(str(dest), _ChainSel(chain))
        return dest
    except Exception:
        return None


# ── API FoldDisco (avec retry backoff sur 429/5xx) ───────────────────────────
def _wait_after(r, delay):
    """Respecte Retry-After si présent, sinon le délai de backoff courant."""
    ra = r.headers.get("Retry-After", "")
    return int(ra) if ra.isdigit() else delay


def _get(url):
    """GET avec retry sur 429/5xx."""
    delay = 15
    for k in range(RETRY_TRIES):
        r = requests.get(url, timeout=120)
        if r.status_code == 429 or r.status_code >= 500:
            w = _wait_after(r, delay)
            print(f"    {r.status_code} (GET) — attente {w}s "
                  f"(essai {k+1}/{RETRY_TRIES})", flush=True)
            time.sleep(w)
            delay = min(delay * 2, 120)
            continue
        r.raise_for_status()
        return r
    r.raise_for_status()
    return r


def submit(pdb_path, motif):
    """POST structure + motif -> id du ticket ; rouvre le fichier à chaque essai."""
    data = [("database[]", db) for db in DBS] + [("motif", motif)]
    delay = 15
    for k in range(RETRY_TRIES):
        with open(pdb_path, "rb") as f:
            r = requests.post(SUBMIT, files={"q": f}, data=data, timeout=60)
        if r.status_code == 429 or r.status_code >= 500:
            w = _wait_after(r, delay)
            print(f"    {r.status_code} (submit) — attente {w}s "
                  f"(essai {k+1}/{RETRY_TRIES})", flush=True)
            time.sleep(w)
            delay = min(delay * 2, 120)
            continue
        r.raise_for_status()
        return r.json().get("id")
    r.raise_for_status()


def wait_complete(tid):
    t0 = time.time()
    while time.time() - t0 < POLL_MAX:
        s = _get(f"{BASE}/api/ticket/{tid}").json().get("status")
        if s in ("COMPLETE", "ERROR"):
            return s
        time.sleep(POLL_EVERY)
    return "TIMEOUT"


def _target_chain(targetresidues):
    ch = [t[0] for t in str(targetresidues).split(",") if t and t[0] != "_"]
    return max(set(ch), key=ch.count) if ch else ""


def _parse_target(dbname, target):
    """(kind, target_id) : PDB -> ('pdb','5yu8') ; AlphaFold -> ('afdb','P45591')."""
    b = os.path.basename(target)
    if dbname.startswith("pdb"):
        return "pdb", b.split(".")[0].lower()
    m = re.search(r"AF-([A-Za-z0-9]+)-F", b)
    return "afdb", (m.group(1) if m else b)


def fetch_hits(tid, query_size, top):
    """Hits des DEUX bases, dédupliqués par cible (PDB: id+chaîne ; AFDB: UniProt)."""
    d = _get(f"{BASE}/api/result/folddisco/{tid}").json()
    out = []
    for res in d.get("results", []):
        dbname = res.get("db", "")
        best = {}
        for a in res.get("alignments", []):
            kind, tgt = _parse_target(dbname, a["target"])
            chain = _target_chain(a.get("targetresidues", "")) if kind == "pdb" else ""
            key = (tgt, chain)
            sc = float(a.get("idfscore", 0))
            if key not in best or sc > best[key]["idfscore"]:
                n = int(a.get("nodecount", 0))
                best[key] = {
                    "db": kind,
                    "target_id": tgt,
                    "target_chain": chain,
                    "nodecount": n,
                    "coverage": round(n / query_size, 3) if query_size else None,
                    "idfscore": round(sc, 1),
                    "rmsd": round(float(a.get("rmsd", 0)), 3),
                }
        out.extend(sorted(best.values(), key=lambda x: -x["idfscore"])[:top])
    return out


# ── construction des motifs par (ABP, cluster) ───────────────────────────────
def _res_num(s):
    m = re.search(r"[\d.]+", str(s))
    return float(m.group()) if m else float("inf")


def build_motifs():
    """[(abp, cluster, src_pdb, src_uni, qsize, motif, iid, pdb, abp_chain)] par (ABP, cluster).

    N'extrait PAS les chaînes ici (fait paresseusement dans la boucle, après --limit)."""
    t = pd.read_csv(OUT / "abp_site_table.csv")
    res = pd.read_csv(ROOT / "data/filtered/details/3.interface_residues.csv")
    uni = (pd.read_csv(OUT / "abp_master.csv")
           .dropna(subset=["uniprot"]).set_index("abp_title")["uniprot"].to_dict())
    t["chain"] = t["pdb"] + "_" + t["abp_chain"].astype(str)
    t["res_num"] = t["resolution"].map(_res_num)
    out = []
    for (abp, cl), g in t.groupby(["abp_title", "actin_site_cluster"]):
        r = g.sort_values("res_num").iloc[0]          # meilleure résolution
        sub = res[(res.interaction_id == int(r.interaction_id))
                  & (res.chain == r.chain)]
        contacts = sorted(set(pd.to_numeric(
            sub.residue_number_structure, errors="coerce").dropna().astype(int)))
        if len(contacts) < 3:
            continue
        motif = ",".join(f"{r.abp_chain}{c}" for c in contacts)
        out.append((abp, cl, str(r.pdb).lower(), uni.get(abp), len(contacts),
                    motif, int(r.interaction_id), r.pdb, r.abp_chain))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None,
                    help="ne traiter que les N premières paires (test)")
    ap.add_argument("--top", type=int, default=1000,
                    help="nb de hits gardés par (ABP, cluster) ET par base")
    ap.add_argument("--force", action="store_true",
                    help="refaire les paires déjà présentes dans le CSV")
    args = ap.parse_args()

    jobs = build_motifs()
    if args.limit:
        jobs = jobs[:args.limit]
    print(f"{len(jobs)} paires (ABP, cluster) à traiter")

    done = set()
    prev_rows = []
    if OUT_CSV.exists() and not args.force:
        prev = pd.read_csv(OUT_CSV)
        prev_rows = prev.to_dict("records")
        done = set(zip(prev["query_abp"], prev["query_cluster"]))

    all_rows = list(prev_rows)
    for i, (abp, cl, src_pdb, src_uni, qsize, motif, iid, pdb, ab_chain) in \
            enumerate(jobs, 1):
        if (abp, cl) in done:
            print(f"[{i}/{len(jobs)}] {abp} @ {cl} — déjà fait, saut")
            continue
        print(f"[{i}/{len(jobs)}] {abp} @ {cl} ({qsize} rés.) — soumission…",
              flush=True)
        try:
            cpath = chain_pdb(iid, pdb, ab_chain)     # extraction paresseuse
            if cpath is None:
                print("    chaîne introuvable — ignoré", flush=True)
                continue
            tid = submit(cpath, motif)
            status = wait_complete(tid)
            if status != "COMPLETE":
                print(f"    statut {status} — ignoré", flush=True)
                continue
            hits = fetch_hits(tid, qsize, args.top)
            for h in hits:
                src = src_pdb if h["db"] == "pdb" else src_uni
                all_rows.append({
                    "query_abp": abp, "query_cluster": cl, "query_size": qsize,
                    "is_source": h["target_id"] == src, **h})
            _npdb = sum(h["db"] == "pdb" for h in hits)
            print(f"    {len(hits)} hits ({_npdb} PDB + {len(hits)-_npdb} AFDB)",
                  flush=True)
            pd.DataFrame(all_rows).to_csv(OUT_CSV, index=False)   # incrémental
        except Exception as e:
            print(f"    ERREUR: {e}", flush=True)
        time.sleep(PAUSE_BETWEEN)

    df = pd.DataFrame(all_rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nécrit : {OUT_CSV}  ({len(df)} lignes, "
          f"{df.groupby(['query_abp','query_cluster']).ngroups if len(df) else 0} "
          f"paires ABP×cluster)")


if __name__ == "__main__":
    main()
