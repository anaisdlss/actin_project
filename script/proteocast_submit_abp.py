#!/usr/bin/env python
"""
Soumission automatique de ProteoCast pour chaque ABP (proteocast.ijm.fr).

Le site n'a pas d'API officielle, mais on reproduit le formulaire (comme pour
PPI3D dans get_summary_results.py). Flux validé :
  1. GET  /results/search/                     -> csrfmiddlewaretoken + cookie
  2. POST /results/upload_file/                 -> JSON {job_id}
       inputOption=uniprot (Retrieve from AlphaFoldDB) + uniprot_id=<acc>
       structureSource=alphafold + uniprotId=<acc>   (scores informés 3D)
  3. GET  /results/check_job_status/?job_id=..  -> {status: finished|error|warning}
  4. GET  /results/download_folder/job<id>/     -> ZIP du dossier résultat
  5. extraction dans data/proteocast/abp/<slug>/ (ce que proteocast_view attend)

Reprise-able : saute les ABP dont <slug>/4.query_ProteoCast.csv existe déjà.

Usage :
  pixi run python script/proteocast_submit_abp.py --only Cofilin_1   # 1 ABP (test)
  pixi run python script/proteocast_submit_abp.py                    # tous les manquants
  pixi run python script/proteocast_submit_abp.py --limit 3
"""
import argparse
import io
import time
import zipfile
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data/proteocast/abp_inputs/manifest.csv"
ABP_DIR = ROOT / "data/proteocast/abp"
# ABP dont la soumission a échoué définitivement (UniProt invalide, fusion, etc.)
# → on ne les re-soumet pas à chaque run (sauf --retry-failed).
FAILED = ABP_DIR / "_failed_slugs.txt"


def _load_failed():
    if FAILED.exists():
        return {l.strip() for l in FAILED.read_text().splitlines() if l.strip()}
    return set()


def _save_failed(slugs):
    FAILED.parent.mkdir(parents=True, exist_ok=True)
    FAILED.write_text("\n".join(sorted(slugs)) + ("\n" if slugs else ""))

BASE = "https://proteocast.ijm.fr"
SEARCH = f"{BASE}/results/search/"
UPLOAD = f"{BASE}/results/upload_file/"
STATUS = f"{BASE}/results/check_job_status/"
DL = f"{BASE}/results/download_folder/"

POLL_EVERY = 10        # s entre deux vérifications de statut
POLL_MAX = 3600        # s max d'attente par job (grosses protéines = longues)
UA = "Mozilla/5.0"


def is_done(slug):
    return (ABP_DIR / slug / "4.query_ProteoCast.csv").exists()


def submit(session, uniprot):
    """POST le formulaire -> job_id (ou None)."""
    import re
    html = session.get(SEARCH, timeout=60).text
    m = re.search(r'name="csrfmiddlewaretoken"\s+value="([^"]+)"', html)
    csrf = m.group(1) if m else None
    data = {
        "csrfmiddlewaretoken": csrf,
        "inputOption": "uniprot", "uniprot_id": uniprot,
        "structureSource": "alphafold", "uniprotId": uniprot,
    }
    r = session.post(UPLOAD, data=data, headers={"Referer": SEARCH}, timeout=120)
    r.raise_for_status()
    return r.json().get("job_id")


def wait(session, job_id):
    """Poll le statut ; renvoie (status_final, message)."""
    t0 = time.time()
    while time.time() - t0 < POLL_MAX:
        try:
            j = session.get(f"{STATUS}?job_id={job_id}", timeout=60).json()
        except Exception:
            time.sleep(POLL_EVERY)
            continue
        s = j.get("status", "")
        if s in ("finished", "error", "warning"):
            return s, j.get("message") or j.get("warning_message", "")
        time.sleep(POLL_EVERY)
    return "timeout", ""


def download_extract(session, job_id, slug):
    """Télécharge le ZIP du job et l'extrait dans data/proteocast/abp/<slug>/."""
    r = session.get(f"{DL}job{job_id}/", timeout=300)
    r.raise_for_status()
    dest = ABP_DIR / slug
    dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(r.content)) as z:
        z.extractall(dest)
    return is_done(slug)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="slug d'un seul ABP à traiter")
    ap.add_argument("--limit", type=int, help="nb max d'ABP à traiter")
    ap.add_argument("--force", action="store_true", help="refaire même si présent")
    ap.add_argument("--retry-failed", action="store_true",
                    help="ré-essayer les ABP marqués en échec")
    args = ap.parse_args()

    man = pd.read_csv(MANIFEST).dropna(subset=["uniprot"])
    if args.only:
        man = man[man["slug"] == args.only]
    failed = set() if (args.force or args.retry_failed) else _load_failed()
    todo = [r for _, r in man.iterrows()
            if args.force or (not is_done(r["slug"]) and r["slug"] not in failed)]
    if args.limit:
        todo = todo[:args.limit]
    print(f"{len(todo)} ABP à soumettre (sur {len(man)} ; "
          f"{len(failed)} en échec ignorés)")

    session = requests.Session()
    session.headers.update({"User-Agent": UA})
    ok = 0
    new_failed = set(failed)
    for i, r in enumerate(todo, 1):
        slug, uni = r["slug"], r["uniprot"]
        print(f"[{i}/{len(todo)}] {r['abp_title']} ({uni}) — soumission…", flush=True)
        _err = None
        try:
            jid = submit(session, uni)
            if not jid:
                _err = "pas de job_id"
            else:
                print(f"    job {jid} — attente…", flush=True)
                status, msg = wait(session, jid)
                if status != "finished":
                    _err = f"{status} : {msg[:120]}"
                elif download_extract(session, jid, slug):
                    ok += 1
                    new_failed.discard(slug)
                    print(f"    OK -> data/proteocast/abp/{slug}/", flush=True)
                else:
                    _err = "ZIP téléchargé mais fichier clé absent"
        except Exception as e:
            _err = f"ERREUR: {e}"
        if _err:
            print(f"    {_err} — marqué en échec", flush=True)
            new_failed.add(slug)
        time.sleep(3)
    _save_failed(new_failed)
    print(f"\nterminé : {ok}/{len(todo)} ABP calculés")


if __name__ == "__main__":
    main()
