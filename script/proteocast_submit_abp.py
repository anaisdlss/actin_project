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


def poll_once(session, job_id):
    """Vérifie le statut UNE fois (non bloquant).
    Renvoie (status, message) : status ∈ finished/error/warning/running/unknown."""
    try:
        j = session.get(f"{STATUS}?job_id={job_id}", timeout=60).json()
    except Exception:
        return "unknown", ""
    s = j.get("status", "running")
    if s in ("finished", "error", "warning"):
        return s, j.get("message") or j.get("warning_message", "")
    return "running", ""


def download_extract(session, job_id, slug):
    """Télécharge le ZIP du job et l'extrait dans data/proteocast/abp/<slug>/.

    Robuste à une structure imbriquée : certains ZIP contiennent les fichiers
    dans un sous-dossier (ex. jobXXXX/…). Si le fichier clé n'est pas à la racine
    après extraction, on le cherche récursivement et on APLATIT le contenu."""
    import shutil
    r = session.get(f"{DL}job{job_id}/", timeout=300)
    r.raise_for_status()
    dest = ABP_DIR / slug
    dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(r.content)) as z:
        z.extractall(dest)
    if is_done(slug):
        return True
    # Pas à la racine → chercher le CSV clé dans les sous-dossiers et aplatir.
    hits = list(dest.rglob("4.query_ProteoCast.csv"))
    if hits:
        _src = hits[0].parent
        for _p in _src.iterdir():
            _t = dest / _p.name
            if not _t.exists():
                shutil.move(str(_p), str(_t))
        return is_done(slug)
    # Vraiment absent : ProteoCast n'a pas produit de résultat pour cet ABP.
    # On logue le contenu du ZIP pour diagnostic.
    try:
        with zipfile.ZipFile(io.BytesIO(r.content)) as z:
            _names = z.namelist()
        print(f"    (ZIP {slug} : {len(_names)} entrées, ex. "
              f"{_names[:5]}) — aucun 4.query_ProteoCast.csv", flush=True)
    except Exception:
        pass
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="slug d'un seul ABP à traiter")
    ap.add_argument("--limit", type=int, help="nb max d'ABP à traiter")
    ap.add_argument("--force", action="store_true", help="refaire même si présent")
    ap.add_argument("--retry-failed", action="store_true",
                    help="ré-essayer les ABP marqués en échec")
    ap.add_argument("--jobs", type=int, default=4,
                    help="nb de jobs en vol en parallèle (défaut 4)")
    args = ap.parse_args()

    man = pd.read_csv(MANIFEST).dropna(subset=["uniprot"])
    if args.only:
        man = man[man["slug"] == args.only]
    failed = set() if (args.force or args.retry_failed) else _load_failed()
    todo = [r for _, r in man.iterrows()
            if args.force or (not is_done(r["slug"]) and r["slug"] not in failed)]
    if args.limit:
        todo = todo[:args.limit]
    window = max(1, args.jobs)
    print(f"{len(todo)} ABP à soumettre (sur {len(man)} ; "
          f"{len(failed)} en échec ignorés) — {window} en parallèle")

    session = requests.Session()
    session.headers.update({"User-Agent": UA})
    ok = 0
    done_n = 0
    total = len(todo)
    new_failed = set(failed)
    pending = list(todo)          # ABP pas encore soumis
    inflight = {}                 # slug -> {jid, uni, title, t0}

    def _fail(slug, msg):
        print(f"    {slug} : {msg} — marqué en échec", flush=True)
        new_failed.add(slug)

    while pending or inflight:
        # 1) remplir la fenêtre : soumettre jusqu'à `window` jobs en vol
        while pending and len(inflight) < window:
            r = pending.pop(0)
            slug, uni = r["slug"], r["uniprot"]
            try:
                jid = submit(session, uni)
            except Exception as e:
                _fail(slug, f"ERREUR soumission: {e}")
                continue
            if not jid:
                _fail(slug, "pas de job_id")
                continue
            inflight[slug] = {"jid": jid, "uni": uni,
                              "title": r["abp_title"], "t0": time.time()}
            print(f"[soumis {len(inflight)}/{window} en vol] "
                  f"{r['abp_title']} ({uni}) — job {jid}", flush=True)
            time.sleep(1)         # court délai entre deux soumissions

        # 2) sonder chaque job en vol UNE fois (non bloquant)
        for slug in list(inflight):
            info = inflight[slug]
            status, msg = poll_once(session, info["jid"])
            if status == "running":
                if time.time() - info["t0"] > POLL_MAX:
                    del inflight[slug]
                    _fail(slug, "timeout")
                    done_n += 1
                continue
            del inflight[slug]
            done_n += 1
            if status != "finished":
                _fail(slug, f"{status} : {msg[:120]}")
                continue
            try:
                got = download_extract(session, info["jid"], slug)
            except Exception as e:
                _fail(slug, f"ERREUR téléchargement: {e}")
                continue
            if got:
                ok += 1
                new_failed.discard(slug)
                print(f"[{done_n}/{total}] OK -> data/proteocast/abp/{slug}/",
                      flush=True)
            else:
                _fail(slug, "ZIP téléchargé mais fichier clé absent")

        # 3) battement de cœur : montrer ce qui tourne + temps écoulé, puis
        #    patienter avant le prochain tour (sinon l'UI a l'air figée).
        if inflight:
            _hb = ", ".join(
                f"{v['title'][:22]} {int(time.time() - v['t0'])}s"
                for v in inflight.values())
            _left = len(pending) + len(inflight)
            print(f"    … {len(inflight)} en cours ({_left} restants) : {_hb}",
                  flush=True)
            time.sleep(POLL_EVERY)

    _save_failed(new_failed)
    print(f"\nterminé : {ok}/{total} ABP calculés")


if __name__ == "__main__":
    main()
