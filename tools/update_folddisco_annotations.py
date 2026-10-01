"""Fetch a bounded, reproducible GO cache from the official UniProt REST API.

Default: up to 200 not-yet-cached exact AlphaFold/UniProt hit accessions, with
best hits distributed across query motifs. The app never contacts this API.
Run again to extend coverage; --refresh explicitly refreshes cached entries.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "script"))
from folddisco_annotations import (CACHE_PATH, exact_accession, load_cache,
                                  parse_uniprot_entry, prioritized_accessions)


def save_json(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(content, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--extra-discovery", type=Path, action="append", default=[])
    parser.add_argument("--limit", type=int, default=200)
    parser.add_argument("--batch-size", type=int, default=25)
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--accessions", nargs="*", help="Exact accessions to request, instead of priority selection")
    args = parser.parse_args()
    if args.limit < 1 or not 1 <= args.batch_size <= 50:
        parser.error("Use a positive limit and a batch size from 1 to 50.")
    paths = [args.root / "data/exports/abp_site_domain/folddisco_discovery.csv"] + args.extra_discovery
    sources = []
    frames = []
    for path in paths:
        payload = path.read_bytes()
        frames.append(pd.read_csv(path))
        sources.append({"name": path.name, "sha256": hashlib.sha256(payload).hexdigest()})
    data = pd.concat(frames, ignore_index=True)
    priority = prioritized_accessions(data)
    if args.accessions:
        requested = [exact_accession("afdb", accession) for accession in args.accessions]
        if any(accession is None for accession in requested):
            parser.error("Every accession must be an exact canonical UniProt accession.")
        priority = list(dict.fromkeys(requested))
    cache = load_cache(args.root)
    selected = [accession for accession in priority if args.refresh or accession not in cache["entries"]][:args.limit]
    cache["source"] = "UniProtKB REST API; GO terms and evidence by exact primary accession"
    cache["scope"] = "Partial candidate cache; no PDB-to-UniProt mapping and no inference by protein name."
    cache["selection"] = {"policy": "best idfscore per query motif first, then subsequent ranks; source hits excluded",
                          "eligible_exact_accessions": len(priority), "limit_this_run": args.limit,
                          "source_discoveries": sources}
    for offset in range(0, len(selected), args.batch_size):
        batch = selected[offset:offset + args.batch_size]
        query = " OR ".join(f"accession:{accession}" for accession in batch)
        url = "https://rest.uniprot.org/uniprotkb/search?" + urlencode(
            {"query": f"({query})", "format": "json", "size": 500,
             "fields": "accession,go,organism_name,reviewed"})
        request = Request(url, headers={"Accept": "application/json", "User-Agent": "actin-abp-research/1.0"})
        # Never record a failed network request as an absence of GO annotation.
        with urlopen(request, timeout=45) as response:
            payload = response.read()
            release = response.headers.get("X-UniProt-Release", "")
            release_date = response.headers.get("X-UniProt-Release-Date", "")
            total = int(response.headers.get("X-Total-Results", "0"))
        document = json.loads(payload)
        if total > len(document.get("results", [])):
            raise RuntimeError("UniProt response was paginated unexpectedly; cache was not updated for this batch.")
        fetched = datetime.now(timezone.utc).isoformat()
        digest = hashlib.sha256(payload).hexdigest()
        parsed = {item["primaryAccession"]: parse_uniprot_entry(item, fetched, url, digest)
                  for item in document.get("results", []) if item.get("primaryAccession") in batch}
        for accession in batch:
            cache["entries"][accession] = parsed.get(accession, {
                "status": "No exact UniProt entry returned", "primary_accession": None,
                "terms": [], "fetched_utc": fetched, "source_url": url, "response_sha256": digest})
        record = {"url": url, "requested_accessions": batch, "fetched_utc": fetched,
                  "response_sha256": digest, "uniprot_release": release,
                  "uniprot_release_date": release_date}
        cache.setdefault("requests", []).append(record)
        raw_path = args.root / "data/annotations/folddisco_go_raw" / f"{digest}.json"
        raw_path.parent.mkdir(parents=True, exist_ok=True)
        raw_path.write_bytes(payload)
        save_json(args.root / CACHE_PATH, cache)
        print(f"Cached {offset + len(batch)}/{len(selected)} requested accessions", flush=True)
    save_json(args.root / CACHE_PATH, cache)
    annotated = sum(entry.get("status") == "annotated" for entry in cache["entries"].values())
    print(f"GO annotations available for {annotated} of {len(cache['entries'])} requested entries; "
          f"{len(priority)} eligible exact accessions in current selection.")


if __name__ == "__main__":
    main()
