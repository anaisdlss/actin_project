#!/usr/bin/env python
"""Refresh only RCSB metadata for retained PDBs, or rebuild offline from its cache.

Examples:
  python tools/update_structure_annotations.py           # fetch missing IDs only
  python tools/update_structure_annotations.py --refresh # refresh current IDs
  python tools/update_structure_annotations.py --offline # reproducible derivation
No raw tables, interaction filters or structural coordinates are modified.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "script"))
from structure_annotations import SCHEMA_VERSION, build_annotations, retained_entries

ENDPOINT = "https://data.rcsb.org/graphql"
FIELDS = """rcsb_id struct{title} polymer_entities{
  rcsb_id entity_poly{rcsb_mutation_count rcsb_entity_polymer_type}
  rcsb_polymer_entity{pdbx_description pdbx_mutation rcsb_macromolecular_names_combined{name provenance_source}}
  rcsb_polymer_entity_container_identifiers{entity_id auth_asym_ids}
  rcsb_entity_source_organism{ncbi_scientific_name ncbi_taxonomy_id}
} nonpolymer_entities{rcsb_id nonpolymer_comp{chem_comp{id name}}}"""


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_json(path, value):
    # Atomic replacement keeps the previous usable cache if interrupted.
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(path)


def fetch_batch(ids):
    query = "{entries(entry_ids:" + json.dumps(ids) + "){" + FIELDS + "}}"
    url = ENDPOINT + "?" + urllib.parse.urlencode({"query": query})
    request = urllib.request.Request(url, headers={"User-Agent": "ActinABP-research-metadata-audit/1.0"})
    with urllib.request.urlopen(request, timeout=40) as response:
        payload = json.load(response)
    if payload.get("errors"):
        raise ValueError("RCSB GraphQL errors: " + json.dumps(payload["errors"]))
    return {entry["rcsb_id"].upper(): entry for entry in payload.get("data", {}).get("entries") or [] if entry}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--reuse-cache", type=Path, help="Import cached metadata from another project snapshot, scoped to current retained IDs")
    args = parser.parse_args(argv)
    if args.offline and args.refresh:
        parser.error("--offline and --refresh cannot be combined")
    root = args.root.resolve()
    base = root / "data/annotations"
    base.mkdir(parents=True, exist_ok=True)
    cache_path = base / "rcsb_structure_metadata.json"
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {"records": {}}
    if args.reuse_cache:
        for key, value in json.loads(args.reuse_cache.read_text()).get("records", {}).items():
            cache["records"].setdefault(key, value)
    ids = sorted(retained_entries(root).pdb_id.unique())
    cache["records"] = {k: v for k, v in cache["records"].items() if k in ids}
    cache.update({"schema_version": SCHEMA_VERSION, "provider": "RCSB PDB Data API", "endpoint": ENDPOINT,
                  "graphql_fields": FIELDS, "scope": "Current filtered_all_data.csv retained PDB IDs"})
    requested = ids if args.refresh else [p for p in ids if not cache["records"].get(p, {}).get("entry")]
    failed = []
    if not args.offline:
        for start in range(0, len(requested), 25):
            batch = requested[start:start + 25]
            try:
                result = fetch_batch(batch)
                fetched = datetime.now(timezone.utc).isoformat()
                for pdb in batch:
                    if pdb in result:
                        cache["records"][pdb] = {"entry": result[pdb], "fetched_utc": fetched, "endpoint": ENDPOINT}
                    else:
                        failed.append(pdb)
                save_json(cache_path, cache)
                print(f"RCSB metadata: {start + len(batch)}/{len(requested)} requested IDs checked", flush=True)
            except Exception as exc:
                failed.extend(batch)
                print(f"Metadata request failed; previous records preserved: {exc}", file=sys.stderr, flush=True)
            time.sleep(.15)
    save_json(cache_path, cache)
    structures, entities = build_annotations(root, cache)
    structures.to_csv(base / "structure_annotations.csv", index=False)
    entities.to_csv(base / "entity_annotations.csv", index=False)
    source_paths = [root / "data/filtered/filtered_all_data.csv", root / "data/filtered/proteins_per_pdb.csv"]
    source_paths = [p for p in source_paths if p.exists()]
    missing = [p for p in ids if not cache["records"].get(p, {}).get("entry")]
    manifest = {
        "schema_version": SCHEMA_VERSION, "generated_utc": datetime.now(timezone.utc).isoformat(),
        "retained_pdb_count": len(ids), "metadata_available": len(ids) - len(missing), "metadata_missing": missing,
        "failed_refresh_ids": sorted(set(failed)), "offline_rebuild": args.offline,
        "sources": [{"file": str(p.relative_to(root)), "sha256": digest(p)} for p in source_paths],
        "source_cache": {"file": str(cache_path.relative_to(root)), "sha256": digest(cache_path)},
        "code_sha256": {name: digest(ROOT / name) for name in ("script/structure_annotations.py", "tools/update_structure_annotations.py")},
        "outputs": {name: digest(base / name) for name in ("structure_annotations.csv", "entity_annotations.csv")},
        "official_documentation": ["https://data.rcsb.org/", "https://www.rcsb.org/docs/search-and-browse/advanced-search/attribute-details"],
        "interpretation": {
            "mutation": "RCSB engineered-mutation count and deposited pdbx_mutation; absence is unknown unless explicit count/text says none. Not pathogenicity or proof of wild type.",
            "toxins": "Name-screened candidates only, including bacterial effector names; not exhaustive and never excluded automatically.",
            "names": "Original PPI3D names preserved; aliases retain RCSB provenance. Multiple organism sources or explicit fusion/chimera text trigger review, not identity correction.",
            "scope": "Metadata can cover additional molecules in the PDB entry beyond retained interacting chains. Filters in the app affect the annotation table only.",
        },
        "reproduce": "python tools/update_structure_annotations.py --offline",
        "refresh": "python tools/update_structure_annotations.py --refresh",
    }
    save_json(base / "manifest.json", manifest)
    print(f"Annotations: {len(structures)} PDBs, {len(entities)} polymer entities; {len(missing)} metadata missing; {len(failed)} failed refresh IDs")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
