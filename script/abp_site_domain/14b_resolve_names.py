#!/usr/bin/env python
"""
Résolution des NOMS de protéines pour les hits de découverte (script 14).

Les réponses FoldDisco ne donnent qu'un identifiant :
  - AlphaFold : accession UniProt (ex. P45591)  -> nom + organisme via UniProt REST
  - PDB       : (id, chaîne) (ex. 5yu8 / J)      -> nom de la chaîne via RCSB GraphQL

Sortie (cache, ré-utilisé et complété) :
  data/exports/abp_site_domain/folddisco_names.csv
    (db, key, target_id, target_chain, protein_name, organism)
  key = target_id (afdb) | f"{target_id}_{chain}" (pdb) — clé de jointure.

Idempotent : ne résout que les identifiants pas encore dans le cache.
Usage : pixi run python script/abp_site_domain/14b_resolve_names.py
"""
import time
from io import StringIO
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/exports/abp_site_domain"
DISCO = OUT / "folddisco_discovery.csv"
NAMES = OUT / "folddisco_names.csv"

UNIPROT = "https://rest.uniprot.org/uniprotkb/accessions"
RCSB_GQL = "https://data.rcsb.org/graphql"


def _key(db, tid, chain):
    return tid if db == "afdb" else f"{tid}_{chain}"


def resolve_uniprot(accs):
    """{accession: (nom, organisme)} par lots (UniProt REST, format TSV)."""
    out = {}
    for i in range(0, len(accs), 100):
        chunk = accs[i:i + 100]
        for attempt in range(5):
            r = requests.get(UNIPROT, params={
                "accessions": ",".join(chunk),
                "fields": "accession,protein_name,organism_name",
                "format": "tsv"}, timeout=60)
            if r.status_code == 200:
                break
            time.sleep(5 * (attempt + 1))
        if r.status_code != 200:
            print(f"  UniProt {r.status_code} sur lot {i}")
            continue
        df = pd.read_csv(StringIO(r.text), sep="\t")
        for _, row in df.iterrows():
            out[row["Entry"]] = (str(row.get("Protein names", "")),
                                 str(row.get("Organism", "")))
        print(f"  UniProt {i + len(chunk)}/{len(accs)}", flush=True)
        time.sleep(0.3)
    return out


def resolve_pdb(entries):
    """{(pdb_id, auth_chain): (nom, organisme)} via RCSB GraphQL, par lots."""
    out = {}
    ids = sorted({e.upper() for e in entries})
    for i in range(0, len(ids), 50):
        chunk = ids[i:i + 50]
        q = ('{entries(entry_ids:["' + '","'.join(chunk) + '"]){rcsb_id '
             'polymer_entities{rcsb_polymer_entity{pdbx_description} '
             'rcsb_polymer_entity_container_identifiers{auth_asym_ids} '
             'rcsb_entity_source_organism{scientific_name}}}}')
        for attempt in range(5):
            r = requests.post(RCSB_GQL, json={"query": q}, timeout=60)
            if r.status_code == 200:
                break
            time.sleep(5 * (attempt + 1))
        if r.status_code != 200:
            print(f"  RCSB {r.status_code} sur lot {i}")
            continue
        for ent in (r.json().get("data", {}).get("entries") or []):
            pid = ent["rcsb_id"].lower()
            for pe in (ent.get("polymer_entities") or []):
                desc = (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description", "")
                orgs = pe.get("rcsb_entity_source_organism") or []
                org = orgs[0].get("scientific_name", "") if orgs else ""
                chains = (pe.get("rcsb_polymer_entity_container_identifiers") or {}
                          ).get("auth_asym_ids") or []
                for c in chains:
                    out[(pid, c)] = (str(desc), str(org))
        print(f"  RCSB {i + len(chunk)}/{len(ids)}", flush=True)
        time.sleep(0.3)
    return out


def main():
    df = pd.read_csv(DISCO)
    cache = pd.read_csv(NAMES) if NAMES.exists() else pd.DataFrame(
        columns=["db", "key", "target_id", "target_chain",
                 "protein_name", "organism"])
    known = set(cache["key"])

    df["key"] = [_key(d, t, c) for d, t, c in
                 zip(df["db"], df["target_id"], df["target_chain"].fillna(""))]
    todo = df[~df["key"].isin(known)].drop_duplicates("key")
    afdb = todo[todo["db"] == "afdb"]
    pdb = todo[todo["db"] == "pdb"]
    print(f"à résoudre : {len(afdb)} AlphaFold + "
          f"{pdb['target_id'].nunique()} entrées PDB")

    rows = []
    if len(afdb):
        umap = resolve_uniprot(sorted(afdb["target_id"].unique()))
        for _, r in afdb.iterrows():
            nm, org = umap.get(r["target_id"], ("", ""))
            rows.append({"db": "afdb", "key": r["key"], "target_id": r["target_id"],
                         "target_chain": "", "protein_name": nm, "organism": org})
    if len(pdb):
        pmap = resolve_pdb(pdb["target_id"].unique())
        for _, r in pdb.iterrows():
            nm, org = pmap.get((r["target_id"], r["target_chain"]), ("", ""))
            rows.append({"db": "pdb", "key": r["key"], "target_id": r["target_id"],
                         "target_chain": r["target_chain"],
                         "protein_name": nm, "organism": org})

    if rows:
        cache = pd.concat([cache, pd.DataFrame(rows)], ignore_index=True)
        cache.drop_duplicates("key", keep="last").to_csv(NAMES, index=False)
    print(f"\ncache écrit : {NAMES}  ({len(cache)} noms)")


if __name__ == "__main__":
    main()
