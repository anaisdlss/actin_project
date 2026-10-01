"""Offline GO annotation of FoldDisco hits by exact UniProt accession only."""
import json
import re
from pathlib import Path

import pandas as pd


CACHE_PATH = Path("data/annotations/folddisco_go.json")
ACCESSION = re.compile(r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\Z")
ASPECTS = {"F": "molecular_function", "P": "biological_process", "C": "cellular_component"}
GO_COLUMNS = ["go_uniprot_accession", "go_molecular_function", "go_biological_process",
              "go_cellular_component", "go_annotation_status", "go_evidence",
              "go_fetched_utc", "go_source_url"]


def exact_accession(db, target_id):
    """PDB chains and protein names do not establish a UniProt correspondence."""
    value = str(target_id).strip().upper()
    return value if str(db).lower() == "afdb" and ACCESSION.fullmatch(value) else None


def parse_uniprot_entry(entry, fetched_utc, source_url, response_sha256):
    """Preserve GO aspect, term ID, evidence and exact primary entry identity."""
    accession = entry.get("primaryAccession", "")
    if not ACCESSION.fullmatch(accession):
        raise ValueError("UniProt response has no valid primary accession.")
    terms = []
    for item in entry.get("uniProtKBCrossReferences", []):
        if item.get("database") != "GO":
            continue
        properties = {p.get("key"): p.get("value", "") for p in item.get("properties", [])}
        term = properties.get("GoTerm", "")
        aspect, separator, label = term.partition(":")
        if aspect in ASPECTS and separator and re.fullmatch(r"GO:\d{7}", item.get("id", "")):
            terms.append(dict(id=item["id"], aspect=aspect, term=label,
                              evidence=properties.get("GoEvidenceType", "not supplied")))
    return {"primary_accession": accession,
            "status": "annotated" if terms else "no GO terms returned",
            "terms": terms, "entry_type": entry.get("entryType", ""),
            "organism": entry.get("organism", {}).get("scientificName", ""),
            "fetched_utc": fetched_utc, "source_url": source_url,
            "response_sha256": response_sha256}


def load_cache(root="."):
    path = Path(root) / CACHE_PATH
    if not path.exists():
        return {"schema_version": 1, "entries": {}, "requests": []}
    cache = json.loads(path.read_text())
    if cache.get("schema_version") != 1 or not isinstance(cache.get("entries"), dict):
        raise ValueError("Unsupported FoldDisco GO cache format.")
    return cache


def annotate_hits(frame, root="."):
    """Enrich rows without changing hit names, order, scores or multiplicity.

    No web call occurs in the app. A missing cache entry is not evidence of an
    absent biological function. Returned primary IDs must match exactly; names,
    secondary-ID redirects and PDB-chain associations are never guessed.
    """
    cache = load_cache(root)
    annotations = []
    for db, target in zip(frame.get("db", pd.Series(index=frame.index, dtype=str)),
                          frame.get("target_id", pd.Series(index=frame.index, dtype=str))):
        accession = exact_accession(db, target)
        row = {column: "" for column in GO_COLUMNS}
        row["go_uniprot_accession"] = accession or ""
        if accession is None:
            row["go_annotation_status"] = ("PDB-to-UniProt mapping unresolved" if str(db).lower() == "pdb"
                                             else "No exact UniProt accession")
        elif accession not in cache["entries"]:
            row["go_annotation_status"] = "Not fetched (partial GO cache)"
        else:
            entry = cache["entries"][accession]
            status = entry.get("status", "No exact UniProt entry returned")
            if entry.get("primary_accession") != accession:
                status = "No exact UniProt entry returned"
            else:
                for aspect, column in ASPECTS.items():
                    row["go_" + column] = "; ".join(dict.fromkeys(
                        f"{term['term']} [{term['id']}]" for term in entry.get("terms", [])
                        if term["aspect"] == aspect))
                row["go_evidence"] = "; ".join(dict.fromkeys(
                    f"{term['id']}: {term['evidence']}" for term in entry.get("terms", [])))
            row.update(go_annotation_status=status, go_fetched_utc=entry.get("fetched_utc", ""),
                       go_source_url=entry.get("source_url", ""))
        annotations.append(row)
    result = frame.copy()
    for column in GO_COLUMNS:
        result[column] = [row[column] for row in annotations]
    return result


def prioritized_accessions(frame):
    """Take each query motif's best-scoring hit first, then second hits, etc."""
    if frame.empty:
        return []
    valid = frame.copy()
    valid["_accession"] = [exact_accession(db, target) for db, target in zip(valid.db, valid.target_id)]
    valid = valid[valid._accession.notna()]
    if "is_source" in valid:
        valid = valid[~valid.is_source.astype(str).str.lower().eq("true")]
    valid = valid.sort_values("idfscore", ascending=False, kind="stable")
    groups = [column for column in ("query_abp", "query_cluster") if column in valid]
    if groups:
        valid = valid.drop_duplicates(groups + ["_accession"])
        valid["_rank"] = valid.groupby(groups, dropna=False).cumcount()
        valid = valid.sort_values(["_rank", "idfscore"], ascending=[True, False], kind="stable")
    return list(dict.fromkeys(valid._accession))
