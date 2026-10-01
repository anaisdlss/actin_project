"""Traceable, descriptive annotations of retained PDB entries.

No annotation in this module removes an interaction, renames a source protein,
or interprets an engineered mutation as pathogenic. Network access is confined
to tools/update_structure_annotations.py; the app uses the local snapshot.
"""
from pathlib import Path
import json
import re

import pandas as pd


SCHEMA_VERSION = 1
MISSING = {"", ".", "?", "unknown", "nan", "null", "n/a", "not available"}
NO_MUTATION = {"no", "none", "false", "0", "no mutation", "no mutations"}
TOXIN_PATTERNS = (
    ("toxin name", r"\b(?:toxin|exotoxin|enterotoxin)\b"),
    ("phalloidin name", r"\bphalloidin\b"),
    ("jasplakinolide name", r"\bjasplakinolide\b"),
    ("latrunculin name", r"\blatrunculin\b"),
    ("bacterial effector name (review)", r"\b(?:exoy|vop[lvf]|sipa|tccc\d*)\b"),
)


def text(value):
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return str(value).strip()


def mutation_annotation(entity):
    """Use deposited mutation text and RCSB engineered-mutation count.

    Missing fields remain unknown. A deposited positive annotation wins over
    count=0, with a separate conflict flag, rather than being silently dropped.
    """
    details = entity.get("rcsb_polymer_entity") or {}
    poly = entity.get("entity_poly") or {}
    raw = text(details.get("pdbx_mutation"))
    value = poly.get("rcsb_mutation_count")
    count = None
    if value is not None and not isinstance(value, bool):
        try:
            number = float(value)
            if number >= 0 and number.is_integer():
                count = int(number)
        except (ValueError, TypeError):
            pass
    lower = raw.lower()
    positive_text = lower not in MISSING | NO_MUTATION
    negative_text = lower in NO_MUTATION
    positive = positive_text or (count is not None and count > 0)
    negative = negative_text or count == 0
    return {
        "mutation_status": "annotated" if positive else "none_reported" if negative else "unknown",
        "mutation_text": raw,
        "engineered_mutation_count": count,
        "mutation_field_conflict": bool((positive_text and count == 0) or (negative_text and count is not None and count > 0)),
    }


def toxin_name_candidates(names):
    """Conservative name screening, not a toxin ontology or a safety label."""
    evidence = []
    for name, source_field in names:
        for label, pattern in TOXIN_PATTERNS:
            if re.search(pattern, text(name), flags=re.I):
                evidence.append(f"{label}: {source_field} = {name}")
    return sorted(set(evidence))


def retained_entries(root):
    df = pd.read_csv(Path(root) / "data/filtered/filtered_all_data.csv", low_memory=False)
    df["pdb_id"] = df.pdb_id.astype(str).str.strip().str.upper()
    return df


def local_proteins(root, retained):
    """Preserve author-chain case and original PPI3D protein names."""
    path = Path(root) / "data/filtered/proteins_per_pdb.csv"
    if path.exists():
        out = pd.read_csv(path).rename(columns={"chain": "chain_id", "protein": "source_name"})
        out["source_field"] = "proteins_per_pdb.csv:protein"
    else:
        parts = []
        for side in (1, 2):
            part = retained[["pdb_id", f"subunit_{side}", f"subunit_{side}_title"]].copy()
            part.columns = ["pdb_id", "chain_id", "source_name"]
            part["source_field"] = f"filtered_all_data.csv:subunit_{side}_title"
            parts.append(part)
        out = pd.concat(parts, ignore_index=True).drop_duplicates()
    out["pdb_id"] = out.pdb_id.astype(str).str.upper()
    out["author_chain"] = out.chain_id.map(lambda x: text(x).split("_", 1)[-1])
    return out[out.pdb_id.isin(set(retained.pdb_id))]


def build_annotations(root, cache):
    retained = retained_entries(root)
    proteins = local_proteins(root, retained)
    records = cache.get("records", {})
    structures, entities = [], []
    for pdb in sorted(retained.pdb_id.unique()):
        source = records.get(pdb, {})
        entry = source.get("entry") or {}
        fetched = source.get("fetched_utc", "")
        local = proteins[proteins.pdb_id.eq(pdb)]
        local_titles = sorted(set(retained.loc[retained.pdb_id.eq(pdb), "pdb_annotation"].dropna().map(text)))
        title = text((entry.get("struct") or {}).get("title"))
        polymer_records = []
        toxin_evidence, naming_evidence = [], []
        for entity in entry.get("polymer_entities") or []:
            detail = entity.get("rcsb_polymer_entity") or {}
            ids = entity.get("rcsb_polymer_entity_container_identifiers") or {}
            chains = ids.get("auth_asym_ids") or []
            matched = local[local.author_chain.isin(chains)]
            aliases = detail.get("rcsb_macromolecular_names_combined") or []
            description = text(detail.get("pdbx_description"))
            organisms = entity.get("rcsb_entity_source_organism") or []
            source_names = sorted(set(matched.source_name.dropna().map(text)))
            names = [(description, "rcsb_polymer_entity.pdbx_description")]
            names += [(a.get("name"), "rcsb_macromolecular_names_combined:" + text(a.get("provenance_source"))) for a in aliases]
            names += [(r.source_name, r.source_field) for r in matched.itertuples()]
            candidates = toxin_name_candidates(names)
            toxin_evidence += [f"{entity.get('rcsb_id', pdb)}: {x}" for x in candidates]
            reasons = []
            if re.search(r"\b(?:chimera|chimeric|fusion)\b", description, flags=re.I):
                reasons.append("RCSB description explicitly mentions chimera/fusion")
            taxa = {o.get("ncbi_taxonomy_id") for o in organisms if o.get("ncbi_taxonomy_id") is not None}
            if len(taxa) > 1:
                reasons.append("Multiple source-organism taxonomies for one polymer entity; review construct")
            if any(n.casefold() != description.casefold() for n in source_names):
                reasons.append("PPI3D and RCSB deposited names differ; original names retained")
            naming_evidence += [f"{entity.get('rcsb_id', pdb)}: {r}" for r in reasons]
            annotation = mutation_annotation(entity)
            row = {"pdb_id": pdb, "entity_id": text(ids.get("entity_id")),
                   "rcsb_entity_id": entity.get("rcsb_id", ""),
                   "polymer_type": (entity.get("entity_poly") or {}).get("rcsb_entity_polymer_type", ""),
                   "author_chains": "; ".join(chains), "source_protein_names": " | ".join(source_names),
                   "rcsb_deposited_name": description, "source_name_fields": " | ".join(sorted(set(matched.source_field))),
                   "rcsb_aliases_with_provenance": json.dumps(aliases, ensure_ascii=False),
                   "source_organisms": json.dumps(organisms, ensure_ascii=False), **annotation,
                   "toxin_effector_name_candidate": bool(candidates), "toxin_candidate_evidence": " | ".join(candidates),
                   "name_or_construct_review": bool(reasons), "name_review_evidence": " | ".join(reasons),
                   "annotation_url": f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{ids.get('entity_id', '')}",
                   "fetched_utc": fetched}
            entities.append(row)
            if row["polymer_type"] == "Protein":
                polymer_records.append(row)
        for entity in entry.get("nonpolymer_entities") or []:
            chem = (entity.get("nonpolymer_comp") or {}).get("chem_comp") or {}
            toxin_evidence += [f"{entity.get('rcsb_id', pdb)}: {x}" for x in toxin_name_candidates([(chem.get("name"), "nonpolymer_comp.chem_comp.name")])]
        # Local source names also work offline, but do not establish a negative.
        toxin_evidence += toxin_name_candidates([(r.source_name, r.source_field) for r in local.itertuples()])
        states = [r["mutation_status"] for r in polymer_records]
        status = "annotated" if "annotated" in states else "none_reported" if states and all(s == "none_reported" for s in states) else "unknown"
        mutation_evidence = [f"entity {r['entity_id']} ({r['rcsb_deposited_name']}): count={r['engineered_mutation_count']}; deposited={r['mutation_text'] or '[absent]'}" for r in polymer_records if r["mutation_status"] == "annotated"]
        title_hint = bool(re.search(r"\b(?:mutant|mutation|mutated)\b", " ".join([title, *local_titles]), flags=re.I))
        structures.append({"pdb_id": pdb, "source_structure_title": " | ".join(local_titles), "rcsb_structure_title": title,
                           "mutation_status": status, "mutation_evidence": " | ".join(mutation_evidence),
                           "title_mentions_mutation": title_hint,
                           "mutation_annotation_review": any(r["mutation_field_conflict"] for r in polymer_records) or (title_hint and status != "annotated"),
                           "clinical_interpretation": "Not assessed; engineered mutation is not a pathogenicity classification",
                           "toxin_effector_name_candidate": bool(toxin_evidence), "toxin_candidate_evidence": " | ".join(sorted(set(toxin_evidence))),
                           "name_or_construct_review": bool(naming_evidence), "name_review_evidence": " | ".join(sorted(set(naming_evidence))),
                           "metadata_status": "available" if entry else "unknown",
                           "source_url": f"https://www.rcsb.org/structure/{pdb}", "fetched_utc": fetched})
    entity_columns = ["pdb_id", "entity_id", "rcsb_entity_id", "polymer_type", "author_chains",
                      "source_protein_names", "rcsb_deposited_name", "source_name_fields",
                      "rcsb_aliases_with_provenance", "source_organisms", "mutation_status", "mutation_text",
                      "engineered_mutation_count", "mutation_field_conflict", "toxin_effector_name_candidate",
                      "toxin_candidate_evidence", "name_or_construct_review", "name_review_evidence",
                      "annotation_url", "fetched_utc"]
    return pd.DataFrame(structures), pd.DataFrame(entities, columns=entity_columns)


def render_structure_annotations(root=".", key_prefix="structure_audit"):
    """Drop-in, offline Streamlit summary-table audit; filters affect this table only."""
    import streamlit as st

    base = Path(root) / "data/annotations"
    path = base / "structure_annotations.csv"
    with st.expander("Structure annotations: mutations, toxin candidates and source names"):
        st.caption("Annotations describe the current retained structures. Engineered mutations are not clinical classifications. "
                   "Toxin/effector candidates are name-screening flags for review, not automatic exclusions. "
                   "The filters below change only this audit table; all analyses keep their existing dataset.")
        if not path.exists():
            st.info("The structure annotation snapshot is not installed.")
            return
        frame = pd.read_csv(path, keep_default_na=False)
        # Prevent an old annotation snapshot from reintroducing filtered-out IDs.
        current = set(retained_entries(root).pdb_id)
        frame = frame[frame.pdb_id.isin(current)]
        missing = current - set(frame.pdb_id)
        if missing:
            st.info(f"{len(missing)} retained structures have no annotation record in this snapshot.")
        dates = sorted({text(value)[:10] for value in frame.fetched_utc if text(value)})
        if dates:
            st.caption(f"Official RCSB metadata retrieved: {dates[0]}" + (f" to {dates[-1]}" if len(dates) > 1 else "") + ". Metadata may cover molecules beyond the retained interacting chains.")
        st.caption("Mutation status: annotated = a positive deposited annotation/count; none_reported = all protein entities "
                   "explicitly report zero/none; unknown = insufficient metadata. None_reported does not prove a wild-type sequence. "
                   "Original PPI3D names and dated RCSB names/aliases remain separate.")
        selected = st.selectbox("Structure annotation filter", ["All retained structures", "Mutation annotated", "Mutation status unknown", "Toxin/effector name candidates", "Source-name or construct review", "Mutation metadata to review"], key=f"{key_prefix}_filter")
        filters = {"Mutation annotated": frame.mutation_status.eq("annotated"),
                   "Mutation status unknown": frame.mutation_status.eq("unknown"),
                   "Toxin/effector name candidates": frame.toxin_effector_name_candidate.astype(str).str.lower().eq("true"),
                   "Source-name or construct review": frame.name_or_construct_review.astype(str).str.lower().eq("true"),
                   "Mutation metadata to review": frame.mutation_annotation_review.astype(str).str.lower().eq("true")}
        shown = frame[filters[selected]] if selected in filters else frame
        st.caption(f"{len(shown)} of {len(frame)} annotated retained structures shown.")
        st.dataframe(shown, hide_index=True, width="stretch")
        st.download_button("Download displayed structure annotations", shown.to_csv(index=False).encode(), file_name="structure_annotations.csv", mime="text/csv", key=f"{key_prefix}_csv")
        entity_path = base / "entity_annotations.csv"
        if entity_path.exists():
            details = pd.read_csv(entity_path, keep_default_na=False)
            details = details[details.pdb_id.isin(set(shown.pdb_id))]
            st.markdown("**Polymer entities and source aliases**")
            st.dataframe(details, hide_index=True, width="stretch")
            st.download_button("Download entity annotations and aliases", details.to_csv(index=False).encode(), file_name="entity_annotations.csv", mime="text/csv", key=f"{key_prefix}_entities")
        manifest = base / "manifest.json"
        if manifest.exists():
            st.download_button("Download annotation provenance", manifest.read_bytes(), file_name="structure_annotation_manifest.json", mime="application/json", key=f"{key_prefix}_manifest")
