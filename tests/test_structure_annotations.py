import json
from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "script"))
from structure_annotations import build_annotations, mutation_annotation, toxin_name_candidates


def entity(mutation=None, count=None):
    return {"rcsb_polymer_entity": {"pdbx_mutation": mutation}, "entity_poly": {"rcsb_mutation_count": count}}


class StructureAnnotationTests(unittest.TestCase):
    def test_mutation_unknown_is_not_negative(self):
        for record in ({}, entity(), entity("?"), entity(".")):
            self.assertEqual(mutation_annotation(record)["mutation_status"], "unknown")
        self.assertEqual(mutation_annotation(entity(count=0))["mutation_status"], "none_reported")
        for value in (False, "false", "NO", "none"):
            self.assertEqual(mutation_annotation(entity(value))["mutation_status"], "none_reported")

    def test_mutation_positive_and_conflict_retained(self):
        result = mutation_annotation(entity("R183W", 0))
        self.assertEqual(result["mutation_status"], "annotated")
        self.assertTrue(result["mutation_field_conflict"])
        self.assertEqual(result["mutation_text"], "R183W")
        self.assertEqual(mutation_annotation(entity(count=2))["mutation_status"], "annotated")

    def test_name_screening_has_evidence_and_is_not_general_bacterial_flag(self):
        evidence = toxin_name_candidates([("Adenylate cyclase ExoY", "deposited_name")])
        self.assertTrue(any("deposited_name = Adenylate cyclase ExoY" in s for s in evidence))
        self.assertEqual(toxin_name_candidates([("Bacterial actin", "name"), ("Cofilin", "name")]), [])

    def test_offline_scope_alias_provenance_chain_case_and_missing(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory) / "data/filtered"
            base.mkdir(parents=True)
            pd.DataFrame([{"pdb_id": "1abc", "pdb_annotation": "original title"}, {"pdb_id": "2abc", "pdb_annotation": "second"}]).to_csv(base / "filtered_all_data.csv", index=False)
            pd.DataFrame([{"pdb_id": "1abc", "chain": "1abc_A", "protein": "Original actin"}, {"pdb_id": "1abc", "chain": "1abc_a", "protein": "Original partner"}]).to_csv(base / "proteins_per_pdb.csv", index=False)
            item = entity("NO", 0)
            item.update({"rcsb_id": "1ABC_1", "rcsb_polymer_entity_container_identifiers": {"entity_id": "1", "auth_asym_ids": ["A"]}})
            item["entity_poly"]["rcsb_entity_polymer_type"] = "Protein"
            item["rcsb_polymer_entity"].update({"pdbx_description": "Actin", "rcsb_macromolecular_names_combined": [{"name": "ACTB", "provenance_source": "UniProt"}]})
            cache = {"records": {"1ABC": {"fetched_utc": "2026-10-01T00:00:00Z", "entry": {"struct": {"title": "RCSB title"}, "polymer_entities": [item]}}, "3ABC": {"entry": {"struct": {"title": "out of scope"}}}}}
            structures, entities = build_annotations(directory, cache)
            self.assertEqual(set(structures.pdb_id), {"1ABC", "2ABC"})
            self.assertEqual(structures.set_index("pdb_id").loc["2ABC", "mutation_status"], "unknown")
            row = entities.iloc[0]
            self.assertEqual(row.source_protein_names, "Original actin")
            self.assertEqual(json.loads(row.rcsb_aliases_with_provenance)[0]["provenance_source"], "UniProt")
            self.assertEqual(structures.iloc[0].source_structure_title, "original title")
            self.assertEqual(row.fetched_utc, "2026-10-01T00:00:00Z")


if __name__ == "__main__":
    unittest.main()
