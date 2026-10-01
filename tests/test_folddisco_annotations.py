"""GO terms must remain tied to exact source identifiers, with unknowns explicit."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "script"))
from folddisco_annotations import (CACHE_PATH, annotate_hits, exact_accession,
                                  parse_uniprot_entry, prioritized_accessions)


class FoldDiscoAnnotationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.entry = {
            "primaryAccession": "P23528", "entryType": "reviewed",
            "uniProtKBCrossReferences": [
                {"database": "GO", "id": "GO:0003779", "properties": [
                    {"key": "GoTerm", "value": "F:actin binding"},
                    {"key": "GoEvidenceType", "value": "IDA:UniProt"}]},
                {"database": "GO", "id": "GO:0030042", "properties": [
                    {"key": "GoTerm", "value": "P:actin filament depolymerization"},
                    {"key": "GoEvidenceType", "value": "IEA:UniProt"}]}]}

    def save(self, entries):
        path = self.root / CACHE_PATH
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"schema_version": 1, "entries": entries}))

    def test_accession_matching_does_not_resolve_protein_names_or_pdbs(self):
        self.assertEqual(exact_accession("afdb", "P23528"), "P23528")
        self.assertEqual(exact_accession("afdb", "A0A0K0IYL8"), "A0A0K0IYL8")
        for database, target in (("pdb", "P23528"), ("afdb", "cofilin"),
                                 ("afdb", "P23528-2"), ("pdb", "3j8a")):
            self.assertIsNone(exact_accession(database, target))

    def test_official_go_aspects_ids_and_evidence_are_preserved(self):
        result = parse_uniprot_entry(self.entry, "2026-10-01", "https://rest.uniprot.org/example", "sha")
        self.assertEqual(result["status"], "annotated")
        self.assertEqual(result["terms"][0], dict(id="GO:0003779", aspect="F", term="actin binding", evidence="IDA:UniProt"))
        self.assertEqual(result["terms"][1]["evidence"], "IEA:UniProt")
        self.assertEqual(result["response_sha256"], "sha")

    def test_annotation_keeps_rows_and_cannot_transfer_by_identical_names(self):
        entry = parse_uniprot_entry(self.entry, "today", "url", "sha")
        self.save({"P23528": entry})
        frame = pd.DataFrame({"db": ["afdb", "afdb", "pdb"],
                              "target_id": ["P23528", "P60981", "3j8a"],
                              "name": ["Cofilin"] * 3, "idfscore": [5, 4, 3]}, index=[9, 2, 7])
        result = annotate_hits(frame, self.root)
        self.assertEqual(result.index.tolist(), [9, 2, 7])
        pd.testing.assert_frame_equal(result[frame.columns], frame)
        self.assertIn("GO:0003779", result.iloc[0].go_molecular_function)
        self.assertIn("IDA:UniProt", result.iloc[0].go_evidence)
        self.assertEqual(result.iloc[1].go_molecular_function, "")
        self.assertIn("Not fetched", result.iloc[1].go_annotation_status)
        self.assertIn("unresolved", result.iloc[2].go_annotation_status)

    def test_redirected_primary_id_is_never_assigned_to_requested_hit(self):
        entry = parse_uniprot_entry(self.entry, "today", "url", "sha")
        self.save({"P60981": entry})
        frame = pd.DataFrame({"db": ["afdb"], "target_id": ["P60981"]})
        result = annotate_hits(frame, self.root).iloc[0]
        self.assertEqual(result.go_annotation_status, "No exact UniProt entry returned")
        self.assertEqual(result.go_molecular_function, "")

    def test_no_terms_and_not_yet_fetched_have_distinct_statuses(self):
        entry = parse_uniprot_entry({"primaryAccession": "P23528"}, "today", "url", "sha")
        self.save({"P23528": entry})
        result = annotate_hits(pd.DataFrame({"db": ["afdb", "afdb"],
                                             "target_id": ["P23528", "P60981"]}), self.root)
        self.assertEqual(result.go_annotation_status.tolist(),
                         ["no GO terms returned", "Not fetched (partial GO cache)"])

    def test_priority_distributes_first_hits_across_query_motifs(self):
        rows = [dict(db="afdb", target_id=accession, query_abp=abp, query_cluster="site",
                     is_source=source, idfscore=score)
                for accession, abp, score, source in [
                    ("P23528", "A", 100, False), ("P60981", "A", 90, False),
                    ("P60709", "B", 10, False), ("P68032", "B", 200, True)]]
        self.assertEqual(prioritized_accessions(pd.DataFrame(rows)), ["P23528", "P60709", "P60981"])


if __name__ == "__main__":
    unittest.main()
