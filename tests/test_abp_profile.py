"""Scientific contracts for ABP sensitivity and sequence/contact coordinates."""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import streamlit

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "script"))
from abp_profile import (AMINO_ACIDS, load_abp_scores, score_sequence,
                         annotate_profile, profile_figure, nonempty_alignment)
import proteocast_view


class ABPProfileTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.folder = self.base / "Synthetic"
        self.folder.mkdir()
        self.csv = self.folder / "4.query_ProteoCast.csv"
        self.query = self.folder / "1.query.fasta"
        self.sequence = "MAG"
        self.rows = pd.DataFrame([
            dict(Mutation=f"{ref}{pos}{alt}", Residue=pos,
                 Variant_score=0.0 if alt == ref else -20.0)
            for pos, ref in enumerate(self.sequence, 1) for alt in AMINO_ACIDS])
        self.write(self.rows)

    def write(self, rows):
        rows.to_csv(self.csv, index=False)

    def test_complete_mean_includes_self_substitution_and_full_length(self):
        _, frame = load_abp_scores(self.csv)
        self.assertEqual(frame.position.tolist(), [1, 2, 3])
        np.testing.assert_allclose(frame.mutational_sensitivity, 19.0)
        self.assertEqual(score_sequence(self.csv), self.sequence)
        self.assertTrue(frame.finite_score_count.eq(20).all())
        self.assertIn("reconstructed", frame.sequence_source.iloc[0])

    def test_missing_substitution_and_nonfinite_scores_are_not_imputed(self):
        for value in (None, np.nan, np.inf, "invalid"):
            rows = self.rows.copy()
            if value is None:
                rows = rows.drop(index=0)
            else:
                rows["Variant_score"] = rows.Variant_score.astype(object)
                rows.at[0, "Variant_score"] = value
            with self.subTest(value=value):
                self.write(rows)
                _, frame = load_abp_scores(self.csv)
                self.assertTrue(pd.isna(frame.mutational_sensitivity.iloc[0]))
                self.assertEqual(frame.finite_score_count.iloc[0], 19)
                np.testing.assert_allclose(frame.mutational_sensitivity.iloc[1:], 19.0)
                with self.assertRaises(ValueError):
                    score_sequence(self.csv)

    def test_query_preserves_positions_absent_from_scores(self):
        self.query.write_text(">query\nMAGW\n")
        self.write(self.rows[self.rows.Residue != 2])
        _, frame = load_abp_scores(self.csv)
        self.assertEqual(frame.position.tolist(), [1, 2, 3, 4])
        self.assertEqual("".join(frame.reference_aa), "MAGW")
        self.assertTrue(frame.loc[frame.position.isin([2, 4]), "mutational_sensitivity"].isna().all())

    def test_duplicate_reference_and_query_contradictions_are_rejected(self):
        duplicate = pd.concat([self.rows, self.rows.iloc[:1]])
        conflicting = self.rows.copy()
        conflicting.at[0, "Mutation"] = "A1A"
        bad_position = self.rows.copy()
        bad_position.at[0, "Residue"] = 2
        for rows in (duplicate, conflicting, bad_position):
            self.write(rows)
            with self.assertRaises(ValueError):
                load_abp_scores(self.csv)
        self.write(self.rows)
        self.query.write_text(">query\nMAW\n")
        with self.assertRaises(ValueError):
            load_abp_scores(self.csv)

    def test_flat_result_layout_checks_its_own_query(self):
        flat = self.base / "Synthetic.csv"
        self.rows.to_csv(flat, index=False)
        self.query.write_text(">query\nMAG\n")
        _, frame = load_abp_scores(flat)
        self.assertEqual(frame.sequence_source.iloc[0], "submitted ProteoCast query FASTA")
        self.query.write_text(">query\nMAW\n")
        with self.assertRaises(ValueError):
            load_abp_scores(flat)

    def test_query_fallback_refuses_incomplete_or_conflicting_results(self):
        with patch.object(proteocast_view, "_ABP_DIR", str(self.base)):
            self.assertEqual(proteocast_view._query_seq("Synthetic"), "MAG")
            self.write(self.rows.drop(index=0))
            self.assertIsNone(proteocast_view._query_seq("Synthetic"))
            self.write(self.rows)
            self.query.write_text(">query\nMAW\n")
            self.assertIsNone(proteocast_view._query_seq("Synthetic"))

    def test_missing_contact_asa_and_surface_rsa_stay_unknown(self):
        _, frame = load_abp_scores(self.csv)
        domains = [dict(name="Domain", db="pfam", spans=[(1, 2)])]
        annotated = annotate_profile(frame, iface_asa={1: 37.5, 2: np.inf},
                                     domains=domains, rsa={1: 0.3, 2: 0.1}, surface_only=True)
        self.assertEqual(annotated.displayed.tolist(), [True, False, False])
        self.assertEqual(annotated.buried_ASA_percent_max.iloc[0], 37.5)
        self.assertTrue(annotated.buried_ASA_percent_max.iloc[1:].isna().all())
        self.assertEqual(annotated.domains.tolist(), ["Domain (pfam)", "Domain (pfam)", ""])

    def test_invalid_new_scores_do_not_reuse_old_cached_contacts(self):
        source = self.base / "all.csv"
        contacts = self.base / "contacts.csv"
        pd.DataFrame([dict(subunit_2="test_A", s2_sequence="MAG",
                           subunit_2_title="Synthetic")]).to_csv(source, index=False)
        pd.DataFrame([dict(chain="test_A", residue_number_sequence=2,
                           buried_ASA_percent="37.5%")]).to_csv(contacts, index=False)
        with patch.object(proteocast_view, "_ABP_DIR", str(self.base)), \
             patch.object(proteocast_view, "_ALL", str(source)), \
             patch.object(proteocast_view, "_IFACE", str(contacts)):
            self.assertEqual(proteocast_view.abp_interface_asa_on_query("Synthetic", "Synthetic", 0), {2: 37.5})
            self.write(self.rows.drop(index=0))
            self.assertEqual(proteocast_view.abp_interface_asa_on_query("Synthetic", "Synthetic", 0), {})
            self.assertEqual(proteocast_view.abp_rsa_on_query("Synthetic", "Synthetic", 0), {})

    def test_figure_keeps_full_sequence_and_exact_track_coordinates(self):
        _, frame = load_abp_scores(self.csv)
        domains = [dict(name="Domain", db="pfam", spans=[(1, 3)])]
        frame = annotate_profile(frame, iface_asa={1: 2.5}, domains=domains)
        fig = profile_figure(frame, domains=domains)
        self.assertEqual(fig.layout.hovermode, "x unified")
        self.assertEqual({trace.xaxis for trace in fig.data}, {"x3"})
        self.assertEqual(tuple(fig.layout.xaxis3.range), (0.5, 3.5))
        for trace in fig.data:
            self.assertEqual(list(trace.x), [1, 2, 3])
        self.assertEqual(list(fig.data[2].y), ["Domain (pfam)"] * 3)
        self.assertTrue(np.isnan(fig.data[1].y[1]))
        self.assertFalse(fig.data[0].connectgaps)

    def test_only_nonempty_multisequence_alignment_is_offered(self):
        path = self.folder / "2.alignment.fasta"
        for text in ("", ">query\nMAG\n", ">query\nMAG\n>empty\n"):
            path.write_text(text)
            self.assertIsNone(nonempty_alignment(self.folder))
        path.write_text(">query\nMAG-\n>homolog\nMA-G\n")
        self.assertEqual(nonempty_alignment(self.folder), path)


if __name__ == "__main__":
    unittest.main()
