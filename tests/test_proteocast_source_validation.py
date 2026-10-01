"""Source-contract tests independent of the installed scientific dataset."""
import sys
import tempfile
import unittest
from pathlib import Path

import streamlit  # Import before adding script/ (which contains streamlit.py).
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from scientific_analysis import load_actin_scores


AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'


class ProteoCastSourceValidationTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        base = self.root / 'data/proteocast/actin'
        base.mkdir(parents=True)
        self.reference = 'M' + (AMINO_ACIDS * 19)[:374]
        fasta = '>synthetic_375_residue_reference\n' + self.reference + '\n'
        (self.root / 'data/P60709_ref.fasta').write_text(fasta)
        (base / '1.query.fasta').write_text(fasta)
        self.csv = base / '4.query_ProteoCast.csv'
        self.rows = pd.DataFrame([
            {
                'Mutation': f'{reference}{position}{alternate}',
                'Residue': position,
                # Self substitutions deliberately differ from all 19 changes.
                'Variant_score': 0.0 if alternate == reference else -20.0,
            }
            for position, reference in enumerate(self.reference, 1)
            for alternate in AMINO_ACIDS
        ])

    def write_source(self, rows):
        rows.to_csv(self.csv, index=False)

    def test_complete_finite_grid_loads_all_reference_positions(self):
        self.write_source(self.rows)
        raw, scores = load_actin_scores(self.root)
        self.assertEqual(len(raw), 375 * 20)
        self.assertEqual(scores.position.tolist(), list(range(1, 376)))
        self.assertEqual(''.join(scores.aa), self.reference)
        self.assertTrue(raw.groupby('position').size().eq(20).all())
        self.assertTrue(np.isfinite(scores.mean_variant_score).all())

    def test_mean_includes_the_unchanged_amino_acid_score(self):
        self.write_source(self.rows)
        _, scores = load_actin_scores(self.root)
        # (19 * -20 + 0) / 20 = -19, not the -20 mean of changes only.
        np.testing.assert_allclose(scores.mean_variant_score, -19.0)
        np.testing.assert_allclose(scores.sensitivity, 19.0)

    def test_missing_mutation_is_rejected(self):
        self.write_source(self.rows.drop(index=0))
        with self.assertRaises(ValueError):
            load_actin_scores(self.root)

    def test_invalid_or_nonfinite_scores_are_rejected(self):
        for value in ('invalid', np.nan, np.inf, -np.inf):
            with self.subTest(score=value):
                rows = self.rows.copy()
                rows['Variant_score'] = rows.Variant_score.astype(object)
                rows.at[0, 'Variant_score'] = value
                self.write_source(rows)
                with self.assertRaises(ValueError):
                    load_actin_scores(self.root)


if __name__ == '__main__':
    unittest.main()
