import sys
import unittest
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "script"))
from numbering_view import chain_coordinates


class CorrespondenceTest(unittest.TestCase):
    def test_preserves_chain_case_insertions_and_distinct_pdb_numbers(self):
        df = pd.DataFrame([
            ["9y52_A", "10A", 10, "A", 237],
            ["9y52_A", "10A", 10, "A", 237],
            ["9y52_A", "11", 11, "F", 380],
            ["9y52_a", "11", 11, "F", 380],
        ], columns=["chain", "residue_number_structure", "residue_number_sequence",
                    "residue_name", "residue_number_canon_mafft"])
        result = chain_coordinates(df, "9y52_A")
        self.assertEqual(len(result), 2)
        self.assertEqual(result.iloc[0]["residue_number_structure"], "10A")
        self.assertTrue(pd.isna(result.iloc[0]["P60709 position"]))
        self.assertEqual(result.iloc[1]["P60709 position"], 375)
