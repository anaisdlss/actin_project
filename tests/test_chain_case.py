import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IFACE = ROOT / "data/filtered/details/3.interface_residues.csv"


@unittest.skipUnless(IFACE.exists(), "public dataset not present")
class ChainCaseTest(unittest.TestCase):
    """PDB chain identifiers are case-sensitive: in 9y52, chain A is actin and
    chain a is cofilin. Cofilin residues must never be counted on the actin side."""

    def test_actin_side_excludes_lowercase_partner_chain(self):
        import os
        sys.path.append(str(ROOT / "script"))   # appended: keep the streamlit package first
        cwd = os.getcwd()
        os.chdir(ROOT)
        try:
            from residue_passport import build_passport, pp_mtimes
            pp = build_passport(pp_mtimes())
        finally:
            os.chdir(cwd)
        import pandas as pd
        iface = pd.read_csv(IFACE)
        interactions = pd.read_csv(ROOT / "data/filtered/details/1.interactions.csv")
        matching = interactions[(interactions["chain_A_id"] == "9y52_A") &
                                (interactions["chain_B_id"] == "9y52_a")]
        self.assertFalse(matching.empty, "The chain-case regression requires PDB 9y52 A/a")
        interaction_id = matching.iloc[0]["interaction_id"]
        cofilin = iface[(iface["interaction_id"] == interaction_id) & (iface["chain"] == "9y52_a")]
        self.assertFalse(cofilin.empty)
        actin_rows = pp["res_long"][pp["res_long"]["interaction_id"] == interaction_id]
        actin = iface[(iface["interaction_id"] == interaction_id) & (iface["chain"] == "9y52_A")]
        self.assertEqual(len(actin_rows), len(actin))


if __name__ == "__main__":
    unittest.main()
