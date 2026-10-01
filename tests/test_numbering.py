import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# loaded by path: putting script/ on sys.path would shadow the streamlit package
_spec = importlib.util.spec_from_file_location("numbering", ROOT / "script/numbering.py")
numbering = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(numbering)
ALN = ROOT / "data/alignments/cluster_6685.aln"
REF = ROOT / "data/P60709_ref.fasta"
CONS = ROOT / "data/proteocast/conservation_vs_asa_per_position.csv"


class NumberingTest(unittest.TestCase):
    def test_round_trip_covers_whole_actin(self):
        for pos in range(1, numbering.ACTIN_LENGTH + 1):
            self.assertEqual(numbering.to_uniprot(numbering.to_canon(pos)), pos)

    def test_insertions_have_no_residue(self):
        for col in (1, 2, 4, 5, 237, 381):
            self.assertIsNone(numbering.to_uniprot(col))
        self.assertEqual(numbering.label(237), "ins. (col. 237)")
        self.assertEqual(numbering.to_uniprot(380), 375)

    @unittest.skipUnless(ALN.exists() and REF.exists(), "public dataset not present")
    def test_rule_matches_p60709_row_of_alignment(self):
        from Bio import AlignIO
        ref = "".join(l.strip() for l in REF.read_text().splitlines() if not l.startswith(">"))
        rows = [r for r in AlignIO.read(str(ALN), "fasta") if str(r.seq).replace("-", "") == ref]
        self.assertTrue(rows, "no P60709 sequence in the actin alignment")
        pos = 0
        for col, aa in enumerate(str(rows[0].seq), start=1):
            if aa == "-":
                self.assertIsNone(numbering.to_uniprot(col), col)
            else:
                pos += 1
                self.assertEqual(numbering.to_uniprot(col), pos, col)

    @unittest.skipUnless(CONS.exists(), "public dataset not present")
    def test_rule_matches_proteocast_table(self):
        import pandas as pd
        cons = pd.read_csv(CONS).dropna(subset=["canon", "Residue"])
        cons = cons[cons["canon"] >= 6]   # N-terminal acidic residues are not homologous
        self.assertTrue((cons["canon"].map(numbering.to_uniprot) == cons["Residue"]).all())


if __name__ == "__main__":
    unittest.main()
