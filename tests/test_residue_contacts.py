import sys
import unittest
from pathlib import Path
import streamlit
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "script"))
from residue_contacts import orient_actin_contacts

class BidirectionalContactsTests(unittest.TestCase):
    def test_homo_has_both_sides_and_chain_case_distinguishes_abp(self):
        contacts = pd.DataFrame([
            [1,"p_A", "10A", "D", 6, "p_B", "20", "F", 380, 1.2, None, 20, 40],
            [2,"p_a", "30", "K", 70, "p_A", "10A", "D", 6, 2., "hydrogen", 50, 30],
        ], columns=["interaction_id","chain_A_id","residue_A_structure","residue_A_name",
                    "residue_A_canon_mafft","chain_B_id","residue_B_structure","residue_B_name",
                    "residue_B_canon_mafft","contact_area","contact_type","asa_pct_A","asa_pct_B"])
        proteins = pd.DataFrame({"chain":["p_A","p_B","p_a"],"protein":["Actin","Actin","ABP"],
                                 "is_actin":[True,True,False]})
        result = orient_actin_contacts(contacts,proteins)
        self.assertEqual(len(result),3)
        homo = result[result.interaction_id.eq(1)]
        self.assertEqual(set(homo["P60709 position"]),{2,375})
        self.assertEqual(set(homo["Interaction type"]),{"Actin–actin"})
        hetero = result[result.interaction_id.eq(2)].iloc[0]
        self.assertEqual(hetero["Interaction type"],"Actin–ABP")
        self.assertEqual(hetero["Actin PDB residue"],"10A")
        self.assertEqual(hetero["Partner buried ASA (%)"],50)
        self.assertEqual(hetero["Actin buried ASA (%)"],30)
