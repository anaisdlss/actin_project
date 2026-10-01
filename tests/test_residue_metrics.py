import sys
import unittest
from pathlib import Path
import pandas as pd
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "script"))
from residue_metrics import complete_actin_positions, surface_masks, pdb_residue_frequencies
import numbering

class ResidueMetricsTests(unittest.TestCase):
    def test_absent_rsa_is_not_exposed_or_buried(self):
        exposed, buried = surface_masks([1,2,3,4,5,6], {1:.2,2:.1,3:None,4:np.nan,5:-1})
        self.assertEqual(exposed.tolist(), [True,False,False,False,False,False])
        self.assertEqual(buried.tolist(), [False,True,False,False,False,False])

    def test_all_375_positions_available_without_fabricating_values(self):
        result = complete_actin_positions(pd.DataFrame({"canon":[380,237],"rsa":[.8,.5]}))
        mapped = result["canon"].map(numbering.to_uniprot).dropna()
        self.assertEqual(set(mapped), set(range(1,376)))
        self.assertEqual(result.loc[result.canon.eq(380),"rsa"].iloc[0], .8)
        self.assertTrue(pd.isna(result.loc[result.canon.eq(3),"rsa"].iloc[0]))

    def test_pdbs_not_chains_define_frequency_and_mixed_pdb_counts_twice(self):
        observations = pd.DataFrame({"pdb_id":["a","a","a","b",None],
                                     "residue_name":["D","D","E","D","E"]})
        result,total = pdb_residue_frequencies(observations)
        self.assertEqual(total,2)
        self.assertEqual(result.set_index("residue_name")["pct"].to_dict(), {"D":100.,"E":50.})

    def test_chain_must_belong_to_this_interaction_not_another(self):
        from residue_metrics import select_interaction_chains
        residues = pd.DataFrame({"interaction_id":[1,1,2,2], "chain":["A","B","B","C"],
                                 "asa":[10,90,20,80]})
        selected = select_interaction_chains(residues,{1:"A",2:"B"})
        self.assertEqual(selected["asa"].tolist(),[10,20])
