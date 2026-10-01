import unittest
import sys
from pathlib import Path
import streamlit
import pandas as pd
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'script'))
from footprint_comparison import footprint_records,jaccard

class FootprintTests(unittest.TestCase):
    def test_jaccard_uses_union_and_empty_union_is_undefined(self):
        self.assertAlmostEqual(jaccard({1,2},{2,3}),1/3)
        self.assertEqual(jaccard({1},set()),0)
        self.assertTrue(np.isnan(jaccard(set(),set())))

    def test_homo_sides_use_their_own_cluster_and_exact_chain(self):
        meta=pd.DataFrame({'subunit_1':['p_A'],'subunit_2':['p_B'],'s1_actine':[True],'s2_actine':[True],
                           'subunit_1_title':['Actin'],'subunit_2_title':['Actin'],
                           's1_binding_site_cluster_data_70':['6685_1'],'s2_binding_site_cluster_data_70':['6685_2']})
        interactions=pd.DataFrame({'interaction_id':[1],'chain_A_id':['p_A'],'chain_B_id':['p_B']})
        residues=pd.DataFrame({'interaction_id':[1,1,1],'chain':['p_A','p_B','p_a'],
                               'residue_number_canon_mafft':[6,380,20],'buried_ASA_percent':[10,20,90]})
        result=footprint_records(meta,interactions,residues)
        self.assertEqual(result.set_index('group').position.to_dict(),{'6685_1':2,'6685_2':375})
