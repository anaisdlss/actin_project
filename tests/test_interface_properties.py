"""Numerical checks independent of the UI and actual structural sample."""
import sys
import unittest
from pathlib import Path
import streamlit
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from interface_properties import chemistry_summary, chemistry_records, dominant_classes, CLASSES


class ChemistryTests(unittest.TestCase):
    def fixture(self):
        # One PDB heavily duplicated relative to another must not dominate.
        rows = [['A','p1',10,20,'Hydrophobic',9,90] for _ in range(20)]
        rows += [['A','p2',10,20,'Polar',1,10], ['B','p3',10,20,'Acidic (D/E)',1,20],
                 ['B','p3',11,0,'Polar',1,20], ['B','p3',12,10,'Polar',np.nan,20]]
        return pd.DataFrame(rows,columns=['Partner','PDB','P60709 position','Actin buried ASA (%)',
                                          'Partner class','Pair contact area (Å²)','Partner buried ASA (%)']).assign(
                                              interaction_id=1, **{'Actin chain':'A','Partner chain':'B'})

    def test_equal_pdb_then_equal_abp_weight_and_missingness(self):
        profile, abp, used = chemistry_summary(self.fixture())
        p=profile.set_index('P60709 position')
        self.assertAlmostEqual(p.loc[10,'Hydrophobic'],.25)
        self.assertAlmostEqual(p.loc[10,'Polar'],.25)
        self.assertAlmostEqual(p.loc[10,'Acidic (D/E)'],.5)
        self.assertEqual(p.loc[10,'ABP source names'],2)
        self.assertTrue(p.loc[11,list(CLASSES)].isna().all())
        self.assertTrue(p.loc[12,list(CLASSES)].isna().all())
        self.assertEqual(len(profile),375)
        self.assertEqual(abp.set_index('Partner').loc['A','PDBs with usable contacts'],2)

    def test_pair_area_and_asa_weights_are_distinct(self):
        rows=self.fixture().iloc[[0,20]].copy()
        rows['PDB']='p1'
        profile,_,_=chemistry_summary(rows)
        self.assertAlmostEqual(profile.iloc[9]['Hydrophobic'],.9)
        profile,_,_=chemistry_summary(rows,'Residue-pair count')
        self.assertAlmostEqual(profile.iloc[9]['Hydrophobic'],.5)
        rows.loc[rows['Partner class'].eq('Hydrophobic'),'Partner buried ASA (%)']=30
        profile,_,_=chemistry_summary(rows,'Partner buried ASA (%)')
        self.assertAlmostEqual(profile.iloc[9]['Hydrophobic'],.75)

    def test_equal_observed_interfaces_before_pdb(self):
        rows=self.fixture().iloc[[0,20]].copy()
        rows['PDB']='p1'
        rows['interaction_id']=[1,2]
        # Nine times the area in one interface must not give it nine times the weight.
        profile,_,_=chemistry_summary(rows)
        self.assertAlmostEqual(profile.iloc[9]['Hydrophobic'],.5)
        self.assertAlmostEqual(profile.iloc[9]['Polar'],.5)

    def test_strict_threshold_and_empty_output(self):
        p,abp,used=chemistry_summary(self.fixture(),cutoff=20)
        self.assertTrue(used.empty)
        self.assertTrue(abp.empty)
        self.assertEqual(len(p),375)
        self.assertTrue(p['Dominant partner class'].eq('No usable contact').all())

    def test_ties_are_not_arbitrary_class(self):
        f=pd.DataFrame([[.5,.5],[np.nan,np.nan],[1,0]],columns=['a','b'])
        self.assertEqual(list(dominant_classes(f)),['Tied classes','No usable contact','a'])

    def test_chain_direction_case_and_duplicate_measurements(self):
        contacts=pd.DataFrame([[1,'p_a','20','K',10,'p_A','10A','D',6,2.,'hydrogen',50,30]],
            columns=['interaction_id','chain_A_id','residue_A_structure','residue_A_name',
                     'residue_A_canon_mafft','chain_B_id','residue_B_structure','residue_B_name',
                     'residue_B_canon_mafft','contact_area','contact_type','asa_pct_A','asa_pct_B'])
        proteins=pd.DataFrame({'chain':['p_A','p_a'],'protein':['Actin','ABP'],'is_actin':[True,False]})
        meta=pd.DataFrame({'subunit_1':['p_a'],'subunit_2':['p_A'],'s1_actine':[False],'s2_actine':[True],
                           's1_binding_site_cluster_data_70':['other'],'s2_binding_site_cluster_data_70':['6685_1']})
        ints=contacts[['interaction_id','chain_A_id','chain_B_id']]
        result=chemistry_records(pd.concat([contacts,contacts]),proteins,meta,ints)
        self.assertEqual(len(result),1)
        self.assertEqual(result.iloc[0]['Partner class'],'Basic (K/R/H)')
        self.assertEqual(result.iloc[0]['Binding site'],'6685_1')
        self.assertEqual(result.iloc[0]['P60709 position'],2)
        self.assertEqual(result.iloc[0]['Actin buried ASA (%)'],30)
        inconsistent=pd.concat([contacts,contacts.assign(contact_area=3.)])
        with self.assertRaises(ValueError): chemistry_records(inconsistent,proteins,meta,ints)


if __name__=='__main__': unittest.main()
