"""Regression checks for exact site scope, AA mixtures and structural sampling."""
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "script"))
from msa_contact_stats import oriented_pairs, scoped_contacts, conditional_area_profile


def metadata(a, b, site_a, site_b, actin_a=True, actin_b=False):
    return dict(subunit_1=a, subunit_2=b, s1_binding_site_cluster_data_70=site_a,
                s2_binding_site_cluster_data_70=site_b, s1_actine=actin_a, s2_actine=actin_b,
                s1_sequence="MAG", s2_sequence="KR", subunit_1_title="Actin",
                subunit_2_title="Partner")


def contact(iid, a, b, area="10", pos_a=3, pos_b=4):
    return dict(interaction_id=iid, chain_A_id=a, chain_B_id=b,
                residue_A_structure="12A", residue_B_structure="17B",
                residue_A_canon_mafft=pos_a, residue_B_canon_mafft=pos_b,
                residue_A_name="G", residue_B_name="K", contact_area=area,
                asa_pct_A=30, asa_pct_B=40, contact_type="")


class MSAContactStatisticsTests(unittest.TestCase):
    def setUp(self):
        self.metadata = pd.DataFrame([
            metadata("1abc_A", "1abc_p", "S0", "P"),
            metadata("1abc_B", "1abc_p", "OTHER", "P"),
            metadata("2ABC_p", "2ABC_a", "P", "S0", False, True),
            metadata("3abc_A", "3abc_a", "S3", "S4", True, True)])
        self.interactions = pd.DataFrame([
            dict(interaction_id=i, chain_A_id=a, chain_B_id=b)
            for i, a, b in [(1,"1abc_A","1abc_p"),(2,"1abc_B","1abc_p"),
                             (3,"2abc_a","2abc_p"),(4,"3abc_A","3abc_a")]])
        self.contacts = pd.DataFrame([
            contact(1,"1ABC_A","1abc_p"), contact(2,"1abc_B","1abc_p",area="999"),
            contact(3,"2abc_p","2abc_a"), contact(4,"3abc_A","3abc_a"),
            contact(999,"1abc_A","1abc_p",area="999")])

    def test_site_scope_uses_interaction_and_both_chains_not_partner_alone(self):
        pairs = oriented_pairs(self.metadata, self.interactions)
        rows = scoped_contacts(self.contacts, pairs, site="S0")
        self.assertEqual(set(rows.interaction_id), {1, 3})
        self.assertNotIn(2, rows.interaction_id.tolist())
        self.assertNotIn(999, rows.interaction_id.tolist())
        reverse = rows[rows.interaction_id == 3].iloc[0]
        self.assertEqual(reverse.chain_A_id, "2abc_a")
        self.assertEqual(reverse.residue_A_structure, "17B")
        self.assertEqual(reverse.canon_a, 4)
        self.assertEqual(reverse.asa_pct_A, 40)

    def test_homo_second_side_and_case_sensitive_chain_identifiers_are_kept(self):
        pairs = oriented_pairs(self.metadata, self.interactions)
        rows = scoped_contacts(self.contacts, pairs, site="S4")
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows.iloc[0].actin_chain, "3abc_a")
        self.assertEqual(rows.iloc[0].partner_chain, "3abc_A")
        self.assertTrue(rows.iloc[0].partner_is_actin)
        both = scoped_contacts(self.contacts, pairs)
        self.assertEqual(len(both[both.interaction_id == 4]), 2)

    def test_censored_invalid_and_nonpositive_areas_are_unavailable(self):
        pairs = oriented_pairs(self.metadata, self.interactions)
        for value in ("<0.1", "bad", np.nan, np.inf, 0, -1):
            with self.subTest(value=value):
                rows = scoped_contacts(pd.DataFrame([contact(1,"1abc_A","1abc_p",value)]), pairs, site="S0")
                self.assertEqual(len(rows), 1)
                self.assertTrue(pd.isna(rows.area_f.iloc[0]))

    def test_chain_then_pdb_means_do_not_zero_fill_absent_positions(self):
        rows = pd.DataFrame([
            dict(seq_low="seq", label="Partner", pdb=pdb, chain_B_id=chain,
                 chain_A_id=chain, canon_b=position, canon_a=position,
                 aa_b=aa, aa_a=aa, area_f=area)
            for pdb, chain, position, aa, area in [
                ("1abc","1abc_A",3,"K",4), ("1abc","1abc_A",3,"R",6),
                ("1abc","1abc_B",3,"K",30), ("2abc","2abc_C",3,"R",100),
                ("1abc","1abc_A",4,"G",20)]])
        for side in ("a", "b"):
            profile, aa = conditional_area_profile(rows, side)
            by_position = profile.set_index("canon_" + side)
            self.assertEqual(by_position.loc[3,"corr_mean"], 60)
            self.assertEqual(by_position.loc[4,"corr_mean"], 20)
            self.assertEqual(by_position.loc[3,"observed_pdbs"], 2)
            self.assertEqual(by_position.loc[3,"observed_chains"], 3)
            self.assertEqual(by_position.loc[4,"observed_pdbs"], 1)
            self.assertEqual(by_position.loc[3,"pct"], 75)
            mixture = aa[aa["canon_" + side].eq(3)].set_index("aa_" + side).corr_mean
            self.assertEqual(mixture.to_dict(), {"K":8.5,"R":51.5})
            self.assertEqual(mixture.sum(), by_position.loc[3,"corr_mean"])


if __name__ == "__main__":
    unittest.main()
