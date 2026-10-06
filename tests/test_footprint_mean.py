import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
import streamlit

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from footprint_comparison import mean_footprint


class FootprintMeanTests(unittest.TestCase):
    def sources(self, rows):
        return pd.DataFrame(rows, columns=['kind', 'group', 'interaction_id', 'chain', 'pdb_id'])

    def residues(self, rows):
        return pd.DataFrame(rows, columns=['interaction_id', 'chain', 'residue_number_canon_mafft', 'buried_ASA_percent'])

    def test_noncontacts_count_zero_but_unresolved_and_invalid_are_missing(self):
        sources = self.sources([('abp', 'A', i, f'{i}_A', str(i)) for i in range(1, 5)])
        residues = self.residues([(1, '1_A', 3, 40), (1, '1_A', 6, 'unknown')])
        resolved = {f'{i}_A': {1, 2} for i in range(1, 5)}
        result = mean_footprint(sources, residues, resolved, 'abp', ['A'])
        self.assertEqual(result.loc[1, 'mean_ASA_percent'], 10)
        self.assertEqual(result.loc[1, 'PDB_count'], 4)
        self.assertEqual(result.loc[2, 'mean_ASA_percent'], 0)
        self.assertEqual(result.loc[2, 'PDB_count'], 3)
        self.assertTrue(np.isnan(result.loc[3, 'mean_ASA_percent']))
        self.assertEqual(result.loc[3, 'PDB_count'], 0)

    def test_pdb_weight_is_not_increased_by_repeated_chains_or_interactions(self):
        sources = self.sources([('abp', 'A', 1, 'x_A', 'x'), ('abp', 'A', 2, 'x_B', 'x'),
                                ('abp', 'A', 3, 'y_A', 'y'), ('abp', 'A', 4, 'x_A', 'x')])
        residues = self.residues([(1, 'x_A', 3, 40), (2, 'x_B', 3, 40), (4, 'x_A', 3, 40)])
        result = mean_footprint(sources, residues, {c: {1} for c in sources.chain}, 'abp', ['A'])
        self.assertEqual(result.loc[1, 'mean_ASA_percent'], 20)
        self.assertEqual(result.loc[1, 'PDB_count'], 2)
        self.assertEqual(result.loc[1, 'chain_count'], 3)

    def test_selected_sites_form_chain_union_without_diluting_by_other_sites(self):
        sources = self.sources([('homo', 's1', 1, 'x_A', 'x'), ('homo', 's2', 2, 'x_A', 'x'),
                                ('homo', 's1', 3, 'y_A', 'y')])
        residues = self.residues([(1, 'x_A', 3, 40), (2, 'x_A', 3, 20), (2, 'x_A', 6, 50)])
        result = mean_footprint(sources, residues, {'x_A': {1, 2}, 'y_A': {1, 2}}, 'homo', ['s1', 's2'])
        self.assertEqual(result.loc[1, 'mean_ASA_percent'], 20)
        self.assertEqual(result.loc[2, 'mean_ASA_percent'], 25)


if __name__ == '__main__':
    unittest.main()
