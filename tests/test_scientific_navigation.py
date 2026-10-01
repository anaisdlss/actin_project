"""Regression checks for the first reorganization, without external services."""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import pandas as pd
import streamlit  # Load the package before adding script/ (which contains streamlit.py).
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from app_sections import select_cluster_types
from residue_explorer import other_abps_by_position
import make_slim_deploy


class ScientificNavigationTests(unittest.TestCase):
    def test_mixed_sites_are_available_in_both_interface_views(self):
        table = pd.DataFrame({'patch': ['h', 'm', 'a'],
                              'interaction_type': ['homo', 'mixed', 'hetero']})
        self.assertEqual(set(select_cluster_types(table, 'homo').patch), {'h', 'm'})
        self.assertEqual(set(select_cluster_types(table, 'hetero').patch), {'m', 'a'})
        self.assertEqual(len(select_cluster_types(table, 'all')), 3)
        self.assertEqual(len(table), 3)

    def test_pair_details_find_partners_outside_the_selected_pair(self):
        contacts = pd.DataFrame({'canon': [12, 12, 12, 12, 13, 13],
                                 'abp': ['A', 'B', 'C', 'C', 'A', 'B']})
        others = other_abps_by_position(contacts, {'A', 'B'})
        self.assertEqual(others.loc[12], 'C')
        self.assertEqual(others.loc[13], '')

    def test_missing_sources_preserve_an_existing_deploy(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            old = root / 'deploy' / 'data'
            old.mkdir(parents=True)
            sentinel = old / 'preserve.csv'
            sentinel.write_text('existing results\n')
            with patch.object(make_slim_deploy, 'DST', old), patch.object(
                    make_slim_deploy, 'validate_source',
                    side_effect=lambda: make_slim_deploy_validate(root / 'missing')):
                with self.assertRaises(SystemExit):
                    make_slim_deploy.main()
            self.assertEqual(sentinel.read_text(), 'existing results\n')


make_slim_deploy_validate = make_slim_deploy.validate_source
if __name__ == '__main__':
    unittest.main()
