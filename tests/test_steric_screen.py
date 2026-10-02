import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from steric_screen import rigid_proximity


class RigidProximityTests(unittest.TestCase):
    def test_partner_follows_anchor_fit_and_residues_are_counted_once(self):
        anchor = {1: np.array([0., 0., 0.]), 2: np.array([1., 0., 0.]),
                  3: np.array([0., 1., 0.]), 4: np.array([0., 0., 1.])}
        rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        translation = np.array([20., -10., 3.])
        reference = {p: xyz @ rotation + translation for p, xyz in anchor.items()}
        partner = np.array([[2., 0., 0.], [2.2, 0., 0.], [8., 0., 0.]])
        neighbours = np.array([[3., 0., 0.]]) @ rotation + translation
        result = rigid_proximity(anchor, reference, partner, ['10', '10', '11A'], neighbours, minimum=4)
        self.assertAlmostEqual(result['anchor_CA_RMSD_A'], 0.)
        self.assertAlmostEqual(result['minimum_neighbor_distance_A'], .8)
        self.assertEqual(result['ABP_atoms_below_1p5_A'], 2)
        self.assertEqual(result['ABP_residues_below_1p5_A'], 1)
        self.assertEqual(result['ABP_resolved_residues'], 2)

    def test_inadequate_anchor_and_missing_neighbours_are_rejected(self):
        anchor = {1: [0., 0., 0.]}
        with self.assertRaises(ValueError):
            rigid_proximity(anchor, anchor, [[1., 0., 0.]], ['2'], [[3., 0., 0.]])
        with self.assertRaises(ValueError):
            rigid_proximity(anchor, anchor, [[1., 0., 0.]], ['2'], [], minimum=1)
