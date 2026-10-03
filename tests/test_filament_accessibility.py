from pathlib import Path
import sys
import unittest

import numpy as np
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.SASA import ShrakeRupley

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "script"))
from filament_accessibility import exact_reference_map, protein_heavy_model, residue_profile


def residue(name, number, coordinates):
    r = Residue((" " if name != "HIC" else "H_HIC", number, " "), name, "")
    for i, (label, element, coord) in enumerate(coordinates):
        r.add(Atom(label, np.array(coord, float), 0, 1, " ", label, i, element=element))
    return r


class FilamentAccessibilityTests(unittest.TestCase):
    def test_exact_sequence_mapping_does_not_assume_author_numbers(self):
        self.assertEqual(exact_reference_map("ACDE", "MACDE"), {0: 2, 1: 3, 2: 4, 3: 5})
        with self.assertRaises(ValueError):
            exact_reference_map("ACDF", "MACDE")
        with self.assertRaises(ValueError):
            exact_reference_map("AA", "AAA")

    def test_neighbor_occlusion_and_heavy_atom_selection(self):
        model = Model(0)
        for chain_id, x in [("I", 0), ("J", 3)]:
            chain = Chain(chain_id)
            chain.add(residue("GLY", 100, [("CA", "C", (x, 0, 0)), ("H", "H", (x, 0, 1))]))
            model.add(chain)
        isolated = protein_heavy_model(model, ["I"])
        pair = protein_heavy_model(model, ["I", "J"])
        self.assertEqual(sum(1 for _ in isolated.get_atoms()), 1)
        atoms = list(pair.get_atoms())
        self.assertEqual(len({a.full_id for a in atoms}), len(atoms))
        self.assertTrue(all(a.full_id == a.get_full_id() for a in atoms))
        sr = ShrakeRupley(n_points=100)
        sr.compute(isolated, level="R")
        sr.compute(pair, level="R")
        self.assertAlmostEqual(isolated["I"][100].sasa, 4 * np.pi * (1.7 + 1.4) ** 2)
        self.assertLess(pair["I"][100].sasa, isolated["I"][100].sasa)

    def test_missing_modified_and_incomplete_rsa_remain_nan(self):
        standard = residue("GLY", 104, [(name, element, (i, 0, 0)) for i, (name, element) in enumerate([("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")])])
        modified = residue("HIC", 105, [("CA", "C", (0, 0, 0))])
        incomplete = residue("ALA", 106, [("CA", "C", (0, 0, 0))])
        areas = {"isolated": {standard.id: 52.0, modified.id: 10.0, incomplete.id: 20.0},
                 "actin_fragment": {standard.id: 26.0, modified.id: 5.0, incomplete.id: 10.0},
                 "with_abp": {standard.id: 13.0, modified.id: 2.0, incomplete.id: 8.0}}
        frame = residue_profile("MGHA", [standard, modified, incomplete], {0: 2, 1: 3, 2: 4}, areas, "TEST", "I")
        self.assertEqual(frame.coordinate_status.tolist(), ["missing_coordinates", "complete_standard_residue", "modified_residue", "incomplete_atoms"])
        self.assertEqual(frame.loc[1, "rsa_isolated"], .5)
        self.assertEqual(frame.loc[1, "delta_rsa_actin"], .25)
        self.assertTrue(frame.loc[[0, 2, 3], "rsa_isolated"].isna().all())
        self.assertEqual(frame.loc[2, "sasa_isolated_A2"], 10.0)

    def test_acetyl_cap_occludes_solvent_without_getting_a_residue_rsa(self):
        model = Model(0)
        chain = Chain('B')
        model.add(chain)
        cap = residue('ACE', 0, [('CH3', 'C', (3, 0, 0))])
        gly = residue('GLY', 1, [('CA', 'C', (0, 0, 0))])
        chain.add(cap)
        chain.add(gly)
        without = protein_heavy_model(model, ['B'])
        capped = protein_heavy_model(model, ['B'], include_acetyl_caps=True)
        self.assertEqual(len(without['B']), 1)
        self.assertEqual(len(capped['B']), 2)
        sr = ShrakeRupley(n_points=100)
        sr.compute(without, level='R')
        sr.compute(capped, level='R')
        self.assertLess(capped['B'][1].sasa, without['B'][1].sasa)

    def test_separately_modelled_modification_excludes_standard_normalization(self):
        gly = residue('GLY', 1, [(name, element, (i, 0, 0)) for i, (name, element) in enumerate(
            [('N', 'N'), ('CA', 'C'), ('C', 'C'), ('O', 'O')])])
        frame = residue_profile('MG', [gly], {0: 2}, {'isolated': {gly.id: 52.0}},
                                'TEST', 'B', modified_positions={2})
        self.assertEqual(frame.loc[1, 'coordinate_status'], 'modified_residue')
        self.assertTrue(np.isnan(frame.loc[1, 'rsa_isolated']))
        self.assertEqual(frame.loc[1, 'sasa_isolated_A2'], 52.0)


if __name__ == "__main__":
    unittest.main()
