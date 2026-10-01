"""Variant displays must distinguish binary presence from substitution counts."""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'script'))
from variant_heatmaps import presence_matrix,presence_figure,substitution_count_trace


class VariantHeatmapTests(unittest.TestCase):
    def test_duplicate_records_cannot_create_a_presence_gradient(self):
        variants=pd.DataFrame({'position':[10,10,10,11],'aa_alt':['A','A','C','A']})
        original=variants.copy(deep=True)
        matrix=presence_matrix(variants)
        self.assertEqual(matrix.shape,(20,375))
        self.assertEqual(matrix.loc['A',10],1)
        self.assertEqual(matrix.loc['C',10],1)
        self.assertEqual(matrix.loc['D',10],0)
        self.assertEqual(set(np.unique(matrix.to_numpy())),{0,1})
        self.assertEqual(matrix[10].sum(),2)
        pd.testing.assert_frame_equal(variants,original)
        figure=presence_figure(matrix,'ACTA1','pathogenic')
        self.assertFalse(figure.data[0].showscale)
        self.assertEqual([trace.name for trace in figure.data[1:]],['Not recorded (0)','Recorded (1)'])
        self.assertTrue(all(trace.showlegend for trace in figure.data[1:]))

    def test_count_gradient_keeps_absolute_integer_counts(self):
        counts=pd.DataFrame([[0,1,2],[0,2,3]],index=['ACTA1','All genes (sum)'],columns=[1,2,3])
        trace=substitution_count_trace(counts)
        np.testing.assert_array_equal(trace.z,counts.values)
        self.assertEqual(trace.zmin,0)
        self.assertEqual(trace.zmax,3)
        self.assertEqual(list(trace.colorbar.tickvals),[0,1,2,3])
        self.assertIn('Distinct substitutions',trace.hovertemplate)
        self.assertEqual(trace.colorbar.title.text,'Substitutions<br>per position')

    def test_empty_category_stays_absent_without_an_invented_score(self):
        matrix=presence_matrix(pd.DataFrame(columns=['position','aa_alt']))
        self.assertEqual(matrix.shape,(20,375))
        self.assertEqual(matrix.to_numpy().sum(),0)
        figure=presence_figure(matrix,'ACTA1','benign')
        self.assertTrue(np.all(figure.data[0].customdata=='No record in selected category'))
        trace=substitution_count_trace(pd.DataFrame([[0]*375]))
        self.assertEqual(trace.zmax,1)
        self.assertEqual(np.asarray(trace.z).sum(),0)


if __name__=='__main__':unittest.main()
