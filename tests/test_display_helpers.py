"""Display regressions: numeric IDs, true colours, and structure-score provenance."""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import streamlit
import numpy as np
import pandas as pd
import plotly.graph_objects as go
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from display_helpers import ordered_interactions, cell_color_hover, viewer_html
from proteocast_view import score_structure, structure_for
import proteocast_view
from abp_profile import AMINO_ACIDS


def atom(residue, position, confidence=90):
    return f'ATOM  {position:5d}  CA  {residue:3s} A{position:4d}    {0.:8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{confidence:6.2f}           C'


class DisplayTests(unittest.TestCase):
    def test_numeric_order_is_stable_and_does_not_erase_unusual_identifiers(self):
        frame = pd.DataFrame({'interaction_id': ['1000','100','2','100',None], 'value': [1,2,3,4,5]})
        result = ordered_interactions(frame)
        self.assertEqual(result.value.tolist(), [3,2,4,1,5])
        self.assertEqual(str(result.interaction_id.dtype), 'Int64')
        result = ordered_interactions(pd.DataFrame({'interaction_id':['unknown','10','2']}))
        self.assertEqual(result.interaction_id.tolist(), ['2','10','unknown'])
        self.assertEqual(ordered_interactions(pd.DataFrame({'interaction_id':pd.Series([10,None,2],dtype='Int64')})).interaction_id.iloc[0],2)

    def test_hover_uses_cell_colour_and_preserves_measurements_and_metadata(self):
        fig=go.Figure(go.Heatmap(z=[[0,50,100]], customdata=[['a','b','c']],
            colorscale=[[0,'#ffffff'],[1,'#ff0000']],zmin=0,zmax=100,
            hovertemplate='%{customdata}: %{z}<extra></extra>'))
        cell_color_hover(fig)
        self.assertIn('#ffffff',fig.data[0].hovertext[0][0])
        self.assertIn('#ff8080',fig.data[0].hovertext[0][1])
        self.assertIn('#ff0000',fig.data[0].hovertext[0][2])
        np.testing.assert_equal(fig.data[0].z,[[0,50,100]])
        self.assertEqual(fig.data[0].customdata[0][1],'b')
        before=fig.to_json();cell_color_hover(fig);self.assertEqual(before,fig.to_json())

    def test_binary_colours_follow_shared_reversed_scale(self):
        fig=go.Figure(go.Heatmap(z=[[0,1]],coloraxis='coloraxis'))
        fig.update_layout(coloraxis=dict(cmin=0,cmax=1,colorscale=[[0,'#ffffff'],[1,'#0072b2']],reversescale=True))
        cell_color_hover(fig)
        self.assertIn('#0072b2',fig.data[0].hovertext[0][0]);self.assertIn('#ffffff',fig.data[0].hovertext[0][1])

    def test_responsive_observer_targets_actual_viewer_element(self):
        import py3Dmol,re
        html=viewer_html(py3Dmol.view(width='100%',height=470))
        element=re.search(r'id="(3dmolviewer_[^"]+)"',html).group(1)
        self.assertIn(f"getElementById('{element}')",html)
        self.assertIn('ResizeObserver',html)
        self.assertIn('clientWidth === lastWidth',html)

    def test_structure_scores_require_complete_grid_and_matching_reference(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'scores.csv'
            rows=pd.DataFrame([dict(Mutation=f'M1{aa}',Residue=1,Variant_score=-2.) for aa in AMINO_ACIDS])
            rows.to_csv(path,index=False)
            text,low,high=score_structure(atom('MET',1),path)
            self.assertEqual(float(text[60:66]),2.)
            self.assertEqual((low,high),(2.,2.))
            with self.assertRaises(ValueError):score_structure(atom('GLY',1),path)
            rows.iloc[:-1].to_csv(path,index=False)
            with self.assertRaises(ValueError):score_structure(atom('MET',1),path)

    def test_plain_alphafold_bfactor_is_confidence_not_sensitivity(self):
        structure_for.clear()
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)/'Synthetic';folder.mkdir()
            (folder/'10.AF-X-F1-model_v6.pdb').write_text(atom('MET',1,confidence=91.))
            with patch.object(proteocast_view,'_ABP_DIR',tmp), patch.object(proteocast_view,'_fetch_af_pdb') as fetch:
                result=structure_for('Synthetic','X',1)
                self.assertFalse(result['metric']);self.assertEqual((result['vmin'],result['vmax']),(0,100))
                self.assertIn('pLDDT',result['by']);fetch.assert_not_called()
                # A newly imported score file must invalidate the cached confidence view.
                pd.DataFrame([dict(Mutation=f'M1{aa}',Residue=1,Variant_score=-2.)
                              for aa in AMINO_ACIDS]).to_csv(folder/'4.query_ProteoCast.csv',index=False)
                updated=structure_for('Synthetic','X',1)
                self.assertTrue(updated['metric'])
                self.assertEqual(float(updated['pdb'][60:66]),2.)
        structure_for.clear()

    def test_region_bounds_remain_available_for_small_proteins(self):
        from proteocast_abp import _abp_actin_focus
        with patch.object(proteocast_view,'abp_interface_asa_on_query',return_value={5:30,60:20}), patch.object(proteocast_view,'_query_seq',return_value='A'*70):
            self.assertEqual(_abp_actin_focus('small','small'),(1,70,70))
        with patch.object(proteocast_view,'abp_interface_asa_on_query',return_value={}):
            self.assertIsNone(_abp_actin_focus('missing','missing'))

if __name__ == '__main__':unittest.main()
