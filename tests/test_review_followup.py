"""Scientific and remote-job regression checks for the mentor follow-up."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import Mock
import numpy as np
import pandas as pd
import streamlit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'script'))
from actin_contact_surface import contact_categories, atom_colors, CATEGORIES
from folddisco_jobs import combine_sites, refresh, submit, result_table, alignment_rows
from interface_geometry import compare_pair


class FollowupTests(unittest.TestCase):
    def test_surface_colors_use_mapped_atoms_and_leave_unmapped_gray(self):
        text = 'ATOM      4\nATOM      7\nHETATM    9\n'
        colors = atom_colors(text, {47: [4, 7]}, {47: 'Hetero — ABP only'})
        self.assertEqual(colors[4], CATEGORIES['Hetero — ABP only'])
        self.assertEqual(colors[7], colors[4])
        self.assertEqual(colors[9], CATEGORIES['No observed contact'])

    def test_oversized_motif_is_not_submitted_or_truncated(self):
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d)
            positions = list(map(str, range(1, 45)))
            (folder/'job.json').write_text(json.dumps({'state': 'prepared', 'positions': positions}))
            session = Mock()
            record = submit(folder, session)
            self.assertEqual(record['state'], 'unsupported')
            self.assertEqual(record['positions'], positions)
            session.post.assert_not_called()

    def test_types_need_positive_observations_and_union_both_partner_types(self):
        r=pd.DataFrame({'kind':['homo','abp','homo','abp','abp','homo'],
                        'position':[1,2,3,3,4,5],'asa':[10,5,2,4,0,np.nan]})
        categories=contact_categories(r)
        self.assertEqual(categories[1],'Homo — actin only')
        self.assertEqual(categories[2],'Hetero — ABP only')
        self.assertEqual(categories[3],'Mixed — both')
        self.assertEqual(categories[4],'No observed contact')
        self.assertEqual(categories[5],'No observed contact')
        self.assertEqual(len(categories),375)

    def test_combination_uses_one_chain_and_all_requested_sites(self):
        sites=pd.DataFrame([
            dict(abp_title='ABP',actin_site_cluster=c,pdb=p,abp_chain=ch,interaction_id=i,resolution=res)
            for c,p,ch,i,res in [('x','1abc','A',1,'3 Å'),('y','1abc','A',2,'3 Å'),('x','2abc','B',3,'1 Å')]])
        residues=pd.DataFrame({'interaction_id':[1,1,2,2,3], 'chain':['1abc_A']*4+['2abc_B'],
                               'residue_number_structure':[1,2,2,4,99]})
        combined=combine_sites(sites,residues,'ABP',['x','y'])
        self.assertEqual(combined['source_pdb'],'1abc')
        self.assertEqual(combined['contact_positions'],'1, 2, 4')
        self.assertEqual(combined['interaction_ids'],[1,2])
        sites.loc[sites.actin_site_cluster.eq('y'),'abp_chain']='a'
        with self.assertRaisesRegex(ValueError,'No single observed'):combine_sites(sites,residues,'ABP',['x','y'])

    def test_neighbour_is_not_independently_refitted(self):
        a={1:np.array([0.,0,0]),2:np.array([1.,0,0]),3:np.array([0.,1,0]),4:np.array([0.,0,1])}
        b={p:v+np.array([0,0,8]) for p,v in a.items()}
        shifted_a={p:v+np.array([12,8,2]) for p,v in a.items()}
        shifted_b={p:v+np.array([12,8,2]) for p,v in b.items()}
        self.assertAlmostEqual(compare_pair(shifted_a,shifted_b,a,b,minimum=4)['neighbor_CA_RMSD_A'],0)
        displaced={p:v+np.array([0,3,0]) for p,v in shifted_b.items()}
        self.assertAlmostEqual(compare_pair(shifted_a,displaced,a,b,minimum=4)['neighbor_CA_RMSD_A'],3)
        self.assertIsNone(compare_pair(a,b,a,b,minimum=300))

    def test_empty_completed_search_is_distinct_from_failure_and_not_resubmitted(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            record={'state':'pending','ticket':'existing','positions':['1','2','3']}
            (folder/'job.json').write_text(json.dumps(record))
            session=Mock()
            session.get.side_effect=[Mock(json=lambda:{'status':'COMPLETE'}),Mock(json=lambda:{'results':[{'db':'pdb','alignments':None}]})]
            result=refresh(folder,session)
            self.assertEqual(result['state'],'complete');self.assertEqual(result['hit_rows'],0)
            self.assertTrue(result_table(folder).empty)
            submit(folder,session);session.post.assert_not_called()

    def test_nested_api_alignments_and_transient_network_error(self):
        hit={'target':'1abc','nodecount':3}
        self.assertEqual(alignment_rows([[hit], [hit]]),[hit,hit])
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);(folder/'job.json').write_text(json.dumps({'state':'pending','ticket':'keep','positions':['1','2','3']}))
            session=Mock();session.get.return_value.json.side_effect=ValueError('not JSON')
            record=refresh(folder,session)
            self.assertEqual(record['state'],'pending');self.assertEqual(record['ticket'],'keep')

if __name__=='__main__':unittest.main()
