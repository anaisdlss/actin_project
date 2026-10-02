import json
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

import streamlit
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'script'))
from folddisco_audit import (motif_catalog,coverage_tables,normalize_discovery,
                            residue_token,validate_motif,read_chain_ca,prepared_query_json)


class FoldDiscoAuditTests(unittest.TestCase):
    def test_discovery_cache_refreshes_results_and_resolved_names(self):
        import folddisco_view
        with tempfile.TemporaryDirectory() as directory:
            result_path=Path(directory)/'results.csv';name_path=Path(directory)/'names.csv'
            result=pd.DataFrame([dict(query_abp='A',query_cluster='s1',db='afdb',target_id='P12345',
                                     target_chain='',is_source=False,idfscore=10,nodecount=4,coverage=1,rmsd=1)])
            names=pd.DataFrame([dict(key='P12345',protein_name='First',organism='Test organism')])
            result.to_csv(result_path,index=False);names.to_csv(name_path,index=False)
            with patch.object(folddisco_view,'_DISCO_CSV',str(result_path)),patch.object(folddisco_view,'_DISCO_NAMES',str(name_path)):
                folddisco_view._load_discovery.clear()
                first=folddisco_view._load_discovery(folddisco_view._discovery_stamp())
                self.assertEqual(first.iloc[0]['name'],'First')
                names.loc[0,'protein_name']='Changed protein name';names.to_csv(name_path,index=False)
                second=folddisco_view._load_discovery(folddisco_view._discovery_stamp())
                self.assertEqual(second.iloc[0]['name'],'Changed protein name')
                result.loc[0,'idfscore']=250;result.to_csv(result_path,index=False)
                third=folddisco_view._load_discovery(folddisco_view._discovery_stamp())
                self.assertEqual(third.iloc[0].idfscore,250)
                folddisco_view._load_discovery.clear()

    def test_catalog_uses_exact_interaction_and_case_sensitive_chain(self):
        sites=pd.DataFrame([
            ('Protein','s1',1,'1abc','A','4.0 Å'),
            ('Protein','s1',2,'1abc','a','2.0 Å'),
        ],columns=['abp_title','actin_site_cluster','interaction_id','pdb','abp_chain','resolution'])
        residues=pd.DataFrame([(1,'1abc_A','7'),(2,'1abc_A','99'),(2,'1abc_a','10'),
                               (2,'1abc_a','10A'),(2,'1abc_a','11.0'),(2,'1abc_a','11.2')],
                              columns=['interaction_id','chain','residue_number_structure'])
        row=motif_catalog(sites,residues).iloc[0]
        self.assertEqual(row.source_chain,'a')
        self.assertEqual(row.interaction_id,2)
        self.assertEqual(row.contact_positions,'10, 10A, 11')
        self.assertEqual(row.unparsed_contact_positions,'11.2')
        self.assertEqual(row.reconstructed_size,3)

    def test_normalization_reports_best_hit_fallback_and_keeps_ratios_above_one(self):
        rows=[('a','s','pdb',True,10),('a','s','pdb',False,20),
              ('a','s','afdb',False,8),('a','s','afdb',False,4),
              ('b','s','pdb',True,0),('b','s','pdb',False,np.nan)]
        frame=pd.DataFrame(rows,columns=['query_abp','query_cluster','db','is_source','idfscore'])
        result=normalize_discovery(frame)
        self.assertEqual(result.loc[1,'score_norm'],2)
        self.assertEqual(result.loc[1,'normalization_denominator'],10)
        self.assertIn('not verified self-match',result.loc[1,'normalization_basis'])
        self.assertEqual(result.loc[2,'score_norm'],1)
        self.assertEqual(result.loc[3,'score_norm'],.5)
        self.assertEqual(result.loc[2,'normalization_basis'],'best available hit')
        self.assertTrue(result.loc[4:,'score_norm'].isna().all())

    def test_preparation_validates_coordinates_and_preserves_insertion_codes(self):
        accepted=validate_motif('10, 10, 10A, 11',['10','10A','11'],'a')
        self.assertEqual(accepted['positions'],['10','10A','11'])
        self.assertEqual(accepted['errors'],[])
        self.assertEqual(accepted['duplicate_tokens_removed'],1)
        self.assertFalse(accepted['text_export_supported'])
        self.assertIsNone(accepted['query'])
        plain=validate_motif('-1, 2, 3',['-1','2','3'],'a')
        self.assertIsNone(plain['query'])
        invalid=validate_motif('1,2,99,2.5',['1','2','3'],'A')
        self.assertEqual(len(invalid['errors']),2)
        self.assertIsNone(invalid['query'])
        self.assertIsNone(residue_token('nan'))

    def test_numeric_chain_export_is_unambiguous_and_retains_coordinates(self):
        from folddisco_jobs import prepare
        from folddisco_status import query_problem, saved_searches, status_inventory
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'source.pdb'
            records=[]
            for serial,(chain,position) in enumerate([('A',99),('9',33),('9',34),('9',37)],1):
                records.append(f'ATOM  {serial:5d}  CA  ALA {chain}{position:4d}    {float(serial):8.3f}{0.:8.3f}{0.:8.3f}  1.00 20.00           C  \n')
            path.write_text(''.join(records)+'END\n')
            row=dict(source_chain='9',source_pdb='1abc',query_abp='Adducin',query_cluster='s1')
            folder,record=prepare(row,'33,34,37',path,root=Path(directory)/'jobs')
            self.assertEqual(record['motif'],'A33,A34,A37')
            self.assertEqual(record['submitted_chain'],'A')
            self.assertEqual(record['source']['source_chain'],'9')
            coords=read_chain_ca(folder/'query.pdb','A')
            self.assertEqual(coords.position.tolist(),['33','34','37'])
            self.assertEqual(coords.x.tolist(),[2.,3.,4.])
            self.assertIsNone(query_problem(record))
            legacy={**record,'state':'complete','hit_rows':0};legacy.pop('submitted_chain')
            self.assertIn('not a negative search',query_problem(legacy))
            (folder/'job.json').write_text(json.dumps(legacy))
            loaded=saved_searches(Path(directory)/'jobs')
            inventory=status_inventory(pd.DataFrame([dict(query_abp='Adducin',query_cluster='s1')]),pd.DataFrame(),loaded)
            self.assertEqual(inventory.iloc[0]['Invalid old queries'],1)
            self.assertTrue(pd.isna(inventory.iloc[0]['Tracked alignments']))

    def test_structure_reading_and_export_do_not_lose_chain_case_or_insertion_code(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'query.pdb'
            records=[]
            for serial,(chain,position,insertion) in enumerate([('A',99,' '),('a',10,' '),('a',10,'A'),('a',11,' ')],1):
                records.append(f'ATOM  {serial:5d}  CA  ALA {chain}{position:4d}{insertion}   {float(serial):8.3f}{0.:8.3f}{0.:8.3f}  1.00 20.00           C  \n')
            path.write_text(''.join(records)+'END\n')
            coords=read_chain_ca(path,'a')
            self.assertEqual(coords.position.tolist(),['10','10A','11'])
            checked=validate_motif('10,10A,11',coords.position,'a')
            row=dict(query_abp='Protein',query_cluster='s1',source_pdb='1abc',source_chain='a',interaction_id=2,contact_positions='10, 10A, 11')
            exported=json.loads(prepared_query_json(row,checked,path))
            self.assertEqual(exported['edited_positions'],['10','10A','11'])
            self.assertIn('not submitted',exported['status'])
            self.assertEqual(len(exported['structure_sha256']),64)

    def test_coverage_distinguishes_missing_rows_and_different_query_sizes(self):
        catalog=pd.DataFrame([
            ('A','s1','1abc','A',1,2.,'1, 2, 3',3,''),
            ('A','s2','1abc','A',2,2.,'4, 5, 6',3,''),
        ],columns=['query_abp','query_cluster','source_pdb','source_chain','interaction_id','resolution_A','contact_positions','reconstructed_size','unparsed_contact_positions'])
        discovery=pd.DataFrame([('A','s1','pdb',4)],columns=['query_abp','query_cluster','db','query_size'])
        with tempfile.TemporaryDirectory() as directory:
            summary,details=coverage_tables(catalog,discovery,pd.DataFrame(),Path(directory))
        self.assertEqual(summary.iloc[0].current_site_motifs,2)
        self.assertEqual(summary.iloc[0].sites_with_saved_PDB_hits,1)
        self.assertEqual(summary.iloc[0].motifs_without_saved_results,1)
        self.assertEqual(summary.iloc[0].motifs_with_local_structure,0)
        self.assertIn('Different query size',details.loc[details.query_cluster.eq('s1'),'query_size_check'].iloc[0])
        self.assertEqual(details.loc[details.query_cluster.eq('s2'),'query_size_check'].iloc[0],'No saved result rows')


if __name__=='__main__':unittest.main()
