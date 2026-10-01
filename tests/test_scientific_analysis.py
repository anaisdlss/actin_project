import sys,unittest
from pathlib import Path
import streamlit
import pandas as pd
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'script'))
from scientific_analysis import bh_adjust,isoform_map,interface_evidence,load_variants,load_actin_scores,disease_associations,canonical_conservation,clinvar_identity_status
ROOT=Path(__file__).resolve().parents[1]

class ScientificTests(unittest.TestCase):
    def test_clinvar_gene_must_match_title_even_when_search_gene_matches(self):
        row={'gene':'ACTA2','pos':4,'aa_ref':'I','aa_alt':'M'}
        wrong={'title':'NM_000043.6(FAS):c.12C>G (p.Ile4Met)','genes':[{'symbol':'ACTA2'}]}
        self.assertEqual(clinvar_identity_status(row,wrong),'different gene in source title')
        self.assertEqual(clinvar_identity_status(row,{'title':'NM_001613.4(ACTA2):c.12C>G (p.Ile4Met)'}),'verified')
        self.assertEqual(clinvar_identity_status(row,None),'not verified')

    def test_bh_preserves_missing_and_order(self):
        np.testing.assert_allclose(bh_adjust([.04,np.nan,.01,.03]),[.04,np.nan,.03,.04],equal_nan=True)

    def test_isoform_insertion_not_forced_onto_reference(self):
        self.assertEqual(isoform_map('MACDEFGHIK','MADEFGHIK'),{1:1,2:2,4:3,5:4,6:5,7:6,8:7,9:8,10:9})

    def test_evidence_deduplicates_pdb_and_includes_s2_site(self):
        data=pd.DataFrame({'pdb_id':['p','p','q','p'],'s1_actine':[True]*4,'s2_actine':[True,True,True,False],
            's1_binding_site_cluster_data_70':['a','a','a','x'],'s2_binding_site_cluster_data_70':['b','b','c','y'],
            'cluster_data_70':['h','h','h','het'],'subunit_1_title':['Actin']*4,'subunit_2_title':['Actin','Actin','Actin','Cofilin-1']})
        sites,c70,occ,matches=interface_evidence(data)
        s=sites.set_index('site');self.assertEqual(s.loc['a','PDB_count'],2);self.assertEqual(s.loc['b','PDB_count'],1)
        self.assertEqual(s.loc['b','context_present'],1);self.assertEqual(s.loc['c','other_present'],1)
        self.assertEqual(len(matches),1);self.assertEqual(c70.iloc[0].PDB_count,2)

    def test_disease_counts_positions_not_variants(self):
        v=pd.DataFrame({'gene':['ACTB']*4,'mapping_valid':[True]*4,'source':['clinvar']*4,
            'classif_cat':['pathogenic']*4,'position':[1,1,2,3],'maladies':['A','A','B','B']})
        fp=pd.DataFrame({'kind':['abp']*2,'group':['x']*2,'position':[1,2],'asa':[10,20]})
        r=disease_associations(v,fp,'ACTB','A').iloc[0]
        self.assertEqual(r.disease_positions_in_footprint,1);self.assertEqual(r.disease_positions_total,1)
        self.assertEqual(r.other_PL_positions_total,2)

    @unittest.skipUnless((ROOT/'data/human_variants/manifest.json').exists(),'variant source snapshot absent')
    def test_variant_reference_and_all_six_genes(self):
        v,m=load_variants(ROOT)
        self.assertEqual(v.gene.nunique(),6)
        self.assertFalse(v.loc[~v.reference_matches,'mapping_valid'].any())
        self.assertTrue(v.loc[v.mapping_valid,'position'].between(1,375).all())
        self.assertTrue(m[m.gene.eq('ACTB')].gene_position.eq(m[m.gene.eq('ACTB')].position).all())
        self.assertTrue(v.loc[v.source.eq('clinvar'),'clinvar_id'].notna().all())
        population=v[v.source.eq('gnomad')]
        self.assertTrue(population.population_dataset.notna().all())
        self.assertFalse(population.loc[population.ac.eq(0),'analysis_eligible'].any())

    @unittest.skipUnless((ROOT/'data/proteocast/actin/4.query_ProteoCast.csv').exists(),'ProteoCast source absent')
    def test_proteocast_all_positions_and_m1(self):
        raw,s=load_actin_scores(ROOT)
        self.assertEqual(len(s),375);self.assertEqual(s.iloc[0].aa,'M')
        self.assertEqual(raw.groupby('position').size().min(),20)
        self.assertAlmostEqual(s.iloc[0].sensitivity,-raw[raw.position.eq(1)].Variant_score.mean())
        canonical=canonical_conservation(ROOT)
        self.assertEqual(canonical.loc[canonical.Residue.eq(1),'canon'].iloc[0],3)
        self.assertFalse(canonical.loc[canonical.Residue.eq(1),'at_interface'].iloc[0])
        self.assertAlmostEqual(canonical.iloc[0].conservation,s.iloc[0].sensitivity)
