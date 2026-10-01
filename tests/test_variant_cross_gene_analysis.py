"""Checks for cross-gene cohort selection, footprint counts and ASA weighting."""
import sys
import unittest
from pathlib import Path

import streamlit
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from scientific_analysis import (abp_position_asa, cross_gene_observations,
                                 cross_gene_footprint_analysis, conflicting_variant_counts)


def variant(gene, position, reference, alternate, source, category='', **extra):
    return dict(gene=gene, position=position, pos=position, aa_ref=reference,
                aa_alt=alternate, source=source, classif_cat=category,
                analysis_eligible=True, mapping_valid=True, **extra)


class CrossGeneVariantTests(unittest.TestCase):
    def setUp(self):
        changes=[(10,'A','V'),(10,'A','G'),(20,'C','D'),(30,'D','E'),
                 (40,'E','F'),(50,'F','G'),(60,'G','H')]
        rows=[]
        for n,(position,reference,alternate) in enumerate(changes):
            rows.append(variant('ACTB',position,reference,alternate,'clinvar','pathogenic',
                                clinvar_id=n+1,classification='Pathogenic'))
            rows.append(variant('ACTG1',position,reference,alternate,'gnomad',
                                variant_id=f'population-{n}',population_dataset='exome'))
        # Multiple submissions/alleles do not change substitution or position counts.
        rows.append(dict(rows[0],clinvar_id=101))
        rows.append(dict(rows[1],variant_id='population-duplicate',population_dataset='genome'))
        # A VUS still counts as a ClinVar annotation in the population gene.
        rows.append(variant('ACTG1',50,'F','G','clinvar','uncertain',clinvar_id=201))
        rows.append(variant('ACTG1',10,'A','C','clinvar','uncertain',clinvar_id=202))
        self.variants=pd.DataFrame(rows)
        self.variants.loc[self.variants.source.eq('gnomad')&self.variants.position.eq(40),'analysis_eligible']=False
        self.variants.loc[self.variants.source.eq('gnomad')&self.variants.position.eq(60),'mapping_valid']=False
        # Duplicate records within one interface chain contribute only its maximum ASA.
        self.records=pd.DataFrame([
            ('abp','X','s1',1,'chain-A',10,10),
            ('abp','X','s1',1,'chain-A',10,20),
            ('abp','X','s1',1,'chain-A',10,20),
            ('abp','X','s2',2,'chain-B',10,60),
            ('abp','X','s3',3,'chain-C',20,100),
            ('abp','X','s4',4,'chain-D',90,90),
            ('abp','Y','s5',5,'chain-E',30,0),
            ('abp','Z','s6',6,'chain-F',99,30),
            ('homo','h','s7',7,'chain-G',30,50),
        ],columns=['kind','group','site','interaction_id','chain','position','asa'])

    def test_substitution_join_excludes_annotations_and_ineligible_records(self):
        before=self.variants.copy(deep=True)
        cross=cross_gene_observations(self.variants)
        self.assertEqual(set(zip(cross.position,cross.aa_alt)),{(10,'V'),(10,'G'),(20,'D'),(30,'E')})
        self.assertEqual(len(cross),4)
        self.assertTrue(cross.gene_PL.eq('ACTB').all())
        self.assertTrue(cross.gene_population.eq('ACTG1').all())
        first=cross[cross.position.eq(10)&cross.aa_alt.eq('V')].iloc[0]
        self.assertIn('101',first.ClinVar_ids)
        self.assertEqual(first.gnomAD_datasets,'exome; genome')
        pd.testing.assert_frame_equal(self.variants,before)

    def test_counts_and_asa_weight_positions_not_substitutions(self):
        stats=abp_position_asa(self.records)
        position=stats[stats.ABP.eq('X')&stats.position.eq(10)].iloc[0]
        self.assertEqual(position.mean_ASA_percent,40)
        self.assertEqual(position.observed_interface_chains,2)
        cross,pairs,summary,details=cross_gene_footprint_analysis(self.variants,self.records)
        pair=pairs.iloc[0]
        self.assertEqual(pair.shared_substitutions,4)
        self.assertEqual(pair.cross_gene_positions_total,3)
        self.assertEqual(pair.positions_with_any_ABP_contact,2)
        row=summary[summary.ABP.eq('X')].iloc[0]
        self.assertEqual(row.footprint_positions,3)
        self.assertEqual(row.shared_substitutions_in_footprint,3)
        self.assertEqual(row.cross_gene_positions_in_footprint,2)
        self.assertAlmostEqual(row.fraction_cross_gene_positions_contacted,2/3)
        self.assertAlmostEqual(row.fraction_footprint_positions_cross_gene,2/3)
        self.assertEqual(row.mean_ASA_at_cross_gene_positions,70)
        self.assertEqual(row.positions_contributing_to_ASA,2)
        self.assertEqual(row.observed_chain_position_records,3)
        self.assertEqual(row.contacted_P60709_positions,'10; 20')
        unobserved=details[details.position.eq(30)].iloc[0]
        self.assertFalse(unobserved.positive_ABP_contact_observed)
        self.assertTrue(pd.isna(unobserved.mean_ASA_percent))
        for name in ['Y','Z']:
            row=summary[summary.ABP.eq(name)].iloc[0]
            self.assertEqual(row.cross_gene_positions_in_footprint,0)
            self.assertTrue(pd.isna(row.mean_ASA_at_cross_gene_positions))
        self.assertTrue(pd.isna(summary.loc[summary.ABP.eq('Y'),'fraction_footprint_positions_cross_gene'].iloc[0]))

    def test_no_contacts_preserves_observations_and_missing_asa(self):
        cross,pairs,summary,details=cross_gene_footprint_analysis(self.variants,self.records.iloc[:0])
        self.assertEqual(len(cross),4)
        self.assertEqual(pairs.iloc[0].positions_with_any_ABP_contact,0)
        self.assertTrue(summary.empty)
        self.assertEqual(len(details),4)
        self.assertFalse(details.positive_ABP_contact_observed.any())
        self.assertTrue(details.mean_ASA_percent.isna().all())
        self.assertTrue(details.cross_gene_positions_total.eq(3).all())

    def test_empty_cross_gene_cohort_has_export_columns(self):
        only_clinvar=self.variants[self.variants.source.eq('clinvar')]
        cross,pairs,summary,details=cross_gene_footprint_analysis(only_clinvar,self.records)
        self.assertTrue(all(frame.empty for frame in [cross,pairs,summary,details]))
        self.assertIn('cross_gene_positions_total',pairs)
        self.assertIn('mean_ASA_at_cross_gene_positions',summary)
        self.assertIn('positive_ABP_contact_observed',details)

    def test_conflict_counts_keep_vus_unchanged_and_deduplicate_substitutions(self):
        rows=[variant('ACTB',10,'A','V','clinvar','conflicting'),
              variant('ACTB',10,'A','V','clinvar','conflicting'),
              variant('ACTB',10,'A','G','clinvar','conflicting'),
              variant('ACTB',11,'C','D','clinvar','uncertain'),
              variant('ACTB',12,'D','E','gnomad','conflicting'),
              variant('ACTB',13,'E','F','clinvar','conflicting')]
        frame=pd.DataFrame(rows);frame.loc[5,'analysis_eligible']=False
        before=frame.copy(deep=True)
        counts=conflicting_variant_counts(frame)
        self.assertEqual(counts.shape,(6,375))
        self.assertEqual(counts.loc['ACTB',10],2)
        self.assertEqual(int(counts.to_numpy().sum()),2)
        pd.testing.assert_frame_equal(frame,before)


if __name__=='__main__':
    unittest.main()
