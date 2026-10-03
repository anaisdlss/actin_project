"""Offline rebuilding must preserve source evidence and reject invalid model grids."""
import importlib.util
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'script'))


def module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT/'tools'/f'{name}.py')
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


class SourceReconstructionTests(unittest.TestCase):
    def test_missing_annotations_do_not_become_a_shared_fold_or_domain(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            out = root/'data/exports/abp_site_domain'
            out.mkdir(parents=True)
            pd.DataFrame({'abp_title':['A','B'], 'actin_site_cluster':['site','site']}).to_csv(out/'abp_site_table.csv',index=False)
            pd.DataFrame({'abp_title':['A','B'], 'pdb':['1abc','1abc'], 'abp_chain':['A','B']}).to_csv(out/'abp_representatives.csv',index=False)
            pd.DataFrame({'abp_title':['A','B'], 'uniprot':['Q00001','Q00002'],
                          'pfam_domains':['Domain',None], 'interpro_domains':['Domain',None]}).to_csv(out/'abp_interpro.csv',index=False)
            (out/'fold_cluster_cluster.tsv').write_text('A__1abc_A\tA__1abc_A\n')
            relative = 'script/abp_site_domain/04_crosstab.py'
            with contextlib.redirect_stdout(io.StringIO()):
                exec(compile((ROOT/relative).read_text(), relative, 'exec'),
                     {'__file__':str(root/relative), '__name__':'__main__'})
            result = pd.read_csv(out/'site_cluster_summary.csv').iloc[0]
            self.assertEqual(result.n_folds, 1)
            self.assertFalse(result.mono_fold)
            self.assertFalse(result.mono_domaine)
            self.assertFalse(result.fold_annotation_complete)
            self.assertFalse(result.domain_annotation_complete)
            import matplotlib.pyplot as plt
            plt.close('all')

    def test_chemistry_uses_actin_from_the_selected_interaction_in_either_orientation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            out = root/'data/exports/abp_site_domain'
            detail = root/'data/filtered/details'
            out.mkdir(parents=True)
            detail.mkdir(parents=True)
            pd.DataFrame([dict(abp_title='Example', pdb='1abc', abp_chain='X',
                               abp_subunit='1abc_X', interaction_id=2)]).to_csv(out/'abp_representatives.csv', index=False)
            pd.DataFrame([
                dict(subunit_1='1abc_A', subunit_2='1abc_X', s1_actine=True, s2_actine=False),
                dict(subunit_1='1abc_X', subunit_2='1abc_B', s1_actine=False, s2_actine=True),
            ]).to_csv(root/'data/filtered/filtered_all_data.csv', index=False)
            pd.DataFrame([
                dict(interaction_id=1, chain_A_id='1abc_A', chain_B_id='1abc_X', interface_area=10,
                     num_contacts=2, num_hbonds=1, num_salt_bridges=0),
                dict(interaction_id=2, chain_A_id='1abc_X', chain_B_id='1abc_B', interface_area=20,
                     num_contacts=3, num_hbonds=2, num_salt_bridges=1),
            ]).to_csv(detail/'1.interactions.csv', index=False)
            pd.DataFrame([
                dict(interaction_id=1, chain='1abc_A', residue_name='LYS'),
                dict(interaction_id=2, chain='1abc_B', residue_name='GLU'),
                dict(interaction_id=2, chain='1abc_X', residue_name='ARG'),
            ]).to_csv(detail/'3.interface_residues.csv', index=False)
            relative = 'script/abp_site_domain/16_interface_chemistry.py'
            with contextlib.redirect_stdout(io.StringIO()):
                exec(compile((ROOT/relative).read_text(), relative, 'exec'),
                     {'__file__': str(root/relative), '__name__': '__main__'})
            result = pd.read_csv(out/'interface_chemistry.csv').iloc[0]
            self.assertEqual(result.charge_actine_patch, -1)
            self.assertEqual(result.charge_nette, 1)
            self.assertEqual(result.interface_area, 20)

    def test_offline_annotations_preserve_original_cache_bytes_and_make_no_request(self):
        annotations = module('update_structure_annotations')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            base = root/'data/annotations'
            base.mkdir(parents=True)
            cache = base/'rcsb_structure_metadata.json'
            original = '{"records":{"1ABC":{"entry":{"rcsb_id":"1ABC"}}},"historical_note":"retain"}'
            cache.write_text(original)
            with patch.object(annotations, 'fetch_batch', side_effect=AssertionError('Offline calculation contacted network')), \
                 patch.object(annotations, 'retained_entries', return_value=pd.DataFrame({'pdb_id':['1ABC']})), \
                 patch.object(annotations, 'build_annotations', return_value=(pd.DataFrame({'pdb':['1ABC']}), pd.DataFrame({'entity':[1]}))):
                self.assertEqual(annotations.main(['--root', str(root), '--offline']), 0)
            self.assertEqual(cache.read_text(), original)
            receipt = json.loads((base/'manifest.json').read_text())
            self.assertTrue(receipt['offline_rebuild'])
            self.assertEqual(receipt['source_cache']['sha256'], annotations.digest(cache))

    def test_availability_validates_query_and_does_not_invent_a_job(self):
        audit = module('audit_proteocast_availability')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root/'data/proteocast/abp_inputs/manifest.csv'
            manifest.parent.mkdir(parents=True)
            pd.DataFrame([{'abp_title':'Example','slug':'Example','uniprot':'Q00000'}]).to_csv(manifest,index=False)
            folder = root/'data/proteocast/abp/Example'
            folder.mkdir(parents=True)
            scores = folder/'4.query_ProteoCast.csv'
            grid = pd.DataFrame([{'Mutation':f'A1{aa}','Residue':1,'Variant_score':-1.0}
                                 for aa in 'ACDEFGHIKLMNPQRSTVWY'])
            grid.to_csv(scores,index=False)
            frame = audit.availability(root).iloc[0]
            self.assertTrue(frame.scores_available)
            self.assertFalse(frame.submitted_query_checked)
            self.assertEqual(frame.job_id, '')
            (folder/'1.query.fasta').write_text('>query\nG\n')
            frame = audit.availability(root).iloc[0]
            self.assertFalse(frame.scores_available)
            self.assertIn('contradict', frame.diagnostic)
            (folder/'1.query.fasta').write_text('>query\nA\n')
            grid.iloc[:-1].to_csv(scores,index=False)
            frame = audit.availability(root).iloc[0]
            self.assertFalse(frame.scores_available)
            self.assertIn('Incomplete', frame.diagnostic)

    def test_missing_offline_sources_are_blocked_before_calculation(self):
        from local_rebuild import registry, status
        with tempfile.TemporaryDirectory() as directory:
            tasks = {t.name:t for t in registry()}
            for name in ('structure_annotations','abp_cross_tables','abp_families','proteocast_availability'):
                state, reason = status(Path(directory), tasks[name])
                self.assertEqual(state, 'blocked')
                self.assertIn('Missing required source', reason)

    def test_unknown_rsa_archive_is_not_certified_by_presence(self):
        from local_rebuild import inventory
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            file = root/'data/proteocast/abp/Example/rsa_values.csv'
            file.parent.mkdir(parents=True)
            file.write_text('residue,rsa\n1,87\n')
            row = inventory(root, [])[0]
            self.assertEqual(row['status'], 'external bundle; RSA recipe unresolved')
            self.assertIn('no replacement invented', row['automatic_action'])


if __name__ == '__main__':
    unittest.main()
