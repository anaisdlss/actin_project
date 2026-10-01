"""Result availability and diagnostics must survive retries and app reruns."""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import pandas as pd
import streamlit

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from proteocast_results import result_file, record_job, job_status, missing_result_reason
import proteocast_view


class ProteoCastResultTests(unittest.TestCase):
    def test_empty_download_does_not_count_as_a_result(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp) / 'Cofilin_1'
            folder.mkdir()
            (folder / '2.aliAF-P23528-F1-msa_v6.fasta').touch()
            (folder / '4.query_ProteoCast.csv').touch()
            self.assertIsNone(result_file(tmp, 'Cofilin_1'))
            self.assertIn('alignment (MSA) is empty', missing_result_reason(tmp, 'Cofilin_1'))

    def test_failure_keeps_job_identity_and_can_be_replaced_by_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            record_job(tmp, 'Cofilin_1', 'submitted', job_id='example', uniprot='P23528')
            record_job(tmp, 'Cofilin_1', 'failed', 'Server download failed: HTTP 403')
            saved = job_status(tmp, 'Cofilin_1')
            self.assertEqual(saved['job_id'], 'example')
            self.assertEqual(saved['uniprot'], 'P23528')
            self.assertIn('HTTP 403', missing_result_reason(tmp, 'Cofilin_1'))
            record_job(tmp, 'Cofilin_1', 'submitted', job_id='retry')
            (Path(tmp) / 'Cofilin_1.csv').write_text('Mutation,Variant_score\nM1A,-1\n')
            record_job(tmp, 'Cofilin_1', 'finished')
            self.assertEqual(missing_result_reason(tmp, 'Cofilin_1'), '')
            self.assertEqual(job_status(tmp, 'Cofilin_1')['job_id'], 'retry')

    def test_manifest_and_new_results_are_read_on_each_rerun(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / 'manifest.csv'
            pd.DataFrame([dict(abp_title='A', slug='A')]).to_csv(manifest, index=False)
            with patch.object(proteocast_view, '_MANIFEST', str(manifest)), \
                 patch.object(proteocast_view, '_ABP_DIR', str(root)):
                self.assertFalse(proteocast_view.load_status(0).iloc[0].fait)
                (root / 'A.csv').write_text('Mutation,Variant_score\nM1A,-1\n')
                pd.DataFrame([dict(abp_title='A', slug='A'),
                              dict(abp_title='B', slug='B')]).to_csv(manifest, index=False)
                latest = proteocast_view.load_status(0)
                self.assertEqual(len(latest), 2)
                self.assertEqual(latest.fait.tolist(), [True, False])
                self.assertEqual(latest.iloc[0].diagnostic, '')
                self.assertIn('No UniProt identifier', latest.iloc[1].diagnostic)


if __name__ == '__main__':
    unittest.main()
