import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'script'))
from local_rebuild import Task, run_task, status, run, REPORT, inventory, build_lock


class LocalRebuildTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        (self.root/'input.csv').write_text('value\n1\n')
        (self.root/'calculate.py').write_text(
            "from pathlib import Path\nPath('result.csv').write_text(Path('input.csv').read_text())\n")
        self.task = Task('example', [sys.executable, 'calculate.py'], ['input.csv'],
                         ['result.csv'], ['calculate.py'], packages=())

    def test_source_and_code_contents_invalidate_even_when_mtime_unchanged(self):
        run_task(self.root, self.task, emit=lambda _: None)
        self.assertEqual(status(self.root, self.task)[0], 'current')
        original = (self.root/'input.csv').stat()
        (self.root/'input.csv').write_text('value\n2\n')
        os.utime(self.root/'input.csv', ns=(original.st_atime_ns, original.st_mtime_ns))
        self.assertEqual(status(self.root, self.task)[0], 'needs rebuild')
        run_task(self.root, self.task, emit=lambda _: None)
        (self.root/'calculate.py').write_text((self.root/'calculate.py').read_text()+'# changed\n')
        self.assertEqual(status(self.root, self.task)[0], 'needs rebuild')

    def test_failed_run_restores_last_good_outputs_and_never_certifies_success(self):
        run_task(self.root, self.task, emit=lambda _: None)
        old = (self.root/'result.csv').read_bytes()
        (self.root/'calculate.py').write_text("from pathlib import Path\nPath('result.csv').write_text('partial')\nraise RuntimeError('network stopped')\n")
        with patch('local_rebuild.registry', return_value=[self.task]):
            with self.assertRaises(RuntimeError):
                run(self.root, emit=lambda _: None)
        self.assertEqual((self.root/'result.csv').read_bytes(), old)
        self.assertEqual(json.loads((self.root/REPORT/'last_run.json').read_text())['state'], 'failed')
        self.assertNotEqual(status(self.root, self.task)[0], 'current')

    def test_changed_or_missing_output_is_not_current(self):
        run_task(self.root, self.task, emit=lambda _: None)
        (self.root/'result.csv').write_text('tampered')
        self.assertEqual(status(self.root, self.task)[0], 'needs rebuild')
        (self.root/'result.csv').unlink()
        self.assertNotEqual(status(self.root, self.task)[0], 'current')

    def test_unknown_csv_origin_is_not_invented(self):
        (self.root/'data').mkdir()
        (self.root/'data/orphan.csv').write_text('rsa\n0.8\n')
        records = inventory(self.root, [])
        self.assertEqual(records[0]['status'], 'unresolved historical provenance')
        self.assertEqual(records[0]['producer'], '')

    def test_concurrent_builds_cannot_write_the_same_outputs(self):
        with build_lock(self.root):
            with self.assertRaises(RuntimeError):
                with build_lock(self.root):
                    pass


class RsaSourceTests(unittest.TestCase):
    def test_absent_documented_calculation_never_uses_legacy_values(self):
        from rsa_source import load_rsa
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/'data/proteocast').mkdir(parents=True)
            (root/'data/proteocast/conservation_vs_asa_per_position.csv').write_text('canon,rsa\n33,0.017667\n')
            frame, message = load_rsa(root)
            self.assertTrue(frame.empty)
            self.assertIn('rebuild', message)

    def test_current_dataset_uses_documented_rsa_and_preserves_missing_positions(self):
        from rsa_source import load_rsa
        from scientific_analysis import load_actin_scores
        root = Path(__file__).resolve().parents[1]
        if not (root/'data/proteocast/actin/4.query_ProteoCast.csv').exists():
            self.skipTest('local research sources absent')
        frame, state = load_rsa(root)
        self.assertEqual(state, 'current')
        _, scores = load_actin_scores(root)
        self.assertTrue(scores.rsa.equals(scores.rsa_isolated))
        a29 = scores.set_index('position').loc[29]
        self.assertEqual(a29.rsa, frame.set_index('position').loc[29, 'rsa_isolated'])
        self.assertTrue(scores.set_index('position').loc[1, ['rsa_isolated','rsa_actin_fragment','rsa_with_abp']].isna().all())


class PipelineSafetyTests(unittest.TestCase):
    def module(self, name, relative):
        path = Path(__file__).resolve().parents[1]/relative
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def test_an_existing_family_table_does_not_hide_a_failed_network_step(self):
        module = self.module('legacy_pipeline_test', 'script/abp_site_domain/run_all.py')
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            (out/'familles.csv').write_text('old output')
            def step(name):
                if name == '03_interpro.py':
                    module._failed.append(name)
            with patch.object(module, 'OUT', out), patch.object(module, '_py', side_effect=step), \
                 patch.object(module, '_sh'), patch.object(module.subprocess, 'run'):
                with self.assertRaises(SystemExit) as error:
                    module.main()
            self.assertNotEqual(error.exception.code, 0)
            self.assertEqual(json.loads((out/'pipeline_incomplete.json').read_text())['state'], 'failed')
            self.assertEqual((out/'familles.csv').read_text(), 'old output')

    def test_changed_alignment_is_blocked_instead_of_guessing_new_numbers(self):
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
        module = self.module('preflight_test', 'tools/rebuild_local.py')
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/'data/alignments').mkdir(parents=True)
            (root/'data/P60709_ref.fasta').write_text('>ref\nM\n')
            (root/'data/alignments/cluster_6685.aln').write_text('>ref\n---M\n')
            with patch('check_dataset.audit', return_value={'checks': []}):
                with self.assertRaisesRegex(ValueError, 'mapping changed'):
                    module.preflight(root)


if __name__ == '__main__':
    unittest.main()
