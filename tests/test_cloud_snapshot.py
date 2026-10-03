import importlib.util
from pathlib import Path
import subprocess
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('sync_cloud_dataset', Path(__file__).resolve().parents[1]/'tools/sync_cloud_dataset.py')
cloud = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cloud)


class CloudSnapshotTests(unittest.TestCase):
    def test_duplicate_copy_requires_identical_original(self):
        with tempfile.TemporaryDirectory() as temporary:
            original=Path(temporary)/'structure.pdb'; duplicate=Path(temporary)/'structure 2.pdb'
            duplicate.write_text('coordinates')
            self.assertFalse(cloud.duplicate_copy(duplicate))
            original.write_text('coordinates')
            self.assertTrue(cloud.duplicate_copy(duplicate))
            duplicate.write_text('different coordinates')
            self.assertFalse(cloud.duplicate_copy(duplicate))

    def test_displayed_data_and_alignments_are_kept_but_indexes_are_not(self):
        for path in ('data/raw/pdb_entry_results.csv', 'data/proteocast/abp/example/2.msa.fasta',
                     'data/exports/folddisco_controls/hits.csv', 'data/exports/folddisco_jobs/123/query.pdb',
                     'data/filtered/details/structures_files/assembly/1abc.cif'):
            self.assertIsNone(cloud.omitted(path), path)
        for path in ('data/exports/folddisco_controls/raw/123.tsv',
                     'data/exports/abp_site_domain/folddisco_index/idx.offset',
                     'data/exports/abp_site_domain/_sweep_pad20/example.pdb'):
            self.assertIsNotNone(cloud.omitted(path), path)

    def test_cleanup_only_removes_tracked_work_files_and_keeps_unknown_results(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); source = root/'source'; destination = root/'public'
            destination.mkdir(); (source/'data').mkdir(parents=True)
            subprocess.run(['git','init','-q',str(destination)],check=True)
            tracked = 'data/exports/abp_site_domain/folddisco_index/idx.offset'
            untracked = 'data/exports/abp_site_domain/folddisco_index/user-note.txt'
            unknown = 'data/exports/previous_result.csv'
            for name in (tracked, untracked, unknown):
                path=destination/name; path.parent.mkdir(parents=True,exist_ok=True); path.write_text('keep')
            subprocess.run(['git','-C',str(destination),'add',tracked,unknown],check=True)
            _, removed=cloud.plan(source,destination)
            self.assertEqual(removed,[tracked])

    def test_conflicting_untracked_source_copy_stops_the_sync(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); source=root/'source'; destination=root/'public'
            for folder in (source,destination):
                (folder/'data').mkdir(parents=True); (folder/'data/example.csv').write_text(folder.name)
            subprocess.run(['git','init','-q',str(destination)],check=True)
            with self.assertRaisesRegex(ValueError,'untracked'):
                cloud.plan(source,destination)


if __name__=='__main__':
    unittest.main()
