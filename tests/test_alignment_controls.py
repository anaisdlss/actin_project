import sys
import unittest
from pathlib import Path
from unittest.mock import Mock, patch
import streamlit
sys.path.append(str(Path(__file__).resolve().parents[1] / 'script'))
import msa_analysis


class AlignmentControlsTest(unittest.TestCase):
    def test_public_build_cannot_run_even_when_mafft_is_installed(self):
        with patch.object(msa_analysis, '_Path') as path, patch('shutil.which', return_value='/bin/mafft'), patch.object(msa_analysis.st, 'caption'):
            path.return_value.exists.return_value = True
            button = Mock()
            msa_analysis._mafft_controls(button, 'public')
            self.assertTrue(button.button.call_args.kwargs['disabled'])

    def test_local_run_requires_executable(self):
        for executable in [None, '/bin/mafft']:
            with self.subTest(executable=executable), patch.object(msa_analysis, '_Path') as path, patch('shutil.which', return_value=executable), patch.object(msa_analysis.st, 'caption'):
                path.return_value.exists.return_value = False
                button = Mock()
                msa_analysis._mafft_controls(button, 'local')
                self.assertEqual(button.button.call_args.kwargs['disabled'], executable is None)
