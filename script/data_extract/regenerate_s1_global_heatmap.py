#!/usr/bin/env python
"""Compatibility entry point: regenerate using the app's corrected calculations."""
from pathlib import Path
import subprocess, sys
if __name__ == '__main__':
    root = Path(__file__).resolve().parents[2]
    subprocess.run([sys.executable, str(root/'tools/regenerate_scientific_figures.py'), '--scope', 'global'], cwd=root, check=True)
