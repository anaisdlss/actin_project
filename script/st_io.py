"""Shared Streamlit IO helpers (CSV/PDB loading), extracted to keep streamlit.py light."""
import csv
import os
import pandas as pd
import streamlit as st


@st.cache_data
def _load_pdb_file(path, mtime=None):
    with open(path) as f:
        return f.read()


def extract_chain(pdb_data: str, chain_letter: str) -> str:
    lines = [l for l in pdb_data.splitlines()
             if l.startswith('ATOM')
             and len(l) > 21 and l[21] == chain_letter]
    return "\n".join(lines)


@st.cache_data
def load_csv(path, mtime=None):
    with open(path, newline="") as f:
        sample = f.read(10000)
    try:
        sep = csv.Sniffer().sniff(sample, delimiters=";,\t").delimiter
    except csv.Error:
        sep = ","
    return pd.read_csv(path, sep=sep)


def read_csv(path):
    """Charge un CSV en invalidant le cache si le fichier a changé."""
    mtime = os.path.getmtime(path) if os.path.exists(path) else None
    return load_csv(path, mtime=mtime)
