"""Actin residue numbering shown to users: UniProt P60709 (human beta-actin).

Internal tables keep the MAFFT column of the actin alignment
(`data/alignments/cluster_6685.aln`): `residue_*_canon_mafft`, `canon`.
Those columns are alignment coordinates, not residue numbers. This module
converts them for display only, so joins between tables stay unchanged.

The conversion is read from the P60709 row of that alignment (7pdz_I):
  column 3            -> M1
  columns 6-236       -> column - 4
  columns 238-380     -> column - 5
  columns 1-2, 4-5    -> no P60709 residue (N-terminal insertion of alpha-actins)
  column 237          -> no P60709 residue (insertion in some actins)
tests/test_numbering.py checks this rule against the alignment when present.
"""
import pandas as pd

ACTIN_LENGTH = 375
AXIS_TITLE = "Actin residue (UniProt P60709)"
SHORT_LABEL = "P60709 position"


def to_uniprot(col):
    """MAFFT column of the actin alignment -> P60709 residue number (or None)."""
    try:
        c = int(col)
    except (TypeError, ValueError):
        return None
    if c == 3:
        return 1
    if 6 <= c <= 236:
        return c - 4
    if 238 <= c <= 380:
        return c - 5
    return None


def to_canon(pos):
    """P60709 residue number -> MAFFT column (inverse of to_uniprot)."""
    try:
        p = int(pos)
    except (TypeError, ValueError):
        return None
    if p == 1:
        return 3
    if 2 <= p <= 232:
        return p + 4
    if 233 <= p <= ACTIN_LENGTH:
        return p + 5
    return None


def label(col):
    """Text for a MAFFT column: its P60709 number, or an explicit insertion tag."""
    pos = to_uniprot(col)
    if pos is not None:
        return str(pos)
    try:
        return f"ins. (col. {int(col)})"
    except (TypeError, ValueError):
        return "?"


def uniprot_series(values):
    """Vectorised to_uniprot for a pandas Series (nullable integers)."""
    return pd.Series(values).map(to_uniprot).astype("Int64")
