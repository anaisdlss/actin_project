"""PDB menu entries derived from the current retained interaction dataset."""
import pandas as pd


def _pdb_ids(values):
    return values.astype("string").str.strip().str.upper().replace("", pd.NA)


def build_pdb_catalog(retained, summary=None, interactions=None):
    """Use retained PDB IDs as the scope; metadata only supplies display titles.

    A missing title or local coordinate file does not invalidate retained data.
    Earlier screening lists and leftover detail files cannot add menu entries.
    """
    if retained is None or "pdb_id" not in retained:
        return pd.DataFrame(columns=["pdb_id", "title"])
    ids = sorted(_pdb_ids(retained["pdb_id"]).dropna().unique())
    titles = {}
    for table, id_col, title_col in (
        (summary, "PDB ID", "Structure title"),
        (interactions, "pdb_id", "structure_title"),
        (retained, "pdb_id", "pdb_annotation"),
    ):
        if table is None or id_col not in table or title_col not in table:
            continue
        for pdb_id, title in zip(_pdb_ids(table[id_col]), table[title_col]):
            if pd.isna(pdb_id) or pd.isna(title):
                continue
            title = str(title).strip()
            if title and title.lower() not in {"nan", "none", "null", "<na>"}:
                titles.setdefault(pdb_id, title)
    return pd.DataFrame({"pdb_id": ids,
                         "title": [titles.get(p, "Title unavailable") for p in ids]})
