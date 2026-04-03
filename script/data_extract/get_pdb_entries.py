import os
import json
import pandas as pd
from io import StringIO
from bs4 import BeautifulSoup

from script.utils.utils import (
    create_session,
    dataset_is_valid,
    save_dataset_metadata
)

# -----------------------------
# CONFIG
# -----------------------------

DATA_DIR = "data/raw"
SUMMARY_PATH = os.path.join(DATA_DIR, "ppi3d_actin_summary.csv")
MAIN_META_PATH = os.path.join(DATA_DIR, "metadata.json")

BASE_PAGE = "https://bioinformatics.lt/ppi3d/site/interactions_in_pdb_entry/"
BASE_AJAX = "https://bioinformatics.lt/ppi3d/site/load_interactions_in_pdb_entry_table/"

REQUEST_TIMEOUT = 30

session = create_session()
session.headers.update({
    "X-Requested-With": "XMLHttpRequest"
})


# -----------------------------
# HELPERS
# -----------------------------

def get_current_job_id():
    if os.path.exists(MAIN_META_PATH):
        with open(MAIN_META_PATH) as f:
            return json.load(f).get("job_id")
    return None


def fetch_pdb_table(pdb: str):
    pdb = str(pdb).strip()

    page_pdb = pdb.lower()
    ajax_pdb = pdb.upper()

    page_url = BASE_PAGE + page_pdb
    ajax_url = BASE_AJAX + ajax_pdb

    # 1) ouvrir la page normale
    page_response = session.get(page_url, timeout=REQUEST_TIMEOUT)
    page_response.raise_for_status()

    # 2) appeler la requête AJAX comme le navigateur
    ajax_response = session.get(
        ajax_url,
        headers={
            "Referer": page_url,
            "X-Requested-With": "XMLHttpRequest"
        },
        timeout=REQUEST_TIMEOUT
    )
    ajax_response.raise_for_status()

    html = ajax_response.text
    soup = BeautifulSoup(html, "lxml")

    table_html = soup.find("table", {"id": "detailed_table"})
    if table_html is None:
        table_html = soup.find("table")

    if table_html is None:
        return None

    table = pd.read_html(StringIO(str(table_html)), flavor="lxml")[0]
    table["pdb_id"] = ajax_pdb

    return table


# -----------------------------
# MAIN
# -----------------------------

def main():

    if not os.path.exists(SUMMARY_PATH):
        raise FileNotFoundError(f"Summary file not found: {SUMMARY_PATH}")

    summary_dir = os.path.dirname(SUMMARY_PATH)
    output_path = os.path.join(summary_dir, "pdb_entry_results.csv")
    meta_path = os.path.join(summary_dir, "metadata_pdb_entries.json")
    output_files = ["pdb_entry_results.csv"]

    print("Using summary:", SUMMARY_PATH)
    print("Output:", output_path)

    current_job_id = get_current_job_id()

    if dataset_is_valid(
        meta_path,
        current_job_id,
        summary_dir,
        output_files
    ):
        print("PDB entries dataset unchanged")
        print("Nothing to do")
        return

    summary = pd.read_csv(SUMMARY_PATH, sep=";")

    if "PDB ID" not in summary.columns:
        raise ValueError("Column 'PDB ID' missing in summary file")

    pdb_ids = (
        summary["PDB ID"]
        .dropna()
        .astype(str)
        .str.strip()
        .str.upper()
        .unique()
    )

    print("Unique PDB:", len(pdb_ids))

    results = []

    for i, pdb in enumerate(pdb_ids):
        print(f"\n{i+1}/{len(pdb_ids)} {pdb}")

        try:
            table = fetch_pdb_table(pdb)

            if table is None:
                print(f"{pdb} -> no interactions table found")
                continue

            print("Interactions:", len(table))
            results.append(table)

        except Exception as e:
            print("Error:", pdb, e)

    if results:
        df = pd.concat(results, ignore_index=True)
        df.to_csv(output_path, sep=";", index=False)

        print("\nSaved:", output_path)

        save_dataset_metadata(
            meta_path,
            current_job_id,
            summary_dir,
            output_files
        )

        print("PDB entries metadata saved")

    else:
        print("No results")


if __name__ == "__main__":
    main()
