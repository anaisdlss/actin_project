import os
import requests
import pandas as pd
from io import StringIO
from bs4 import BeautifulSoup

DATA_DIR = "data"

BASE_TABLE = "https://bioinformatics.lt/ppi3d/site/load_interactions_in_pdb_entry_table/"

session = requests.Session()
session.headers.update({
    "User-Agent": "Mozilla/5.0",
    "X-Requested-With": "XMLHttpRequest"
})


def main():

    mode = input(
        "Use cluster summary instead of default summary? (y/[n]): ").lower()

    if mode == "y":

        SUMMARY_PATH = os.path.join(
            DATA_DIR,
            "clusters",
            "clusters_summary.csv"
        )

        OUTPUT = os.path.join(
            DATA_DIR,
            "clusters",
            "pdb_entry_results.csv"
        )

    else:

        SUMMARY_PATH = os.path.join(
            DATA_DIR,
            "ppi3d_actin_summary.csv"
        )

        OUTPUT = os.path.join(
            DATA_DIR,
            "pdb_entry_results.csv"
        )

    summary = pd.read_csv(SUMMARY_PATH, sep=";")

    pdb_ids = summary["PDB ID"].dropna().unique()

    print("Unique PDB:", len(pdb_ids))

    results = []

    for i, pdb in enumerate(pdb_ids):

        pdb = pdb.upper()

        print(f"\n{i+1}/{len(pdb_ids)} {pdb}")

        url = BASE_TABLE + pdb

        try:

            r = session.get(url, timeout=30)

            html = r.text

            soup = BeautifulSoup(html, "lxml")

            table_html = soup.find("table", {"id": "detailed_table"})

            if table_html is None:
                print(f"{pdb} → no interactions")
                continue

            table = pd.read_html(StringIO(str(table_html)))[0]

            table["pdb_id"] = pdb

            print("Interactions:", len(table))

            results.append(table)

        except Exception as e:

            print("Error:", pdb, e)

    if results:

        df = pd.concat(results, ignore_index=True)

        df.to_csv(OUTPUT, index=False)

        print("\nSaved:", OUTPUT)

    else:

        print("No results")


if __name__ == "__main__":
    main()
