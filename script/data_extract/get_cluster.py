import os
import requests
import pandas as pd
from io import StringIO
from bs4 import BeautifulSoup
from concurrent.futures import ThreadPoolExecutor
import threading
from script.utils.utils import (
    create_session,
    dataset_is_valid,
    save_dataset_metadata,
    load_metadata
)
import json

# -----------------------------
# CONFIG
# -----------------------------

DATA_DIR = "data/raw"
SUMMARY_PATH = os.path.join(DATA_DIR, "ppi3d_actin_summary.csv")
CLUSTER_DIR = os.path.join(DATA_DIR, "clusters")

CLUSTER_META_PATH = os.path.join(CLUSTER_DIR, "metadata.json")

OUTPUT_FILES = [
    "clusters_summary.csv"
]

REQUEST_TIMEOUT = 30

os.makedirs(CLUSTER_DIR, exist_ok=True)

session = create_session()

interaction_counter = 0
counter_lock = threading.Lock()
total_interactions = 0


def download_cluster(row):

    cluster_size = row.get("No. of members in cluster", 1)

    if cluster_size <= 1:
        cluster_size = 0

    global interaction_counter

    with counter_lock:
        interaction_counter += cluster_size
        print(
            f"Downloaded {interaction_counter}/{total_interactions} interactions")
    url = row["cluster_url"]

    if pd.isna(url):
        return None

    try:
        r = session.get(url, timeout=REQUEST_TIMEOUT)
        r.raise_for_status()

        html = r.text

        tables = pd.read_html(StringIO(html), flavor="lxml")

        if len(tables) == 0:
            return None

        table = tables[0]

        soup = BeautifulSoup(html, "lxml")

        links = []
        for a in soup.select("#detailed_table a[href*='interaction_details']"):
            links.append(a["href"])

        table["detail_url"] = links[:len(table)]
        table["cluster_id"] = row["cluster_id"]

        return table

    except Exception as e:
        print("Cluster error:", e)
        return None


# -----------------------------
# MAIN
# -----------------------------

def main():

    current_job_id = None

    if os.path.exists(os.path.join(DATA_DIR, "metadata.json")):

        with open(os.path.join(DATA_DIR, "metadata.json")) as f:
            current_job_id = json.load(f).get("job_id")

    if dataset_is_valid(
        CLUSTER_META_PATH,
        current_job_id,
        CLUSTER_DIR,
        OUTPUT_FILES
    ):

        print("Clusters dataset unchanged")
        print("Nothing to do")
        return

    global total_interactions

    summary = pd.read_csv(SUMMARY_PATH, sep=";")

    total_interactions = summary["No. of members in cluster"]
    total_interactions = total_interactions[total_interactions > 1].sum()

    print("Total interactions dans clusters to download:", total_interactions)

    if "cluster_url" not in summary.columns:
        raise ValueError("cluster_url column missing in summary file")

    clusters = summary[["Link to details", "cluster_url",
                        "No. of members in cluster"]].copy()

    clusters["cluster_id"] = clusters["Link to details"].str.extract(
        r"(\d+)").astype(int)

    print("Number of clusters:", len(clusters))

    with ThreadPoolExecutor(max_workers=10) as executor:
        results = list(executor.map(download_cluster,
                                    clusters.to_dict("records")))

    clusters_data = [r for r in results if r is not None]

    # -----------------------------
    # SAVE DATA
    # -----------------------------

    if clusters_data:

        clusters_df = pd.concat(clusters_data, ignore_index=True)

        output = os.path.join(CLUSTER_DIR, "clusters_summary.csv")

        clusters_df.to_csv(output, sep=";", index=False)

        print("Clusters saved:", output)

        save_dataset_metadata(
            CLUSTER_META_PATH,
            current_job_id,
            CLUSTER_DIR,
            OUTPUT_FILES
        )

        print("Clusters metadata saved")

    else:

        print("No clusters extracted")


if __name__ == "__main__":
    main()
