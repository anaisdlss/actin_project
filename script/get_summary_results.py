import time
import requests
import pandas as pd
from io import StringIO
from bs4 import BeautifulSoup
import os
import json
import re

session = requests.Session()

BASE_URL = "https://bioinformatics.lt/ppi3d"
SUBMIT_URL = f"{BASE_URL}/site/submit_data"

NAME_PROT = "actin"
UNIPROT_ID = "P60709"

DATA_DIR = "data"
META_PATH = os.path.join(DATA_DIR, "metadata.json")

os.makedirs(DATA_DIR, exist_ok=True)


# -----------------------------
# récupérer date update PPI3D
# -----------------------------

def get_ppi3d_update():

    html = requests.get(BASE_URL).text
    soup = BeautifulSoup(html, "html.parser")

    text = soup.get_text()

    for line in text.split("\n"):

        if "Protein Data Bank structures" in line:

            match = re.search(r"\((.*?)\)", line)

            if match:
                return match.group(1)

    return "unknown"


# -----------------------------
# lire metadata
# -----------------------------

def load_metadata():

    if os.path.exists(META_PATH):

        with open(META_PATH) as f:
            return json.load(f)

    return None


# -----------------------------
# sauvegarder metadata
# -----------------------------

def save_metadata(update_date, job_id, results_url):

    metadata = {
        "uniprot": UNIPROT_ID,
        "ppi3d_last_update": update_date,
        "job_id": job_id,
        "results_url": results_url
    }

    with open(META_PATH, "w") as f:
        json.dump(metadata, f, indent=2)


# -----------------------------
# vérifier update
# -----------------------------

current_update = get_ppi3d_update()
metadata = load_metadata()

if metadata:

    if metadata["ppi3d_last_update"] == current_update:

        print("No database update detected")
        job_id = metadata["job_id"]

    else:

        print("Database updated → new job required")
        metadata = None

# -----------------------------
# soumettre nouveau job
# -----------------------------

if metadata is None:

    print("Submitting job...")

    data = {
        "sequence": "",
        "uniprot_ac": UNIPROT_ID,
        "job_name": "",
        "email_value": "",
        "two_sequences_query": "",
        "nucleic_acid_query": "",
        "sequence_search_method": "blast",
        "psiblast_num_iterations": "2",
        "psiblast_inclusion_ethresh": "0.002",
        "e_value": "0.01",
        "word_size": "3",
        "matrix": "BLOSUM62",
        "gap_costs": "11 1",
        "composition_based_statistics": "2",
        "mask_lowercase_letters": "0",
        "filter_low_complexity": "0",
        "use_soft_masking": "0",
        "pseudocount": "0",
        "submit": "Submit"
    }

    response = session.post(SUBMIT_URL, data=data, allow_redirects=True)

    final_url = response.url

    job_id = final_url.split("/")[-1]

    print("New Job ID:", job_id)

    results_url = f"{BASE_URL}/site/results/{job_id}"
    save_metadata(current_update, job_id, results_url)

else:

    print("Using frozen job:", job_id)
    results_url = f"{BASE_URL}/site/results/{job_id}"
    save_metadata(current_update, job_id, results_url)


# -----------------------------
# attendre BLAST
# -----------------------------

results_url = f"{BASE_URL}/site/results/{job_id}"

print("Waiting for BLAST results...")

while True:

    html = session.get(results_url).text

    if "Processing your query" not in html:
        break

    print("Job running...")
    time.sleep(10)

print("Job finished")


# -----------------------------
# detailed page
# -----------------------------

detailed_url = (
    f"{BASE_URL}/site/detailed/single_sequence/"
    f"{job_id}/1/null/1/0/1/0"
)

html = session.get(detailed_url).text

print("Detailed page downloaded")

soup = BeautifulSoup(html, "html.parser")


# -----------------------------
# table interactions
# -----------------------------

tables = pd.read_html(StringIO(html), flavor="bs4")
summary_df = tables[0]

print("Number of interactions:", len(summary_df))


# -----------------------------
# extraire liens DETAILS et CLUSTERS
# -----------------------------

rows = soup.select("#detailed_table tbody tr")

detail_links = []
cluster_links = []

for i, row in enumerate(rows):

    # ------------------
    # detail link
    # ------------------

    detail = row.find("a", href=lambda x: x and "interaction_details" in x)

    if detail:
        link = detail["href"]
        if not link.startswith("http"):
            link = BASE_URL + link
        detail_links.append(link)
    else:
        detail_links.append(None)

    # ------------------
    # cluster link
    # ------------------

    cluster = row.find("a", href=lambda x: x and "detailed_cluster" in x)

    if cluster:
        link = cluster["href"]
        if not link.startswith("http"):
            link = BASE_URL + link
        cluster_links.append(link)
    else:
        cluster_links.append(None)


# -----------------------------
# ajouter les colonnes
# -----------------------------

if len(detail_links) == len(summary_df):

    summary_df["detail_url"] = detail_links
    summary_df["cluster_url"] = cluster_links

else:

    print("Warning: mismatch between table rows and detail links")


# -----------------------------
# sauvegarder summary
# -----------------------------

summary_path = os.path.join(
    DATA_DIR,
    f"ppi3d_{NAME_PROT}_summary.csv"
)

summary_df.to_csv(summary_path, index=False, sep=";")

print("Summary saved:", summary_path)
print("\nSummary extraction finished.")
