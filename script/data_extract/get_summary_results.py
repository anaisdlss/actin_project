import time
import requests
import pandas as pd
from io import StringIO
from bs4 import BeautifulSoup
import os
import json
import re
import hashlib
import sys

session = requests.Session()

BASE_URL = "https://bioinformatics.lt/ppi3d"
SUBMIT_URL = f"{BASE_URL}/site/submit_data"

NAME_PROT = "actin"
UNIPROT_ID = "P60709"

DATA_DIR = "data/raw"
META_PATH = os.path.join(DATA_DIR, "metadata.json")

summary_path = os.path.join(
    DATA_DIR,
    f"ppi3d_{NAME_PROT}_summary.csv"
)

os.makedirs(DATA_DIR, exist_ok=True)

# ---------
# hash
# ---------


def file_hash(path):

    with open(path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


# -----------------------------
# récupérer date update PPI3D
# -----------------------------

def get_ppi3d_update():
    try:
        html = requests.get(BASE_URL).text
    except Exception:
        return "unknown"
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

def save_metadata(
    update_date,
    job_id,
    summary_path,
    source="submitted_job"
):

    metadata = {
        "uniprot": UNIPROT_ID,
        "ppi3d_last_update": update_date,
        "job_id": job_id,
        "summary_hash": file_hash(summary_path),
        "source": source
    }

    with open(META_PATH, "w") as f:
        json.dump(metadata, f, indent=2)


def extract_job_id_from_url(url):
    if not isinstance(url, str) or not url:
        return None

    match = re.search(r"/single_sequence/([^/]+)/", url)

    if match:
        return match.group(1)

    match = re.search(r"/results/([^/]+)", url)

    if match:
        return match.group(1)

    return None


def bootstrap_metadata_from_existing_summary(update_date):
    if not os.path.exists(summary_path):
        return None

    try:
        df = pd.read_csv(summary_path, sep=";")
    except Exception:
        return None

    for column in ["detail_url", "cluster_url"]:
        if column not in df.columns:
            continue

        for value in df[column].dropna().astype(str):
            job_id = extract_job_id_from_url(value)

            if job_id:
                save_metadata(
                    update_date,
                    job_id,
                    summary_path,
                    source="bootstrapped_existing_summary"
                )
                return job_id

    return None


# -----------------------------
# vérifier update
# -----------------------------

current_update = get_ppi3d_update()
metadata = load_metadata()

print("\n--- CHECK DATASET STATUS ---")

if metadata:
    print("Metadata found")
    print("job id:", metadata["job_id"])
    print("Stored DB update:", metadata["ppi3d_last_update"])
    print("Current DB update:", current_update)

if metadata and os.path.exists(summary_path):

    current_hash = file_hash(summary_path)

    print("Stored hash :", metadata.get("summary_hash"))
    print("Current hash:", current_hash)

    if metadata["ppi3d_last_update"] != current_update:

        print("\nDatabase update detected → regenerate dataset")
        metadata = None

    else:

        if metadata.get("summary_hash") != current_hash:
            print("\nDataset file changed but DB update is unchanged")
            print("Refreshing metadata hash and keeping the same job id")
            save_metadata(
                current_update,
                metadata["job_id"],
                summary_path,
                source=metadata.get("source", "submitted_job")
            )
        else:
            print("\nDataset valid")

        print("No changes detected at PPI3D level → using existing dataset")
        sys.exit()

else:

    if os.path.exists(summary_path):
        print("Summary found without metadata → rebuilding metadata")
        job_id = bootstrap_metadata_from_existing_summary(current_update)

        if job_id:
            print("Recovered Job ID from existing summary:", job_id)
            print("Metadata recreated:", META_PATH)
            print("Using existing dataset")
            sys.exit()

    print("No metadata or dataset missing → generating dataset")

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

else:

    print("Using frozen job:", job_id)


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
# summary page
# -----------------------------

clustered_summary_url = (
    f"{BASE_URL}/site/detailed/single_sequence/"
    f"{job_id}/1/null/1/0/0/0"
)

html = session.get(clustered_summary_url).text

soup = BeautifulSoup(html, "html.parser")

settings_form = soup.find(
    "form",
    action=lambda x: x and "change_visualisation_settings" in x
)

if settings_form:

    settings_url = settings_form.get("action")
    target_path = None

    if not settings_url.startswith("http"):
        settings_url = BASE_URL + settings_url

    form_data = {}

    for field in settings_form.find_all(["input", "select"]):
        name = field.get("name")

        if not name:
            continue

        if field.name == "select":
            selected = field.find("option", selected=True)

            if selected is None:
                selected = field.find("option")

            form_data[name] = selected.get("value", "") if selected else ""
            continue

        if field.get("type") in {"checkbox", "radio"}:
            if field.has_attr("checked"):
                form_data[name] = field.get("value", "on")
            continue

        if field.get("type") == "submit":
            continue

        form_data[name] = field.get("value", "")

    target_path = form_data.get("to")
    form_data["clustering"] = "none"

    response = session.post(settings_url, data=form_data, allow_redirects=True)

    if target_path:
        if not target_path.startswith("http"):
            target_url = f"{BASE_URL}/{target_path.lstrip('/')}"
        else:
            target_url = target_path
    else:
        target_url = response.url

    html = session.get(target_url).text

    print("Summary page switched to no clustering")
    print("No clustering URL:", target_url)

else:

    print("Warning: clustering settings form not found, using default summary page")

print("Summary page downloaded")
soup = BeautifulSoup(html, "html.parser")


# -----------------------------
# table interactions
# -----------------------------

if "detailed_table" not in html:
    with open("debug_summary.html", "w") as f:
        f.write(html)
    raise ValueError(
        "No detailed interaction table found after switching to no clustering. "
        "HTML saved to debug_summary.html"
    )

tables = pd.read_html(StringIO(html), flavor="lxml")
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

summary_df.to_csv(summary_path, index=False, sep=";")

save_metadata(current_update, job_id, summary_path)

print("Summary saved:", summary_path)
print("\nSummary extraction finished.")
