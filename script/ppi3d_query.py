import time
import requests
import pandas as pd
from io import StringIO
from bs4 import BeautifulSoup
import os

session = requests.Session()

submit_url = "https://bioinformatics.lt/ppi3d/site/submit_data"

DATA_DIR = "data"
os.makedirs(DATA_DIR, exist_ok=True)
NAME_PROT = "actin"  # human
UNIPROT_ID = "P60709"

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

# envoyer la requête
response = session.post(submit_url, data=data, allow_redirects=True)

html = response.text

print("Page récupérée :", len(html), "caractères")

print(response.url)
# print(html)

job_id = response.url.split("/")[-1]

results_url = f"https://bioinformatics.lt/ppi3d/site/results/{job_id}"

while True:
    try:
        response = session.get(results_url, timeout=30)
        html = response.text
    except requests.exceptions.RequestException as e:
        print("Erreur réseau :", e)
        print("On attend 10s puis on réessaie...")
        time.sleep(10)
        continue
    if "Processing your query" not in html:
        print("Job terminé")
        break
    print("Job running...")
    time.sleep(10)


detailed_url = f"https://bioinformatics.lt/ppi3d/site/detailed/single_sequence/{
    job_id}/1/null/1/0/1/0"

while True:

    try:
        html = session.get(detailed_url, timeout=30).text

    except requests.exceptions.RequestException:
        print("Erreur réseau, on réessaie...")
        time.sleep(10)
        continue

    # tant que la table n'est pas là
    if "Structure title" not in html:
        print("Detailed page not ready...")
        time.sleep(10)
        continue

    break

print("Detailed page ready")
print("Length page detailed :", len(html))


tables = pd.read_html(StringIO(html))
print("Number of tables :", len(tables))

df = tables[0]

print("Number of protéine-protéine results:", len(df))
print(df.head())

# -----------------------------------
# récupérer les liens des pages DETAILS
# -----------------------------------

soup = BeautifulSoup(html, "html.parser")

detail_links = []

for a in soup.find_all("a", href=True):

    if "interaction_details" in a["href"]:

        link = a["href"]

        if not link.startswith("http"):
            link = "https://bioinformatics.lt" + link

        detail_links.append(link)

print("Number of detail pages :", len(detail_links))
print("Example links :", detail_links[:5])

# -----------------------------------
# télécharger les pages details
# -----------------------------------

details_tables = []

for i, url in enumerate(detail_links):

    print(f"Downloading {i+1}/{len(detail_links)}")

    try:

        r = session.get(url, timeout=30)

        detail_html = r.text

        tables = pd.read_html(StringIO(detail_html))

        if len(tables) > 0:

            table = tables[0]

            table["source_url"] = url

            details_tables.append(table)

    except Exception as e:

        print("Error:", e)

    time.sleep(1)

# -----------------------------------
# construire dataset final
# -----------------------------------

print("\nBuilding final dataset...")

if len(details_tables) > 0:

    details_df = pd.concat(details_tables, ignore_index=True)

    csv_details = os.path.join(
        DATA_DIR, f"ppi3d_{NAME_PROT}_details_dataset.csv")

    details_df.to_csv(csv_details, index=False)

    print(f"Details dataset saved : {csv_details}")
    print(f"Number of detail rows : {len(details_df)}")

else:
    print("No detail tables found.")

# -----------------------------------
# sauvegarder les résultats principaux
# -----------------------------------

print("\nSaving summary results...")

csv_results = os.path.join(DATA_DIR, f"ppi3d_{NAME_PROT}_results.csv")

df.to_csv(csv_results, index=False)

print(f"Results saved : {csv_results}")
print(f"Number of summary rows : {len(df)}")

print("\nPipeline finished successfully.")
