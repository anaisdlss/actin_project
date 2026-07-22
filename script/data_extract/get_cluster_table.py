import gzip
import json
import os
import re
import sys
import time

from script.utils.utils import create_session


PPI3D_HOME_URL = "https://bioinformatics.lt/ppi3d/start"
BASE_URL = "https://bioinformatics.lt/ppi3d/clusters"
SUBMIT_URL = f"{BASE_URL}/submit_interfaces_request"
OUTPUT_DIR = "data/raw"
META_PATH = os.path.join(OUTPUT_DIR, "metadata_all_data.json")


def get_ppi3d_update(session):
    # IMPORTANT : extraire le TEXTE via BeautifulSoup (le regex sur le HTML brut
    # échouait car des balises séparent « structures » de la date → 'unknown',
    # ce qui faisait sauter à tort le re-téléchargement). Même logique que
    # get_summary_results.get_ppi3d_update.
    from bs4 import BeautifulSoup
    for url in [PPI3D_HOME_URL, BASE_URL, "https://bioinformatics.lt/ppi3d"]:
        try:
            response = session.get(url, timeout=60)
            response.raise_for_status()
        except Exception:
            continue
        text = BeautifulSoup(response.text, "html.parser").get_text()
        for line in text.split("\n"):
            if "Protein Data Bank structures" in line:
                match = re.search(r"\((.*?)\)", line)
                if match:
                    return match.group(1).strip()

    return "unknown"


def extract_release_date_to_from_form(html):
    patterns = [
        r'name="PDB_data\[release_date_to\]"\s+value="([^"]+)"',
        r'id="release_date_to_input"\s+value="([^"]+)"',
    ]

    for pattern in patterns:
        match = re.search(pattern, html)
        if match:
            return match.group(1).strip()

    raise ValueError("release_date_to not found in clusters form page")


def extract_result_url_from_text(text):
    match = re.search(r"/ppi3d/clusters/data_request/([A-Za-z0-9]+)", text)
    if not match:
        return None
    return f"https://bioinformatics.lt/ppi3d/clusters/data_request/{match.group(1)}"


def extract_data_link(html):
    match = re.search(
        r'href="([^"]*?/ppi3d/downloads/data_requests/[^"]+?\.csv\.gz)"',
        html
    )

    if not match:
        raise ValueError("Data table link not found in the response page")

    link = match.group(1)

    if link.startswith("http"):
        return link

    if link.startswith("/"):
        return "https://bioinformatics.lt" + link

    return "https://bioinformatics.lt/" + link.lstrip("/")


def extract_request_id(url):
    return url.rstrip("/").split("/")[-1]


def save_text(path, content):
    with open(path, "w") as f:
        f.write(content)


def load_metadata():
    if not os.path.exists(META_PATH):
        return None

    with open(META_PATH) as f:
        return json.load(f)


def save_metadata(update_date, request_id, result_url, data_link, csv_path):
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    metadata = {
        "ppi3d_last_update": update_date,
        "request_id": request_id,
        "result_url": result_url,
        "data_link": data_link,
        "csv_path": csv_path,
    }

    with open(META_PATH, "w") as f:
        json.dump(metadata, f, indent=2)


def slim_to_actin(csv_path, summary_path):
    """Réduit all_data.csv aux seules interactions où l'actine est en S1 OU S2.

    all_data brut = tout PPI3D (~1 M lignes / 132k PDB), dont l'analyse actine
    n'utilise que ~0,3 %. On garde uniquement les lignes dont subunit_1_title OU
    subunit_2_title est une actine réelle (« Result protein » du summary avec
    Expect value == 0.0, même définition que graphe_filter). C'est un sur-ensemble
    EXACT de ce que le filtrage downstream conserve → filtered_all_data reste
    identique (validé). Garde-fou : au moindre doute (summary/colonnes absents,
    aucune actine) on ne touche à rien — mieux vaut garder le full que risquer
    une perte.
    """
    import pandas as pd
    if not os.path.exists(summary_path):
        print(f"  slim ignoré : summary absent ({summary_path})")
        return
    try:
        sm = pd.read_csv(summary_path, sep=";")
        if "Expect value" not in sm.columns or "Result protein" not in sm.columns:
            print("  slim ignoré : colonnes summary inattendues")
            return
        sm["Expect value"] = pd.to_numeric(sm["Expect value"], errors="coerce")
        real_actin = set(sm.loc[sm["Expect value"] == 0.0, "Result protein"]
                         .dropna().unique())
        if not real_actin:
            print("  slim ignoré : aucune actine (Expect==0) dans le summary")
            return
        df = pd.read_csv(csv_path, low_memory=False)
        if "subunit_1_title" not in df.columns or "subunit_2_title" not in df.columns:
            print("  slim ignoré : colonnes subunit_*_title absentes")
            return
        keep = (df["subunit_1_title"].isin(real_actin)
                | df["subunit_2_title"].isin(real_actin))
        n0 = len(df)
        df[keep].to_csv(csv_path, index=False)
        print(f"  slim actine S1/S2 : {n0} -> {int(keep.sum())} lignes "
              f"(sur-ensemble exact du set filtré)")
    except Exception as e:
        print(f"  slim ignoré (erreur, all_data conservé intact) : {e}")


def main():
    session = create_session()
    metadata = load_metadata()
    manual_result_url = None

    try:
        current_update = get_ppi3d_update(session)
    except Exception as e:
        print(f"Warning: could not reach PPI3D server ({e})")
        csv_path = metadata.get("csv_path") if metadata else None
        if csv_path and os.path.exists(csv_path):
            print("Server unavailable but existing cluster table found — using it.")
            return
        else:
            raise RuntimeError(
                "Server unavailable and no existing cluster table found.") from e

    if len(sys.argv) > 1:
        manual_result_url = sys.argv[1].strip()

    if metadata and not manual_result_url:
        csv_path = metadata.get("csv_path")

        if (
            metadata.get("ppi3d_last_update") == current_update
            and csv_path
            and os.path.exists(csv_path)
        ):
            print("PPI3D update unchanged")
            print("Using existing cluster table")
            print("Decompressed file:", csv_path)
            return

    if manual_result_url:
        result_url = manual_result_url
        print("Using provided cluster table request")
    else:
        form_response = session.get(BASE_URL, timeout=60)
        form_response.raise_for_status()
        release_date_to = extract_release_date_to_from_form(form_response.text)

        payload = [
            ("interaction_types[protein_protein_interactions]", "0"),
            ("interaction_types[protein_protein_interactions]", "1"),
            ("interaction_types[protein_peptide_interactions]", "0"),
            ("interaction_types[protein_nucleic_interactions]", "0"),
            ("interaction_types[domain_nucleic_interactions]", "0"),
            ("interaction_types[intra_chain_domain_interactions]", "0"),
            ("PDB_data[resolution]", "4"),
            ("PDB_data[release_date_from]", "1973-01-01"),
            ("PDB_data[release_date_to]", release_date_to),
            ("complex[min_number_of_subunits]", "1"),
            ("complex[max_number_of_subunits]", "10000"),
            ("complex[min_number_of_protein_subunits]", "1"),
            ("complex[max_number_of_protein_subunits]", "10000"),
            ("complex[min_number_of_different_subunits]", "1"),
            ("complex[max_number_of_different_subunits]", "10000"),
            ("complex[min_number_of_residues]", "1"),
            ("complex[max_number_of_residues]", "1000000000"),
            ("complex[min_number_of_residues_visible]", "1"),
            ("complex[max_number_of_residues_visible]", "1000000000"),
            ("subunits[min_number_of_residues]", "1"),
            ("subunits[max_number_of_residues]", "100000"),
            ("subunits[min_number_of_residues_visible]", "1"),
            ("subunits[max_number_of_residues_visible]", "100000"),
            ("interface[min_area]", "100"),
            ("interface[max_area]", "100000"),
            ("interface[min_number_of_contacts]", "1"),
            ("interface[max_number_of_contacts]", "10000"),
            ("interface[allow_ligands]", "0"),
            ("interface[allow_ligands]", "1"),
            ("interface[homo]", "0"),
            ("interface[homo]", "1"),
            ("interface[hetero]", "0"),
            ("interface[hetero]", "1"),
            ("clustering", "none"),
            ("include_clustering_data_in_table", "0"),
            ("include_clustering_data_in_table", "1"),
            ("submit", "Submit"),
        ]

        print("Submitting PPI3D cluster table request...")
        response = session.post(SUBMIT_URL, data=payload, timeout=60)
        response.raise_for_status()

        result_url = ""

        if "/ppi3d/clusters/data_request/" in response.url:
            result_url = response.url

        if not result_url:
            for previous_response in response.history:
                location = previous_response.headers.get("Location", "")
                if "/ppi3d/clusters/data_request/" in location:
                    if location.startswith("http"):
                        result_url = location
                    else:
                        result_url = "https://bioinformatics.lt" + location
                    break

        if not result_url:
            result_url = extract_result_url_from_text(response.text)

        if not result_url:
            debug_path = os.path.join(
                OUTPUT_DIR, "cluster_table_submit_debug.html")
            os.makedirs(OUTPUT_DIR, exist_ok=True)
            save_text(debug_path, response.text)
            raise ValueError(
                "Result page URL not found after form submission. "
                f"Final URL: {response.url}. Debug HTML saved to {debug_path}"
            )

    print("Request page:", result_url)

    data_link = None
    max_wait = 600
    poll_interval = 15
    elapsed = 0

    while elapsed < max_wait:
        result_response = session.get(result_url, timeout=60)
        result_response.raise_for_status()
        try:
            data_link = extract_data_link(result_response.text)
            break
        except ValueError:
            if "queued" in result_response.text or "Preparing" in result_response.text or "loading.gif" in result_response.text:
                print(f"Data still being prepared, waiting {poll_interval}s... ({elapsed}s elapsed)")
                time.sleep(poll_interval)
                elapsed += poll_interval
            else:
                debug_path = os.path.join(OUTPUT_DIR, "cluster_table_result_debug.html")
                os.makedirs(OUTPUT_DIR, exist_ok=True)
                save_text(debug_path, result_response.text)
                raise ValueError(
                    "Data table link not found in the result page. "
                    f"Debug HTML saved to {debug_path}"
                )

    if data_link is None:
        debug_path = os.path.join(OUTPUT_DIR, "cluster_table_result_debug.html")
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        save_text(debug_path, result_response.text)
        raise TimeoutError(
            f"Data table not ready after {max_wait}s. Debug HTML saved to {debug_path}"
        )

    print("Downloading:", data_link)

    for attempt in range(5):
        download_response = session.get(data_link, timeout=180)
        download_response.raise_for_status()
        try:
            csv_content = gzip.decompress(download_response.content).decode("utf-8")
            break
        except (EOFError, gzip.BadGzipFile, OSError) as e:
            print(f"Download incomplete ({e}), retry {attempt + 1}/5 in 10s...")
            time.sleep(10)
    else:
        raise RuntimeError("Failed to download a complete gzip file after 5 attempts")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    csv_path = os.path.join(OUTPUT_DIR, "all_data.csv")
    save_text(csv_path, csv_content)
    # Ne garder du raw que la partie utile (actine en S1 ou S2) : ~0,3 % des lignes,
    # le reste (tout PPI3D non-actine) n'est jamais utilisé par le filtrage.
    slim_to_actin(csv_path, os.path.join(OUTPUT_DIR, "ppi3d_actin_summary.csv"))

    request_id = extract_request_id(result_url)
    save_metadata(current_update, request_id, result_url, data_link, csv_path)

    print("Saved decompressed file:", csv_path)
    print("Saved metadata:", META_PATH)


if __name__ == "__main__":
    main()
