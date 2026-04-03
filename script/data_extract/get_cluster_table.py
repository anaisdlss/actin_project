import gzip
import json
import os
import re
import sys

from script.utils.utils import create_session


PPI3D_HOME_URL = "https://bioinformatics.lt/ppi3d/start"
BASE_URL = "https://bioinformatics.lt/ppi3d/clusters"
SUBMIT_URL = f"{BASE_URL}/submit_interfaces_request"
OUTPUT_DIR = "data/raw"
META_PATH = os.path.join(OUTPUT_DIR, "metadata_all_data.json")


def get_ppi3d_update(session):
    for url in [PPI3D_HOME_URL, BASE_URL]:
        response = session.get(url, timeout=60)
        response.raise_for_status()

        match = re.search(
            r"Protein Data Bank structures\s*\(([^)]+)\)",
            response.text
        )
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

    result_response = session.get(result_url, timeout=60)
    result_response.raise_for_status()

    try:
        data_link = extract_data_link(result_response.text)
    except ValueError:
        debug_path = os.path.join(
            OUTPUT_DIR, "cluster_table_result_debug.html")
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        save_text(debug_path, result_response.text)
        raise ValueError(
            "Data table link not found in the result page. "
            f"Debug HTML saved to {debug_path}"
        )

    print("Downloading:", data_link)

    download_response = session.get(data_link, timeout=120)
    download_response.raise_for_status()

    csv_content = gzip.decompress(download_response.content).decode("utf-8")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    csv_path = os.path.join(OUTPUT_DIR, "all_data.csv")
    save_text(csv_path, csv_content)

    request_id = extract_request_id(result_url)
    save_metadata(current_update, request_id, result_url, data_link, csv_path)

    print("Saved decompressed file:", csv_path)
    print("Saved metadata:", META_PATH)


if __name__ == "__main__":
    main()
