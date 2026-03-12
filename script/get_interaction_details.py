import os
import time
import gzip
import re
import requests
import pandas as pd
from bs4 import BeautifulSoup
from io import StringIO


# -----------------------------
# UTILS
# -----------------------------

def create_directories():
    """Create all necessary folders."""
    os.makedirs(PAIRWISE_DIR, exist_ok=True)
    os.makedirs(ASSEMBLY_DIR, exist_ok=True)
    os.makedirs(PYMOL_DIR, exist_ok=True)


def download_text(session, url):
    try:
        r = session.get(url, timeout=REQUEST_TIMEOUT)
        r.raise_for_status()
        return r.text
    except Exception as e:
        print(f"Download error: {url} -> {e}")
        return None


def download_binary(session, url):
    try:
        r = session.get(url, timeout=REQUEST_TIMEOUT)
        r.raise_for_status()
        return r.content
    except Exception as e:
        print(f"Download error: {url} -> {e}")
        return None


def save_file(content, folder, filename, binary=False):
    path = os.path.join(folder, filename)

    mode = "wb" if binary else "w"

    with open(path, mode) as f:
        f.write(content)

    return path


def save_dataframe(data, filename):
    df = pd.DataFrame(data)
    path = os.path.join(DETAIL_DIR, filename)
    df.to_csv(path, index=False)
    print(f"{filename} saved:", path)

# -----------------------------
# CONFIG
# -----------------------------


DATA_DIR = "data"
DETAIL_DIR = os.path.join(DATA_DIR, "details")

REQUEST_TIMEOUT = 30
PAUSE = 0.05

session = requests.Session()
session.headers.update({
    "User-Agent": "ppi3d-actin-scraper/1.0"
})
adapter = requests.adapters.HTTPAdapter(pool_connections=50, pool_maxsize=50)
session.mount("https://", adapter)


def main():

    # -----------------------------
    # CHOOSE INPUT
    # -----------------------------

    global DETAIL_DIR, PAIRWISE_DIR, ASSEMBLY_DIR, PYMOL_DIR

    mode = input(
        "Use clusters summary instead of default summary? (y/[n]): "
    ).strip().lower()

    if mode == "y":
        SUMMARY_PATH = os.path.join(
            DATA_DIR, "clusters", "clusters_summary.csv")
        DETAIL_DIR = os.path.join(DATA_DIR, "clusters", "details")
        is_cluster = True
    else:
        SUMMARY_PATH = os.path.join(DATA_DIR, "ppi3d_actin_summary.csv")
        DETAIL_DIR = os.path.join(DATA_DIR, "details")
        is_cluster = False

    PAIRWISE_DIR = os.path.join(DETAIL_DIR, "structures_files/pairwise")
    ASSEMBLY_DIR = os.path.join(DETAIL_DIR, "structures_files/assembly")
    PYMOL_DIR = os.path.join(DETAIL_DIR, "structures_files/pymol")

    create_directories()

    # -----------------------------
    # LOAD SUMMARY
    # -----------------------------

    summary = pd.read_csv(SUMMARY_PATH, sep=";")
    print("Using summary file:", SUMMARY_PATH)

    detail_urls = summary["detail_url"].tolist()

    # -----------------------------
    # PARSE PAGE DETAILS
    # -----------------------------

    interactions = []
    proteins = []
    interface_residues = []
    residue_contacts = []
    ligands = []
    metadata_alignements = []
    alignment_sequences = []
    structures = []

    for i, url in enumerate(detail_urls):

        print(f"Downloading {i+1}/{len(detail_urls)}")

        try:

            r = session.get(url, timeout=REQUEST_TIMEOUT)
            r.raise_for_status()
            html = r.text

            soup = BeautifulSoup(html, "lxml")
            tables = pd.read_html(StringIO(html))

            # récupérer le texte de la page
            text = soup.get_text(separator="\n", strip=True)
            lines = text.split("\n")

            if is_cluster:
                interaction_id = summary.iloc[i]["Link to details"].split(" ")[
                    0].strip()
                cluster_id = summary.iloc[i]["cluster_id"]
                interaction = {
                    "interaction_id": f"{cluster_id}.{interaction_id}"
                }
            else:
                interaction = {"interaction_id": i+1}

            current_id = interaction["interaction_id"]
            safe_id = str(current_id).replace(".", "_")

            for idx, line in enumerate(lines):
                line = line.strip()

            # -----------------------------
            # table 1 interaction
            # -----------------------------

                if line.startswith("PDB ID"):
                    interaction["pdb_id"] = lines[idx + 1].strip()

                if line.startswith("Title") and "structure_title" not in interaction:
                    interaction["structure_title"] = line.replace(
                        "Title:", "").strip()

                if line.startswith("Release date"):
                    interaction["release_date"] = line.split(":")[-1].strip()

                if line.startswith("Resolution"):
                    interaction["resolution"] = line.split(":")[-1].strip()

                if line.startswith("Chain A"):
                    interaction["chain_A_id"] = line.split(":")[-1].strip()

                if line.startswith("Chain B"):
                    interaction["chain_B_id"] = line.split(":")[-1].strip()

                if "Buried interface area" in line:
                    area = line.split(":")[-1].strip().split()[0]
                    interaction["interface_area"] = f"{area} Å²"

                if "Number of inter-residue contacts" in line:
                    interaction["num_contacts"] = line.split(":")[-1].strip()

                if "Number of H-bonds" in line:
                    interaction["num_hbonds"] = line.split(":")[-1].strip()

                if "Number of salt bridges" in line:
                    interaction["num_salt_bridges"] = line.split(
                        ":")[-1].strip()

            # -----------------------------
            # table 2 proteins
            # -----------------------------

                if line.startswith("Chain A") or line.startswith("Chain B"):

                    protein = {}

                    protein["interaction_id"] = current_id
                    protein["chain_id"] = line.split(":")[-1].strip()

                    protein["protein_name"] = lines[idx +
                                                    1].replace("Title:", "").strip()
                    protein["organism"] = lines[idx +
                                                2].replace("Source organism:", "").strip()

                    residues = lines[idx + 3].split(":")[-1].strip()

                    match = re.search(r"\d+", residues)

                    if match:
                        protein["num_residues"] = int(match.group())
                    else:
                        protein["num_residues"] = None

                    proteins.append(protein)

            interaction["source_url"] = url
            interactions.append(interaction)

            # -----------------------------
            # table 3 interface residues
            # -----------------------------

            interface_tables = [
                t for t in tables
                if "Residue no. in structure" in " ".join(t.columns.astype(str))
            ]

            chains = []

            for tag in soup.find_all("h4"):
                if "Interface residues in" in tag.text:
                    chains.append(tag.text.split("in")[-1].strip())

            for chain, table in zip(chains, interface_tables):

                for _, row in table.iterrows():

                    residue = {}

                    residue["interaction_id"] = current_id
                    residue["chain"] = chain

                    residue["residue_number_structure"] = row["Residue no. in structure"]
                    residue["residue_number_sequence"] = row["Residue no. in sequence"]
                    residue["residue_name"] = row["Residue name"]

                    residue["buried_ASA_Å²"] = row["Buried ASA, Å2"]
                    residue["buried_ASA_percent"] = row["Buried ASA, %"]

                    interface_residues.append(residue)

            # -----------------------------
            # table 4 residue contacts
            # -----------------------------

            for table in tables:

                cols = " ".join(table.columns.astype(str))

                if "Residue no. in chain A structure" in cols and "Residue in chain B" in cols:

                    for _, row in table.iterrows():

                        contact = {}

                        contact["interaction_id"] = current_id

                        contact["chain_A_id"] = interaction["chain_A_id"]
                        contact["residue_A_structure"] = row["Residue no. in chain A structure"]
                        contact["residue_A_sequence"] = row["Residue no. in chain A sequence"]
                        contact["residue_A_name"] = row["Residue in chain A"]

                        contact["chain_B_id"] = interaction["chain_B_id"]
                        contact["residue_B_structure"] = row["Residue no. in chain B structure"]
                        contact["residue_B_sequence"] = row["Residue no. in chain B sequence"]
                        contact["residue_B_name"] = row["Residue in chain B"]

                        contact["contact_area"] = row["Contact area, Å2"]

                        # parfois vide
                        contact["contact_type"] = row.get("Contact type", None)

                        residue_contacts.append(contact)

            # -----------------------------
            # table 5 ligands at interface
            # -----------------------------

            for table in tables:

                cols = " ".join(table.columns.astype(str))

                if "Ligand name" in cols:

                    for _, row in table.iterrows():

                        ligand = {}

                        ligand["interaction_id"] = current_id
                        ligand["chain"] = row["Chain"]
                        ligand["protein"] = row["Protein"]
                        ligand["residue_number"] = row["Residue no."]
                        ligand["ligand_name"] = row["Ligand name"]

                        ligand["contact_area_same_domain"] = row["Contact area with same domain"]
                        ligand["contact_area_other_domain"] = row["Contact area with other domain"]

                        ligands.append(ligand)

            # -----------------------------
            # table 6 sequence alignment metadata
            # -----------------------------

            metadata = {"interaction_id": current_id}

            in_interface_block = False

            for line in lines:

                line = line.strip()

                if line.startswith("Query protein"):
                    protein_line = line.split("Query protein:")[1].strip()
                    metadata["query_protein"] = protein_line.split("|")[1]

                if line.startswith("Result domain"):
                    metadata["template_chain"] = line.split(
                        ":")[1].strip().split(";")[0]

                if "Expectation value" in line:
                    metadata["e_value"] = line.split(
                        "=")[1].split(",")[0].strip()
                    metadata["score"] = line.split(
                        "Score =")[-1].split("bits")[0].strip()

                if line.startswith("Interface alignment data"):
                    in_interface_block = True
                    continue

                if "Interface residues in alignment" in line:
                    metadata["interface_residue_percent"] = line.split(":")[1].split("%")[
                        0].strip()
                    metadata["interface_residue_count"] = line.split(
                        "(")[1].split("/")[0]

                if "Identities =" in line:

                    parts = line.split(",")

                    if in_interface_block:

                        metadata["interface_identity_percent"] = parts[0].split("=")[1].split("%")[
                            0].strip()
                        metadata["interface_positive_percent"] = parts[1].split("=")[1].split("%")[
                            0].strip()
                        metadata["interface_gap_percent"] = parts[2].split("=")[1].split("%")[
                            0].strip()

                    else:

                        metadata["identity_percent"] = parts[0].split("=")[1].split("%")[
                            0].strip()
                        metadata["positive_percent"] = parts[1].split("=")[1].split("%")[
                            0].strip()
                        metadata["gap_percent"] = parts[2].split(
                            "=")[1].split("%")[0].strip()

            metadata_alignements.append(metadata)

            # -----------------------------
            # table 7 alignment
            # -----------------------------

            alignment_seq = {"interaction_id": current_id}

            # ---- récupérer FASTA ----

            for a in soup.find_all("a", href=True):

                if "sequence_alignments" in a["href"] and a["href"].endswith(".fasta"):

                    fasta_text = download_text(session, a["href"])
                    if fasta_text is None:
                        continue
                    fasta_lines = [l.strip()
                                   for l in fasta_text.split("\n") if l.strip()]

                    query_header = fasta_lines[0]
                    template_header = fasta_lines[2]

                    # sequences
                    alignment_seq["query_sequence"] = fasta_lines[1]
                    alignment_seq["template_sequence"] = fasta_lines[3]
                    alignment_seq["interface_positions"] = fasta_lines[5]

                    # query info
                    alignment_seq["query_id"] = query_header.split("|")[1]

                    m = re.search(r"starts at (\d+)", query_header)
                    if m:
                        alignment_seq["query_start"] = int(m.group(1))

                    m = re.search(r"ends at (\d+)", query_header)
                    if m:
                        alignment_seq["query_end"] = int(m.group(1))

                    # template info
                    alignment_seq["template_id"] = template_header.replace(">", "").split()[
                        0]

                    m = re.search(r"starts at (\d+)", template_header)
                    if m:
                        alignment_seq["template_start"] = int(m.group(1))

                    m = re.search(r"ends at (\d+)", template_header)
                    if m:
                        alignment_seq["template_end"] = int(m.group(1))

                    break

            # ---- récupérer DSSP ----

            aligned_block = soup.find(id="aligned_sequence")

            if aligned_block:

                dssp = ""

                for p in aligned_block.find_all("p"):

                    text = p.get_text()

                    if text.startswith("dssp:"):
                        dssp += text.replace("dssp:", "").strip()

                alignment_seq["secondary_structure"] = dssp

            # ---- ajouter résultat ----

            alignment_sequences.append(alignment_seq)

            # -----------------------------
            # table 8 structures
            # -----------------------------

            structure = {"interaction_id": current_id}

            for a in soup.find_all("a", href=True):

                href = a["href"]

                try:

                    # -----------------------------
                    # pairwise PDB
                    # -----------------------------
                    if "pdb_files" in href and href.endswith(".pdb"):

                        pdb_text = download_text(session, href)

                        if pdb_text:
                            pdb_id = interaction["pdb_id"]
                            chainA = interaction["chain_A_id"].split("_")[-1]
                            chainB = interaction["chain_B_id"].split("_")[-1]

                            filename = (
                                f"{safe_id}_{pdb_id}_{chainA}_{chainB}.pdb"
                            )

                            file_path = save_file(
                                pdb_text,
                                PAIRWISE_DIR,
                                filename
                            )

                            structure["pairwise_pdb_file"] = file_path

                    # -----------------------------
                    # PyMOL script
                    # -----------------------------
                    if href.endswith(".py"):

                        pymol_script = download_text(session, href)

                        if pymol_script:

                            pdb_id = interaction["pdb_id"]
                            chainA = interaction["chain_A_id"].split("_")[-1]
                            chainB = interaction["chain_B_id"].split("_")[-1]

                            filename = (
                                f"{safe_id}_{pdb_id}_{chainA}_{chainB}.py")

                            script_path = save_file(
                                pymol_script,
                                PYMOL_DIR,
                                filename
                            )

                            structure["pymol_script_file"] = script_path

                    # -----------------------------
                    # biological assembly
                    # -----------------------------

                    if "biounits" in href:

                        data = download_binary(session, href)

                        if data:

                            content = gzip.decompress(data).decode()

                            # mmCIF assembly
                            if "mmCIF" in href or content.startswith("data_"):

                                pdb_id = interaction["pdb_id"]
                                filename = f"{safe_id}_{pdb_id}.cif"

                                cif_path = save_file(
                                    content,
                                    ASSEMBLY_DIR,
                                    filename
                                )

                                structure["biological_assembly_cif_file"] = cif_path

                            # PDB assembly
                            if "pdb" in href:

                                pdb_id = interaction["pdb_id"]

                                filename = f"{safe_id}_{pdb_id}.pdb"

                                pdb_path = save_file(
                                    content,
                                    ASSEMBLY_DIR,
                                    filename
                                )

                                structure["biological_assembly_pdb_file"] = pdb_path
                except Exception as e:

                    print("Structure download error:", e)

            # -----------------------------
            # biological assembly info
            # -----------------------------

            assembly_info = None

            h4 = soup.find(
                "h4", string=lambda x: x and "Biological assembly" in x)

            if h4:
                div = h4.find_next("div")
                if div:
                    assembly_info = div.get_text(strip=True)

            structure["biological_assembly_info"] = assembly_info

            structures.append(structure)

        except Exception as e:
            print("Error:", e)

        # pause pour éviter surcharge serveur
        time.sleep(PAUSE)

    # -----------------------------
    # SAVE DF
    # -----------------------------

    datasets = {
        "1.interactions.csv": interactions,
        "2.proteins.csv": proteins,
        "3.interface_residues.csv": interface_residues,
        "4.inter-residue_contacts.csv": residue_contacts,
        "5.ligands.csv": ligands,
        "6.meta_alignement.csv": metadata_alignements,
        "7.alignment_sequences.csv": alignment_sequences,
    }

    for filename, data in datasets.items():
        save_dataframe(data, filename)

    structures_df = pd.DataFrame(structures)
    structures_df = structures_df.fillna("")
    structures_df = structures_df.reindex(columns=[
        "interaction_id",
        "pairwise_pdb_file",
        "pymol_script_file",
        "biological_assembly_cif_file",
        "biological_assembly_pdb_file",
        "biological_assembly_info"
    ])

    structures_output = os.path.join(DETAIL_DIR, "8.structures.csv")
    structures_df.to_csv(structures_output, index=False)
    print("Structures saved:", structures_output)


if __name__ == "__main__":
    main()
