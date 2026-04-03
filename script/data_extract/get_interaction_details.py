from concurrent.futures import ThreadPoolExecutor
from bs4 import BeautifulSoup
import csv
import pandas as pd
import re
import gzip
import json
import os

from script.utils.utils import (
    create_session,
    file_hash,
    load_metadata,
    save_metadata,
    compute_outputs_hashes,
)

# -----------------------------
# CONFIG
# -----------------------------

DATA_DIR = "data/raw"
SUMMARY_META_PATH = os.path.join(DATA_DIR, "metadata.json")

REQUEST_TIMEOUT = 30
MAX_WORKERS = 10

OUTPUT_FILES = [
    "1.interactions.csv",
    "2.proteins.csv",
    "3.interface_residues.csv",
    "4.inter-residue_contacts.csv",
    "5.ligands.csv",
    "6.meta_alignement.csv",
    "7.alignment_sequences.csv",
    "8.structures.csv"
]

INTERACTIONS_COLUMNS = [
    "interaction_id",
    "pdb_id",
    "structure_title",
    "release_date",
    "resolution",
    "chain_A_id",
    "chain_B_id",
    "interface_area",
    "num_contacts",
    "num_hbonds",
    "num_salt_bridges",
    "source_url"
]

ALIGNMENT_SEQUENCE_COLUMNS = [
    "interaction_id",
    "query_sequence",
    "template_sequence",
    "interface_positions",
    "query_id",
    "query_start",
    "query_end",
    "template_id",
    "template_start",
    "template_end",
    "secondary_structure"
]

STRUCTURES_COLUMNS = [
    "interaction_id",
    "pdb_id",
    "pairwise_pdb_file",
    "pymol_script_file",
    "biological_assembly_cif_file",
    "biological_assembly_pdb_file",
    "biological_assembly_info"
]

INTERACTION_COUNT_COLUMNS = [
    "num_contacts",
    "num_hbonds",
    "num_salt_bridges"
]

# -----------------------------
# UTILS
# -----------------------------


def load_failed_ids(failed_path):
    if not os.path.exists(failed_path):
        return set()

    df = pd.read_csv(failed_path, dtype={"interaction_id": str})

    if "interaction_id" not in df.columns:
        return set()

    return set(df["interaction_id"].astype(str))


def remove_file_if_exists(path):
    if os.path.exists(path):
        os.remove(path)


def create_directories(pairwise_dir, assembly_dir, pymol_dir):
    os.makedirs(pairwise_dir, exist_ok=True)
    os.makedirs(assembly_dir, exist_ok=True)
    os.makedirs(pymol_dir, exist_ok=True)


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
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, filename)

    if os.path.exists(path):
        return path

    mode = "wb" if binary else "w"
    with open(path, mode) as f:
        f.write(content)

    return path


def get_current_job_id():
    if os.path.exists(SUMMARY_META_PATH):
        with open(SUMMARY_META_PATH) as f:
            return json.load(f).get("job_id")
    return None


def details_dataset_is_valid(
    metadata_path,
    current_job_id,
    summary_path,
    mode,
    directory,
    files
):
    metadata = load_metadata(metadata_path)

    if metadata is None:
        return False

    if not os.path.exists(summary_path):
        return False

    current_hash = compute_outputs_hashes(directory, files)

    if metadata.get("job_id") != current_job_id:
        return False

    if metadata.get("summary_hash") != file_hash(summary_path):
        return False

    if metadata.get("mode") != mode:
        return False

    if metadata.get("files_hash", {}) != current_hash:
        return False

    if len(current_hash) != len(files):
        return False

    return True


def save_details_metadata(
    metadata_path,
    job_id,
    summary_path,
    mode,
    directory,
    files
):
    metadata = {
        "job_id": job_id,
        "summary_hash": file_hash(summary_path),
        "mode": mode,
        "files_hash": compute_outputs_hashes(directory, files)
    }

    save_metadata(metadata_path, metadata)


def details_outputs_exist(directory, files):
    return all(os.path.exists(os.path.join(directory, filename)) for filename in files)


def append_rows(rows, path):
    if not rows:
        return

    df = pd.DataFrame(rows)

    if os.path.basename(path) == "1.interactions.csv":
        df = df.reindex(columns=INTERACTIONS_COLUMNS)

    if os.path.basename(path) == "7.alignment_sequences.csv":
        df = df.reindex(columns=ALIGNMENT_SEQUENCE_COLUMNS)

    if os.path.basename(path) == "8.structures.csv":
        if os.path.exists(path):
            repair_structures_csv(path)
        df = df.reindex(columns=STRUCTURES_COLUMNS)

    write_header = not os.path.exists(path)
    df.to_csv(path, mode="a", header=write_header, index=False)


def load_done_ids(progress_path):
    if not os.path.exists(progress_path):
        return set()

    df = pd.read_csv(progress_path, dtype={"interaction_id": str})

    if "interaction_id" not in df.columns:
        return set()

    return set(df["interaction_id"].astype(str))


def repair_structures_csv(path):
    if not os.path.exists(path):
        return

    rows = []

    with open(path, newline="") as f:
        reader = csv.reader(f)
        next(reader, None)

        for row in reader:
            if not row:
                continue

            structure = {column: None for column in STRUCTURES_COLUMNS}
            structure["interaction_id"] = row[0] if len(row) > 0 else None
            structure["pdb_id"] = row[1] if len(row) > 1 else None

            if len(row) > 2:
                structure["biological_assembly_info"] = row[-1]

            for value in row[2:-1]:
                if value.endswith(".py"):
                    structure["pymol_script_file"] = value
                elif "/pairwise/" in value and value.endswith(".pdb"):
                    structure["pairwise_pdb_file"] = value
                elif value.endswith(".cif"):
                    structure["biological_assembly_cif_file"] = value
                elif "/assembly/" in value and value.endswith(".pdb"):
                    structure["biological_assembly_pdb_file"] = value

            rows.append(structure)

    pd.DataFrame(rows, columns=STRUCTURES_COLUMNS).to_csv(path, index=False)


def find_incomplete_structure_ids(path):
    if not os.path.exists(path):
        return set()

    repair_structures_csv(path)
    df = pd.read_csv(path, dtype=str)

    if df.empty or "interaction_id" not in df.columns:
        return set()

    pairwise_missing = df["pairwise_pdb_file"].fillna("").str.strip() == ""
    pymol_missing = df["pymol_script_file"].fillna("").str.strip() == ""
    assembly_cif_missing = (
        df["biological_assembly_cif_file"].fillna("").str.strip() == ""
    )
    assembly_pdb_missing = (
        df["biological_assembly_pdb_file"].fillna("").str.strip() == ""
    )

    incomplete = df[
        pairwise_missing |
        pymol_missing |
        (assembly_cif_missing & assembly_pdb_missing)
    ]

    return set(incomplete["interaction_id"].astype(str))


def purge_interaction_ids_from_csv(path, interaction_ids):
    if not interaction_ids or not os.path.exists(path):
        return

    if os.path.basename(path) == "8.structures.csv":
        repair_structures_csv(path)

    df = pd.read_csv(path, dtype=str)

    if df.empty or "interaction_id" not in df.columns:
        return

    filtered_df = df[~df["interaction_id"].astype(str).isin(interaction_ids)]
    filtered_df.to_csv(path, index=False)


def requeue_interaction_ids(detail_dir, progress_path, failed_path, interaction_ids, reason):
    if not interaction_ids:
        return set()

    print(f"{reason}: {len(interaction_ids)}")

    for filename in OUTPUT_FILES:
        purge_interaction_ids_from_csv(
            os.path.join(detail_dir, filename),
            interaction_ids
        )

    purge_interaction_ids_from_csv(progress_path, interaction_ids)
    purge_interaction_ids_from_csv(failed_path, interaction_ids)

    return interaction_ids


def requeue_incomplete_structure_ids(detail_dir, progress_path, failed_path):
    structures_path = os.path.join(detail_dir, "8.structures.csv")
    incomplete_ids = find_incomplete_structure_ids(structures_path)

    return requeue_interaction_ids(
        detail_dir,
        progress_path,
        failed_path,
        incomplete_ids,
        "Requeueing interactions with incomplete required structure files"
    )


def find_invalid_interaction_metric_ids(path):
    if not os.path.exists(path):
        return set()

    df = pd.read_csv(path, dtype={"interaction_id": str})

    if df.empty or "interaction_id" not in df.columns:
        return set()

    invalid_mask = pd.Series(False, index=df.index)

    for column in INTERACTION_COUNT_COLUMNS:
        if column not in df.columns:
            continue

        values = df[column].fillna("").astype(str).str.strip()
        invalid_mask |= (values != "") & ~values.str.fullmatch(r"\d+")

    if "source_url" in df.columns:
        source_values = df["source_url"].fillna("").astype(str).str.strip()
        invalid_mask |= ~source_values.str.startswith("http")

    return set(df.loc[invalid_mask, "interaction_id"].astype(str))


def requeue_invalid_interaction_metric_ids(detail_dir, progress_path, failed_path):
    interactions_path = os.path.join(detail_dir, "1.interactions.csv")
    invalid_ids = find_invalid_interaction_metric_ids(interactions_path)

    return requeue_interaction_ids(
        detail_dir,
        progress_path,
        failed_path,
        invalid_ids,
        "Requeueing interactions with invalid interaction counts"
    )


def find_invalid_alignment_sequence_ids(path):
    if not os.path.exists(path):
        return set()

    df = pd.read_csv(path, dtype={"interaction_id": str})

    if df.empty or "interaction_id" not in df.columns:
        return set()

    required_columns = [
        "query_sequence",
        "template_sequence",
        "interface_positions",
        "query_id",
        "template_id",
    ]

    invalid_mask = pd.Series(False, index=df.index)

    for column in required_columns:
        if column not in df.columns:
            return set(df["interaction_id"].astype(str))

        values = df[column].fillna("").astype(str).str.strip()
        invalid_mask |= values == ""

    return set(df.loc[invalid_mask, "interaction_id"].astype(str))


def requeue_invalid_alignment_sequence_ids(detail_dir, progress_path, failed_path):
    path = os.path.join(detail_dir, "7.alignment_sequences.csv")
    invalid_ids = find_invalid_alignment_sequence_ids(path)

    return requeue_interaction_ids(
        detail_dir,
        progress_path,
        failed_path,
        invalid_ids,
        "Requeueing interactions with invalid alignment sequences"
    )


def find_legacy_cluster_interaction_ids(path):
    if not os.path.exists(path):
        return set()

    df = pd.read_csv(path, dtype={"interaction_id": str})

    if df.empty or "interaction_id" not in df.columns:
        return set()

    interaction_ids = df["interaction_id"].fillna("").astype(str).str.strip()
    legacy_ids = interaction_ids[interaction_ids.str.contains(
        r"\.", regex=True)]

    return set(legacy_ids)


def requeue_legacy_cluster_interaction_ids(detail_dir, progress_path, failed_path):
    interactions_path = os.path.join(detail_dir, "1.interactions.csv")
    legacy_ids = find_legacy_cluster_interaction_ids(interactions_path)

    return requeue_interaction_ids(
        detail_dir,
        progress_path,
        failed_path,
        legacy_ids,
        "Requeueing legacy cluster interaction ids"
    )


def normalize_text(value):
    return " ".join(str(value).split())


def extract_metric_value(soup, lines, label):
    normalized_label = normalize_text(label).lower()
    tag_names = ["p", "div", "li", "td", "span", "strong", "b"]

    for tag in soup.find_all(tag_names):
        tag_text = normalize_text(tag.get_text(" ", strip=True))
        if normalized_label not in tag_text.lower():
            continue

        match = re.search(
            rf"{re.escape(label)}\s*:?\s*(\d+)\b",
            tag_text,
            flags=re.IGNORECASE
        )
        if match:
            return match.group(1)

    for idx, line in enumerate(lines):
        normalized_line = normalize_text(line)
        if normalized_label not in normalized_line.lower():
            continue

        match = re.search(
            rf"{re.escape(label)}\s*:?\s*(\d+)\b",
            normalized_line,
            flags=re.IGNORECASE
        )
        if match:
            return match.group(1)

        for offset in range(1, 4):
            if idx + offset >= len(lines):
                break

            candidate = normalize_text(lines[idx + offset])
            if re.fullmatch(r"\d+", candidate):
                return candidate

    return None


def normalize_csv(path, sort_columns=None):
    if not os.path.exists(path):
        return

    if os.path.basename(path) == "1.interactions.csv":
        df = pd.read_csv(path, dtype={"interaction_id": str})

        if df.empty:
            return

        df = df.reindex(columns=INTERACTIONS_COLUMNS)
        df = df.drop_duplicates(subset=["interaction_id"])

        if sort_columns:
            existing_sort_columns = [
                c for c in sort_columns if c in df.columns]
            if existing_sort_columns:
                df = df.sort_values(existing_sort_columns)

        df.to_csv(path, index=False)
        return

    if os.path.basename(path) == "7.alignment_sequences.csv":
        df = pd.read_csv(path, dtype={"interaction_id": str})

        if df.empty:
            return

        df = df.reindex(columns=ALIGNMENT_SEQUENCE_COLUMNS)
        df = df.drop_duplicates()

        if sort_columns:
            existing_sort_columns = [c for c in sort_columns if c in df.columns]
            if existing_sort_columns:
                df = df.sort_values(existing_sort_columns)

        df.to_csv(path, index=False)
        return

    if os.path.basename(path) == "8.structures.csv":
        repair_structures_csv(path)

    df = pd.read_csv(path, dtype={"interaction_id": str})

    if df.empty:
        return

    df = df.drop_duplicates()

    if sort_columns:
        existing_sort_columns = [c for c in sort_columns if c in df.columns]
        if existing_sort_columns:
            df = df.sort_values(existing_sort_columns)

    df.to_csv(path, index=False)


def normalize_outputs(detail_dir):
    normalize_csv(
        os.path.join(detail_dir, "1.interactions.csv"),
        ["interaction_id"]
    )
    normalize_csv(
        os.path.join(detail_dir, "2.proteins.csv"),
        ["interaction_id", "chain_id"]
    )
    normalize_csv(
        os.path.join(detail_dir, "3.interface_residues.csv"),
        ["interaction_id", "chain", "residue_number_structure"]
    )
    normalize_csv(
        os.path.join(detail_dir, "4.inter-residue_contacts.csv"),
        ["interaction_id", "chain_A_id", "residue_A_structure",
            "chain_B_id", "residue_B_structure"]
    )
    normalize_csv(
        os.path.join(detail_dir, "5.ligands.csv"),
        ["interaction_id", "chain", "residue_number"]
    )
    normalize_csv(
        os.path.join(detail_dir, "6.meta_alignement.csv"),
        ["interaction_id"]
    )
    normalize_csv(
        os.path.join(detail_dir, "7.alignment_sequences.csv"),
        ["interaction_id"]
    )
    normalize_csv(
        os.path.join(detail_dir, "8.structures.csv"),
        ["interaction_id"]
    )

# -----------------------------
# WORKER
# -----------------------------


def process_interaction(task):
    (
        global_index,
        pass_index,
        row,
        interaction_id,
        is_cluster,
        pairwise_dir,
        assembly_dir,
        pymol_dir,
        pass_total
    ) = task

    session = create_session()
    url = row["detail_url"]

    print(f"Downloading {pass_index}/{pass_total}")

    try:
        r = session.get(url, timeout=REQUEST_TIMEOUT)
        r.raise_for_status()
        html = r.text

        soup = BeautifulSoup(html, "lxml")
        tables = soup.find_all("table")

        links = soup.find_all("a", href=True)
        h4_tags = soup.find_all("h4")

        text = soup.get_text(separator="\n", strip=True)
        lines = text.split("\n")

        interaction = {"interaction_id": interaction_id}
        current_id = str(interaction_id)
        safe_id = current_id.replace(".", "_")

        interactions = []
        proteins = []
        interface_residues = []
        residue_contacts = []
        ligands = []
        metadata_alignements = []
        alignment_sequences = []
        structures = []

        # -----------------------------
        # table 1 interaction + table 2 proteins
        # -----------------------------

        for idx, line in enumerate(lines):
            line = line.strip()

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

            if line.startswith("Chain A") or line.startswith("Chain B"):
                protein = {}
                protein["interaction_id"] = current_id
                protein["chain_id"] = line.split(":")[-1].strip()
                protein["protein_name"] = None
                protein["organism"] = None
                protein["num_residues"] = None

                next_idx = idx + 1
                while next_idx < len(lines):
                    next_line = lines[next_idx].strip()

                    if next_line.startswith("Chain A") or next_line.startswith("Chain B"):
                        break

                    if next_line.startswith("Title:"):
                        protein["protein_name"] = next_line.replace(
                            "Title:", ""
                        ).strip()
                    elif next_line.startswith("Source organism:"):
                        protein["organism"] = next_line.replace(
                            "Source organism:", ""
                        ).strip()
                    elif next_line.startswith("Number of residues"):
                        residues = next_line.split(":")[-1].strip()
                        match = re.search(r"\d+", residues)
                        protein["num_residues"] = int(
                            match.group()) if match else None

                    next_idx += 1

                proteins.append(protein)

        metric_labels = {
            "num_contacts": "Number of inter-residue contacts at the interface",
            "num_hbonds": "Number of H-bonds",
            "num_salt_bridges": "Number of salt bridges",
        }

        for field, label in metric_labels.items():
            value = extract_metric_value(soup, lines, label)
            if value is not None:
                interaction[field] = value

        interaction["source_url"] = url
        interactions.append(interaction)

        # -----------------------------
        # table 3 interface residues
        # -----------------------------

        interface_tables = []

        for table in tables:
            headers = [th.get_text(strip=True) for th in table.find_all("th")]
            if "Residue no. in structure" in headers and "Buried ASA, %" in headers:
                interface_tables.append((headers, table))

        chains = []
        for tag in h4_tags:
            if "Interface residues in" in tag.text:
                chains.append(tag.text.split("in")[-1].strip())

        for chain, (headers, table) in zip(chains, interface_tables):
            rows = table.find_all("tr")

            for r in rows:
                cols = [c.get_text(strip=True)
                        for c in r.find_all(["th", "td"])]

                if not cols:
                    continue

                if len(cols) == len(headers) - 1 and headers[0].startswith("No."):
                    row_table = dict(zip(headers[1:], cols))
                else:
                    if len(cols) != len(headers):
                        continue

                    row_table = dict(zip(headers, cols))

                if (
                    "Residue no. in structure" not in row_table
                    or "Residue no. in sequence" not in row_table
                    or "Residue name" not in row_table
                    or "Buried ASA, Å2" not in row_table
                    or "Buried ASA, %" not in row_table
                ):
                    continue

                residue = {}
                residue["interaction_id"] = current_id
                residue["chain"] = chain
                residue["residue_number_structure"] = row_table["Residue no. in structure"]
                residue["residue_number_sequence"] = row_table["Residue no. in sequence"]
                residue["residue_name"] = row_table["Residue name"]
                residue["buried_ASA_Å²"] = row_table["Buried ASA, Å2"]
                residue["buried_ASA_percent"] = row_table["Buried ASA, %"]

                interface_residues.append(residue)

        # -----------------------------
        # table 4 residue contacts
        # -----------------------------

        for table in tables:
            headers = [th.get_text(strip=True) for th in table.find_all("th")]
            header_text = " | ".join(headers)

            if "Residue no. in chain A structure" in header_text and "Residue in chain B" in header_text:
                rows = table.find_all("tr")

                for r in rows:
                    cols = [c.get_text(strip=True) for c in r.find_all("td")]

                    if len(cols) != len(headers):
                        continue

                    row_table = dict(zip(headers, cols))

                    contact = {}
                    contact["interaction_id"] = current_id
                    contact["chain_A_id"] = interaction["chain_A_id"]
                    contact["residue_A_structure"] = row_table["Residue no. in chain A structure"]
                    contact["residue_A_sequence"] = row_table["Residue no. in chain A sequence"]
                    contact["residue_A_name"] = row_table["Residue in chain A"]
                    contact["chain_B_id"] = interaction["chain_B_id"]
                    contact["residue_B_structure"] = row_table["Residue no. in chain B structure"]
                    contact["residue_B_sequence"] = row_table["Residue no. in chain B sequence"]
                    contact["residue_B_name"] = row_table["Residue in chain B"]
                    contact["contact_area"] = row_table["Contact area, Å2"]
                    contact["contact_type"] = row_table.get(
                        "Contact type", None)

                    residue_contacts.append(contact)

        # -----------------------------
        # table 5 ligands
        # -----------------------------

        ligand_tables = []

        for table in tables:
            headers = [th.get_text(strip=True) for th in table.find_all("th")]
            header_text = " | ".join(headers)

            if (
                "Chain" in headers
                and "Protein" in headers
                and "Residue no." in headers
                and "Ligand name" in headers
                and "Contact area with same domain" in header_text
                and "Contact area with other domain" in header_text
            ):
                ligand_tables.append((headers, table))

        for headers, table in ligand_tables:
            rows = table.find_all("tr")

            for r in rows:
                cols = [c.get_text(strip=True) for c in r.find_all(["th", "td"])]

                if not cols or all(value.lower() == "filter" for value in cols):
                    continue

                if len(cols) != len(headers):
                    continue

                row_table = dict(zip(headers, cols))

                required_keys = [
                    "Chain",
                    "Protein",
                    "Residue no.",
                    "Ligand name",
                    "Contact area with same domain",
                    "Contact area with other domain",
                ]

                if any(key not in row_table for key in required_keys):
                    continue

                if row_table["Chain"] == "Chain":
                    continue

                ligand = {}
                ligand["interaction_id"] = current_id
                ligand["chain"] = row_table["Chain"]
                ligand["protein"] = row_table["Protein"]
                ligand["residue_number"] = row_table["Residue no."]
                ligand["ligand_name"] = row_table["Ligand name"]
                ligand["contact_area_same_domain"] = row_table["Contact area with same domain"]
                ligand["contact_area_other_domain"] = row_table["Contact area with other domain"]

                ligands.append(ligand)

        # -----------------------------
        # table 6 metadata alignement
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
                metadata["e_value"] = line.split("=")[1].split(",")[0].strip()
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

        for a in links:
            href = a["href"]

            if "sequence_alignments" in href and href.endswith(".fasta"):
                fasta_text = download_text(session, href)
                if fasta_text is None:
                    continue

                fasta_lines = [l.strip()
                               for l in fasta_text.split("\n") if l.strip()]
                if len(fasta_lines) < 6:
                    continue

                query_header = fasta_lines[0]
                template_header = fasta_lines[2]

                alignment_seq["query_sequence"] = fasta_lines[1]
                alignment_seq["template_sequence"] = fasta_lines[3]
                alignment_seq["interface_positions"] = fasta_lines[5]
                alignment_seq["query_id"] = query_header.split("|")[1]

                m = re.search(r"starts at (\d+)", query_header)
                if m:
                    alignment_seq["query_start"] = int(m.group(1))

                m = re.search(r"ends at (\d+)", query_header)
                if m:
                    alignment_seq["query_end"] = int(m.group(1))

                alignment_seq["template_id"] = template_header.replace(">", "").split()[
                    0]

                m = re.search(r"starts at (\d+)", template_header)
                if m:
                    alignment_seq["template_start"] = int(m.group(1))

                m = re.search(r"ends at (\d+)", template_header)
                if m:
                    alignment_seq["template_end"] = int(m.group(1))

                break

        aligned_block = soup.find(id="aligned_sequence")
        if aligned_block:
            dssp_parts = []
            for p in aligned_block.find_all("p"):
                text_p = p.get_text("", strip=False).replace("\xa0", " ")
                if text_p.startswith("dssp:"):
                    dssp_parts.append(text_p[len("dssp:"):].rstrip("\r\n"))
            alignment_seq["secondary_structure"] = "".join(dssp_parts)

        alignment_sequences.append(alignment_seq)

        # -----------------------------
        # table 8 structures
        # -----------------------------

        structure = {
            "interaction_id": current_id,
            "pdb_id": interaction.get("pdb_id", "")
        }

        for a in links:
            href = a["href"]

            try:
                if "pdb_files" in href and href.endswith(".pdb"):
                    pdb_text = download_text(session, href)
                    if pdb_text is None:
                        print(
                            f"Required structure file missing: pairwise pdb -> {href}")
                        continue

                    pdb_id = interaction["pdb_id"]
                    chainA = interaction["chain_A_id"].split("_")[-1]
                    chainB = interaction["chain_B_id"].split("_")[-1]
                    filename = f"{safe_id}_{pdb_id}_{chainA}_{chainB}.pdb"
                    file_path = save_file(pdb_text, pairwise_dir, filename)
                    structure["pairwise_pdb_file"] = file_path

                if href.endswith(".py"):
                    pymol_script = download_text(session, href)
                    if pymol_script is None:
                        print(
                            f"Required structure file missing: pymol script -> {href}")
                        continue

                    pdb_id = interaction["pdb_id"]
                    chainA = interaction["chain_A_id"].split("_")[-1]
                    chainB = interaction["chain_B_id"].split("_")[-1]
                    filename = f"{safe_id}_{pdb_id}_{chainA}_{chainB}.py"
                    script_path = save_file(pymol_script, pymol_dir, filename)
                    structure["pymol_script_file"] = script_path

                if "biounits" in href:
                    data = download_binary(session, href)
                    if data is None:
                        print(
                            f"Biological assembly file missing: {href}"
                        )
                        continue

                    content = gzip.decompress(data).decode()

                    if "mmCIF" in href or content.startswith("data_"):
                        pdb_id = interaction["pdb_id"]
                        filename = f"{safe_id}_{pdb_id}.cif"
                        cif_path = save_file(content, assembly_dir, filename)
                        structure["biological_assembly_cif_file"] = cif_path

                    if "pdb" in href:
                        pdb_id = interaction["pdb_id"]
                        filename = f"{safe_id}_{pdb_id}.pdb"
                        pdb_path = save_file(content, assembly_dir, filename)
                        structure["biological_assembly_pdb_file"] = pdb_path

            except Exception as e:
                print(f"Structure download skipped: {href} -> {e}")

        assembly_info = None
        h4 = soup.find("h4", string=lambda x: x and "Biological assembly" in x)
        if h4:
            div = h4.find_next("div")
            if div:
                assembly_info = div.get_text(strip=True)

        structure["biological_assembly_info"] = assembly_info
        structures.append(structure)

        missing_required_files = []

        if not structure.get("pairwise_pdb_file"):
            missing_required_files.append("pairwise_pdb_file")

        if not structure.get("pymol_script_file"):
            missing_required_files.append("pymol_script_file")

        if not structure.get("biological_assembly_cif_file") and not structure.get(
            "biological_assembly_pdb_file"
        ):
            missing_required_files.append(
                "biological_assembly_cif_file|biological_assembly_pdb_file"
            )

        if missing_required_files:
            print(
                f"Required structure files missing for interaction {current_id}: "
                + ", ".join(missing_required_files)
            )
            return {
                "interaction_id": current_id,
                "detail_url": url,
                "status": "failed"
            }

        return {
            "interaction_id": current_id,
            "detail_url": url,
            "status": "done",
            "interactions": interactions,
            "proteins": proteins,
            "interface_residues": interface_residues,
            "residue_contacts": residue_contacts,
            "ligands": ligands,
            "metadata_alignements": metadata_alignements,
            "alignment_sequences": alignment_sequences,
            "structures": structures,
        }

    except Exception as e:
        print(
            f"Error on interaction {interaction_id} (row {global_index + 1}): {e}")
        return {
            "interaction_id": str(interaction_id),
            "detail_url": url,
            "status": "failed"
        }

# -----------------------------
# MAIN
# -----------------------------


def main():
    mode = input(
        "Use filtered summary (f), clusters summary (y), or default summary ([n]): "
    ).strip().lower()

    if mode == "y":
        summary_path = os.path.join(
            DATA_DIR, "clusters", "clusters_summary.csv")
        detail_dir = os.path.join(DATA_DIR, "clusters", "details")
        is_cluster = True
        mode_name = "cluster"
    elif mode == "f":
        summary_path = os.path.join("data", "filtered", "filtered_summary.csv")
        detail_dir = os.path.join("data", "filtered", "details")
        is_cluster = False
        mode_name = "filtered"
    else:
        summary_path = os.path.join(DATA_DIR, "ppi3d_actin_summary.csv")
        detail_dir = os.path.join(DATA_DIR, "details")
        is_cluster = False
        mode_name = "summary"

    detail_meta_path = os.path.join(detail_dir, "metadata.json")
    progress_path = os.path.join(detail_dir, "progress.csv")
    failed_path = os.path.join(detail_dir, "failed.csv")

    pairwise_dir = os.path.join(detail_dir, "structures_files", "pairwise")
    assembly_dir = os.path.join(detail_dir, "structures_files", "assembly")
    pymol_dir = os.path.join(detail_dir, "structures_files", "pymol")

    interactions_path = os.path.join(detail_dir, "1.interactions.csv")
    proteins_path = os.path.join(detail_dir, "2.proteins.csv")
    interface_residues_path = os.path.join(
        detail_dir, "3.interface_residues.csv")
    residue_contacts_path = os.path.join(
        detail_dir, "4.inter-residue_contacts.csv")
    ligands_path = os.path.join(detail_dir, "5.ligands.csv")
    metadata_alignements_path = os.path.join(
        detail_dir, "6.meta_alignement.csv")
    alignment_sequences_path = os.path.join(
        detail_dir, "7.alignment_sequences.csv")
    structures_path = os.path.join(detail_dir, "8.structures.csv")

    os.makedirs(detail_dir, exist_ok=True)
    create_directories(pairwise_dir, assembly_dir, pymol_dir)

    current_job_id = get_current_job_id()
    existing_meta = load_metadata(detail_meta_path)

    if details_dataset_is_valid(
        detail_meta_path,
        current_job_id,
        summary_path,
        mode_name,
        detail_dir,
        OUTPUT_FILES
    ):
        print("Details dataset unchanged")
        print("Nothing to do")
        return

    if (
        existing_meta is not None
        and os.path.exists(summary_path)
        and existing_meta.get("job_id") == current_job_id
        and existing_meta.get("summary_hash") == file_hash(summary_path)
        and existing_meta.get("mode") == mode_name
        and details_outputs_exist(detail_dir, OUTPUT_FILES)
    ):
        print("Details files already present")
        print("Refreshing metadata hashes without re-downloading")
        save_details_metadata(
            detail_meta_path,
            current_job_id,
            summary_path,
            mode_name,
            detail_dir,
            OUTPUT_FILES
        )
        print("Details metadata refreshed")
        print("Nothing to do")
        return

    if is_cluster:
        requeue_legacy_cluster_interaction_ids(
            detail_dir, progress_path, failed_path)
    requeue_incomplete_structure_ids(detail_dir, progress_path, failed_path)
    requeue_invalid_interaction_metric_ids(
        detail_dir, progress_path, failed_path)
    requeue_invalid_alignment_sequence_ids(
        detail_dir, progress_path, failed_path)

    sep = "," if mode_name == "filtered" else ";"
    summary = pd.read_csv(summary_path, sep=sep)
    print("Using summary file:", summary_path)

    if "detail_url" not in summary.columns:
        raise ValueError("detail_url column missing in summary file")

    total = len(summary)
    records = summary.to_dict("records")

    current_pass = 1

    while True:
        done_ids = load_done_ids(progress_path)
        remove_file_if_exists(failed_path)

        pending_records = []

        for i, row in enumerate(records):
            if is_cluster:
                raw_interaction_id = row["Link to details"].split(" ")[
                    0].strip()
                interaction_id = raw_interaction_id
            else:
                interaction_id = str(i + 1)

            if interaction_id in done_ids:
                continue

            pending_records.append(
                (
                    i,
                    row,
                    interaction_id,
                )
            )

        tasks = []
        remaining_total = len(pending_records)

        for pass_index, (global_index, row, interaction_id) in enumerate(
            pending_records,
            start=1
        ):
            tasks.append(
                (
                    global_index,
                    pass_index,
                    row,
                    interaction_id,
                    is_cluster,
                    pairwise_dir,
                    assembly_dir,
                    pymol_dir,
                    remaining_total
                )
            )

        print(
            f"\n===== PASS {current_pass} | REMAINING {remaining_total} =====")
        print("Already done:", len(done_ids))
        print("Remaining interactions to process:", remaining_total)

        if not tasks:
            print("No remaining interactions.")
            break

        progress_rows = []
        failed_rows = []

        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            for result in executor.map(process_interaction, tasks):
                if result is None:
                    continue

                if result["status"] == "failed":
                    failed_rows.append({
                        "interaction_id": result["interaction_id"],
                        "detail_url": result["detail_url"],
                        "status": "failed"
                    })

                    if len(failed_rows) >= 50:
                        append_rows(failed_rows, failed_path)
                        failed_rows = []

                    continue

                append_rows(result["interactions"], interactions_path)
                append_rows(result["proteins"], proteins_path)
                append_rows(result["interface_residues"],
                            interface_residues_path)
                append_rows(result["residue_contacts"], residue_contacts_path)
                append_rows(result["ligands"], ligands_path)
                append_rows(result["metadata_alignements"],
                            metadata_alignements_path)
                append_rows(result["alignment_sequences"],
                            alignment_sequences_path)
                append_rows(result["structures"], structures_path)

                progress_rows.append({
                    "interaction_id": result["interaction_id"],
                    "detail_url": result["detail_url"],
                    "status": "done"
                })

                if len(progress_rows) >= 50:
                    append_rows(progress_rows, progress_path)
                    progress_rows = []

        if progress_rows:
            append_rows(progress_rows, progress_path)

        if failed_rows:
            append_rows(failed_rows, failed_path)

        failed_ids = load_failed_ids(failed_path)
        print("Failed after this pass:", len(failed_ids))

        if len(failed_ids) == 0:
            print("All interactions processed successfully.")
            break
        current_pass += 1

    normalize_outputs(detail_dir)

    save_details_metadata(
        detail_meta_path,
        current_job_id,
        summary_path,
        mode_name,
        detail_dir,
        OUTPUT_FILES
    )

    remove_file_if_exists(progress_path)

    print("Details metadata saved")


if __name__ == "__main__":
    main()
