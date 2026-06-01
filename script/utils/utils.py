import os
import json
import pandas as pd
import hashlib
import requests

REQUEST_TIMEOUT = 30


def create_session():
    session = requests.Session()
    session.headers.update({
        "User-Agent": "ppi3d-actin-scraper/1.0"
    })

    adapter = requests.adapters.HTTPAdapter(
        pool_connections=50,
        pool_maxsize=50
    )

    session.mount("https://", adapter)

    return session


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

    mode = "wb" if binary else "w"

    with open(path, mode) as f:
        f.write(content)

    return path


def save_dataframe(data, path):

    df = pd.DataFrame(data)
    df.to_csv(path, index=False)


def file_hash(path):

    with open(path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def load_metadata(path):

    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)

    return None


def save_metadata(path, data):

    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def compute_outputs_hashes(directory, files):

    hashes = {}

    for f in files:

        path = os.path.join(directory, f)

        if os.path.exists(path):
            hashes[f] = file_hash(path)

    return hashes


def dataset_is_valid(metadata_path, current_job_id, directory, files):

    metadata = load_metadata(metadata_path)

    if metadata is None:
        return False

    stored_job = metadata.get("job_id")
    stored_hash = metadata.get("files_hash", {})

    current_hash = compute_outputs_hashes(directory, files)

    if stored_job != current_job_id:
        return False

    if stored_hash != current_hash:
        return False

    if len(current_hash) != len(files):
        return False

    return True


def save_dataset_metadata(metadata_path, job_id, directory, files):

    metadata = {
        "job_id": job_id,
        "files_hash": compute_outputs_hashes(directory, files)
    }

    save_metadata(metadata_path, metadata)
