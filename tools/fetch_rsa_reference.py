"""Download the explicitly selected human ACTB reference from RCSB.

Run once before the offline rebuild, or explicitly refresh the deposited files.
The calculation checks these saved files; it never silently downloads a new model.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import requests

ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / 'data/reference_structures'
URLS = {
    '8dnh.cif': 'https://files.rcsb.org/download/8DNH.cif',
    '8dnh_entity.json': 'https://data.rcsb.org/rest/v1/core/polymer_entity/8DNH/1',
}


def validate_entity(entity):
    identifiers = entity['rcsb_polymer_entity_container_identifiers']
    organisms = entity['rcsb_entity_source_organism']
    if (identifiers['entry_id'] != '8DNH' or identifiers['entity_id'] != '1'
            or identifiers['uniprot_ids'] != ['P60709']
            or {o['ncbi_taxonomy_id'] for o in organisms} != {9606}
            or set(identifiers['auth_asym_ids']) != set('ABCD')):
        raise ValueError('The selected reference must be human ACTB (P60709), 8DNH entity 1.')


def main():
    downloaded = {}
    for name, url in URLS.items():
        response = requests.get(url, timeout=90)
        response.raise_for_status()
        downloaded[name] = response.content
    validate_entity(json.loads(downloaded['8dnh_entity.json']))
    if not downloaded['8dnh.cif'].startswith(b'data_8DNH'):
        raise ValueError('Unexpected coordinate response')
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    for name, content in downloaded.items():
        (DIRECTORY / name).write_bytes(content)
    manifest = {'retrieved_utc': datetime.now(timezone.utc).isoformat(),
                'pdb_id': '8DNH', 'entity_id': '1', 'organism': 'Homo sapiens',
                'uniprot': 'P60709', 'files': {
                    name: {'url': URLS[name], 'sha256': hashlib.sha256(content).hexdigest()}
                    for name, content in downloaded.items()},
                'reproduce': 'python tools/fetch_rsa_reference.py'}
    (DIRECTORY / '8dnh_source.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print('Downloaded human ACTB reference 8DNH and its RCSB metadata.')


if __name__ == '__main__':
    main()
