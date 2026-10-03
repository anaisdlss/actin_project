"""Only the documented, reproducible structural RSA is used by current analyses."""
from pathlib import Path
import hashlib
import json
import pandas as pd

TABLE = 'reports/scientific_audit/filament_accessibility_8dnh_B.csv'
MANIFEST = 'reports/scientific_audit/filament_accessibility_manifest.json'
CONTEXTS = {'actin_fragment': 'In the actin fragment', 'isolated': 'Chain alone'}
LABEL = 'Human beta-actin, 8DNH chain B, isolated in the same conformation'


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def source_files(root=Path('.')):
    names = [TABLE, MANIFEST, 'reports/scientific_audit/filament_accessibility_convergence.csv', 'data/P60709_ref.fasta',
             'data/reference_structures/8dnh.cif', 'data/reference_structures/8dnh_entity.json',
             'data/reference_structures/8dnh_source.json', 'tools/fetch_rsa_reference.py', 'script/filament_accessibility.py',
             'tools/calculate_filament_accessibility.py', 'script/rsa_source.py']
    return [root / name for name in names]


def load_rsa(root=Path('.')):
    """Return a checked profile and status. Never fall back to the orphan CSV."""
    root = Path(root)
    try:
        manifest = json.loads((root / MANIFEST).read_text())
        if (manifest['pdb_id'], manifest['target_author_chain']) != ('8DNH', 'B'):
            raise ValueError('unexpected reference structure')
        if manifest['taxonomy_id'] != 9606 or manifest['uniprot'] != 'P60709':
            raise ValueError('expected human ACTB reference')
        checks = dict(manifest['code_sha256'])
        checks.update(manifest['outputs_sha256'])
        for key in ('reference', 'structure', 'entity_metadata', 'download_receipt'):
            item = manifest[key]
            checks[item['file']] = item['sha256']
        for name, expected in checks.items():
            if sha256(root / name) != expected:
                raise ValueError(f'changed input, code or output: {name}')
        frame = pd.read_csv(root / TABLE)
        if frame.position.duplicated().any() or set(frame.position) != set(range(1, 376)):
            raise ValueError('incomplete or duplicate reference positions')
        reference = ''.join(s.strip() for s in (root / 'data/P60709_ref.fasta').read_text().splitlines() if not s.startswith('>'))
        frame = frame.sort_values('position')
        if ''.join(frame.aa) != reference:
            raise ValueError('residue identities differ from P60709')
        return frame, 'current'
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return pd.DataFrame(), f'RSA needs a local rebuild: {exc}'
