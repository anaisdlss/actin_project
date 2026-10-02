"""Reproducible FoldDisco submissions for single or combined observed motifs."""
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path
import hashlib
import json
import math
import re

import pandas as pd
import requests

from folddisco_audit import residue_token, residue_order, read_chain_ca, validate_motif

BASE = 'https://search.foldseek.com'
DATABASES = ['pdb_folddisco', 'afdb-proteome_folddisco']
JOBS = Path('data/exports/folddisco_jobs')


def combine_sites(sites, residues, abp, clusters):
    """Choose one observed PDB chain containing every requested site; never splice structures."""
    clusters = sorted(set(clusters))
    sub = sites[sites.abp_title.eq(abp) & sites.actin_site_cluster.isin(clusters)].copy()
    if not clusters or sub.empty:
        raise ValueError('Select at least one observed binding site.')
    sub['resolution_A'] = pd.to_numeric(sub.resolution.astype(str).str.extract(r'(\d+(?:\.\d+)?)', expand=False), errors='coerce')
    candidates = []
    for (pdb, chain), rows in sub.groupby(['pdb', 'abp_chain']):
        if set(rows.actin_site_cluster) != set(clusters):
            continue
        ids = sorted(set(rows.interaction_id.astype(int)))
        source = residues[residues.interaction_id.isin(ids) & residues.chain.eq(f'{str(pdb).lower()}_{chain}')]
        positions = [residue_token(v) for v in source.residue_number_structure]
        if not positions or any(v is None for v in positions):
            continue
        # Every site must have residue observations in this same chain.
        if any(not source.interaction_id.isin(rows.loc[rows.actin_site_cluster.eq(c), 'interaction_id']).any() for c in clusters):
            continue
        positions = sorted(set(positions), key=residue_order)
        resolution = rows.resolution_A.min()
        candidates.append(dict(query_abp=abp, query_cluster='+'.join(clusters), clusters=clusters,
            source_pdb=str(pdb).lower(), source_chain=str(chain), interaction_id=ids[0], interaction_ids=ids,
            resolution_A=float(resolution) if pd.notna(resolution) else None,
            contact_positions=', '.join(positions), reconstructed_size=len(positions), unparsed_contact_positions=''))
    if not candidates:
        raise ValueError('No single observed PDB chain covers all selected sites. Choose another combination; coordinates from different structures are not merged.')
    return min(candidates, key=lambda r: (r['resolution_A'] if r['resolution_A'] is not None else math.inf, r['source_pdb'], r['source_chain']))


def chain_pdb_text(path, chain):
    from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
    from Bio.PDB.Polypeptide import protein_letters_3to1_extended
    class Chain(Select):
        def accept_chain(self, c): return c.id == chain
        def accept_residue(self, r): return r.resname in protein_letters_3to1_extended
        def accept_atom(self, a): return a.element not in {'H', 'D'} and (not a.is_disordered() or a.get_altloc() in {' ', 'A'})
    if len(chain) != 1:
        raise ValueError('The FoldDisco PDB export requires a one-character chain identifier.')
    parser = MMCIFParser(QUIET=True) if Path(path).suffix == '.cif' else PDBParser(QUIET=True)
    model = next(iter(parser.get_structure('query', str(path))))
    if chain not in model: raise ValueError('Selected chain is absent from the structure.')
    io = PDBIO(); io.set_structure(model); target = StringIO(); io.save(target, Chain())
    return target.getvalue()


def prepare(row, positions, path, root=JOBS):
    text = chain_pdb_text(path, row['source_chain'])
    coords = read_chain_ca(path, row['source_chain'])
    checked = validate_motif(positions, coords.position, row['source_chain'])
    if checked['errors']: raise ValueError('; '.join(checked['errors']))
    if not checked['query']: raise ValueError('This query contains unsupported insertion codes or chain identifiers.')
    sha = hashlib.sha256(text.encode()).hexdigest()
    request = dict(structure_sha256=sha, motif=checked['query'], databases=DATABASES)
    identity = hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest()
    folder = Path(root) / identity
    folder.mkdir(parents=True, exist_ok=True)
    file = folder / 'job.json'
    if file.exists(): return folder, json.loads(file.read_text())
    source = {k: None if isinstance(v, float) and not math.isfinite(v) else v for k, v in row.items()}
    record = dict(**request, request_id=identity, source=source, positions=checked['positions'],
                  state='prepared', endpoint=BASE, created_utc=datetime.now(timezone.utc).isoformat())
    (folder/'query.pdb').write_text(text)
    save(folder, record)
    return folder, record


def save(folder, record):
    record['updated_utc'] = datetime.now(timezone.utc).isoformat()
    temp = Path(folder)/'job.tmp';temp.write_text(json.dumps(record, indent=2, allow_nan=False));temp.replace(Path(folder)/'job.json')


def submit(folder, session=None):
    folder = Path(folder); record = json.loads((folder/'job.json').read_text())
    if record['state'] != 'prepared': return record
    if len(record['positions']) > 32:
        record.update(state='unsupported', message='The public FoldDisco server accepts at most 32 motif residues. Edit the selection or use a local FoldDisco index; no automatic truncation is performed.')
        save(folder, record); return record
    session = session or requests.Session()
    record.update(state='submitting', message='Submission started; do not submit the same request twice.')
    save(folder, record)
    try:
        with (folder/'query.pdb').open('rb') as f:
            response = session.post(f'{BASE}/api/ticket/folddisco', files={'q': f},
                data=[('database[]', db) for db in DATABASES]+[('motif', record['motif'])], timeout=45)
        response.raise_for_status(); payload = response.json()
        ticket = payload.get('id')
        if not ticket or not re.fullmatch(r'[A-Za-z0-9_-]+', str(ticket)):
            raise ValueError(f'No valid ticket in server response: {str(payload)[:250]}')
        record.update(state='pending', ticket=str(ticket), message='Submitted to FoldDisco.')
    except (requests.RequestException, ValueError) as exc:
        record.update(state='submission_uncertain', message=str(exc))
    save(folder, record); return record



def retry(folder, session=None):
    """Explicit retry of an unresolved submission; keep the previous attempt."""
    folder = Path(folder); record = json.loads((folder/'job.json').read_text())
    if record['state'] not in {'failed', 'submission_uncertain'}:
        return record
    record.setdefault('previous_attempts', []).append({k: record.get(k) for k in ('state', 'ticket', 'message', 'updated_utc')})
    record['state'] = 'prepared'
    save(folder, record)
    return submit(folder, session)


def alignment_rows(value):
    """The API may group alignments by hit or return a flat legacy list."""
    if isinstance(value, dict):
        if 'target' in value:
            return [value]
        return [row for group in value.values() for row in alignment_rows(group)]
    if isinstance(value, list):
        return [row for group in value for row in alignment_rows(group)]
    if value is None:
        return []
    raise ValueError('Unrecognized alignment result format.')


def refresh(folder, session=None):
    folder = Path(folder); record = json.loads((folder/'job.json').read_text())
    if record['state'] not in {'pending', 'results_pending'}: return record
    session = session or requests.Session()
    try:
        response = session.get(f"{BASE}/api/ticket/{record['ticket']}", timeout=30)
        response.raise_for_status(); payload = response.json(); status = payload.get('status')
        record.update(server_status=status, message=f'Server status: {status}')
        if status in {'ERROR', 'UNKNOWN'}: record.update(state='failed')
        elif status == 'COMPLETE':
            record['state'] = 'results_pending'
            result = session.get(f"{BASE}/api/result/folddisco/{record['ticket']}", timeout=45)
            result.raise_for_status(); body = result.json()
            if not isinstance(body.get('results'), list): raise ValueError('Missing result list in server response.')
            # Preserve original data, including empty completed searches.
            (folder/'results.json').write_text(json.dumps(body))
            count = sum(len(alignment_rows(r.get('alignments'))) for r in body['results'])
            record.update(state='complete', hit_rows=count, message=f'Completed: {count} returned alignments.')
    except (requests.RequestException, ValueError) as exc:
        record['message'] = f'Status/results retrieval failed; the existing ticket is retained: {exc}'
    save(folder, record); return record


def result_table(folder):
    folder = Path(folder); path = folder/'results.json'
    if not path.exists(): return pd.DataFrame()
    record = json.loads((folder/'job.json').read_text()); rows=[]
    for result in json.loads(path.read_text()).get('results', []):
        for hit in alignment_rows(result.get('alignments')):
            row = {'database': result.get('db'), **hit}
            row['query_size'] = len(record['positions'])
            try: row['coverage'] = float(hit['nodecount']) / len(record['positions'])
            except (KeyError, ValueError, TypeError): row['coverage'] = None
            rows.append(row)
    return pd.DataFrame(rows)
