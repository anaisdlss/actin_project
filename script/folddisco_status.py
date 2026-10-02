"""Saved-search inventory. Never infer a negative result from missing rows."""
import json
import re
from pathlib import Path
import pandas as pd


def query_problem(record):
    chain = str(record.get('source', {}).get('source_chain', ''))
    if chain and not re.fullmatch('[A-Za-z]', chain) and not record.get('submitted_chain'):
        return 'Invalid legacy query: a numeric/non-letter chain was concatenated with residue numbers. This response is not a negative search.'
    if record.get('state') == 'invalid_query':
        return record.get('message', 'Invalid legacy query; not a negative search.')
    return None


def saved_searches(root):
    records = []
    for path in sorted(Path(root).glob('*/job.json')):
        try:
            record = json.loads(path.read_text())
            if not isinstance(record, dict):
                raise ValueError('Search manifest is not an object.')
        except (OSError, ValueError) as exc:
            records.append((path.parent, dict(state='unreadable', message=str(exc), source={})))
            continue
        problem = query_problem(record)
        if problem:
            record = {**record, 'state': 'invalid_query', 'message': problem}
        records.append((path.parent, record))
    return sorted(records, key=lambda pair: pair[1].get('updated_utc', ''), reverse=True)


def search_label(record):
    state = record.get('state', 'unknown')
    if state == 'complete':
        n = record.get('hit_rows')
        return f'{n} returned alignments' if n else 'Completed — no returned alignments'
    return {'invalid_query': 'Invalid old query — excluded', 'unsupported': 'Public-server motif limit',
            'pending': 'Running on server', 'results_pending': 'Retrieving results',
            'submission_uncertain': 'Submission could not be confirmed', 'failed': 'Server reported failure',
            'prepared': 'Prepared — not submitted', 'unreadable': 'Unreadable saved manifest'}.get(state, state)


def status_inventory(catalog, discovery, records):
    keys = set(map(tuple, catalog[['query_abp', 'query_cluster']].itertuples(index=False, name=None)))
    if not discovery.empty:
        keys.update(map(tuple, discovery[['query_abp', 'query_cluster']].itertuples(index=False, name=None)))
    keys.update((r.get('source', {}).get('query_abp'), r.get('source', {}).get('query_cluster')) for _, r in records)
    counts = discovery.groupby(['query_abp', 'query_cluster']).size() if not discovery.empty else {}
    rows = []
    for abp, site in sorted(keys, key=lambda k: (str(k[0]), str(k[1]))):
        matching = [(p, r) for p, r in records if (r.get('source', {}).get('query_abp'), r.get('source', {}).get('query_cluster')) == (abp, site)]
        valid = [(p, r) for p, r in matching if r.get('state') != 'invalid_query']
        current = (valid or matching)
        record = current[0][1] if current else None
        historical = int(counts.get((abp, site), 0))
        rows.append({'ABP': abp or 'Unknown', 'Site': site or 'Unknown',
                     'Saved search': search_label(record) if record else ('Historical hits available' if historical else 'No recorded search'),
                     'Historical rows': historical, 'Tracked alignments': record.get('hit_rows') if record and record.get('state') == 'complete' else None,
                     'Saved searches': len(matching), 'Invalid old queries': len(matching)-len(valid)})
    return pd.DataFrame(rows)
