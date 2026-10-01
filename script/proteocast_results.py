"""Read result availability and preserve ProteoCast job diagnostics locally."""
import json
from datetime import datetime, timezone
from pathlib import Path


def result_file(base, slug):
    base = Path(base)
    for path in (base / slug / '4.query_ProteoCast.csv', base / f'{slug}.csv'):
        if path.is_file() and path.stat().st_size > 0:
            return path
    return None


def job_status(base, slug):
    path = Path(base) / '.job_status' / f'{slug}.json'
    try:
        return json.loads(path.read_text())
    except (FileNotFoundError, ValueError):
        return {}


def record_job(base, slug, state, message='', **details):
    record = job_status(base, slug)
    record.update(slug=slug, state=state, message=message,
                  updated_utc=datetime.now(timezone.utc).isoformat())
    record.update(details)
    path = Path(base) / '.job_status' / f'{slug}.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(record, ensure_ascii=False, indent=2) + '\n')
    temporary.replace(path)


def missing_result_reason(base, slug):
    if result_file(base, slug):
        return ''
    record = job_status(base, slug)
    if record.get('state') == 'failed' and record.get('message'):
        return record['message']
    folder = Path(base) / slug
    alignments = list(folder.glob('2.ali*.fasta'))
    if alignments and all(p.stat().st_size == 0 for p in alignments):
        return ('The downloaded alignment (MSA) is empty; no variant scores were '
                'produced. Check alignment retrieval or provide a valid MSA before retrying.')
    if folder.is_dir() and any(folder.iterdir()):
        return 'Only partial files are present; the variant-score CSV is missing.'
    return 'No ProteoCast scores are available yet.'
