"""Rebuild registered scientific outputs offline, or inspect their provenance."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'script'))
from local_rebuild import run, registry, status, inventory, atomic_json, REPORT


def preflight(root):
    """Do not reuse the internal alignment coordinate rule on a changed MSA."""
    from Bio import AlignIO
    from check_dataset import audit
    import numbering
    report = audit(root)
    atomic_json(root/'reports/dataset_integrity.json', report)
    failures = [check['check'] for check in report['checks'] if check['status'] == 'FAIL']
    if failures:
        raise ValueError('Source integrity checks failed: ' + '; '.join(failures))
    reference = ''.join(s.strip() for s in (root/'data/P60709_ref.fasta').read_text().splitlines() if not s.startswith('>'))
    alignment = AlignIO.read(str(root/'data/alignments/cluster_6685.aln'), 'fasta')
    rows = [r for r in alignment if str(r.seq).replace('-', '') == reference]
    if not rows:
        raise ValueError('P60709 reference is absent from the actin alignment.')
    for row in rows:
        position = 0
        for column, aa in enumerate(str(row.seq), 1):
            if aa != '-':
                position += 1
            expected = None if aa == '-' else position
            if numbering.to_uniprot(column) != expected:
                raise ValueError('The MAFFT reference mapping changed. Update the shared numbering rule before rebuilding; no positions were guessed.')
    print('Source integrity and P60709 alignment coordinates verified.', flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--only', nargs='+', choices=[t.name for t in registry()])
    parser.add_argument('--status', action='store_true')
    parser.add_argument('--inventory', action='store_true')
    parser.add_argument('--force', action='store_true')
    args = parser.parse_args()
    if args.status:
        for task in registry():
            print(task.name, *status(ROOT, task), sep=': ')
    elif args.inventory:
        print(f'Inventoried {len(inventory(ROOT, registry()))} CSV files')
    else:
        try:
            preflight(ROOT)
        except Exception as exc:
            atomic_json(ROOT/REPORT/'last_run.json',
                        {'state': 'failed', 'stage': 'source checks', 'error': str(exc),
                         'started_utc': datetime.now(timezone.utc).isoformat(), 'completed': []})
            raise SystemExit(str(exc))
        run(ROOT, args.only, args.force, emit=lambda line: print(line, flush=True))
