"""Rebuild the availability report from installed scores and saved job evidence.

This validates existing external results; it does not run the ProteoCast model,
verify a UniProt accession against a live database, or assign a missing job ID.
"""
from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'script'))
from abp_profile import load_abp_scores, query_file
from proteocast_results import result_file, missing_result_reason, job_status


def availability(root):
    base = root / 'data/proteocast/abp'
    manifest = pd.read_csv(root / 'data/proteocast/abp_inputs/manifest.csv').fillna('')
    rows = []
    for row in manifest.sort_values('abp_title').to_dict('records'):
        path = result_file(base, row['slug'])
        job = job_status(base, row['slug'])
        valid, count, complete, query = False, 0, 0, False
        diagnostic = missing_result_reason(base, row['slug'])
        if path:
            try:
                _, profile = load_abp_scores(path)
                count = len(profile)
                complete = int(profile.complete_score_grid.sum())
                valid = complete == count and count > 0
                query = profile.sequence_source.eq('submitted ProteoCast query FASTA').all()
                diagnostic = ('Complete score grid; submitted query checked.' if query else
                              'Complete score grid; sequence reconstructed from mutation labels; submitted query unavailable.')
                if not valid:
                    diagnostic = f'Incomplete score grid: {complete}/{count} positions have 20 finite scores.'
            except (ValueError, KeyError, pd.errors.ParserError) as exc:
                diagnostic = 'Invalid score file: ' + str(exc)
        elif not row['uniprot']:
            diagnostic = 'No verified UniProt identifier in the manifest; no automatic submission.'
        rows.append(dict(abp=row['abp_title'], uniprot=row['uniprot'], scores_available=valid,
                         diagnostic=diagnostic, last_attempt_utc=job.get('updated_utc', ''),
                         job_id=job.get('job_id', ''), saved_job_state=job.get('state', ''),
                         score_file=str(path.relative_to(root)) if path else '',
                         positions=count, complete_positions=complete, submitted_query_checked=bool(query),
                         origin='External ProteoCast result; this report checks the installed files only'))
    return pd.DataFrame(rows)


if __name__ == '__main__':
    frame = availability(ROOT)
    output = ROOT / 'reports/scientific_audit/proteocast_availability.csv'
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False)
    print(f'ProteoCast: {int(frame.scores_available.sum())}/{len(frame)} complete score grids; no remote jobs submitted.')
