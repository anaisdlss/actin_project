"""Publish a checked local snapshot into an existing public checkout.

Dry run by default. Only explicitly identified calculation work files are removed;
source structures, displayed results, alignments and public-only imported bundles
are retained. This does not submit jobs or deploy the application.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def omitted(relative):
    parts = Path(relative).parts
    if any(p in {'__pycache__', '.job_status'} for p in parts):
        return 'runtime or local job state'
    if Path(relative).suffix in {'.pyc', '.lock', '.log', '.pse'} or '.DS_Store' in parts:
        return 'runtime file or desktop session'
    if any(p.startswith('_sweep_pad') or p == 'folddisco_index' for p in parts):
        return 'local search index or temporary structural comparison'
    if str(relative).startswith('data/exports/folddisco_controls/') and (
            any(p in {'raw', 'chains'} for p in parts) or Path(relative).name.startswith('index')):
        return 'local FoldDisco working files; result tables retained'
    if str(relative).startswith('data/filtered/details/structures_files/') and any(
            p in {'pairwise', 'pymol', 'filament'} for p in parts):
        return 'local structural work files; assemblies retained'
    return None


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def duplicate_copy(path):
    """Only remove Finder-style duplicates when the original is byte-identical."""
    if ' 2.' not in path.name:
        return False
    original = path.with_name(path.name.replace(' 2.', '.'))
    return original.is_file() and sha(original) == sha(path)


def plan(source, destination):
    selected = {}
    for directory in ('data', 'reports/scientific_audit', 'reports/local_rebuild'):
        for path in (source/directory).rglob('*'):
            if path.is_file() and not path.is_symlink():
                relative = path.relative_to(source).as_posix()
                if not omitted(relative) and not duplicate_copy(path):
                    selected[relative] = path
    for name in ('reports/dataset_integrity.json', 'reports/s1_figures_all.json'):
        if (source/name).is_file():
            selected[name] = source/name
    # Never clean an untracked user file, or remove an unknown extra result.
    tracked = subprocess.check_output(['git', '-C', str(destination), 'ls-files', '-z']).decode().split('\0')
    removed = [r for r in tracked if r.startswith('data/') and (destination/r).is_file()
               and (omitted(r) or duplicate_copy(destination/r))]
    conflicts = []
    tracked_set = set(tracked)
    for relative, path in selected.items():
        target = destination/relative
        if target.is_file() and relative not in tracked_set and sha(path) != sha(target):
            conflicts.append(relative)
    if conflicts:
        raise ValueError('Preserving conflicting untracked destination files: '+', '.join(conflicts))
    return selected, removed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    source, destination = args.source.resolve(), args.destination.resolve()
    if source == destination or source in destination.parents or destination in source.parents:
        parser.error('Use distinct, non-nested source and destination checkouts.')
    if not (destination/'data/.slim_deploy').is_file():
        parser.error('Destination must already be the public slim checkout.')
    sys.path.insert(0, str(source/'tools'))
    sys.path.insert(0, str(source/'script'))
    from check_dataset import audit
    from local_rebuild import registry, status
    failures = [c['check'] for c in audit(source)['checks'] if c['status'] == 'FAIL']
    stale = [t.name for t in registry() if status(source, t)[0] != 'current']
    if failures or stale:
        raise ValueError(f'Incomplete source checks: {failures}; calculations needing rebuild: {stale}')
    selected, removed = plan(source, destination)
    before = sum(p.stat().st_size for p in (destination/'data').rglob('*') if p.is_file())
    removable = sum((destination/r).stat().st_size for r in removed)
    summary = dict(source_files=len(selected), remove_files=len(removed),
                   remove_MB=round(removable/1e6, 1), before_MB=round(before/1e6, 1), applied=args.apply)
    if not args.apply:
        print(json.dumps(summary, indent=2))
        return
    hashes = {}
    for relative, path in sorted(selected.items()):
        target = destination/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        expected = sha(path)
        if not target.is_file() or sha(target) != expected:
            temporary = target.with_name(target.name+'.sync-tmp')
            shutil.copy2(path, temporary)
            if sha(temporary) != expected:
                raise ValueError('Copy verification failed: '+relative)
            temporary.replace(target)
        hashes[relative] = expected
    for relative in removed:
        (destination/relative).unlink()
    (destination/'data/.slim_deploy').write_text('1\n')
    after = sum(p.stat().st_size for p in (destination/'data').rglob('*') if p.is_file())
    summary.update(after_MB=round(after/1e6, 1), created_utc=datetime.now(timezone.utc).isoformat(),
                   source_commit=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip(),
                   source_working_tree_has_changes=bool(subprocess.check_output(
                       ['git','-C',str(source),'status','--porcelain'],text=True).strip()),
                   files_sha256=hashes, removed={r:omitted(r) or 'byte-identical duplicate copy' for r in removed},
                   scope='Current local scientific data and successful calculation receipts. Public-only imported bundles and historical results retained; local executable indexes omitted.')
    output = destination/'reports/cloud_snapshot.json'
    output.write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in {'files_sha256','removed'}},indent=2))


if __name__ == '__main__':
    main()
