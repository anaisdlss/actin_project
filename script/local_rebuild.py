"""Content-checked, resumable offline calculations over a fixed source snapshot.

The registry is deliberately explicit. Unknown historical files are reported,
never assigned an invented producer. A successful build is not scientific validation.
"""
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from contextlib import contextmanager
import csv
import fcntl
import hashlib
import importlib.metadata
import json
import os
import shutil
import subprocess
import sys
import tempfile

REPORT = 'reports/local_rebuild'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')
    temporary.replace(path)


@dataclass
class Task:
    name: str
    command: list
    inputs: list
    outputs: list
    code: list
    dependencies: list = field(default_factory=list)
    packages: tuple = ('numpy', 'pandas', 'biopython', 'scipy')


def registry():
    details = 'data/filtered/details/'
    tables = ['data/filtered/filtered_all_data.csv', details+'1.interactions.csv',
              details+'3.interface_residues.csv']
    reference = ['data/P60709_ref.fasta', 'data/alignments/cluster_6685.aln']
    structures = [details+'structures_files/assembly/*.pdb', details+'structures_files/assembly/*.cif']
    report = 'reports/scientific_audit/'
    profiles = ['variant_mapping', 'variant_exclusions', 'conservation_profiles',
                'conservation_correlations', 'conservation_footprints', 'interface_sites',
                'interface_c70', 'interface_occurrences', 'interface_context',
                'cross_gene_observations', 'cross_gene_pair_counts', 'cross_gene_abp_summaries',
                'cross_gene_abp_details', 'partner_chemistry_positions', 'partner_chemistry_per_abp']
    # Explicit code closure includes shared numerical and display utilities. This
    # conservatively invalidates a result after a utility edit, never misses one.
    common_code = ['script/numbering.py', 'script/residue_contacts.py', 'script/display_helpers.py']
    def task(name, script, inputs, outputs, code, dependencies=(), args=()):
        return Task(name, [sys.executable, '-u', script, *args], inputs, outputs,
                    [script, 'script/local_rebuild.py', *common_code, *code], list(dependencies))
    return [
        task('rsa', 'tools/calculate_filament_accessibility.py',
             ['data/P60709_ref.fasta', 'data/filtered/proteins_per_pdb.csv',
              details+'structures_files/assembly/7pdz.pdb'],
             [report+'filament_accessibility_7pdz_I.csv', report+'filament_accessibility_convergence.csv',
              report+'filament_accessibility_manifest.json'], ['script/filament_accessibility.py']),
        task('abp_catalog', 'script/abp_site_domain/01_build_table.py', tables[:2],
             ['data/exports/abp_site_domain/abp_site_table.csv', 'data/exports/abp_site_domain/abp_representatives.csv'], []),
        task('scientific_tables', 'tools/export_scientific_audit.py',
             [*tables, *reference, details+'4.inter-residue_contacts.csv', 'data/filtered/proteins_per_pdb.csv',
              'data/human_variants/sources/*.csv', 'data/human_variants/ref/*.fasta',
              'data/human_variants/clinvar_identity.json', 'data/proteocast/actin/*',
              report+'filament_accessibility_7pdz_I.csv', report+'filament_accessibility_manifest.json'],
             [*[report+name+'.csv' for name in profiles], report+'manifest.json'],
             ['script/scientific_analysis.py', 'script/footprint_comparison.py',
              'script/interface_properties.py', 'script/rsa_source.py'], ['rsa']),
        task('interface_geometry', 'tools/audit_all_interface_geometry.py',
             [tables[0], *reference, *structures],
             [report+'all_interface_geometry'+suffix for suffix in
              ('.csv', '_coverage.csv', '_summary.csv', '_manifest.json')],
             ['script/interface_geometry.py', 'script/scientific_analysis.py']),
        task('filament_proximity', 'tools/audit_steric_proximity.py',
             [*tables, *reference, *structures, 'data/exports/abp_site_domain/abp_site_table.csv'],
             [report+'rigid_filament_proximity.csv', report+'rigid_filament_proximity_manifest.json'],
             ['script/steric_screen.py', 'script/scientific_analysis.py', 'script/folddisco_audit.py'], ['abp_catalog']),
        task('figures', 'tools/regenerate_scientific_figures.py', [*tables, *reference],
             ['reports/s1_figures_all.json', 'data/visualisations/actin_s1_*.png',
              'data/visualisations/actin_s1_*_p60709.csv', 'data/visualisations/actin_s1_clusters/*.png',
              'data/visualisations/actin_s1_clusters_by_c70/*.png', 'data/filtered/actin_s1_canon_area_by_cluster.csv'],
             ['script/s1_heatmaps.py', 'script/plot_interaction.py']),
        task('folddisco_controls', 'tools/audit_folddisco_controls.py',
             ['data/exports/abp_site_domain/abp_site_table.csv', tables[2], *structures],
             ['data/exports/folddisco_controls/queries.csv', 'data/exports/folddisco_controls/hits.csv',
              'data/exports/folddisco_controls/shared_site_controls.csv', 'data/exports/folddisco_controls/manifest.json',
              report+'folddisco_shared_site_controls.csv', report+'folddisco_local_queries.csv', report+'folddisco_controls_manifest.json'],
             ['script/folddisco_audit.py', 'script/folddisco_jobs.py'], ['abp_catalog']),
    ]


def files(root, patterns, required=False):
    paths = set()
    for pattern in patterns:
        found = [p for p in root.glob(pattern) if p.is_file()]
        if not found and required and not any(c in pattern for c in '*?['):
            raise FileNotFoundError(f'Missing required source: {pattern}')
        paths.update(found)
    return sorted(paths)


def hashes(root, patterns, required=False):
    return {str(p.relative_to(root)): digest(p) for p in files(root, patterns, required)}


def fingerprint(root, task):
    runtime = {'python': sys.version.split()[0]}
    for package in task.packages:
        runtime[package] = importlib.metadata.version(package)
    if task.name == 'folddisco_controls':
        binary = Path(sys.executable).parent / 'folddisco'
        runtime['folddisco_binary_sha256'] = digest(binary)
    return {'inputs': hashes(root, task.inputs, True), 'code': hashes(root, task.code, True),
            'runtime': runtime, 'command': task.command[1:]}


def status(root, task):
    try:
        current = fingerprint(root, task)
        manifest = json.loads((root / REPORT / (task.name+'.json')).read_text())
        if manifest['fingerprint'] != current:
            return 'needs rebuild', 'Sources, code, parameters or software changed'
        if not manifest['outputs'] or manifest['outputs'] != hashes(root, task.outputs, True):
            return 'needs rebuild', 'Output missing or changed'
        return 'current', 'Input, code and output checksums verified'
    except FileNotFoundError as exc:
        return ('needs rebuild', 'No successful tracked build') if REPORT in str(exc) else ('blocked', str(exc))
    except (ValueError, KeyError, importlib.metadata.PackageNotFoundError) as exc:
        return 'needs rebuild', str(exc)


@contextmanager
def build_lock(root):
    path = root / REPORT / 'build.lock'
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('A local rebuild is already running in this project.')
        yield


def run_task(root, task, emit=print, force=False):
    state, reason = status(root, task)
    if state == 'current' and not force:
        emit(f'{task.name}: current — verified, skipped')
        return
    before = fingerprint(root, task)
    previous = files(root, task.outputs)
    # Last good tables are restored on failure; an interrupted process leaves no
    # new success receipt and therefore cannot be reported as current.
    with tempfile.TemporaryDirectory(prefix='actin-rebuild-') as temporary:
        backup = Path(temporary)
        for path in previous:
            dest = backup / path.relative_to(root)
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, dest)
        log = root / REPORT / (task.name+'.log')
        log.parent.mkdir(parents=True, exist_ok=True)
        emit(f'{task.name}: running ({reason})')
        try:
            with log.open('w') as stream:
                process = subprocess.Popen(task.command, cwd=root, stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT, text=True)
                for line in process.stdout:
                    stream.write(line); stream.flush()
                    emit(line.rstrip())
                process.stdout.close()
                if process.wait():
                    raise RuntimeError(f'{task.name} failed; see {log.relative_to(root)}')
            if fingerprint(root, task) != before:
                raise RuntimeError('Sources changed during calculation; results were not accepted.')
            result = hashes(root, task.outputs, True)
            if not result:
                raise RuntimeError('Calculation produced no registered outputs.')
            atomic_json(root / REPORT / (task.name+'.json'),
                        {'finished_utc': datetime.now(timezone.utc).isoformat(),
                         'fingerprint': before, 'outputs': result})
        except BaseException:
            for path in files(root, task.outputs):
                path.unlink()
            for path in previous:
                path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(backup / path.relative_to(root), path)
            raise


def inventory(root, tasks):
    """Inventory every CSV without mistaking a checksum for scientific provenance."""
    produced = {}
    for task in tasks:
        state, _ = status(root, task)
        for path in files(root, task.outputs):
            produced[str(path.relative_to(root))] = (task.name, state)
    imported = {}
    source_manifest = root/'data/human_variants/manifest.json'
    if source_manifest.exists():
        for item in json.loads(source_manifest.read_text()).get('sources', []):
            imported[item['file']] = item
    from source_recipes import recovered_origin
    rows = []
    for base in ('data', 'reports'):
        for path in sorted((root/base).rglob('*.csv')):
            rel = str(path.relative_to(root))
            if rel.startswith(REPORT+'/') or '/figure_backups/' in rel:
                continue
            checksum = digest(path)
            producer, state = produced.get(rel, ('', 'unresolved historical provenance'))
            note = 'Producer and original source must be established; file presence is not validation.'
            if producer:
                note = 'Locally reproducible from the inputs listed in its build receipt; upstream source limitations still apply.'
            elif rel in imported:
                state = ('imported snapshot; release date unknown' if checksum == imported[rel]['sha256']
                         else 'imported snapshot changed since import')
                producer = imported[rel].get('source_project', 'external snapshot')
                note = 'An imported database/model result, not a local measurement. Import time is not database release time.'
            elif rel == 'data/proteocast/conservation_vs_asa_per_position.csv':
                note = 'Legacy RSA origin unresolved; excluded from current RSA values and sensitivity summaries.'
            if state == 'unresolved historical provenance':
                recovered = recovered_origin(root, rel)
                if recovered:
                    state, producer, note = recovered
            rows.append(dict(file=rel, status=state, producer=producer, sha256=checksum, note=note))
    out = root / REPORT / 'csv_inventory.csv'
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=['file','status','producer','sha256','note'])
        writer.writeheader(); writer.writerows(rows)
    return rows


def run(root, selected=None, force=False, emit=print):
    root = Path(root).resolve()
    tasks = registry()
    names = {t.name for t in tasks}
    wanted = set(selected or names)
    if wanted - names:
        raise ValueError(f'Unknown calculations: {wanted-names}')
    for _ in tasks:
        wanted |= {d for t in tasks if t.name in wanted for d in t.dependencies}
    with build_lock(root):
        state = {'started_utc': datetime.now(timezone.utc).isoformat(), 'state': 'running',
                 'requested': sorted(wanted), 'completed': []}
        atomic_json(root/REPORT/'last_run.json', state)
        try:
            for task in tasks:
                if task.name in wanted:
                    run_task(root, task, emit, force)
                    state['completed'].append(task.name)
                    atomic_json(root/REPORT/'last_run.json', state)
            inventory(root, tasks)
            state['state'] = 'complete'
        except Exception as exc:
            state.update(state='failed', error=str(exc))
            raise
        finally:
            atomic_json(root/REPORT/'last_run.json', state)
