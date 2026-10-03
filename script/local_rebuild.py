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
    packages: tuple = ('numpy', 'pandas', 'biopython', 'scipy', 'matplotlib')


def registry():
    details = 'data/filtered/details/'
    tables = ['data/filtered/filtered_all_data.csv', details+'1.interactions.csv',
              details+'3.interface_residues.csv']
    reference = ['data/P60709_ref.fasta', 'data/alignments/cluster_6685.aln']
    structures = [details+'structures_files/assembly/*.pdb', details+'structures_files/assembly/*.cif']
    report = 'reports/scientific_audit/'
    abp = 'data/exports/abp_site_domain/'
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
             ['data/P60709_ref.fasta', 'data/reference_structures/8dnh.cif',
              'data/reference_structures/8dnh_entity.json', 'data/reference_structures/8dnh_source.json'],
             [report+'filament_accessibility_8dnh_B.csv', report+'filament_accessibility_convergence.csv',
              report+'filament_accessibility_manifest.json'], ['script/filament_accessibility.py', 'tools/fetch_rsa_reference.py']),
        task('abp_catalog', 'script/abp_site_domain/01_build_table.py', tables[:2],
             ['data/exports/abp_site_domain/abp_site_table.csv', 'data/exports/abp_site_domain/abp_representatives.csv'], []),
        task('structure_annotations', 'tools/update_structure_annotations.py',
             [tables[0], 'data/filtered/proteins_per_pdb.csv', 'data/annotations/rcsb_structure_metadata.json'],
             ['data/annotations/structure_annotations.csv', 'data/annotations/entity_annotations.csv',
              'data/annotations/manifest.json'], ['script/structure_annotations.py'], args=['--offline']),
        task('proteocast_availability', 'tools/audit_proteocast_availability.py',
             ['data/proteocast/abp_inputs/manifest.csv', 'data/proteocast/abp/*.csv',
              'data/proteocast/abp/*/*', 'data/proteocast/abp/.job_status/*.json'],
             [report+'proteocast_availability.csv'], ['script/abp_profile.py', 'script/proteocast_results.py']),
        task('abp_cross_tables', 'script/abp_site_domain/04_crosstab.py',
             [abp+n for n in ('abp_site_table.csv', 'abp_representatives.csv', 'abp_interpro.csv', 'fold_cluster_cluster.tsv')],
             [abp+n for n in ('abp_master.csv', 'site_cluster_summary.csv', 'figure_site_vs_domaine.png')],
             [], ['abp_catalog']),
        task('abp_families', 'script/abp_site_domain/11_family_regroup.py',
             [abp+n for n in ('abp_interpro.csv', 'abp_master.csv', 'abp_site_table.csv', 'whole_pairs_all.tsv')],
             [abp+n for n in ('familles.csv', 'convergences_inter_familles.csv', 'figure_convergence_familles.png')],
             [], ['abp_cross_tables']),
        task('abp_footprints', 'script/abp_site_domain/13_actin_footprint_overlap.py',
             [*tables, *reference, abp+'familles.csv'],
             [abp+'actin_footprint_overlap.csv', abp+'figure_footprint_overlap.png'], [], ['abp_families']),
        task('abp_chemistry', 'script/abp_site_domain/16_interface_chemistry.py',
             [*tables, abp+'abp_representatives.csv'], [abp+'interface_chemistry.csv'], [], ['abp_catalog']),
        task('abp_determinants', 'script/abp_site_domain/18_site_determinants.py',
             [*tables, *reference, abp+'familles.csv', 'data/proteocast/actin/*',
              report+'filament_accessibility_8dnh_B.csv', report+'filament_accessibility_manifest.json'],
             [abp+'actin_residue_determinants.csv', abp+'figure_site_determinants.png'],
             ['script/scientific_analysis.py', 'script/footprint_comparison.py', 'script/rsa_source.py'],
             ['rsa', 'abp_families']),
        task('representative_geometry', 'tools/check_interface_geometry.py',
             [tables[0], 'data/P60709_ref.fasta', *[details+'structures_files/assembly/'+p+'.pdb' for p in ('3j8a','5yu8','6vao')]],
             [report+'representative_geometry'+s for s in ('.csv','_summary.csv','_manifest.json')],
             ['script/scientific_analysis.py']),
        task('scientific_tables', 'tools/export_scientific_audit.py',
             [*tables, *reference, details+'4.inter-residue_contacts.csv', 'data/filtered/proteins_per_pdb.csv',
              'data/human_variants/sources/*.csv', 'data/human_variants/ref/*.fasta',
              'data/human_variants/clinvar_identity.json', 'data/proteocast/actin/*',
              report+'filament_accessibility_8dnh_B.csv', report+'filament_accessibility_manifest.json'],
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
            if rel in produced:
                action = 'Rebuilt automatically when its sources or code change'
            elif rel == 'data/proteocast/conservation_vs_asa_per_position.csv':
                action = 'Archived, excluded from current RSA; replaced by documented local calculation'
            elif rel.endswith('/rsa_values.csv') and '/proteocast/abp/' in rel:
                action = 'Supplementary imported archive; not read by the current application; no replacement invented'
            elif rel == 'data/proteocast/abp_inputs/manifest.csv':
                action = 'Keep reviewed accession choices; do not regenerate ambiguous identities from a protein name'
            elif state.startswith('producer recipe found'):
                action = 'Existing source pipeline; outside the offline result rebuild (historical execution not certified)'
            elif state.startswith(('external', 'imported')):
                action = 'Preserve external snapshot; local summaries use these identified inputs without recreating the model/database'
            else:
                action = 'Not automatically regenerated: insufficient provenance'
            rows.append(dict(file=rel, status=state, producer=producer, sha256=checksum, note=note,
                             automatic_action=action))
    out = root / REPORT / 'csv_inventory.csv'
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=['file','status','producer','sha256','note','automatic_action'], lineterminator='\n')
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
