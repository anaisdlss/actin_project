"""Read-only integrity audit. Run from any directory; no data are regenerated."""
import json
import hashlib
from pathlib import Path
from datetime import datetime, timezone
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

def audit(root=ROOT):
    checks = []
    def add(name, count, detail='', warning=False):
        checks.append(dict(check=name, status=('WARN' if warning else 'FAIL') if count else 'PASS',
                           count=int(count), detail=detail))
    def read(name):
        return pd.read_csv(root / 'data/filtered' / name, low_memory=False)
    inter = read('details/1.interactions.csv')
    residues = read('details/3.interface_residues.csv')
    contacts = read('details/4.inter-residue_contacts.csv')
    proteins = read('proteins_per_pdb.csv')
    structures = read('filtered_pdb_entry.csv')
    add('Duplicate interaction identifiers', inter.interaction_id.duplicated().sum())
    ids = set(inter.interaction_id)
    for label, table in [('Interface residues',residues), ('Residue contacts',contacts)]:
        add(f'{label}: unknown interaction IDs', (~table.interaction_id.isin(ids)).sum())
    pairs = inter[['interaction_id','chain_A_id','chain_B_id']].drop_duplicates('interaction_id')
    joined = contacts.merge(pairs, on='interaction_id', how='inner', suffixes=('','_expected'))
    for side in ('A','B'):
        col = f'chain_{side}_id'
        add(f'Contact chain {side} differs from interaction source', joined[col].ne(joined[col+'_expected']).sum())
    annotated = set(proteins.chain)
    add('Interaction chains without protein annotation',
        len((set(inter.chain_A_id)|set(inter.chain_B_id))-annotated), warning=True)
    add('Conflicting actin annotations for one chain',
        (proteins.groupby('chain').is_actin.nunique()>1).sum())
    for col, table in [('buried_ASA_percent',residues),('asa_pct_A',contacts),('asa_pct_B',contacts)]:
        values = pd.to_numeric(table[col].astype(str).str.replace('%','',regex=False),errors='coerce')
        add(f'{col}: values outside 0–100', (values.notna() & ~values.between(0,100)).sum())
        add(f'{col}: missing or nonnumeric measurements', values.isna().sum(), warning=True)
    pdb_col = next(c for c in structures if c.lower() == 'pdb_id')
    listed = set(structures[pdb_col].astype(str).str.lower())
    detailed = set(inter.pdb_id.astype(str).str.lower())
    retained = set(read('filtered_all_data.csv').pdb_id.astype(str).str.lower())
    excluded = sorted(listed-retained)
    add('Structure-list PDBs outside retained interaction scope', len(excluded), ', '.join(excluded), warning=True)
    expected = set(read('filtered_summary.csv').interaction_id)
    add('Expected interaction identifiers without details', len(expected-ids))
    add('Retained PDBs without interaction details', len(retained-detailed), ', '.join(sorted(retained-detailed)))
    add('Duplicate source contact rows', contacts.duplicated().sum(), warning=True)
    for folder in ['data/proteocast/abp','data/alignments']:
        path = root/folder
        # Availability only: no freshness or scientific validity inferred.
        add(f'Optional analysis folder absent: {folder}', not path.exists(), warning=True)
    source_hashes = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in [root/'data/filtered'/f for f in [
                         'details/1.interactions.csv','details/3.interface_residues.csv',
                         'details/4.inter-residue_contacts.csv','proteins_per_pdb.csv',
                         'filtered_pdb_entry.csv']]}
    return dict(source_sha256=source_hashes, timestamp_utc=datetime.now(timezone.utc).isoformat(),
                repository=root.name, interaction_count=len(inter),
                pdbs_listed=len(listed), pdbs_with_details=len(detailed), checks=checks,
                limitation='Integrity checks do not validate biological interpretation, RSA provenance, or freshness of derived figures.')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = audit()
    text = json.dumps(report, ensure_ascii=False, indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text+'\n')
    print(text)
    raise SystemExit(any(c['status']=='FAIL' for c in report['checks']))
