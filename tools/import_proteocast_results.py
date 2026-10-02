"""Reuse verified existing scores from another project with matching accession and query."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys

ROOT=Path(__file__).resolve().parents[1]
os.chdir(ROOT);sys.path.insert(0,str(ROOT/'script'))
import pandas as pd
from abp_profile import load_abp_scores
from proteocast_results import result_file


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--from-project',required=True);args=parser.parse_args()
    source=Path(args.from_project).resolve()
    target_manifest=pd.read_csv(ROOT/'data/proteocast/abp_inputs/manifest.csv').fillna('')
    source_manifest=pd.read_csv(source/'data/proteocast/abp_inputs/manifest.csv').fillna('')
    audit=[]
    for row in target_manifest.to_dict('records'):
        dest=ROOT/'data/proteocast/abp'/row['slug']
        if result_file(dest.parent,row['slug']):continue
        if not row['uniprot']:continue
        candidates=source_manifest[source_manifest.uniprot.eq(row['uniprot']) & source_manifest.slug.eq(row['slug'])]
        for candidate in candidates.to_dict('records'):
            file=result_file(source/'data/proteocast/abp',candidate['slug'])
            if file is None:continue
            try:
                raw,profile=load_abp_scores(file)
                if not profile.complete_score_grid.all():raise ValueError('Incomplete score grid')
                target_query=dest/'1.query.fasta'
                if target_query.exists() and target_query.stat().st_size:
                    load_abp_scores(file,query_path=target_query)
                dest.mkdir(parents=True,exist_ok=True)
                shutil.copy2(file,dest/'4.query_ProteoCast.csv')
                source_query=file.parent/'1.query.fasta'
                if source_query.exists() and source_query.stat().st_size and not target_query.exists():shutil.copy2(source_query,target_query)
                audit.append(dict(abp=row['abp_title'],uniprot=row['uniprot'],source_project=source.name,
                    file=str(file.relative_to(source)),sha256=hashlib.sha256(file.read_bytes()).hexdigest(),positions=len(profile),scores=len(raw),state='imported'))
            except ValueError as exc:audit.append(dict(abp=row['abp_title'],state='rejected',reason=str(exc)))
    out=ROOT/'reports/scientific_audit/proteocast_import.json';out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(audit,indent=2));print(pd.Series([r['state'] for r in audit]).value_counts().to_string())

if __name__=='__main__':main()
