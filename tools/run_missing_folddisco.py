"""Run previously missing motifs with explicit state and complete query provenance."""
import argparse
import json
import os
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
os.chdir(ROOT);sys.path.insert(0,str(ROOT/'script'))
import pandas as pd
from folddisco_audit import motif_catalog, structure_path, CATALOG_FILES
from folddisco_jobs import prepare, submit, refresh, JOBS


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--limit',type=int,default=1)
    parser.add_argument('--wait-seconds',type=int,default=300)
    args=parser.parse_args()
    catalog=motif_catalog(*(pd.read_csv(p) for p in CATALOG_FILES))
    saved=pd.read_csv('data/exports/abp_site_domain/folddisco_discovery.csv')
    keys=set(zip(saved.query_abp,saved.query_cluster))
    missing=catalog[[tuple(r) not in keys for r in catalog[['query_abp','query_cluster']].itertuples(index=False,name=None)]]
    audit=[]
    for row in missing.head(args.limit).to_dict('records'):
        print(row['query_abp'],row['query_cluster'],flush=True)
        path=structure_path(row)
        try:
            if path is None:raise ValueError('No local structure')
            folder,record=prepare(row,row['contact_positions'],path)
            record=submit(folder)
            deadline=time.monotonic()+args.wait_seconds
            while record['state'] in {'pending','results_pending'} and time.monotonic()<deadline:
                time.sleep(5);record=refresh(folder)
            outcome=dict(abp=row['query_abp'],site=row['query_cluster'],state=record['state'],message=record.get('message'),request_id=record['request_id'])
        except (ValueError,OSError) as exc:
            outcome=dict(abp=row['query_abp'],site=row['query_cluster'],state='unavailable',message=str(exc))
        audit.append(outcome);print(outcome,flush=True)
        if outcome['state'] in {'submission_uncertain','failed','pending','results_pending'}:
            print('Stopping campaign after an unresolved external result; existing ticket retained.',flush=True);break
    out=ROOT/'reports/scientific_audit/folddisco_new_campaign.json';out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(audit,indent=2))

if __name__=='__main__':main()
