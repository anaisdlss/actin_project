"""Preserve invalid numeric-chain searches, then submit corrected chain exports.

Run from the full project; bounded sequential retries stop on unresolved service
errors. No old query or response is deleted. The >32-residue motif stays unsupported.
"""
import argparse
import json
import os
from pathlib import Path
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
os.chdir(ROOT);sys.path.insert(0,str(ROOT/'script'))
from folddisco_jobs import JOBS,prepare,submit,refresh,save
from folddisco_status import query_problem
from folddisco_audit import structure_path


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--limit',type=int,default=1)
    args=parser.parse_args();out=[];n=0
    for file in sorted(JOBS.glob('*/job.json')):
        old=json.loads(file.read_text())
        if not query_problem(old):continue
        if old['state'] != 'invalid_query':
            old['previous_outcome']={k:old.get(k) for k in ['state','message','hit_rows','updated_utc']}
            old.update(state='invalid_query',message=query_problem(old))
            save(file.parent,old)
        if n>=args.limit:continue
        path=structure_path(old['source'])
        folder,new=prepare(old['source'],','.join(old['positions']),path)
        new['supersedes_invalid_request']=old['request_id'];save(folder,new)
        new=submit(folder);deadline=time.monotonic()+180
        while new['state'] in {'pending','results_pending'} and time.monotonic()<deadline:
            time.sleep(5);new=refresh(folder)
        item=dict(abp=old['source']['query_abp'],site=old['source']['query_cluster'],
                  old_request=old['request_id'],new_request=new['request_id'],
                  state=new['state'],hit_rows=new.get('hit_rows'))
        out.append(item);print(item,flush=True);n+=1
        if new['state'] not in {'complete','unsupported'}:break
    (ROOT/'reports/scientific_audit/folddisco_query_repair.json').write_text(json.dumps(out,indent=2))

if __name__=='__main__':main()
