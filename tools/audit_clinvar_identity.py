"""Verify cached ClinVar records against official ESummary titles (no reclassification).
The old gene search can return multi-gene records whose protein change is on another gene.
"""
from pathlib import Path
import csv,json,urllib.request,time,datetime,hashlib
ROOT=Path(__file__).resolve().parents[1];base=ROOT/'data/human_variants'
ids=set()
for p in (base/'sources').glob('clinvar_*.csv'):
    with p.open() as f:ids.update(r['clinvar_id'] for r in csv.DictReader(f))
ids=sorted(ids,key=int);result={}
cache=base/'clinvar_identity.json'
# Continue a partially completed audit, preserving previously fetched records.
if cache.exists():result=json.loads(cache.read_text()).get('records',{})
missing=[x for x in ids if x not in result]
for i in range(0,len(missing),150):
    subset=missing[i:i+150]
    url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&retmode=json&id='+','.join(subset)
    req=urllib.request.Request(url,headers={'User-Agent':'ActinABP-research-source-audit'})
    with urllib.request.urlopen(req,timeout=60) as response:payload=json.load(response)
    fetched=datetime.datetime.now(datetime.timezone.utc).isoformat()
    for uid in subset:
        record=payload.get('result',{}).get(uid)
        if record is None or record.get('error'):raise ValueError(f'ClinVar record {uid} not returned; retry audit.')
        result[uid]={'title':record.get('title',''),'genes':record.get('genes',[]),'checked_utc':fetched,
                     'current_classification':record.get('germline_classification',{}).get('description',''),
                     'last_evaluated':record.get('germline_classification',{}).get('last_evaluated','')}
    cache.write_text(json.dumps({'provider':'NCBI ClinVar ESummary','purpose':'Identity validation; current classifications are retained for comparison, original annotations remain unchanged.','records':result},indent=2))
    print(f'Checked {len(result)}/{len(ids)} IDs',flush=True);time.sleep(.35)
print('Identity audit complete')
manifest_path=base/'manifest.json'
if manifest_path.exists():
    manifest=json.loads(manifest_path.read_text())
    manifest['identity_audit']={'file':'data/human_variants/clinvar_identity.json',
                                'sha256':hashlib.sha256(cache.read_bytes()).hexdigest(),
                                'records_checked':len(result)}
    manifest_path.write_text(json.dumps(manifest,indent=2))
