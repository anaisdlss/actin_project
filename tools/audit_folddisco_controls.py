"""Search every current motif against observed ABP chains and audit known contacts.

This local, finite positive-control panel is not a PDB/AlphaFold discovery search.
Contact validation uses exact PDB/chain/residue identifiers, not inferred domains.
"""
import hashlib
import io
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from datetime import datetime, timezone
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
os.chdir(ROOT);sys.path.insert(0,str(ROOT/'script'))
from folddisco_audit import CATALOG_FILES,motif_catalog,structure_path,read_chain_ca,validate_motif,residue_token,query_chain
from folddisco_jobs import chain_pdb_text

OUT=ROOT/'data/exports/folddisco_controls'


def fingerprint(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    binary=Path(sys.executable).parent/'folddisco'
    if not binary.exists():raise SystemExit('Install folddisco in the research environment first.')
    sites,residues=(pd.read_csv(p,low_memory=False) for p in CATALOG_FILES)
    catalog=motif_catalog(sites,residues)
    OUT.mkdir(parents=True,exist_ok=True);target_dir=OUT/'chains';target_dir.mkdir(exist_ok=True)
    targets=[];missing=[]
    for (pdb,chain),group in catalog.groupby(['source_pdb','source_chain']):
        row=group.iloc[0].to_dict();path=structure_path(row)
        try:
            if path is None:raise ValueError('No local structure')
            text=chain_pdb_text(path,chain)
            # A fresh index contains only files explicitly listed in this run.
            key=hashlib.sha256(f'{pdb}:{chain}'.encode()).hexdigest()[:16]
            dest=target_dir/f'{key}.pdb';dest.write_text(text)
            targets.append(dict(key=key,pdb=pdb,chain=chain,submitted_chain=query_chain(chain),
                abps=sorted(group.query_abp.unique()),file=str(dest),sha256=fingerprint(dest)))
        except (ValueError,OSError) as exc:missing.append(dict(pdb=pdb,chain=chain,error=str(exc)))
    wanted={Path(t['file']).name for t in targets}
    stale=[p.name for p in target_dir.glob('*.pdb') if p.name not in wanted]
    # Only generated representative-chain files in this tool-owned directory.
    for name in stale:(target_dir/name).unlink()
    index=OUT/'index';cmd=[str(binary),'index','-p',str(target_dir),'-i',str(index),'-t','4']
    subprocess.run(cmd,check=True,capture_output=True,text=True)
    target_map={t['key']:t for t in targets}
    lookup={(t['pdb'],t['chain']):t for t in targets}
    contacts={}
    raw={key:set(filter(None,(residue_token(p) for p in g.residue_number_structure)))
         for key,g in residues.groupby(['interaction_id','chain'])}
    for r in sites.itertuples():
        key=(str(r.pdb).lower(),str(r.abp_chain),r.actin_site_cluster)
        contacts.setdefault(key,set()).update(raw.get((r.interaction_id,f'{r.pdb}_{r.abp_chain}'),set()))
    hits=[];queries=[];controls=[]
    raw_dir=OUT/'raw';raw_dir.mkdir(exist_ok=True)
    for number,row in enumerate(catalog.to_dict('records'),1):
        abp,site=row['query_abp'],row['query_cluster'];source=lookup.get((row['source_pdb'],row['source_chain']))
        status=dict(query_abp=abp,query_cluster=site,query_size=row['reconstructed_size'],source_pdb=row['source_pdb'],source_chain=row['source_chain'])
        if source is None:
            queries.append({**status,'state':'missing_structure'});continue
        coords=read_chain_ca(source['file'],source['submitted_chain'])
        checked=validate_motif(row['contact_positions'],coords.position,source['submitted_chain'])
        if checked['errors'] or not checked['query']:
            queries.append({**status,'state':'unsupported_query','reason':'; '.join(checked['errors']) or 'Unsupported residue identifiers'});continue
        cmd=[str(binary),'query','-p',source['file'],'-q',checked['query'],'-i',str(index),'-t','2',
             '--header','--per-match','--format-output','tid,node_count,idf,rmsd,matching_residues,query_residues']
        result=subprocess.run(cmd,capture_output=True,text=True,timeout=120)
        if result.returncode:raise RuntimeError(f'{abp}/{site}: {result.stderr}')
        qkey=hashlib.sha256(f'{abp}:{site}'.encode()).hexdigest()[:16]
        (raw_dir/f'{qkey}.txt').write_text(result.stdout)
        frame=pd.read_csv(io.StringIO(result.stdout),sep='\t')
        # Keep the highest IDF match per observed target chain, with its own geometry.
        best=frame.sort_values(['idf','rmsd'],ascending=[False,True]).drop_duplicates('tid')
        found=[]
        for r in best.to_dict('records'):
            target=target_map[Path(r['tid']).stem]
            matched={m[1:] for m in str(r['matching_residues']).split(',') if re.fullmatch('[A-Za-z][0-9]+',m)}
            observed=contacts.get((target['pdb'],target['chain'],site),set())
            any_site=set().union(*(v for (p,c,_),v in contacts.items() if (p,c)==(target['pdb'],target['chain'])))
            for name in target['abps']:
                hit={**status,'target_abp':name,'target_pdb':target['pdb'],'target_chain':target['chain'],
                     'matched_residues':int(r['node_count']),'coverage':float(r['node_count'])/len(checked['positions']),
                     'idfscore':r['idf'],'rmsd_A':r['rmsd'],'target_positions':', '.join(sorted(matched,key=int)),
                     'observed_contact_count_same_site':len(observed),'matched_contact_count_same_site':len(matched&observed),
                     'fraction_matched_on_observed_contacts':len(matched&observed)/len(matched) if matched else None,
                     'observed_contact_count_any_site':len(any_site),'matched_contact_count_any_site':len(matched&any_site),
                     'fraction_matched_on_any_observed_contact':len(matched&any_site)/len(matched) if matched else None,
                     'same_source_chain':(target['pdb'],target['chain'])==(source['pdb'],source['chain'])}
                hits.append(hit);found.append(hit)
        self_hits=[h for h in found if h['same_source_chain']]
        queries.append({**status,'state':'complete','target_chains_returned':len(best),
                        'self_coverage':max((h['coverage'] for h in self_hits),default=0)})
        expected=set(catalog.loc[catalog.query_cluster.eq(site),'query_abp'])-{abp}
        for target in sorted(expected):
            candidates=[h for h in found if h['target_abp']==target and h['observed_contact_count_same_site']>0]
            top=max(candidates,key=lambda h:h['idfscore']) if candidates else {}
            controls.append({**status,'target_abp':target,'returned_match':bool(top),
                **{k:top.get(k) for k in ['target_pdb','target_chain','matched_residues','coverage','idfscore','rmsd_A',
                    'observed_contact_count_same_site','matched_contact_count_same_site','fraction_matched_on_observed_contacts']}})
        print(f'{number}/{len(catalog)} {abp} {site}: {len(best)} target chains',flush=True)
    pd.DataFrame(hits).to_csv(OUT/'hits.csv',index=False)
    pd.DataFrame(queries).to_csv(OUT/'queries.csv',index=False)
    pd.DataFrame(controls).to_csv(OUT/'shared_site_controls.csv',index=False)
    manifest=dict(created_utc=datetime.now(timezone.utc).isoformat(),tool_version=subprocess.check_output([str(binary),'version'],text=True).strip(),
        query_count=len(queries),target_chain_count=len(targets),sources={str(p):fingerprint(p) for p in CATALOG_FILES},
        targets=targets,missing_structures=missing,script_sha256=fingerprint(__file__),
        method='Default FoldDisco local query; highest IDF hit per target chain. Same-site controls use the highest IDF among representative target chains observed contacting this site. Fractions compare exact PDB residue identifiers to observed contacts in the source dataset. No homology or biological negative classification.',
        scope='Current best-resolution site motifs against the finite panel of observed representative ABP chains. Not PDB/AlphaFold discovery; not independent of the source dataset. No curated UniProt ABD boundaries or negative-control specificity estimate.')
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2))
    report=ROOT/'reports/scientific_audit';report.mkdir(exist_ok=True)
    pd.DataFrame(controls).to_csv(report/'folddisco_shared_site_controls.csv',index=False)
    pd.DataFrame(queries).to_csv(report/'folddisco_local_queries.csv',index=False)
    (report/'folddisco_controls_manifest.json').write_text(json.dumps({k:v for k,v in manifest.items() if k!='targets'},indent=2))

if __name__=='__main__':main()
