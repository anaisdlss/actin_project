"""Export source-identified variant, conservation and interface research tables."""
from pathlib import Path
import os,sys,json,hashlib,datetime
import streamlit
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];os.chdir(ROOT);sys.path.insert(0,str(ROOT/'script'))
from scientific_analysis import load_variants,load_actin_scores,conservation_summary,interface_evidence
from footprint_comparison import footprint_records,FILES
out=ROOT/'reports/scientific_audit';out.mkdir(parents=True,exist_ok=True)
v,m=load_variants(ROOT);raw,scores=load_actin_scores(ROOT)
records=footprint_records(*(pd.read_csv(p,low_memory=False) for p in FILES))
profiles,cor,summary=conservation_summary(scores,records)
a,c,o,ctx=interface_evidence(pd.read_csv(FILES[0],low_memory=False))
for name,frame in [('variant_mapping',m),('variant_exclusions',v[~v.analysis_eligible]),('conservation_profiles',profiles),('conservation_correlations',cor),('conservation_footprints',summary),('interface_sites',a),('interface_c70',c),('interface_occurrences',o),('interface_context',ctx)]:
 frame.to_csv(out/f'{name}.csv',index=False)
paths=[ROOT/'data/human_variants/clinvar_identity.json']+list(FILES)+list((ROOT/'data/human_variants').rglob('*.csv'))+list((ROOT/'data/human_variants').rglob('*.fasta'))+list((ROOT/'data/proteocast/actin').glob('*'))+[ROOT/'data/proteocast/conservation_vs_asa_per_position.csv']
manifest={'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'variant_records':len(v),'mapped_variant_records':int(v.mapping_valid.sum()),'analysis_eligible_records':int(v.analysis_eligible.sum()),'positions':len(scores),'ASA_threshold':0,'numbering':'P60709',
 'scientific_status':'Exploratory; structure-state labels and clinical effects are not validated by these statistics.',
 'source_sha256':{str(p.relative_to(ROOT) if p.is_absolute() else p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
 'code_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'script/scientific_analysis.py',ROOT/'script/footprint_comparison.py']}}
(out/'manifest.json').write_text(json.dumps(manifest,indent=2))
print(a.head(10).to_string(index=False))
print(f'Exported {out}')
