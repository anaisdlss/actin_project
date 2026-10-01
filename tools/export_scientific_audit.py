"""Export source-identified variant, conservation and interface research tables."""
from pathlib import Path
import os,sys,json,hashlib,datetime
import streamlit
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];os.chdir(ROOT);sys.path.insert(0,str(ROOT/'script'))
from scientific_analysis import (load_variants,load_actin_scores,conservation_summary,interface_evidence,
                                 cross_gene_footprint_analysis)
from footprint_comparison import footprint_records,FILES
from interface_properties import chemistry_records,chemistry_summary,SOURCES as CHEMISTRY_SOURCES
out=ROOT/'reports/scientific_audit';out.mkdir(parents=True,exist_ok=True)
v,m=load_variants(ROOT);raw,scores=load_actin_scores(ROOT)
records=footprint_records(*(pd.read_csv(p,low_memory=False) for p in FILES))
profiles,cor,summary=conservation_summary(scores,records)
a,c,o,ctx=interface_evidence(pd.read_csv(FILES[0],low_memory=False))
cross,cross_pairs,cross_abp,cross_detail=cross_gene_footprint_analysis(v,records)
chem_records=chemistry_records(*(pd.read_csv(p,low_memory=False,dtype={'residue_A_structure':str,'residue_B_structure':str}) for p in CHEMISTRY_SOURCES))
chem_profile,chem_abp,_=chemistry_summary(chem_records)
for name,frame in [('variant_mapping',m),('variant_exclusions',v[~v.analysis_eligible]),('conservation_profiles',profiles),('conservation_correlations',cor),('conservation_footprints',summary),('interface_sites',a),('interface_c70',c),('interface_occurrences',o),('interface_context',ctx)]:
 frame.to_csv(out/f'{name}.csv',index=False)
for name,frame in [('cross_gene_observations',cross),('cross_gene_pair_counts',cross_pairs),
                   ('cross_gene_abp_summaries',cross_abp),('cross_gene_abp_details',cross_detail),
                   ('partner_chemistry_positions',chem_profile),('partner_chemistry_per_abp',chem_abp)]:
 frame.to_csv(out/f'{name}.csv',index=False)
paths=[ROOT/'data/human_variants/clinvar_identity.json']+list(FILES)+list(CHEMISTRY_SOURCES)+list((ROOT/'data/human_variants').rglob('*.csv'))+list((ROOT/'data/human_variants').rglob('*.fasta'))+list((ROOT/'data/proteocast/actin').glob('*'))+[ROOT/'data/proteocast/conservation_vs_asa_per_position.csv']
manifest={'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'variant_records':len(v),'mapped_variant_records':int(v.mapping_valid.sum()),'analysis_eligible_records':int(v.analysis_eligible.sum()),'positions':len(scores),'ASA_threshold':0,'numbering':'P60709',
 'scientific_status':'Exploratory; structure-state labels and clinical effects are not validated by these statistics.',
 'source_sha256':{str(p.relative_to(ROOT) if p.is_absolute() else p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
 'chemistry_weight':'Positive uncensored pair contact area; fractions per interface and position, equal interface means within PDB/ABP, equal PDB then ABP means',
 'code_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'script/scientific_analysis.py',ROOT/'script/footprint_comparison.py',ROOT/'script/interface_properties.py',ROOT/'script/residue_contacts.py',ROOT/'tools/export_scientific_audit.py']}}
(out/'manifest.json').write_text(json.dumps(manifest,indent=2))
print(a.head(10).to_string(index=False))
print(f'Exported {out}')
