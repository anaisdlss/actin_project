"""Reproducible residue-level research summaries; no clinical reclassification."""
from pathlib import Path
import json
import re
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, fisher_exact
import numbering

GENES = {'ACTA1':'P68133','ACTA2':'P62736','ACTB':'P60709','ACTC1':'P68032','ACTG1':'P63261','ACTG2':'P63267'}


def sequence(path):
    return ''.join(x.strip() for x in Path(path).read_text().splitlines() if not x.startswith('>'))


def bh_adjust(values):
    p=np.asarray(values,dtype=float); out=np.full(len(p),np.nan)
    valid=np.flatnonzero(np.isfinite(p)); order=valid[np.argsort(p[valid])]
    if len(order):out[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(order)/np.arange(1,len(order)+1))[::-1])[::-1])
    return out


def isoform_map(src,ref):
    """Global BLOSUM62 alignment; exclude positions with ambiguous optimal mappings."""
    from Bio.Align import PairwiseAligner, substitution_matrices
    al=PairwiseAligner(mode='global',open_gap_score=-10,extend_gap_score=-0.5)
    al.substitution_matrix=substitution_matrices.load('BLOSUM62')
    maps=[]
    for i,a in enumerate(al.align(src,ref)):
        if i>=100:raise ValueError('Too many equally optimal isoform alignments; review mapping.')
        maps.append({s+k+1:r+k+1 for (s,e),(r,_) in zip(*a.aligned) for k in range(e-s)})
    return {p:maps[0][p] for p in maps[0] if all(m.get(p)==maps[0][p] for m in maps)}


def clinvar_identity_status(row, record):
    if record is None:return 'not verified'
    title=record.get('title','')
    gene=re.search(r'\(([A-Za-z0-9_-]+)\):',title)
    if not gene:return 'gene not resolved in source title'
    if gene.group(1)!=row['gene']:return 'different gene in source title'
    mutation=re.search(r'\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)',title)
    if not mutation:return 'protein change not resolved in source title'
    from Bio.SeqUtils import seq1
    ref,pos,alt=mutation.groups()
    if (int(pos),seq1(ref),seq1(alt))!=(int(row['pos']),row['aa_ref'],row['aa_alt']):return 'different protein change in source title'
    return 'verified'


def load_variants(root=Path('.')):
    base=root/'data/human_variants'; ref=sequence(base/'ref/P60709.fasta')
    audit_path=base/'clinvar_identity.json'
    audit=json.loads(audit_path.read_text()).get('records',{}) if audit_path.exists() else {}
    frames=[]; mapping=[]
    for gene,acc in GENES.items():
        seq=sequence(base/f'ref/{acc}.fasta'); mp=isoform_map(seq,ref)
        mapping.extend({'gene':gene,'gene_position':p,'gene_aa':aa,'position':mp.get(p)} for p,aa in enumerate(seq,1))
        for source in ['clinvar','gnomad']:
            df=pd.read_csv(base/f'sources/{source}_{gene}.csv')
            if 'source' in df:
                df['population_dataset' if source=='gnomad' else 'original_source']=df['source']
            df['source']=source
            df['position']=df.pos.map(mp)
            df['reference_matches']=df.apply(lambda r:pd.notna(r.pos) and 1<=r.pos<=len(seq) and seq[int(r.pos)-1]==r.aa_ref,axis=1)
            df['mapping_valid']=df.reference_matches & df.position.notna()
            if source=='clinvar':
                identities=[audit.get(str(int(uid))) for uid in df.clinvar_id]
                df['source_identity_status']=[clinvar_identity_status(row,record) for (_,row),record in zip(df.iterrows(),identities)]
                df['source_title']=[(r or {}).get('title','') for r in identities]
                df['identity_checked_utc']=[(r or {}).get('checked_utc','') for r in identities]
                df['current_ClinVar_classification']=[(r or {}).get('current_classification','') for r in identities]
                df['mapping_valid'] &= df.source_identity_status.eq('verified')
            if source=='gnomad':
                df['source_identity_status']='local transcript-scoped source; reference sequence checked'
                df['classif_cat']='population_unclassified';df['classification']='Not classified by this population source'
                df['population_observed']=pd.to_numeric(df.ac,errors='coerce').gt(0) & pd.to_numeric(df.an,errors='coerce').gt(0)
            df['analysis_eligible']=df.mapping_valid & (df.population_observed if source=='gnomad' else True)
            frames.append(df)
    return pd.concat(frames,ignore_index=True),pd.DataFrame(mapping)


def load_actin_scores(root=Path('.')):
    base=root/'data/proteocast/actin'; ref=sequence(root/'data/P60709_ref.fasta')
    if sequence(base/'1.query.fasta')!=ref:raise ValueError('ProteoCast query differs from P60709; explicit mapping required.')
    raw=pd.read_csv(base/'4.query_ProteoCast.csv')
    parsed=raw.Mutation.str.extract(r'^([A-Z])(\d+)([A-Z])$'); pos=pd.to_numeric(parsed[1],errors='coerce')
    ok=pos.between(1,len(ref)) & pos.eq(raw.Residue) & parsed[0].eq(pos.map(lambda p:ref[int(p)-1] if pd.notna(p) and 1<=p<=len(ref) else None))
    if not ok.all():raise ValueError('ProteoCast mutation identities do not match query.')
    raw['alternate']=parsed[2];raw['position']=pos.astype(int)
    if raw.duplicated(['position','alternate']).any():raise ValueError('Duplicate ProteoCast mutation scores.')
    raw['Variant_score']=pd.to_numeric(raw.Variant_score,errors='coerce')
    if not np.isfinite(raw.Variant_score).all():raise ValueError('Nonfinite ProteoCast scores.')
    expected=set('ACDEFGHIKLMNPQRSTVWY')
    coverage=raw.groupby('position').alternate.agg(set)
    if set(coverage.index)!=set(range(1,376)) or not coverage.map(lambda x:x==expected).all():
        raise ValueError('ProteoCast requires all 20 supplied scores at each of the 375 positions.')
    scores=raw.groupby('position').Variant_score.mean().reindex(range(1,376)).to_frame('mean_variant_score')
    scores['sensitivity']=-scores.mean_variant_score
    scores['aa']=list(ref);scores.index.name='position'
    old=root/'data/proteocast/conservation_vs_asa_per_position.csv'
    scores['rsa']=np.nan
    if old.exists():
        rs=pd.read_csv(old)
        rs['position']=numbering.uniprot_series(rs.canon)
        # RSA uses structural MAFFT mapping, not the ProteoCast sequence index.
        rs=rs.dropna(subset=['position']).set_index('position')
        scores['rsa']=pd.to_numeric(rs.rsa,errors='coerce').reindex(scores.index)
        scores.loc[~scores.rsa.between(0,1),'rsa']=np.nan
    return raw,scores.reset_index()


def conservation_summary(scores,records,cutoff=0):
    records=records[records.asa.gt(cutoff)]
    out=scores.copy().set_index('position')
    homo=set(records.loc[records.kind.eq('homo'),'position']);abp=records[records.kind.eq('abp')]
    out['homo_contact']=out.index.isin(homo).astype(int)
    out['n_abp_names']=abp.groupby('position')['group'].nunique().reindex(out.index,fill_value=0)
    out['abp_contact']=out.n_abp_names.gt(0).astype(int)
    out['category']=np.select([out.homo_contact.eq(1)&out.abp_contact.eq(1),out.homo_contact.eq(1),out.abp_contact.eq(1)],['Both','Actin only','ABP only'],default='No observed contact')
    correlations=[]
    for metric in ['rsa','n_abp_names','homo_contact','abp_contact']:
        sub=out[['sensitivity',metric]].dropna()
        rho,p=spearmanr(sub.sensitivity,sub[metric]) if len(sub)>2 and sub[metric].nunique()>1 and sub.sensitivity.nunique()>1 else (np.nan,np.nan)
        correlations.append(dict(metric=metric,n_positions=len(sub),spearman_rho=rho,p_value=p))
    cor=pd.DataFrame(correlations);cor['q_BH']=bh_adjust(cor.p_value)
    summary=[]
    for (kind,name),g in records.groupby(['kind','group']):
        pos=set(g.position.astype(int));v=out.loc[out.index.isin(pos),'sensitivity'].dropna()
        summary.append(dict(kind=kind,footprint=name,n_positions=len(pos),n_scored=len(v),mean_sensitivity=v.mean(),median_sensitivity=v.median()))
    for (name,site),g in records[records.kind.eq('abp')].groupby(['group','site']):
        pos=set(g.position.astype(int));v=out.loc[out.index.isin(pos),'sensitivity'].dropna()
        summary.append(dict(kind='abp_site',footprint=f'{name} / {site}',n_positions=len(pos),n_scored=len(v),mean_sensitivity=v.mean(),median_sensitivity=v.median()))
    return out.reset_index(),cor,pd.DataFrame(summary)


def interface_evidence(data):
    """One observation per distinct PDB, site and C70; both actin sides included."""
    homo=data[data.s1_actine.astype(str).str.lower().eq('true') & data.s2_actine.astype(str).str.lower().eq('true')].copy()
    parts=[]
    for side in [1,2]:
        p=homo[['pdb_id','cluster_data_70',f's{side}_binding_site_cluster_data_70']].copy()
        p.columns=['pdb_id','c70','site'];parts.append(p)
    occurrences=pd.concat(parts).dropna().drop_duplicates()
    universe=set(homo.pdb_id);partners=[]
    for side in [1,2]:
        flag=data[f's{side}_actine'].astype(str).str.lower().eq('false')
        partners.append(data.loc[flag,['pdb_id',f'subunit_{side}_title']].rename(columns={f'subunit_{side}_title':'partner'}))
    partners=pd.concat(partners).drop_duplicates()
    pattern=r'cofilin|coronin|WDR1|WD.repeat.containing.protein.1|actin.interacting.protein.1'
    matches=partners[partners.partner.str.contains(pattern,case=False,na=False,regex=True)]
    context=set(matches.pdb_id)&universe
    tables=[]
    for level in ['site','c70']:
        rows=[]
        for name,g in occurrences.groupby(level):
            present=set(g.pdb_id);a=len(present&context);b=len(context-present);c=len(present-context);d=len(universe-context-present)
            odds,p=fisher_exact([[a,b],[c,d]]) if context and universe-context else (np.nan,np.nan)
            rows.append({level:name,'PDB_count':len(present),'PDB_denominator':len(universe),'PDB_fraction':len(present)/len(universe),
                         'context_present':a,'context_total':len(context),'other_present':c,'other_total':len(universe-context),
                         'odds_ratio':odds,'p_value':p,'PDB_ids':'; '.join(sorted(present))})
        frame=pd.DataFrame(rows)
        if len(frame):frame['q_BH']=bh_adjust(frame.p_value);frame=frame.sort_values(['PDB_count',level],ascending=[False,True])
        tables.append(frame)
    return *tables,occurrences,matches


def disease_associations(variants,records,gene,disease):
    """Condition mentions among one gene's aggregate P/LP variant records.

    Pooled condition labels do not establish condition-specific assertions.
    """
    v=variants[variants.mapping_valid & variants.gene.eq(gene) & variants.source.eq('clinvar') & variants.classif_cat.isin(['pathogenic','likely_pathogenic'])]
    background=set(v.position.astype(int));cases=set(v.loc[v.maladies.fillna('').str.split(' | ',regex=False).map(lambda x:disease in x),'position'].astype(int))
    rows=[]
    for name,g in records[records.kind.eq('abp') & records.asa.gt(0)].groupby('group'):
        fp=set(g.position.astype(int));a=len(cases&fp);b=len(cases-fp);c=len((background-cases)&fp);d=len(background-cases-fp)
        odds,p=fisher_exact([[a,b],[c,d]]) if cases and background-cases else (np.nan,np.nan)
        rows.append(dict(ABP=name,disease_positions_in_footprint=a,disease_positions_total=len(cases),other_PL_positions_in_footprint=c,other_PL_positions_total=len(background-cases),footprint_positions=len(fp),PL_positions_in_footprint=len(background&fp),PL_fraction_of_footprint=len(background&fp)/len(fp) if fp else np.nan,odds_ratio=odds,p_value=p))
    out=pd.DataFrame(rows)
    if len(out):out['q_BH']=bh_adjust(out.p_value)
    return out


def canonical_conservation(root=Path('.')):
    """Refresh score coordinates while preserving structural RSA in its own mapping."""
    raw,scores=load_actin_scores(root)
    out=scores.rename(columns={'position':'Residue','sensitivity':'conservation','mean_variant_score':'mean_vs'}).copy()
    out['canon']=out.Residue.map(numbering.to_canon)
    classes=raw.groupby('position').Residue_class.agg(lambda s:s.iloc[0] if s.nunique()==1 else 'mixed')
    out['residue_class']=out.Residue.map(classes)
    out['frac_impactful']=out.Residue.map(raw.assign(impact=raw.Variant_class.eq('impactful')).groupby('position').impact.mean())
    from footprint_comparison import footprint_records,FILES
    if all((root/p).exists() for p in FILES):
        records=footprint_records(*(pd.read_csv(root/p,low_memory=False) for p in FILES))
        records=records[records.asa.gt(0)]
        out['at_homo']=out.Residue.isin(records.loc[records.kind.eq('homo'),'position'])
        out['at_hetero']=out.Residue.isin(records.loc[records.kind.eq('abp'),'position'])
        out['at_interface']=out.at_homo | out.at_hetero
    else:
        for col in ['at_homo','at_hetero','at_interface']:out[col]=pd.NA
    return out


def variant_footprint_summary(variants,records,gene):
    v=variants[variants.analysis_eligible & variants.gene.eq(gene)]
    categories={'PL':{'pathogenic','likely_pathogenic'},'BL':{'benign','likely_benign'},'VUS':{'uncertain'},'conflicting':{'conflicting'}}
    sets={name:set(v.loc[v.source.eq('clinvar')&v.classif_cat.isin(labels),'position'].astype(int)) for name,labels in categories.items()}
    cv_keys=set(map(tuple,v.loc[v.source.eq('clinvar'),['position','aa_ref','aa_alt']].values))
    pop=v[v.source.eq('gnomad')]
    sets['population_without_ClinVar']=set(int(r.position) for r in pop.itertuples() if (r.position,r.aa_ref,r.aa_alt) not in cv_keys)
    rows=[]
    for name,g in records[records.kind.eq('abp')&records.asa.gt(0)].groupby('group'):
        footprint=set(g.position.astype(int));row={'gene':gene,'ABP':name,'footprint_positions':len(footprint)}
        for label,positions in sets.items():
            row[label+'_positions']=len(positions&footprint);row[label+'_fraction']=len(positions&footprint)/len(footprint)
        # Weight each interface chain once, then each observed position equally.
        per=g.groupby(['position','interaction_id','chain']).asa.max().groupby('position').mean()
        row['mean_ASA_at_PL_positions']=per.reindex(sorted(sets['PL']&footprint)).mean()
        rows.append(row)
    return pd.DataFrame(rows)
