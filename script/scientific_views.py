"""Views for the mentor's variant, conservation and interface questions."""
from pathlib import Path
import json
import pandas as pd
import numpy as np
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scientific_analysis import (load_variants,load_actin_scores,conservation_summary,
                                 interface_evidence,disease_associations,variant_footprint_summary,GENES)
from footprint_comparison import FILES,load_footprints


def signatures(paths):
    # Streamlit hashes this wrapper, not the imported scientific calculations.
    dependencies=[Path(__file__).with_name('scientific_analysis.py'),
                  Path(__file__).with_name('footprint_comparison.py'),
                  Path(__file__).with_name('numbering.py')]
    return tuple((str(p),p.stat().st_mtime_ns,p.stat().st_size) for p in [*paths,*dependencies] if p.exists())


@st.cache_data(show_spinner=False)
def variants_cached(stamp):return load_variants()


@st.cache_data(show_spinner=False)
def scores_cached(stamp):return load_actin_scores()


def current_scores():
    return scores_cached(signatures(list(Path('data/proteocast/actin').glob('*'))+[Path('data/P60709_ref.fasta'),Path('data/proteocast/conservation_vs_asa_per_position.csv')]))


def current_footprints():
    return load_footprints(tuple(p.stat().st_mtime_ns for p in FILES))


def download(frame,label,name,key):
    st.download_button(label,frame.to_csv(index=False).encode(),file_name=name,mime='text/csv',key=key)


def render_conservation():
    if not Path('data/proteocast/actin/4.query_ProteoCast.csv').exists():
        st.info('Import the actin ProteoCast source results to enable the complete profile.');return
    raw,scores=current_scores()
    st.subheader('Actin mutational sensitivity and interface use')
    st.caption('Sensitivity = minus the mean of the 20 supplied ProteoCast scores at a position (including the unchanged amino acid). This is a model-derived '
               'mutational sensitivity proxy, not a clinical classification or a direct sequence-identity percentage. '
               'The query sequence and all mutation reference letters are checked against P60709.')
    cutoff=st.slider('Buried ASA threshold for conservation footprints (%)',0.0,100.0,0.0,step=1.0,key='cons_asa')
    table,cor,clusters=conservation_summary(scores,current_footprints(),cutoff)
    fig=make_subplots(rows=4,cols=1,shared_xaxes=True,vertical_spacing=.06,
                      subplot_titles=('Mutational sensitivity','RSA in existing structural source','Distinct ABP source names','Observed actin–actin contacts'))
    for row,(col,color) in enumerate([('sensitivity','#0072B2'),('rsa','#777777'),('n_abp_names','#E69F00'),('homo_contact','#884EA0')],1):
        fig.add_trace(go.Scatter(x=table.position,y=table[col],mode='lines',name=col,line=dict(color=color),customdata=table.aa,
                                 hovertemplate='%{customdata}%{x}: %{y}<extra></extra>'),row=row,col=1)
    fig.update_layout(height=700,showlegend=False);fig.update_xaxes(title_text='P60709 position',row=4,col=1)
    st.plotly_chart(fig,use_container_width=True,key='cons_tracks')
    with st.expander('Actin ProteoCast landscape and sequence alignment'):
        pivot=raw.pivot(index='alternate',columns='position',values='Variant_score').reindex(columns=range(1,376))
        st.plotly_chart(go.Figure(go.Heatmap(z=pivot.values,x=pivot.columns,y=pivot.index,colorscale='Blues',colorbar=dict(title='Variant score'))),use_container_width=True,key='cons_landscape')
        msa=Path('data/proteocast/actin/2.aliAF-P60709-F1-msa_v6.fasta')
        if msa.exists():st.download_button('Download actin ProteoCast alignment',msa.read_bytes(),file_name=msa.name,key='cons_msa')
        download(raw,'Download actin variant scores','actin_proteocast_scores.csv','cons_raw')
    st.markdown('**Exploratory correlations and footprint summaries**')
    st.caption('One observation per P60709 position; pairwise missing values are excluded. Spearman correlations '
               'and BH corrections cover the four displayed tests. Residues are structurally dependent: '
               'p-values are exploratory. RSA comes from the existing structural table; its filament/monomer '
               'provenance has not been established. No observed contact does not mean no possible interaction.')
    st.dataframe(cor,hide_index=True,width='stretch')
    fig=go.Figure()
    for cat,color in zip(['No observed contact','Actin only','ABP only','Both'],['#999999','#0072B2','#E69F00','#884EA0']):
        sub=table[table.category.eq(cat)]
        fig.add_trace(go.Box(y=sub.sensitivity,name=f'{cat} (n={sub.sensitivity.notna().sum()})',marker_color=color))
    fig.update_layout(yaxis_title='ProteoCast sensitivity',height=330)
    st.plotly_chart(fig,use_container_width=True,key='cons_categories')
    st.dataframe(clusters,hide_index=True,width='stretch')
    with st.expander('Conservation within a selected ABP binding-site cluster'):
        records=current_footprints()
        names=sorted(records.loc[records.kind.eq('abp'),'group'].unique())
        name=st.selectbox('ABP for cluster conservation',names,key='cons_cluster_abp')
        sites=sorted(records.loc[records.kind.eq('abp')&records.group.eq(name),'site'].dropna().unique())
        if sites:
            site=st.selectbox('Binding site for conservation',sites,key='cons_cluster_site')
            selected=records[records.kind.eq('abp')&records.group.eq(name)&records.site.eq(site)&records.asa.gt(cutoff)]
            positions=set(selected.position.astype(int));profile=table[table.position.isin(positions)]
            fig=go.Figure(go.Scatter(x=table.position,y=table.sensitivity,mode='lines',name='All positions',line=dict(color='#BBBBBB')))
            fig.add_trace(go.Scatter(x=profile.position,y=profile.sensitivity,mode='markers',name=site,marker_color='#0072B2',customdata=profile.aa,hovertemplate='%{customdata}%{x}: %{y}<extra></extra>'))
            fig.update_layout(xaxis_title='P60709 position',yaxis_title='ProteoCast sensitivity',height=320)
            st.plotly_chart(fig,use_container_width=True,key='cons_cluster_profile')
            download(profile,'Download selected cluster conservation','cluster_conservation.csv','cons_cluster_csv')
    download(table,'Download aligned residue profiles','actin_conservation_profiles.csv','cons_profiles_csv')
    download(cor,'Download correlations','actin_conservation_correlations.csv','cons_cor_csv')
    download(clusters,'Download footprint conservation summaries','footprint_conservation.csv','cons_clusters_csv')


def render_variants():
    base=Path('data/human_variants')
    if not (base/'manifest.json').exists():
        st.info('Human variant source snapshot is not installed.');return
    variants,mapping=variants_cached(signatures(list(base.rglob('*.csv'))+list(base.rglob('*.fasta'))+[base/'clinvar_identity.json']))
    raw,scores=current_scores()
    fp=current_footprints(); valid=variants[variants.analysis_eligible].copy()
    st.caption('Local ClinVar and gnomAD source snapshot imported from the existing actin-variants project. '
               'Original download dates are unknown; this is not a live database update. Source records and '
               'identifiers are retained. ClinVar genes and protein changes are checked against official record titles. '
               'Unverified identities, sequence mismatches, unmapped positions and population records with zero observed alleles are excluded from analyses '
               'and retained in the audit. Classifications and conditions use the original snapshot; current '
               'classifications are retained separately for comparison. Presence in gnomAD is not a benign label.')
    manifest=base/'manifest.json'
    st.download_button('Download variant source provenance',manifest.read_bytes(),file_name='variant_sources.json',key='hv_manifest')
    c1,c2,c3=st.columns(3)
    c1.metric('Source records',len(variants));c2.metric('Records eligible for analysis',len(valid));c3.metric('Excluded records',len(variants)-len(valid))
    gene=st.selectbox('Human actin gene',sorted(GENES),key='hv_gene')
    categories=sorted(valid.classif_cat.dropna().unique())
    category=st.selectbox('Variant annotation to display',categories,index=categories.index('pathogenic'),key='hv_category')
    selected=valid[valid.classif_cat.eq(category)].drop_duplicates(['gene','position','aa_ref','aa_alt'])
    st.caption('Heatmaps count distinct protein substitutions in the selected category. Categories can overlap '
               'when different source records describe the same substitution. White = no record in this snapshot/category.')
    per=selected[selected.gene.eq(gene)]
    matrix=pd.crosstab(per.aa_alt,per.position).reindex(index=list('ACDEFGHIKLMNPQRSTVWY'),columns=range(1,376),fill_value=0)
    fig=go.Figure(go.Heatmap(z=matrix.values,x=matrix.columns,y=matrix.index,colorscale=[[0,'#ffffff'],[1,'#0072B2']],zmin=0,zmax=max(1,matrix.values.max()),hovertemplate='P60709 %{x} → %{y}<br>Substitutions: %{z}<extra></extra>'))
    fig.update_layout(title=f'{gene}: {category}',xaxis_title='Aligned P60709 position',yaxis_title='Alternate amino acid',height=420)
    st.plotly_chart(fig,use_container_width=True,key='hv_gene_heatmap')
    counts=selected.groupby(['gene','position']).size().unstack(fill_value=0).reindex(index=sorted(GENES),columns=range(1,376),fill_value=0).fillna(0)
    combined=counts.sum(axis=0).to_frame().T;combined.index=['All genes (gene-specific substitutions)']
    counts=pd.concat([counts,combined])
    fig=make_subplots(rows=3,cols=1,shared_xaxes=True,row_heights=[.6,.2,.2],vertical_spacing=.08)
    fig.add_trace(go.Heatmap(z=counts.values,x=counts.columns,y=counts.index,colorscale='Blues',colorbar=dict(title='Count',len=.45,y=.8)),row=1,col=1)
    fig.add_trace(go.Scatter(x=scores.position,y=scores.sensitivity,name='ProteoCast sensitivity',line=dict(color='#884EA0')),row=2,col=1)
    n=fp[fp.kind.eq('abp')&fp.asa.gt(0)].groupby('position')['group'].nunique().reindex(range(1,376),fill_value=0)
    fig.add_trace(go.Bar(x=n.index,y=n.values,name='ABP names',marker_color='#E69F00'),row=3,col=1)
    fig.update_layout(height=630);fig.update_xaxes(title_text='P60709 position',row=3,col=1)
    st.plotly_chart(fig,use_container_width=True,key='hv_combined')
    detail=valid[valid.gene.eq(gene)].copy()
    detail['P60709_reference']=detail.position.map(scores.set_index('position').aa)
    detail=detail.merge(raw[['position','alternate','Variant_score']],left_on=['position','aa_alt'],right_on=['position','alternate'],how='left')
    detail.loc[detail.aa_ref.ne(detail.P60709_reference),'Variant_score']=np.nan
    detail['Variant_score_context']='P60709 model; not gene-specific or clinical prediction'
    detail['clinvar_url']=detail.clinvar_id.map(lambda x:f'https://www.ncbi.nlm.nih.gov/clinvar/variation/{int(x)}/' if pd.notna(x) else '')
    st.dataframe(detail,hide_index=True,width='stretch')
    download(detail,'Download selected gene records and research scores',f'{gene}_variants.csv','hv_gene_csv')
    with st.expander('Conflicting annotations, mapping audit and cross-gene observations'):
        conflicts=valid[valid.classif_cat.eq('conflicting')]
        st.markdown('**ClinVar records explicitly marked conflicting**')
        st.dataframe(conflicts,hide_index=True,width='stretch')
        conflict_counts=conflicts.drop_duplicates(['gene','position','aa_ref','aa_alt']).groupby(['gene','position']).size().unstack(fill_value=0).reindex(index=sorted(GENES),columns=range(1,376),fill_value=0).fillna(0)
        st.plotly_chart(go.Figure(go.Heatmap(z=conflict_counts.values,x=conflict_counts.columns,y=conflict_counts.index,colorscale='Purples')),use_container_width=True,key='hv_conflicts_heatmap')
        keys=['position','aa_ref','aa_alt']
        patho=valid[valid.source.eq('clinvar')&valid.classif_cat.isin(['pathogenic','likely_pathogenic'])][['gene',*keys]].drop_duplicates()
        population=valid[valid.source.eq('gnomad')][['gene',*keys]].drop_duplicates()
        # Exclude substitutions having ANY ClinVar record in the population gene.
        cv=valid[valid.source.eq('clinvar')][['gene',*keys]].drop_duplicates()
        population=population.merge(cv.assign(in_clinvar=True),on=['gene',*keys],how='left')
        population=population[population.in_clinvar.isna()].drop(columns='in_clinvar')
        cross=patho.merge(population,on=keys,suffixes=('_PL','_population'))
        cross=cross[cross.gene_PL.ne(cross.gene_population)]
        st.caption('Same reference and alternate amino acid at the same aligned position: P/LP in one gene, '
                   'gnomAD observation without a ClinVar record in another. This is a cross-gene observation, '
                   'not a conflict within ClinVar and not evidence of tolerance.')
        abp=fp[fp.kind.eq('abp')&fp.asa.gt(0)].groupby('position')['group'].agg(lambda s:'; '.join(sorted(set(s))))
        cross['ABP_source_names']=cross.position.map(abp).fillna('')
        st.dataframe(cross,hide_index=True,width='stretch')
        download(cross,'Download cross-gene observations','cross_gene_observations.csv','hv_cross_csv')
        download(conflicts,'Download conflicting source records','conflicting_variants.csv','hv_conflicts_csv')
        download(mapping,'Download isoform-to-P60709 mapping','human_actin_mapping.csv','hv_mapping_csv')
        st.dataframe(variants[~variants.analysis_eligible],hide_index=True,width='stretch')
        download(variants,'Download all original annotations with mapping checks','human_variant_audit.csv','hv_audit_csv')
    with st.expander('Variant categories within each ABP footprint'):
        proportions=variant_footprint_summary(variants,fp,gene)
        st.caption('Counts are unique aligned positions in each footprint; fractions divide by that footprint’s '
                   'number of positions. PL = pathogenic/likely pathogenic; BL = benign/likely benign. '
                   'Categories can overlap at a position because substitutions or records differ. '
                   'Population without ClinVar means that the same substitution is absent from the local '
                   'ClinVar snapshot, not benign. ASA averages first across interface chains at each position, '
                   'then across the P/LP positions; these are positive observed contacts only.')
        st.dataframe(proportions,hide_index=True,width='stretch')
        download(proportions,'Download variant proportions and ASA by footprint',f'{gene}_footprint_variants.csv','hv_proportions_csv')
    with st.expander('Disease and ABP footprint associations (exploratory)'):
        pl=valid[valid.gene.eq(gene)&valid.source.eq('clinvar')&valid.classif_cat.isin(['pathogenic','likely_pathogenic'])]
        diseases=sorted({d for text in pl.maladies.dropna() for d in text.split(' | ') if d})
        if diseases:
            disease=st.selectbox('Condition in source annotations',diseases,key='hv_disease')
            assoc=disease_associations(variants,fp,gene,disease)
            st.caption('Two-sided Fisher tests compare positions of aggregate P/LP variants mentioning the selected condition with '
                       'the other P/LP positions of the SAME gene, inside/outside each ABP footprint. One count '
                       'per position, regardless of submissions. Pooled labels do not establish condition-specific '
                       'classifications. BH correction across displayed ABPs. '
                       'Conditions use unharmonized source labels; reporting bias and overlapping footprints '
                       'prevent a causal or clinical interpretation. VUS are never reclassified here.')
            st.dataframe(assoc,hide_index=True,width='stretch')
            download(assoc,'Download disease–footprint analysis','disease_footprint_associations.csv','hv_disease_csv')
    st.markdown('[ClinVar classification definitions](https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/) · '
                '[ClinVar review status](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/)')


@st.cache_data(show_spinner=False)
def evidence_cached(stamp):
    return interface_evidence(pd.read_csv('data/filtered/filtered_all_data.csv',low_memory=False))


def render_interface_evidence():
    p=Path('data/filtered/filtered_all_data.csv')
    if not p.exists():return
    with st.expander('Evidence for prevalent and context-associated actin–actin interfaces'):
        sites,c70,occurrences,matches=evidence_cached(signatures([p]))
        st.caption('Counts use distinct PDB structures and both sides of actin–actin interactions. '
                   'The denominator is all retained PDBs with a homo interaction. Cofilin/coronin/WDR1 '
                   '(including actin-interacting protein 1) are matched by source names anywhere in that PDB. '
                   'Co-occurrence does not establish direct contact or a distorted filament conformation.')
        fig=go.Figure(go.Bar(x=sites.site,y=sites.PDB_count,marker_color='#0072B2'))
        fig.update_layout(xaxis_title='Binding-site cluster',yaxis_title='Distinct PDB structures',height=330)
        st.plotly_chart(fig,use_container_width=True,key='interface_evidence_counts')
        st.dataframe(sites,hide_index=True,width='stretch')
        st.markdown('**Interface clusters (C70)**')
        st.dataframe(c70,hide_index=True,width='stretch')
        st.caption('Fisher tests compare occurrence within the context cohort versus other retained PDBs; '
                   'BH correction is separate for sites and C70 clusters. PDB entries need not be independent '
                   'experiments. Statistical association supports prioritization, not a validated major/minor '
                   'biological label. The proposed 6685_1–4 reference group must also be checked structurally.')
        download(sites,'Download site evidence','actin_interface_site_evidence.csv','ie_sites')
        download(c70,'Download C70 evidence','actin_interface_c70_evidence.csv','ie_c70')
        download(occurrences,'Download PDB–site–C70 correspondence','interface_occurrences.csv','ie_occurrences')
        st.markdown('**Source names defining the context cohort**')
        st.dataframe(matches,hide_index=True,width='stretch')
        geometry=Path('reports/scientific_audit/representative_geometry_summary.csv')
        if geometry.exists():
            st.markdown('**Representative structural check**')
            st.caption('One actin subunit is superimposed by its mapped C-alpha atoms; the neighbor is measured '
                       'under the same transform, without refitting. 3J8A is the F-actin/tropomyosin reference; '
                       '5YU8 and 6VAO are cofilin-decorated references. These examples support distinct pair '
                       'geometries; they do not validate every cluster or prove a clash.')
            st.dataframe(pd.read_csv(geometry),hide_index=True,width='stretch')
            st.markdown('[3J8A](https://www.rcsb.org/structure/3J8A) · '
                        '[5YU8](https://www.rcsb.org/structure/5YU8) · '
                        '[6VAO](https://www.rcsb.org/structure/6VAO)')
            download(pd.read_csv(geometry),'Download representative geometry check','representative_geometry.csv','ie_geometry')
