from display_helpers import plotly_chart
"""Views for the mentor's variant, conservation and interface questions."""
from pathlib import Path
import json
import pandas as pd
import numpy as np
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plot_interaction import position_hover
from scientific_analysis import (load_variants,load_actin_scores,conservation_summary,
                                 interface_evidence,disease_associations,variant_footprint_summary,GENES,
                                 abp_position_asa,cross_gene_footprint_analysis,conflicting_variant_counts)
from footprint_comparison import FILES,load_footprints
from variant_heatmaps import presence_matrix,presence_figure,substitution_count_trace


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
    from rsa_source import source_files
    return scores_cached(signatures(source_files()+list(Path('data/proteocast/actin').glob('*'))+[Path('data/P60709_ref.fasta')]))


def current_footprints():
    return load_footprints(tuple(p.stat().st_mtime_ns for p in FILES))


def download(frame,label,name,key):
    st.download_button(label,frame.to_csv(index=False).encode(),file_name=name,mime='text/csv',key=key)


def render_conservation():
    if not Path('data/proteocast/actin/4.query_ProteoCast.csv').exists():
        st.info('Import the actin ProteoCast source results to enable the complete profile.');return
    raw,scores=current_scores()
    if scores.rsa_status.iloc[0] != "current":
        st.warning(scores.rsa_status.iloc[0] + " Open Documentation → Data management.")
    from residue_passport import build_passport, pp_mtimes
    from residue_explorer import render_residue_conservation
    passport = build_passport(pp_mtimes())
    if passport is not None:
        render_residue_conservation(passport)
    st.subheader('Actin mutational sensitivity and interface use')
    st.caption('Sensitivity = minus the mean of the 20 supplied ProteoCast scores at a position (including the unchanged amino acid). This is a model-derived '
               'mutational sensitivity proxy, not a clinical classification or a direct sequence-identity percentage. '
               'The query sequence and all mutation reference letters are checked against P60709.')
    cutoff=st.slider('Buried ASA threshold for sensitivity footprints (%)',0.0,100.0,0.0,step=1.0,key='cons_asa')
    table,cor,clusters=conservation_summary(scores,current_footprints(),cutoff)
    fig=make_subplots(rows=4,cols=1,shared_xaxes=True,vertical_spacing=.06,
                      subplot_titles=('Mutational sensitivity','RSA · isolated human 8DNH chain B','Distinct ABP source names','Observed actin–actin contacts'))
    for row,(col,label,color) in enumerate([('sensitivity','Mutational sensitivity','#0072B2'),('rsa','RSA','#777777'),('n_abp_names','Distinct ABP source names','#E69F00'),('homo_contact','Actin–actin contact (0/1)','#884EA0')],1):
        fig.add_trace(go.Scatter(x=table.position,y=table[col],mode='lines',name=label,line=dict(color=color),customdata=table.aa,
                                 hovertemplate='%{customdata}%{x} · %{fullData.name}: %{y:.3~g}<extra></extra>'),row=row,col=1)
    fig.update_layout(height=700,showlegend=False);fig.update_xaxes(title_text='P60709 position',row=4,col=1)
    position_hover(fig)
    plotly_chart(fig,use_container_width=True,key='cons_tracks')
    with st.expander('Actin ProteoCast landscape and sequence alignment'):
        pivot=raw.pivot(index='alternate',columns='position',values='Variant_score').reindex(columns=range(1,376))
        plotly_chart(position_hover(go.Figure(go.Heatmap(z=pivot.values,x=pivot.columns,y=pivot.index,colorscale='Blues',colorbar=dict(title='Variant score')))),use_container_width=True,key='cons_landscape')
        msa=Path('data/proteocast/actin/2.aliAF-P60709-F1-msa_v6.fasta')
        if msa.exists():st.download_button('Download actin ProteoCast alignment',msa.read_bytes(),file_name=msa.name,key='cons_msa')
        download(raw,'Download actin variant scores','actin_proteocast_scores.csv','cons_raw')
    st.markdown('**Exploratory correlations and footprint summaries**')
    st.caption('One observation per P60709 position; pairwise missing values are excluded. Spearman correlations '
               'and BH corrections cover the four displayed tests. Residues are structurally dependent: '
               'p-values are exploratory. RSA is calculated locally on isolated human 8DNH chain B in its experimental '
               'filament conformation. It is not an average across structures. No observed contact does not mean no possible interaction.')
    st.dataframe(cor,hide_index=True,width='stretch')
    fig=go.Figure()
    for cat,color in zip(['No observed contact','Actin only','ABP only','Both'],['#999999','#0072B2','#E69F00','#884EA0']):
        sub=table[table.category.eq(cat)]
        fig.add_trace(go.Box(y=sub.sensitivity,name=f'{cat} (n={sub.sensitivity.notna().sum()})',marker_color=color))
    fig.update_layout(yaxis_title='ProteoCast sensitivity',height=330)
    plotly_chart(fig,use_container_width=True,key='cons_categories')
    st.dataframe(clusters,hide_index=True,width='stretch')
    with st.expander('Sensitivity footprints: positive contacts versus selected ASA threshold'):
        _,_,baseline=conservation_summary(scores,current_footprints(),0.0)
        st.caption('The first table uses all positive buried-ASA observations. The second uses the selected '
                   'strict threshold. Both retain the same source sequences and ProteoCast scores.')
        st.markdown('**All positive contacts (0%)**')
        st.dataframe(baseline,hide_index=True,width='stretch')
        st.markdown(f'**ASA strictly above {cutoff:g}%**')
        st.dataframe(clusters,hide_index=True,width='stretch')
        download(baseline,'Download sensitivity without an additional ASA threshold','footprint_conservation_asa0.csv','cons_baseline_csv')
    with st.expander('Sensitivity within a selected ABP binding-site cluster'):
        records=current_footprints()
        names=sorted(records.loc[records.kind.eq('abp'),'group'].unique())
        name=st.selectbox('ABP for cluster sensitivity',names,key='cons_cluster_abp')
        sites=sorted(records.loc[records.kind.eq('abp')&records.group.eq(name),'site'].dropna().unique())
        if sites:
            site=st.selectbox('Binding site for sensitivity',sites,key='cons_cluster_site')
            selected=records[records.kind.eq('abp')&records.group.eq(name)&records.site.eq(site)&records.asa.gt(cutoff)]
            positions=set(selected.position.astype(int));profile=table[table.position.isin(positions)]
            fig=go.Figure(go.Scatter(x=table.position,y=table.sensitivity,mode='lines',name='All positions',line=dict(color='#BBBBBB')))
            fig.add_trace(go.Scatter(x=profile.position,y=profile.sensitivity,mode='markers',name=site,marker_color='#0072B2',customdata=profile.aa,hovertemplate='%{customdata}%{x}: %{y}<extra></extra>'))
            fig.update_layout(xaxis_title='P60709 position',yaxis_title='ProteoCast sensitivity',height=320)
            position_hover(fig, unified=False)
            plotly_chart(fig,use_container_width=True,key='cons_cluster_profile')
            download(profile,'Download selected cluster sensitivity','cluster_conservation.csv','cons_cluster_csv')
    with st.expander('Download sensitivity tables'):
        download(table,'Download aligned residue profiles','actin_conservation_profiles.csv','cons_profiles_csv')
        download(cor,'Download correlations','actin_conservation_correlations.csv','cons_cor_csv')
        download(clusters,'Download footprint sensitivity summaries','footprint_conservation.csv','cons_clusters_csv')


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
    c1,c2,c3=st.columns(3)
    c1.metric('Source records',len(variants));c2.metric('Records eligible for analysis',len(valid));c3.metric('Excluded records',len(variants)-len(valid))
    gene=st.selectbox('Human actin gene',sorted(GENES),key='hv_gene')
    categories=sorted(valid.classif_cat.dropna().unique())
    st.info("ClinVar categories describe evidence about disease association: benign, likely benign, uncertain significance, likely pathogenic and pathogenic. They are not degrees of disease severity. Conflicting annotations and population observations are shown separately; gnomAD presence does not mean benign.")
    category=st.selectbox('Variant annotation to display',categories,index=categories.index('pathogenic'),key='hv_category')
    selected=valid[valid.classif_cat.eq(category)].drop_duplicates(['gene','position','aa_ref','aa_alt'])
    st.caption('For one gene, each alternate-amino-acid/position cell is binary: blue = a recorded substitution; '
               'light gray = no record in this snapshot and selected category. Colour does not grade pathogenicity '
               'or severity. Categories can overlap when source records differ.')
    per=selected[selected.gene.eq(gene)]
    matrix=presence_matrix(per)
    fig=presence_figure(matrix,gene,category)
    position_hover(fig,unified=False)
    plotly_chart(fig,use_container_width=True,key='hv_gene_heatmap')
    counts=selected.groupby(['gene','position']).size().unstack(fill_value=0).reindex(index=sorted(GENES),columns=range(1,376),fill_value=0).fillna(0)
    combined=counts.sum(axis=0).to_frame().T;combined.index=['All genes (sum)']
    counts=pd.concat([counts,combined])
    st.markdown('**Number of distinct substitutions per position**')
    st.caption('Here the gradient represents a count: darker cells contain more different substitutions '
               'in the selected category. The total row sums gene-specific substitutions, so the same change '
               'in two genes contributes twice. This is neither a pathogenicity score nor an allele frequency.')
    fig=make_subplots(rows=3,cols=1,shared_xaxes=True,row_heights=[.6,.2,.2],vertical_spacing=.08,
                      subplot_titles=['Substitution counts by gene','P60709 model sensitivity','Observed ABP source names'])
    fig.add_trace(substitution_count_trace(counts),row=1,col=1)
    fig.add_trace(go.Scatter(x=scores.position,y=scores.sensitivity,name='ProteoCast sensitivity',line=dict(color='#884EA0')),row=2,col=1)
    n=fp[fp.kind.eq('abp')&fp.asa.gt(0)].groupby('position')['group'].nunique().reindex(range(1,376),fill_value=0)
    fig.add_trace(go.Bar(x=n.index,y=n.values,name='ABP names',marker_color='#E69F00'),row=3,col=1)
    fig.update_layout(height=720,margin=dict(l=85,r=135,t=65,b=115),
                      legend=dict(orientation='h',x=0,y=-.16,xanchor='left',yanchor='top'))
    fig.update_xaxes(title_text='P60709 position',row=3,col=1)
    position_hover(fig)
    plotly_chart(fig,use_container_width=True,key='hv_combined')
    detail=valid[valid.gene.eq(gene)].copy()
    detail['P60709_reference']=detail.position.map(scores.set_index('position').aa)
    detail=detail.merge(raw[['position','alternate','Variant_score']],left_on=['position','aa_alt'],right_on=['position','alternate'],how='left')
    detail.loc[detail.aa_ref.ne(detail.P60709_reference),'Variant_score']=np.nan
    detail['Variant_score_context']='P60709 model; not gene-specific or clinical prediction'
    detail['clinvar_url']=detail.clinvar_id.map(lambda x:f'https://www.ncbi.nlm.nih.gov/clinvar/variation/{int(x)}/' if pd.notna(x) else '')
    with st.expander('Source records and downloads'):
        st.dataframe(detail,hide_index=True,width='stretch')
        download(detail,'Download selected gene records and research scores',f'{gene}_variants.csv','hv_gene_csv')
        st.download_button('Download variant source provenance',manifest.read_bytes(),file_name='variant_sources.json',key='hv_manifest')
    with st.expander('Conflicting annotations, mapping audit and cross-gene observations'):
        conflicts=valid[valid.source.eq('clinvar')&valid.classif_cat.eq('conflicting')]
        st.markdown('**ClinVar records explicitly marked conflicting**')
        st.dataframe(conflicts,hide_index=True,width='stretch')
        conflict_counts=conflicting_variant_counts(variants)
        stats=abp_position_asa(fp)
        names=sorted(fp.loc[fp.kind.eq('abp'),'group'].dropna().unique())
        selected_abp=st.selectbox('ABP footprint aligned with conflicting annotations',names,key='hv_conflict_abp') if names else None
        track=scores.set_index('position')[['aa','sensitivity']].reindex(range(1,376))
        selected_stats=stats[stats.ABP.eq(selected_abp)].set_index('position')
        track=track.join(selected_stats[['mean_ASA_percent','observed_interface_chains','observed_interactions']])
        track['positive_contact_observed']=track.mean_ASA_percent.notna()
        fig=make_subplots(rows=4,cols=1,shared_xaxes=True,row_heights=[.50,.10,.20,.20],vertical_spacing=.07,
                          subplot_titles=('Conflicting ClinVar substitutions by gene',
                                          f'Observed footprint: {selected_abp or "no ABP available"}',
                                          'Mean buried ASA at observed positions (%)','P60709 ProteoCast sensitivity'))
        fig.add_trace(go.Heatmap(z=conflict_counts.values,x=conflict_counts.columns,y=conflict_counts.index,
                                colorscale='Purples',zmin=0,zmax=max(1,conflict_counts.values.max()),
                                colorbar=dict(title='Substitutions',len=.35,y=.85),
                                hovertemplate='%{y}<br>P60709 %{x}<br>Conflicting substitutions: %{z}<extra></extra>'),row=1,col=1)
        fig.add_trace(go.Heatmap(z=[track.positive_contact_observed.astype(int).values],x=track.index,y=['Positive ASA observed'],
                                colorscale=[[0,'#F2F2F2'],[1,'#E69F00']],zmin=0,zmax=1,showscale=False,
                                hovertemplate='P60709 %{x}<br>Positive ASA observed: %{z}<extra></extra>'),row=2,col=1)
        fig.add_trace(go.Scatter(x=track.index,y=track.mean_ASA_percent,mode='markers',name='Mean ASA',
                                marker=dict(color='#E69F00',size=5),customdata=track.observed_interface_chains,
                                hovertemplate='P60709 %{x}<br>Mean buried ASA: %{y:.2f}%<br>Interface chains: %{customdata}<extra></extra>'),row=3,col=1)
        fig.add_trace(go.Scatter(x=track.index,y=track.sensitivity,mode='lines',name='ProteoCast sensitivity',
                                line=dict(color='#884EA0')),row=4,col=1)
        fig.update_layout(height=730,showlegend=False);fig.update_xaxes(title_text='P60709 position',row=4,col=1)
        plotly_chart(position_hover(fig),use_container_width=True,key='hv_conflicts_heatmap')
        st.caption('Tracks share aligned P60709 positions. Conflicting means the source annotation; VUS are unchanged. '
                   'The footprint aggregates structural observations across the retained dataset, not structures '
                   'of the selected human gene. A gray footprint cell means no positive ASA observation; it does '
                   'not establish absence of binding. ASA averages distinct interface chains at each position.')
        aligned=track.join(conflict_counts.T.add_prefix('conflicting_substitutions_'))
        aligned.insert(0,'ABP',selected_abp or '')
        download(aligned.reset_index(names='position'),'Download conflicting annotations and aligned ABP tracks',
                 'conflicting_variants_abp_tracks.csv','hv_conflict_tracks_csv')
        cross,pairs,footprints,cross_detail=cross_gene_footprint_analysis(variants,fp)
        st.markdown('**Cross-gene P/LP and population observations**')
        st.caption('Same reference and alternate amino acid at the same aligned position: P/LP in one gene, '
                   'gnomAD observation without a matching substitution in the verified local ClinVar snapshot '
                   'of another gene. Any ClinVar category, including VUS, excludes that population substitution. '
                   'This is not a conflicting ClinVar annotation and is not evidence of benignity or tolerance.')
        if cross.empty:
            st.info('No cross-gene substitutions meet these criteria in the eligible source snapshot.')
        else:
            st.dataframe(pairs,hide_index=True,width='stretch')
            pair_options=list(pairs[['gene_PL','gene_population']].itertuples(index=False,name=None))
            pair=st.selectbox('Gene pair for ABP footprint comparison',pair_options,
                              format_func=lambda p:f'{p[0]} P/LP → {p[1]} population observation',key='hv_cross_pair')
            subset=footprints[footprints.gene_PL.eq(pair[0])&footprints.gene_population.eq(pair[1])]
            st.dataframe(subset,hide_index=True,width='stretch')
            st.caption('Substitutions and positions have separate counts. For each ABP, the two fractions use '
                       'the pair’s distinct positions and the ABP footprint’s distinct positions, respectively. '
                       'ASA takes the maximum per interaction/actin chain/position, averages those chain values '
                       'at each position, then averages contacted positions equally. Repeated substitutions do '
                       'not add weight. Positions without a positive contact have missing ASA, not zero.')
            st.dataframe(cross_detail[cross_detail.gene_PL.eq(pair[0])&cross_detail.gene_population.eq(pair[1])],
                         hide_index=True,width='stretch')
        download(cross,'Download cross-gene observations','cross_gene_observations.csv','hv_cross_csv')
        download(pairs,'Download all gene-pair counts and denominators','cross_gene_pair_counts.csv','hv_cross_pairs_csv')
        download(footprints,'Download all gene-pair and ABP footprint summaries','cross_gene_abp_summaries.csv','hv_cross_abp_csv')
        download(cross_detail,'Download substitutions with ABP positions and ASA evidence','cross_gene_abp_details.csv','hv_cross_details_csv')
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
    st.subheader("Evidence for actin–actin interfaces")
    sites,c70,occurrences,matches=evidence_cached(signatures([p]))
    st.caption('Counts use distinct PDB structures and both sides of actin–actin interactions. '
               'The denominator is all retained PDBs with a homo interaction. Cofilin/coronin/WDR1 '
               '(including actin-interacting protein 1) are matched by source names anywhere in that PDB. '
               'Co-occurrence does not establish direct contact or a distorted filament conformation.')
    fig=go.Figure(go.Bar(x=sites.site,y=sites.PDB_count,marker_color='#0072B2'))
    fig.update_layout(xaxis_title='Binding-site cluster',yaxis_title='Distinct PDB structures',height=330)
    plotly_chart(fig,use_container_width=True,key='interface_evidence_counts')
    st.dataframe(sites,hide_index=True,width='stretch')
    st.markdown('**Interface clusters (C70)**')
    st.dataframe(c70,hide_index=True,width='stretch')
    st.caption('Fisher tests compare occurrence within the context cohort versus other retained PDBs; '
               'BH correction is separate for sites and C70 clusters. PDB entries need not be independent '
               'experiments. Statistical association supports prioritization, not a validated major/minor '
               'biological label. The proposed 6685_1–4 reference group must also be checked structurally.')
    with st.expander('Cohort definitions and evidence downloads'):
        st.markdown('**Source names defining the context cohort**')
        st.dataframe(matches,hide_index=True,width='stretch')
        download(sites,'Download site evidence','actin_interface_site_evidence.csv','ie_sites')
        download(c70,'Download C70 evidence','actin_interface_c70_evidence.csv','ie_c70')
        download(occurrences,'Download PDB–site–C70 correspondence','interface_occurrences.csv','ie_occurrences')
    audit = Path('reports/scientific_audit')
    global_summary = audit/'all_interface_geometry_summary.csv'
    global_manifest = audit/'all_interface_geometry_manifest.json'
    if global_summary.exists() and global_manifest.exists():
        manifest = json.loads(global_manifest.read_text())
        import hashlib
        if manifest.get('source_table_sha256') != hashlib.sha256(p.read_bytes()).hexdigest():
            st.warning('The structural audit predates the current interaction dataset. The following results describe its saved snapshot and must be regenerated before interpreting the updated dataset.')
        st.markdown('**Geometry across retained actin–actin interfaces**')
        st.caption(f"{manifest['measured_pairs']} / {manifest['expected_pairs']} chain pairs measured; "
                   f"{manifest['measured_clusters']} / {manifest['expected_clusters']} interface clusters. "
                   'One actin is aligned to the reference; the neighbor is measured under the same transform. '
                   'The nearest pair geometry is retained separately for 3J8A (with tropomyosin) and 5YU8 (with cofilin). '
                   'Summaries use one median per PDB, so repeated chains do not increase a PDB weight.')
        st.dataframe(pd.read_csv(global_summary), hide_index=True, width='stretch')
        st.caption('Smaller RMSD indicates closer pair geometry. These descriptive measurements do not automatically assign major/minor labels, validate every assembly, or test steric compatibility with an ABP.')
        with st.expander('Geometry measurements, coverage and methods'):
            for name, label in [('all_interface_geometry.csv', 'All pair measurements'),
                                ('all_interface_geometry_coverage.csv', 'Measured and unavailable pairs'),
                                ('all_interface_geometry_summary.csv', 'Summary by cluster and reference')]:
                download(pd.read_csv(audit/name), label, name, f'ie_{name}')
            st.download_button('Reference structures, method and source fingerprints', global_manifest.read_bytes(),
                               file_name=global_manifest.name, key='ie_global_manifest')
    with st.expander('Original three-structure example'):
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
