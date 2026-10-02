from display_helpers import plotly_chart
"""Exploratory overlap of observed footprints; no steric or functional verdict."""
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
import streamlit as st
import numbering
from plot_interaction import position_hover

FILES = [Path('data/filtered/filtered_all_data.csv'),
         Path('data/filtered/details/1.interactions.csv'),
         Path('data/filtered/details/3.interface_residues.csv')]


def footprint_records(all_data, interactions, residues):
    m = all_data.merge(interactions[['interaction_id','chain_A_id','chain_B_id']],
                       left_on=['subunit_1','subunit_2'], right_on=['chain_A_id','chain_B_id'])
    parts=[]
    for side,other in ((1,2),(2,1)):
        actin = m[f's{side}_actine'].astype(str).str.lower().eq('true')
        partner_flag = m[f's{other}_actine'].astype(str).str.lower()
        for kind,flag in [('homo','true'),('abp','false')]:
            rows=m[actin & partner_flag.eq(flag)].copy()
            rows['site'] = rows[f's{side}_binding_site_cluster_data_70']
            rows['group'] = rows['site'] if kind=='homo' else rows[f'subunit_{other}_title']
            rows=rows.rename(columns={f'subunit_{side}':'chain'})
            joined=rows[['interaction_id','chain','group','site']].merge(residues,on=['interaction_id','chain'])
            joined['position'] = numbering.uniprot_series(joined['residue_number_canon_mafft'])
            joined['asa'] = pd.to_numeric(joined['buried_ASA_percent'].astype(str).str.replace('%','',regex=False),errors='coerce')
            joined['kind']=kind
            parts.append(joined[['kind','group','site','interaction_id','chain','position','asa']])
    return pd.concat(parts,ignore_index=True).dropna(subset=['group','position','asa']).drop_duplicates()


def jaccard(a,b):
    union=a|b
    return len(a&b)/len(union) if union else np.nan


@st.cache_data(show_spinner=False)
def load_footprints(mtimes):
    return footprint_records(*(pd.read_csv(p,low_memory=False) for p in FILES))


@st.cache_data(show_spinner=False)
def source_signatures(mtimes):
    return {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in FILES}


def render_footprint_comparison():
    if not all(p.exists() for p in FILES):
        return
    st.caption('Select reference and comparison sites to explore candidate major/minor '
               'interfaces. The initial reference sites are 6685_1–4, proposed in the project '
               'brief; no biological major/minor label is assigned automatically. Both sides '
               'of homo interactions are included. Each footprint is the union of observed '
               'P60709 positions, not an average of structures.')
    records=load_footprints(tuple(p.stat().st_mtime_ns for p in FILES))
    sites=sorted(records.loc[records.kind.eq('homo'),'group'].unique(),key=lambda s:int(str(s).split('_')[-1]))
    if st.button('Compare reference F-actin sites with observed 5YU8 cofilactin sites',key='fp_reference_preset'):
        meta=pd.read_csv(FILES[0],low_memory=False)
        subset=meta[meta.pdb_id.eq('5yu8') & meta.s1_actine & meta.s2_actine]
        observed=set(subset.s1_binding_site_cluster_data_70.dropna())|set(subset.s2_binding_site_cluster_data_70.dropna())
        st.session_state['fp_reference']=[p for p in ['6685_1','6685_2','6685_3','6685_4'] if p in sites]
        st.session_state['fp_comparison']=[p for p in sites if p in observed]
    st.caption('The 5YU8 preset uses sites observed in a published cofilin-decorated filament. '
               'Footprints still aggregate all observations belonging to each selected site; '
               'this is not a classification of all minor interfaces.')
    a=st.multiselect('Reference actin–actin sites',sites,
                     default=[p for p in ['6685_1','6685_2','6685_3','6685_4'] if p in sites],key='fp_reference')
    b=st.multiselect('Comparison actin–actin sites',sites,key='fp_comparison')
    cutoff=st.slider('Buried ASA threshold for this comparison (%)',0.0,100.0,0.0,step=1.0,key='fp_cutoff')
    st.caption('Positions must have buried ASA strictly above the threshold in at least one '
               'record. Missing ASA is excluded. This control applies only to this comparison.')
    settings = dict(reference_sites=a, comparison_sites=b, buried_ASA_strictly_above=cutoff,
                    numbering='UniProt P60709', aggregation='union of observed positions',
                    source_sha256=source_signatures(tuple(p.stat().st_mtime_ns for p in FILES)))
    st.download_button('Download comparison settings and source identifiers (JSON)',
                       json.dumps(settings,indent=2),file_name='footprint_comparison_settings.json',
                       mime='application/json',key='fp_settings_download')
    retained=records[records.asa.gt(cutoff)]
    def positions(kind,groups):
        return set(retained.loc[retained.kind.eq(kind)&retained.group.isin(groups),'position'].astype(int))
    aa,bb=positions('homo',a),positions('homo',b)
    if not a and not b:
        st.info('Select at least one actin–actin site.')
        return
    groups={}
    if a:groups['Reference actin sites']=aa
    if b:groups['Comparison actin sites']=bb
    if set(a)&set(b):st.caption('Some sites occur in both selections; shared positions include their common contribution.')
    if a and b:
        st.dataframe(pd.DataFrame({'Position category':['Reference only','Shared','Comparison only'],
                                  'Residue count':[len(aa-bb),len(aa&bb),len(bb-aa)],
                                  'P60709 positions':[', '.join(map(str,sorted(x))) for x in [aa-bb,aa&bb,bb-aa]]}),
                     hide_index=True,width='stretch')
    import plotly.graph_objects as go
    names=sorted(records.loc[records.kind.eq('abp'),'group'].unique(),key=str.casefold)
    abps={name:positions('abp',[name]) for name in names}
    selected_abps=st.multiselect('ABP footprints aligned with actin–actin sites',names,key='fp_aligned_abps')
    displayed={**groups, **{f'ABP: {name}':abps[name] for name in selected_abps}}
    matrix=[[int(p in siteset) for p in range(1,376)] for siteset in displayed.values()]
    fig=go.Figure(go.Heatmap(z=matrix,x=list(range(1,376)),y=list(displayed),
                             colorscale=[[0,'#f2f2f2'],[.4999,'#f2f2f2'],[.5,'#0072B2'],[1,'#0072B2']],zmin=0,zmax=1,showscale=False,
                             customdata=[['Contact observed' if value else 'No contact above threshold' for value in row] for row in matrix],
                             hovertemplate='%{y}<br>P60709 position %{x}<br>%{customdata}<extra></extra>'))
    fig.update_layout(height=max(230,32*len(displayed)+100),xaxis_title=numbering.AXIS_TITLE,margin=dict(l=5,r=5,t=10,b=40))
    position_hover(fig)
    plotly_chart(fig,use_container_width=True,key='homo_footprint_comparison')
    from display_helpers import color_legend
    color_legend([('Contact observed', '#0072B2'), ('No contact above threshold', '#f2f2f2')])
    with st.expander('Compare positive contacts with the selected ASA threshold'):
        counts=[]
        for label,kind,selected in [('Reference actin sites','homo',a),('Comparison actin sites','homo',b)]+[(f'ABP: {n}','abp',[n]) for n in selected_abps]:
            if not selected: continue
            subset=records[records.kind.eq(kind)&records.group.isin(selected)]
            counts.append({'Group':label,'Positive ASA positions':subset.loc[subset.asa.gt(0),'position'].nunique(),
                           'Positions above selected threshold':subset.loc[subset.asa.gt(cutoff),'position'].nunique(),
                           'Threshold (%)':cutoff})
        st.dataframe(pd.DataFrame(counts),hide_index=True,width='stretch')
    result=pd.DataFrame([{'Partner source name':name,'Footprint residues':len(fp),
                          **{label:jaccard(fp,group) for label,group in groups.items()}}
                         for name,fp in abps.items()])
    st.markdown('**Jaccard overlap with ABP footprints**')
    st.caption('Jaccard = shared residues / union of residues (0–1). An empty union is '
               'undefined. Residue overlap alone demonstrates neither a steric clash nor '
               'competition; geometry and filament state require separate validation.')
    st.dataframe(result,hide_index=True,width='stretch')
    st.download_button('Download ABP–actin overlap (CSV)',result.to_csv(index=False).encode(),
                       file_name='abp_actin_jaccard.csv',mime='text/csv',key='fp_overlap_download')
    if st.checkbox('Show the complete ABP and selected actin-site Jaccard matrix',key='fp_all_matrix'):
        all_groups={**groups,**{f'ABP: {name}':fp for name,fp in abps.items()}}
        matrix=pd.DataFrame([[jaccard(x,y) for y in all_groups.values()] for x in all_groups.values()],index=all_groups,columns=all_groups)
        st.dataframe(matrix,width='stretch')
        st.download_button('Download complete Jaccard matrix (CSV)',matrix.to_csv().encode(),
                           file_name='abp_jaccard_matrix.csv',mime='text/csv',key='fp_matrix_download')
    from steric_screen import render_steric_screen
    render_steric_screen()
