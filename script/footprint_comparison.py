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


def footprint_sources(all_data, interactions):
    """All retained actin-side observations, including those without contact rows."""
    merged = all_data.merge(interactions[['interaction_id', 'chain_A_id', 'chain_B_id']],
                            left_on=['subunit_1', 'subunit_2'], right_on=['chain_A_id', 'chain_B_id'])
    parts = []
    for side, other in ((1, 2), (2, 1)):
        rows = merged[merged[f's{side}_actine'].astype(str).str.lower().eq('true')].copy()
        homo = rows[f's{other}_actine'].astype(str).str.lower().eq('true')
        rows['kind'] = np.where(homo, 'homo', 'abp')
        rows['group'] = np.where(homo, rows[f's{side}_binding_site_cluster_data_70'], rows[f'subunit_{other}_title'])
        rows['chain'] = rows[f'subunit_{side}']
        parts.append(rows[['kind', 'group', 'interaction_id', 'chain', 'pdb_id']])
    return pd.concat(parts, ignore_index=True).dropna(subset=['group']).drop_duplicates()


def mean_footprint(sources, residues, resolved, kind, groups):
    """Union selected contacts per chain; mean chains within PDB, then mean PDBs.

    Only resolved mapped positions can supply non-contact zeros. Explicit invalid
    contact measurements stay missing. Repeated interactions never add PDB weight.
    """
    selected = sources[sources.kind.eq(kind) & sources.group.isin(groups)]
    contacts_all = residues.merge(selected[['interaction_id', 'chain']].drop_duplicates(),
                                  on=['interaction_id', 'chain'])
    contacts_all['position'] = numbering.uniprot_series(contacts_all.residue_number_canon_mafft)
    contacts_all['asa'] = pd.to_numeric(contacts_all.buried_ASA_percent.astype(str).str.replace('%', '', regex=False), errors='coerce')
    contacts_all.loc[~np.isfinite(contacts_all.asa) | ~contacts_all.asa.between(0, 100), 'asa'] = np.nan
    by_chain = {chain: rows.dropna(subset=['position']).groupby('position').asa.max()
                for chain, rows in contacts_all.groupby('chain')}
    values = []
    for chain, observations in selected.groupby('chain'):
        contacts = by_chain.get(chain, pd.Series(dtype=float))
        row = pd.Series(np.nan, index=range(1, 376), dtype=float)
        present = sorted(set(resolved.get(chain, ())) & set(row.index))
        row.loc[present] = 0.0
        # A measured contact remains usable even if its coordinates are unavailable.
        for position, value in contacts.items():
            if position in row.index:
                row.loc[position] = value
        row['pdb_id'] = observations.pdb_id.iloc[0]
        values.append(row)
    if not values:
        return pd.DataFrame({'mean_ASA_percent': np.nan, 'PDB_count': 0, 'chain_count': 0}, index=range(1, 376))
    chains = pd.DataFrame(values)
    pdbs = chains.groupby('pdb_id').mean()
    return pd.DataFrame({'mean_ASA_percent': pdbs.mean(), 'PDB_count': pdbs.count(),
                         'chain_count': chains.drop(columns='pdb_id').count()}).reindex(range(1, 376))


@st.cache_data(show_spinner='Checking resolved actin residues…')
def load_mean_inputs(stamp):
    from Bio.PDB import PDBParser, MMCIFParser
    from Bio.PDB.Polypeptide import protein_letters_3to1_extended
    from scientific_analysis import sequence, isoform_map
    sources = footprint_sources(pd.read_csv(FILES[0], low_memory=False), pd.read_csv(FILES[1]))
    residues = pd.read_csv(FILES[2], low_memory=False)
    reference = sequence(Path('data/P60709_ref.fasta'))
    resolved, sequence_maps = {}, {}
    base = Path('data/filtered/details/structures_files/assembly')
    for pdb, rows in sources.groupby('pdb_id'):
        file = next((base / f'{pdb}{ext}' for ext in ('.pdb', '.cif') if (base / f'{pdb}{ext}').exists()), None)
        if file is None:
            continue
        parser = MMCIFParser(QUIET=True) if file.suffix == '.cif' else PDBParser(QUIET=True)
        try:
            model = next(iter(parser.get_structure(str(pdb), str(file))))
        except (ValueError, KeyError, OSError, StopIteration):
            continue
        needed = set(rows.chain)
        for chain in model:
            label = f'{pdb}_{chain.id}'
            if label not in needed:
                continue
            observed = [r for r in chain if 'CA' in r and r.resname in protein_letters_3to1_extended]
            seq = ''.join(protein_letters_3to1_extended[r.resname] for r in observed)
            if not seq:
                continue
            if seq not in sequence_maps:
                sequence_maps[seq] = set(isoform_map(seq, reference).values())
            resolved[label] = sequence_maps[seq]
    return sources, residues, resolved


def mean_input_stamp():
    paths = FILES + [Path('data/P60709_ref.fasta')]
    paths += sorted(Path('data/filtered/details/structures_files/assembly').glob('*'))
    return tuple((str(p), p.stat().st_mtime_ns, p.stat().st_size) for p in paths if p.is_file())


@st.cache_data(show_spinner=False)
def load_footprints(mtimes):
    return footprint_records(*(pd.read_csv(p,low_memory=False) for p in FILES))


@st.cache_data(show_spinner=False)
def source_signatures(mtimes):
    return {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in FILES}


def render_footprint_comparison():
    if not all(p.exists() for p in FILES):
        return
    st.caption('Compare footprints on actin. The initial selection contains the four principal '
               'actin–actin sites 6685_1–4. Both sides of actin–actin interactions are included.')
    records=load_footprints(tuple(p.stat().st_mtime_ns for p in FILES))
    sites=sorted(records.loc[records.kind.eq('homo'),'group'].unique(),key=lambda s:int(str(s).split('_')[-1]))
    if st.button('Compare reference F-actin sites with observed 5YU8 cofilactin sites',key='fp_reference_preset'):
        meta=pd.read_csv(FILES[0],low_memory=False)
        subset=meta[meta.pdb_id.eq('5yu8') & meta.s1_actine & meta.s2_actine]
        observed=set(subset.s1_binding_site_cluster_data_70.dropna())|set(subset.s2_binding_site_cluster_data_70.dropna())
        st.session_state['fp_reference']=[p for p in ['6685_1','6685_2','6685_3','6685_4'] if p in sites]
        st.session_state['fp_comparison']=[p for p in sites if p in observed]
    st.caption('The 5YU8 preset selects sites observed in a cofilin-decorated filament. '
               'Profiles aggregate the retained observations for these sites, not only 5YU8.')
    a=st.multiselect('Reference actin–actin sites',sites,
                     default=[p for p in ['6685_1','6685_2','6685_3','6685_4'] if p in sites],key='fp_reference')
    b=st.multiselect('Comparison actin–actin sites',sites,key='fp_comparison')
    mode=st.radio('Footprint display', ['Mean buried ASA (%)', 'Contact observed (blue/grey)'],
                  horizontal=True, key='fp_display_mode')
    cutoff=st.slider('Buried ASA threshold for contact overlap (%)',0.0,100.0,0.0,step=1.0,key='fp_cutoff')
    st.caption('Positions must have buried ASA strictly above the threshold in at least one '
               'record. This threshold controls the binary view and overlap tables; it does not filter the mean ASA.')
    settings = dict(reference_sites=a, comparison_sites=b, buried_ASA_strictly_above=cutoff,
                    numbering='UniProt P60709', display=mode,
                    mean_aggregation='maximum selected contact per actin chain; mean chains within PDB; mean PDBs',
                    zeros='resolved mapped residue without a selected contact', missing='unresolved residue or invalid contact measurement',
                    overlap_aggregation='union of observed positions',
                    source_sha256=source_signatures(tuple(p.stat().st_mtime_ns for p in FILES)))
    st.download_button('Download comparison settings and source identifiers (JSON)',
                       json.dumps(settings,indent=2),file_name='footprint_comparison_settings.json',
                       mime='application/json',key='fp_settings_download')
    retained=records[records.asa.gt(cutoff)]
    def positions(kind,groups):
        return set(retained.loc[retained.kind.eq(kind)&retained.group.isin(groups),'position'].astype(int))
    aa,bb=positions('homo',a),positions('homo',b)
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
    if not displayed:
        st.info('Select actin–actin sites or one or more ABPs to compare.')
        return
    if mode == 'Mean buried ASA (%)':
        sources, residues, resolved = load_mean_inputs(mean_input_stamp())
        selections = [(label, 'homo', a if label == 'Reference actin sites' else b) for label in groups]
        selections += [(f'ABP: {name}', 'abp', [name]) for name in selected_abps]
        profiles = {label: mean_footprint(sources, residues, resolved, kind, chosen)
                    for label, kind, chosen in selections}
        matrix = [profile.mean_ASA_percent.to_numpy() for profile in profiles.values()]
        counts = [profile[['PDB_count', 'chain_count']].to_numpy() for profile in profiles.values()]
        fig=go.Figure(go.Heatmap(z=matrix,x=list(range(1,376)),y=list(profiles),
                                colorscale='YlOrRd', zmin=0, zmax=100,
                                colorbar=dict(title='Mean buried ASA (%)'), customdata=counts,
                                hovertemplate='%{y}<br>P60709 position %{x}<br>Mean buried ASA: %{z:.2f}%<br>PDBs: %{customdata[0]}<br>Actin chains: %{customdata[1]}<extra></extra>',
                                hoverongaps=False))
        st.caption('Mean buried ASA uses a fixed 0–100% scale. Resolved residues without a selected contact '
                   'count as zero; unresolved residues and invalid measurements are excluded. Selected contacts '
                   'are combined by their maximum for each actin chain, then averaged within each PDB and '
                   'across PDBs with equal PDB weight. Hover shows the number of contributing PDBs and chains. '
                   'Blank cells have no usable measurement; these percentages are not binding strength.')
        exported = pd.concat([profile.rename_axis('P60709_position').reset_index().assign(Group=label)
                              for label, profile in profiles.items()], ignore_index=True)
        st.download_button('Download mean ASA and observation counts (CSV)', exported.to_csv(index=False).encode(),
                           file_name='footprint_mean_asa.csv', key='fp_mean_csv')
    else:
        matrix=[[int(p in siteset) for p in range(1,376)] for siteset in displayed.values()]
        fig=go.Figure(go.Heatmap(z=matrix,x=list(range(1,376)),y=list(displayed),
                             colorscale=[[0,'#f2f2f2'],[.4999,'#f2f2f2'],[.5,'#0072B2'],[1,'#0072B2']],zmin=0,zmax=1,showscale=False,
                             customdata=[['Contact observed' if value else 'No contact above threshold' for value in row] for row in matrix],
                             hovertemplate='%{y}<br>P60709 position %{x}<br>%{customdata}<extra></extra>'))
    fig.update_layout(height=max(230,32*len(displayed)+100),xaxis_title=numbering.AXIS_TITLE,margin=dict(l=5,r=5,t=10,b=40))
    position_hover(fig)
    plotly_chart(fig,use_container_width=True,key='homo_footprint_comparison')
    from display_helpers import color_legend
    if mode != 'Mean buried ASA (%)':
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
