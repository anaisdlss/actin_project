"""Observed partner chemistry, with explicit weighting and structural sampling."""
from pathlib import Path
import numpy as np
import pandas as pd
import streamlit as st
from residue_contacts import orient_actin_contacts

CLASSES = {
    'Hydrophobic': set('AVILMFWP'), 'Polar': set('STCYNQ'),
    'Basic (K/R/H)': set('KRH'), 'Acidic (D/E)': set('DE'),
    'Glycine': set('G'), 'Other / modified': set(),
}
COLORS = dict(zip(CLASSES, ['#0072B2', '#56B4E9', '#E69F00', '#884EA0', '#999999', '#555555']))
WEIGHTS = {'Residue-pair count': None, 'Pair contact area (Å²)': 'Pair contact area (Å²)',
           'Partner buried ASA (%)': 'Partner buried ASA (%)'}
SOURCES = [Path('data/filtered/details/4.inter-residue_contacts.csv'),
           Path('data/filtered/proteins_per_pdb.csv'),
           Path('data/filtered/filtered_all_data.csv'),
           Path('data/filtered/details/1.interactions.csv')]


def aa_class(aa):
    return next((name for name, letters in CLASSES.items() if str(aa).upper() in letters), 'Other / modified')


def chemistry_records(contacts, proteins, all_data, interactions):
    rows = orient_actin_contacts(contacts, proteins)
    rows = rows[rows['Interaction type'].eq('Actin–ABP')].copy()
    meta = all_data.merge(interactions[['interaction_id','chain_A_id','chain_B_id']],
                          left_on=['subunit_1','subunit_2'], right_on=['chain_A_id','chain_B_id'])
    parts = []
    for side, other in ((1,2),(2,1)):
        selected = meta[meta[f's{side}_actine'].astype(str).str.lower().eq('true') &
                        meta[f's{other}_actine'].astype(str).str.lower().eq('false')]
        parts.append(selected[['interaction_id',f'subunit_{side}',f'subunit_{other}',f's{side}_binding_site_cluster_data_70']].rename(
            columns={f'subunit_{side}':'Actin chain',f'subunit_{other}':'Partner chain',
                     f's{side}_binding_site_cluster_data_70':'Binding site'}))
    rows = rows.merge(pd.concat(parts).drop_duplicates(), on=['interaction_id','Actin chain','Partner chain'], validate='many_to_one')
    rows['Partner class'] = rows['Partner aa'].map(aa_class)
    # A residue pair is one observation even if repeated with a different text label.
    pair_keys = ['interaction_id','Actin chain','Partner chain','Actin PDB residue','Partner PDB residue']
    if rows.duplicated(pair_keys).any():
        measures = ['Pair contact area (Å²)','Partner buried ASA (%)','Actin buried ASA (%)']
        if (rows.groupby(pair_keys)[measures].nunique(dropna=False) > 1).any().any():
            raise ValueError('Conflicting measurements for a repeated residue pair; review the contact source.')
    return rows.drop_duplicates(pair_keys)


def chemistry_summary(records, weighting='Pair contact area (Å²)', cutoff=0.0):
    """Class fractions per interface/position, equal interfaces then PDB means.

    Input values below a censoring limit remain missing in area mode. No zero
    contact is inferred for positions not represented in the source.
    """
    if weighting not in WEIGHTS:
        raise ValueError('Unknown weighting')
    rows = records[records['P60709 position'].notna() & records['Actin buried ASA (%)'].gt(cutoff)].copy()
    col = WEIGHTS[weighting]
    rows['weight'] = 1.0 if col is None else pd.to_numeric(rows[col], errors='coerce')
    rows = rows[np.isfinite(rows.weight) & rows.weight.gt(0)].copy()
    keys = ['Partner','PDB','P60709 position','interaction_id','Actin chain','Partner chain']
    amounts = rows.groupby(keys + ['Partner class'], observed=True).weight.sum().unstack('Partner class', fill_value=0).reindex(columns=CLASSES, fill_value=0)
    fractions = amounts.div(amounts.sum(axis=1), axis=0)
    fractions = fractions.groupby(level=['Partner','PDB','P60709 position']).mean()
    per_abp = fractions.groupby(level=['Partner','P60709 position']).mean()
    global_profile = per_abp.groupby(level='P60709 position').mean().reindex(range(1,376))
    n_pdb = fractions.groupby(level=['Partner','P60709 position']).size().rename('PDBs with usable contacts')
    per_abp = per_abp.join(n_pdb).reset_index()
    global_profile['ABP source names'] = per_abp.groupby('P60709 position').Partner.nunique().reindex(range(1,376), fill_value=0)
    global_profile.index.name = 'P60709 position'
    global_profile['Dominant partner class'] = dominant_classes(global_profile[list(CLASSES)])
    return global_profile.reset_index(), per_abp, rows


def dominant_classes(fractions):
    def winner(row):
        if row.isna().all() or row.sum() <= 0:
            return 'No usable contact'
        maxima = row.index[np.isclose(row, row.max(), equal_nan=False)]
        return str(maxima[0]) if len(maxima) == 1 else 'Tied classes'
    return fractions.apply(winner, axis=1)


def family_composition(records, weighting, cutoff, families=None):
    _, per_abp, _ = chemistry_summary(records, weighting, cutoff)
    # Equal residue positions within each ABP: no weighting by chain multiplicity.
    compositions = per_abp.groupby('Partner')[list(CLASSES)].mean()
    compositions['Family annotation'] = compositions.index.map(lambda p: (families or {}).get(p, 'Unresolved'))
    return compositions.reset_index()


@st.cache_data(show_spinner=False)
def load_chemistry(stamp):
    data = [pd.read_csv(p, low_memory=False, dtype={'residue_A_structure':str,'residue_B_structure':str}) for p in SOURCES]
    return chemistry_records(*data)


def render_surface_chemistry(profile):
    """Project on one real actin chain with checked sequence correspondence."""
    import io
    from Bio.PDB import PDBParser, is_aa
    from Bio.SeqUtils import seq1
    from Bio.Data.PDBData import protein_letters_3to1
    from scientific_analysis import isoform_map, sequence
    import py3Dmol
    pdb = Path('data/filtered/details/structures_files/assembly/7pdz.pdb')
    ref = Path('data/P60709_ref.fasta')
    if not pdb.exists() or not ref.exists():
        st.info('The local 7PDZ reference and P60709 sequence are required for this surface view.'); return
    model = next(PDBParser(QUIET=True).get_structure('reference', io.StringIO(pdb.read_text())).get_models())
    if 'I' not in model:
        st.info('Reference actin chain I is absent.'); return
    residues = [r for r in model['I'] if is_aa(r) and 'CA' in r]
    seq = ''.join(seq1(r.resname, custom_map={'HIC':'H'}) for r in residues)
    reference = sequence(ref)
    try:
        mapped = isoform_map(seq, reference)
    except ValueError:
        st.info('Ambiguous sequence correspondence; surface projection unavailable.'); return
    exact = {i:p for i,p in mapped.items() if seq[i-1] == reference[p-1]}
    if len(exact) < 300 or len(exact)/max(1,len(mapped)) < .95:
        st.info('Reference sequence correspondence failed validation.'); return
    residue_keys = {(r.id[1], r.id[2]) for r in residues}
    atom_lines = [line for line in pdb.read_text().splitlines()
                  if line.startswith(('ATOM','HETATM')) and line[21:22] == 'I'
                  and (int(line[22:26]), line[26:27]) in residue_keys]
    pdb_text = '\n'.join(atom_lines)
    dominant = profile.set_index('P60709 position')['Dominant partner class'].to_dict()
    st.caption('Projection onto experimental 7PDZ chain I, aligned to P60709. Only matching reference '
               'residues are coloured. Grey = no usable contact, tie or unmapped position; it is not an accessibility measurement. '
               'The partner view pools observed ABP contacts; it is not a predicted binding surface.')
    for col, title, complementary in zip(st.columns(2), ['Actin residue classes','Observed partner classes'], [False,True]):
        with col:
            st.markdown(f'**{title}**')
            view=py3Dmol.view(width='100%',height=360); view.addModel(pdb_text,'pdb')
            view.setStyle({}, {'cartoon':{'color':'#dddddd'}})
            residue_colors = {}
            for i,p in exact.items():
                category = dominant.get(p) if complementary else (aa_class(reference[p-1])
                    if residues[i-1].resname in protein_letters_3to1 else 'Other / modified')
                resid = residues[i-1].id
                residue_colors[(resid[1],resid[2])] = COLORS.get(category,'#dddddd')
            # Explicit surface colours, with the complete protein defining geometry.
            # See https://3dmol.org/doc/GLViewer.html#addSurface (atomsel/allsel).
            serials = {}
            for line in atom_lines:
                color = residue_colors.get((int(line[22:26]), line[26:27]), '#dddddd')
                serials.setdefault(color, []).append(int(line[6:11]))
            for color, atoms in serials.items():
                view.addSurface(py3Dmol.SES, {'opacity':1.0,'color':color}, {'serial':atoms}, {})
            view.setBackgroundColor('white'); view.zoomTo()
            st.components.v1.html(view._make_html(),height=365)


def render_interface_properties():
    if not all(p.exists() for p in SOURCES): return
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from scientific_analysis import sequence
    from plot_interaction import position_hover
    st.subheader('Actin surface and complementary partner chemistry')
    records = load_chemistry(tuple((p.stat().st_mtime_ns,p.stat().st_size) for p in SOURCES))
    weighting=st.selectbox('Chemistry weighting', list(WEIGHTS), index=1,key='surface_chem_weight')
    cutoff=st.slider('Buried actin ASA threshold for chemistry (%)',0.0,100.0,0.0,step=1.0,key='surface_chem_cutoff')
    st.caption('Only observed actin–ABP residue pairs are used. Contacts require actin buried ASA strictly '
               'above this threshold; zero keeps all positive contacts. Fractions are normalized within '
               'each interface / actin position, averaged equally across observed interfaces within each PDB / ABP, '
               'then equally across PDBs and across ABP source names. '
               'Missing, censored or nonpositive weights are excluded, not imputed. K/R/H are grouped as basic; '
               'these classes do not predict protonation or electrostatic potential.')
    if weighting == 'Partner buried ASA (%)':
        st.caption('This exploratory weight is the partner residue’s burial in its whole interface, not the area of this residue pair.')
    profile, per_abp, usable = chemistry_summary(records, weighting, cutoff)
    st.caption(f'{len(usable):,} usable residue pairs out of {len(records):,} retained actin–ABP pairs (before threshold and weight checks).')
    fractions=profile.set_index('P60709 position')[list(CLASSES)]
    ref_path=Path('data/P60709_ref.fasta')
    reference=sequence(ref_path) if ref_path.exists() else 'X'*375
    letters=list(reference) if len(reference)==375 else ['?']*375
    classes=[aa_class(aa) for aa in letters]
    palette=[]
    for i, color in enumerate(COLORS.values()):
        palette.extend([[i/len(CLASSES),color],[(i+1)/len(CLASSES),color]])
    fig=make_subplots(rows=2,cols=1,shared_xaxes=True,row_heights=[.15,.85],vertical_spacing=.12,
                      subplot_titles=['Actin residue class (P60709 sequence)','Observed partner class fractions'])
    fig.add_trace(go.Heatmap(x=fractions.index,y=['Actin'],z=[[list(CLASSES).index(c) for c in classes]],
        zmin=-.5,zmax=len(CLASSES)-.5,colorscale=palette,showscale=False,
        customdata=[list(zip(letters,classes))],
        hovertemplate='%{customdata[0]}%{x}: %{customdata[1]}<extra></extra>'),row=1,col=1)
    fig.add_trace(go.Heatmap(x=fractions.index,y=list(CLASSES),z=fractions.T.values,zmin=0,zmax=1,colorscale='Blues',
        customdata=[letters]*len(CLASSES),colorbar=dict(title='Fraction',len=.72,y=.38),hoverongaps=False,
        hovertemplate='%{y}<br>%{customdata}%{x}: %{z:.3f}<extra></extra>'),row=2,col=1)
    fig.update_layout(height=410,xaxis2_title='P60709 position')
    st.plotly_chart(position_hover(fig),width='stretch',key='surface_chem_profile')
    st.caption('Class colours in 3D: '+ ' · '.join(f'{name}: {color}' for name,color in COLORS.items()))
    with st.expander('View actin and partner chemistry in 3D'):
        if st.checkbox('Load chemistry surfaces',key='surface_chem_3d'):
            render_surface_chemistry(profile)
    st.download_button('Download partner chemistry per actin position',profile.to_csv(index=False),file_name='actin_partner_chemistry.csv',key='surface_chem_csv')
    st.download_button('Download chemistry per ABP and position',per_abp.to_csv(index=False),file_name='abp_position_chemistry.csv',key='surface_chem_abp_csv')
    with st.expander('Threshold comparison and calculation details'):
        import hashlib
        import json
        baseline, _, baseline_rows = chemistry_summary(records,weighting,0.0)
        st.dataframe(pd.DataFrame({'Threshold (%)':[0.0,cutoff],
            'Usable residue pairs':[len(baseline_rows),len(usable)],
            'Positions with usable contacts':[int(baseline['ABP source names'].gt(0).sum()),int(profile['ABP source names'].gt(0).sum())]}),hide_index=True)
        st.download_button('Download retained contact measurements',usable.to_csv(index=False),file_name='chemistry_contact_measurements.csv',key='surface_chem_raw')
        manifest={'weighting':weighting,'actin_buried_ASA_strictly_above':cutoff,
                  'aggregation':'class fractions per interaction/actin-chain/partner-chain/position, equal interface means within PDB/ABP, then equal PDB and ABP means',
                  'classes':{k:sorted(v) for k,v in CLASSES.items()},
                  'sources':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in SOURCES}}
        st.download_button('Download chemistry settings and source identifiers',json.dumps(manifest,indent=2),file_name='chemistry_settings.json',key='surface_chem_manifest')
    with st.expander('Compare ABP families sharing an actin binding site'):
        sites=sorted(records['Binding site'].dropna().unique())
        if not sites: return
        site=st.selectbox('Site for chemistry comparison',sites,key='surface_chem_site')
        fam_path=Path('data/exports/abp_site_domain/familles.csv'); families={}
        if fam_path.exists():
            for row in pd.read_csv(fam_path).itertuples():
                families.update({name.strip():row.famille for name in str(row.membres).split(' · ')})
        comp=family_composition(records[records['Binding site'].eq(site)],weighting,cutoff,families)
        st.caption('Each ABP composition averages its contacted positions equally, after the PDB normalization above. '
                   'Family labels are existing source annotations; similar compositions alone do not prove a shared motif.')
        st.dataframe(comp,hide_index=True,width='stretch')
        if not comp.empty:
            fig=go.Figure(go.Heatmap(x=list(CLASSES),y=comp.Partner,z=comp[list(CLASSES)].values,zmin=0,zmax=1,colorscale='Blues',colorbar=dict(title='Fraction')))
            fig.update_layout(height=max(280,28*len(comp)+130))
            st.plotly_chart(fig,width='stretch',key='surface_chem_families')
        st.download_button('Download binding-site chemistry comparison',comp.to_csv(index=False),file_name=f'{site}_partner_chemistry.csv',key='surface_chem_family_csv')
