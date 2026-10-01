"""Exercise native pages and representative controls, without starting calculations.
Run from either repository: python tools/check_app.py.
"""
import json
import os
import sys
from pathlib import Path
from unittest.mock import patch
from streamlit.testing.v1 import AppTest

root = Path(__file__).resolve().parents[1]
os.chdir(root)
sys.path.insert(0, str(root / 'script'))
from app_sections import SECTIONS
import proteocast_view

app = AppTest.from_file(str(root / 'script/streamlit.py'), default_timeout=180)


def check(label):
    app.run()
    errors = [e.message for e in app.exception]
    assert not errors, (label, errors)
    print(f'PASS: {label}', flush=True)


def page(section, view=None):
    app.switch_page(f'views/{section}.py')
    check(section)
    if view and app.radio(key=f'page_view_{section}').value != view:
        app.radio(key=f'page_view_{section}').set_value(view)
        check(f'{section} / {view}')


def choose_view(section, view):
    app.radio(key=f'page_view_{section}').set_value(view)
    check(f'{section} / {view}')


check('startup: documentation only')
assert len(app.get('plotly_chart')) == 0
assert not app.selectbox
assert not app.radio
assert any('What this app is' in m.value for m in app.markdown)
assert any(e.label == 'Data management — reload, download and calculations' for e in app.expander)
if Path('data/.slim_deploy').exists():
    assert 'Run / update' not in [b.label for b in app.button]
    assert not any(b.key == 'pc_run_all_missing' for b in app.button)
else:
    assert 'Run / update' in [b.label for b in app.button]
    assert any(b.key == 'pc_run_all_missing' for b in app.button)

# Exercise every page and each visible view selector with the installed data.
# Offline check: an unavailable remote AlphaFold fallback is simulated only;
# local ProteoCast structures and the rest of the renderers remain unchanged.
with patch.object(proteocast_view, '_fetch_af_pdb', return_value=None):
    for section, _, _ in SECTIONS:
        page(section)
        view_widget = next((r for r in app.radio if r.key == f'page_view_{section}'), None)
        if view_widget:
            for view in list(view_widget.options)[1:]:
                choose_view(section, view)

page('comparative-binding-sites', 'ABP networks')
app.radio(key='net_view').set_value('Cooperation')
check('ABP co-presence network')
assert not any(e.label == 'Show the network graph' for e in app.expander)

page('residue-level', 'Residue explorer')
assert len(app.selectbox(key='actin_ov_selbox').options) == 375
app.selectbox(key='actin_ov_selbox').set_value(51)  # P60709 M47
check('manual residue M47')
page('summary-tables', 'Structures')
assert 'actin_ov_selbox' not in {s.key for s in app.selectbox}
page('residue-level')
assert app.selectbox(key='actin_ov_selbox').value == 51
assert any('M47' in str(m.value) for m in app.markdown)
print('PASS: residue and view retained across native pages', flush=True)
app.selectbox(key='actin_ov_selbox').set_value(3)
check('first residue, including missing measurements')
assert any('M1' in str(m.value) for m in app.markdown)
app.selectbox(key='actin_ov_selbox').set_value(380)
check('last residue F375')
assert any('F375' in str(m.value) for m in app.markdown)

page('actin-actin-interfaces', 'Compare footprints')
app.multiselect(key='fp_comparison').set_value(['6685_17'])
app.slider(key='fp_cutoff').set_value(50.0)
app.checkbox(key='fp_all_matrix').check()
check('footprint groups, ASA threshold and full Jaccard matrix')
app.multiselect(key='fp_aligned_abps').set_value(['Cofilin-1'])
check('aligned ABP footprints')
choose_view('actin-actin-interfaces', 'Binding sites')
assert not any(e.label == 'Show binding-site table' for e in app.expander)
app.selectbox(key='homo_site_link').set_value('6685_3')
app.button(key='homo_site_open').click()
check('mixed binding-site navigation')
assert app.selectbox(key='sel_s1').value == '6685_3'
assert app.radio(key='page_view_comparative-binding-sites').value == 'Binding-site clusters'
# AppTest doesn't update its page hash after an application's st.switch_page;
# synchronize it with the already asserted destination before the next event.
app.switch_page('views/comparative-binding-sites.py')
# The existing graph's hidden bridge must now open the actual ABP page.
bridge = next(b for b in app.button if b.key and b.key.startswith('__ab_') and b.label == 'Cofilin-1')
bridge.click()
check('network ABP node opens protein page')
assert app.selectbox(key='sel_abp_detail').value == 'Cofilin-1'
assert app.header[0].value == 'ABP-actin interfaces'
app.switch_page('views/abp-actin-interfaces.py')
page('summary-tables')
page('abp-actin-interfaces')
assert app.selectbox(key='sel_abp_detail').value == 'Cofilin-1'

# A heatmap selection must open cluster details, including the clicked residue.
page('residue-level', 'Binding-site heatmap')
chart = app.get('plotly_chart')[0]
states = app._tree.get_widget_states()
state = next((w for w in states.widgets if w.id == chart.proto.id), None)
if state is None:
    state = states.widgets.add(id=chart.proto.id)
state.string_value = json.dumps({'selection': {'points': [
    {'x': 3, 'y': '6685_0', 'curve_number': 1, 'point_index': 2}],
    'point_indices': [2], 'box': [], 'lasso': []}})
app._run(states)
assert not app.exception, [e.message for e in app.exception]
assert app.selectbox(key='sel_s1').value == '6685_0'
assert app.selectbox(key='s1posdet_6685_0').value == 7
assert app.header[0].value == 'Comparative analyses of binding sites'
app.switch_page('views/comparative-binding-sites.py')
print('PASS: heatmap opens matching cluster and residue', flush=True)

page('comparative-binding-sites', 'ABP pairs and sequences')
app.radio(key='explo_mode').set_value('ABP pair')
check('comparison of two ABPs')
if Path('data/human_variants/manifest.json').exists():
    page('human-actin-variants')
    app.selectbox(key='hv_gene').set_value('ACTG2')
    app.selectbox(key='hv_category').set_value('uncertain')
    check('human ACTG2 uncertain variants and disease associations')
    app.selectbox(key='hv_gene').set_value('ACTB')
    app.selectbox(key='hv_category').set_value('conflicting')
    check('conflicting variants')
    app.selectbox(key='hv_conflict_abp').set_value('Cofilin-1')
    app.selectbox(key='hv_cross_pair').set_value(('ACTG2', 'ACTA1'))
    check('cross-gene observations and aligned conflicting-variant tracks')
page('actin-conservation', 'Overview')
app.slider(key='cons_asa').set_value(100.0)
check('empty conservation footprints')
page('summary-tables', 'Dataset checks')
app.selectbox(key='structure_audit_filter').set_value('Mutation annotated')
check('cached RCSB mutation annotations')
page('interface-properties', 'Actin surface chemistry')
app.selectbox(key='surface_chem_weight').set_value('Partner buried ASA (%)')
app.checkbox(key='surface_chem_3d').check()
check('partner chemistry and generated 3D surfaces')
app.slider(key='surface_chem_cutoff').set_value(100.0)
check('empty chemistry footprints')
page('actin-conservation', 'Solvent accessibility')
app.slider(key='filament_rsa_cutoff').set_value(1.0)
check('strict filament-accessibility threshold')
print('All page and interaction checks passed. Browser/3D and scientific validation remain separate.')
