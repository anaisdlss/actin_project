"""Exercise the screenshot regressions without remote calculations."""
import os,sys,json
from pathlib import Path
from unittest.mock import patch
import streamlit
from streamlit.testing.v1 import AppTest
root=Path(__file__).resolve().parents[1];os.chdir(root);sys.path.insert(0,str(root/'script'))
import proteocast_view

def check():
    assert not app.exception,[e.message for e in app.exception]

app=AppTest.from_file(str(root/'script/streamlit.py'),default_timeout=180).run()
app.switch_page('views/summary-tables.py').run()
app.radio(key='page_view_summary-tables').set_value('Source tables').run()
app.selectbox(key='source_table').set_value('1.interactions.csv').run();check()
frame=next(d.value for d in app.dataframe if 'interaction_id' in d.value)
assert frame.interaction_id.is_monotonic_increasing
assert frame.interaction_id.dtype.kind in 'iu'
print('PASS source interaction IDs are numeric and ordered',flush=True)
app.switch_page('views/comparative-binding-sites.py').run()
app.radio(key='page_view_comparative-binding-sites').set_value('Binding-site clusters').run()
app.selectbox(key='sel_s1').set_value('6685_6').run();check()
select=app.selectbox(key='abp3d_v2_6685_6')
assert any('Myosin-6' in name for name in select.options)
from abp_3d import _s1_abp_3d_options
choice=next(o for o in _s1_abp_3d_options('6685_6',os.path.getmtime('data/filtered/filtered_all_data.csv')) if 'Myosin-6' in o['label'])
select.set_value(dict(choice,mode='one')).run();check()
assert any('Observed pair:' in c.value for c in app.caption)
html=' '.join(i.proto.srcdoc for i in app.get('iframe'))
assert 'ResizeObserver' in html
select=app.selectbox(key='abp3d_v2_6685_6')
select.set_value({'label':'All available partners, aligned on actin','mode':'all'}).run();check()
assert any('available pairs aligned' in c.value for c in app.caption)
print('PASS individual myosin pair and aligned partners',flush=True)
with patch.object(proteocast_view,'_fetch_af_pdb',return_value=None):
    app.switch_page('views/abp-conservation.py').run()
    app.radio(key='page_view_abp-conservation').set_value('3D structure').run()
    protein=app.selectbox(key='conservation_abp')
    has_confidence_only='Actin-interacting protein 1' in protein.options
    protein.set_value('Actin-interacting protein 1' if has_confidence_only else 'Actin-related protein 2').run();check()
    assert app.header[0].value=='ABP mutational sensitivity'
    assert any('Zoom on the observed actin-contact region' in t.label for t in app.toggle)
    expected='AlphaFold confidence' if has_confidence_only else 'Mutational sensitivity (−mean'
    assert any(expected in m.value for m in app.markdown)
print('PASS explicit metric legend and persistent region-zoom control',flush=True)
app.switch_page('views/human-actin-variants.py').run();check()
for chart in app.get('plotly_chart'):
    figure=json.loads(chart.proto.spec);config=json.loads(chart.proto.config)
    assert config['displayModeBar'] is True
    for trace in figure.get('data',[]):
        if trace['type']=='heatmap':assert '<span style="color:' in str(trace.get('hovertext'))
print('PASS visible zoom tools and cell-colour tooltip payloads',flush=True)
