"""Exercise saved FoldDisco searches, combined motifs and residue navigation. No submission.
Run from either repository: python tools/check_followup.py.
"""
import os,sys
from pathlib import Path
from streamlit.testing.v1 import AppTest
root=Path(__file__).resolve().parents[1];os.chdir(root);sys.path.insert(0,str(root/'script'))
app=AppTest.from_file(str(root/'script/streamlit.py'),default_timeout=180)
app.run();app.switch_page('views/homolog-search.py').run()
assert not app.exception,[e.message for e in app.exception]
app.selectbox(key='homolog_abp').set_value('Inverted formin-2').run()
app.selectbox(key='disco_cl_Inverted formin-2').set_value('6685_213').run()
assert not app.exception,[e.message for e in app.exception]
assert any('Completed: 2000' in str(m.value) for m in app.markdown)
assert any(b.proto.label=='Download new search results' for b in app.get('download_button'))
print('PASS tracked new results displayed with the exact matching motif',flush=True)
app.selectbox(key='disco_cl_Inverted formin-2').set_value('6685_0').run()
app.multiselect(key='fd_combination_Inverted formin-2_6685_0').set_value(['6685_200']).run()
assert not app.exception,[e.message for e in app.exception]
assert any('9AZ4' in str(m.value) for m in app.markdown)
assert any('41 distinct contact positions' in str(m.value) for m in app.markdown)
if not (root/'data/.slim_deploy').exists():assert any('32 residues' in i.value for i in app.info)
print('PASS coherent combined motif and oversized-query limit',flush=True)
app.selectbox(key='homolog_abp').set_value('Adducin 1').run()
site=app.selectbox(key='disco_cl_Adducin 1').value
assert any('2000 returned alignments' in str(m.value) for m in app.markdown)
assert any(b.proto.label=='Download matching prepared chain (PDB)' for b in app.get('download_button'))
positions=app.text_area(key=f'fd_query_positions_Adducin 1_{site}')
positions.set_value(', '.join(positions.value.split(',')[:4])).run()
assert not app.exception,[e.message for e in app.exception]
assert any('2000 returned alignments' in str(m.value) for m in app.markdown)
assert any(b.proto.label=='Download saved search results (CSV)' for b in app.get('download_button'))
print('PASS corrected numeric chain results remain accessible after motif editing',flush=True)
app.selectbox(key='homolog_abp').set_value('Beta-adducin').run()
app.selectbox(key='disco_cl_Beta-adducin').set_value('6685_155').run()
assert not app.exception,[e.message for e in app.exception]
assert any(b.proto.label=='Download local matches and contact checks' for b in app.get('download_button'))
assert any('44 residues' in c.value for c in app.caption)
print('PASS local results for the 44-residue motif are accessible',flush=True)
app.switch_page('views/actin-actin-interfaces.py').run()
app.radio(key='page_view_actin-actin-interfaces').set_value('Compare footprints').run()
app.selectbox(key='steric_abp').set_value('Vinculin').run()
assert not app.exception,[e.message for e in app.exception]
assert any(b.proto.label=='Download all rigid-placement measurements' for b in app.get('download_button'))
print('PASS rigid-placement measurements accessible from footprint comparison',flush=True)
app.switch_page('views/human-actin-variants.py').run()
for category in ['pathogenic','likely_pathogenic']:
    app.selectbox(key='hv_category').set_value(category).run()
    assert not app.exception,[e.message for e in app.exception]
    assert any('degrees of disease severity' in i.value for i in app.info)
print('PASS both requested clinical categories',flush=True)
app.switch_page('views/residue-level.py').run()
app.selectbox(key='actin_ov_selbox').set_value(51).run()
app.button(key='residue_conservation_51').click().run()
assert not app.exception,[e.message for e in app.exception]
assert app.header[0].value=='Actin mutational sensitivity'
assert app.selectbox(key='actin_ov_selbox').value==51
print('PASS selected residue opens its conservation detail',flush=True)
