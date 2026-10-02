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
app.switch_page('views/residue-level.py').run()
app.selectbox(key='actin_ov_selbox').set_value(51).run()
app.button(key='residue_conservation_51').click().run()
assert not app.exception,[e.message for e in app.exception]
assert app.header[0].value=='Actin binding sites conservation'
assert app.selectbox(key='actin_ov_selbox').value==51
print('PASS selected residue opens its conservation detail',flush=True)
