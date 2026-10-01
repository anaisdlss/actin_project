"""Run representative Streamlit paths without launching external calculations.
Usage: python tools/check_app.py (from the repository root).
"""
import os
import sys
from pathlib import Path
from streamlit.testing.v1 import AppTest

root = Path(__file__).resolve().parents[1]
os.chdir(root)
sys.path.insert(0, str(root / 'script'))
app = AppTest.from_file(str(root / 'script/streamlit.py'), default_timeout=180)

def check(label):
    app.run()
    errors = [e.message for e in app.exception]
    assert not errors, (label, errors)
    print(f'PASS: {label}', flush=True)

check('startup')
from network_viz import _BIP_REQUIRED_FILES
if all(Path(p).exists() for p in _BIP_REQUIRED_FILES):
    assert not any('Run the pipeline to generate the network data' in x.value for x in app.info)

app.multiselect(key='fp_comparison').set_value(['6685_17'])
app.slider(key='fp_cutoff').set_value(50.0)
app.checkbox(key='fp_all_matrix').check()
check('footprint groups, ASA threshold and full Jaccard matrix')
selector = app.selectbox(key='actin_ov_selbox')
assert len(selector.options) == 375, len(selector.options)
selector.set_value(3)  # MAFFT column corresponding to P60709 M1
check('first residue, including missing measurements')
assert any('M1' in str(m.value) for m in app.markdown)
app.selectbox(key='actin_ov_selbox').set_value(380)
check('last residue F375')
assert any('F375' in str(m.value) for m in app.markdown)
app.radio(key='explo_mode').set_value('ABP pair')
check('comparison of two ABPs')
app.selectbox(key='homo_site_link').set_value('6685_3')
app.button(key='homo_site_open').click()
check('mixed binding-site navigation')
assert app.selectbox(key='sel_s1').value == '6685_3'
print('All representative paths passed. Browser/3D and scientific validation remain separate.')
