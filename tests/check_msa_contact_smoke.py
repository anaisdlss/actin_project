"""Offline Streamlit check of hetero, mixed and homo legacy contact panels."""
import os
import sys
from pathlib import Path
from streamlit.testing.v1 import AppTest

ROOT=Path(__file__).resolve().parents[1]
os.chdir(ROOT)
sys.path.insert(0,str(ROOT/'script'))
APP='''
import streamlit as st
from msa_analysis import _s1_get_ch_maps, _msa_contact_analysis
site=st.selectbox('Smoke site',['6685_0','6685_3','6685_2'],key='smoke_site')
chains,titles=_s1_get_ch_maps(site)
_msa_contact_analysis(None,'s1c_'+site,_ch2seq=chains,_ch2title=titles,tabs=['A','B','C','D','E'])
'''
app=AppTest.from_string(APP,default_timeout=180)
for site in ['6685_0','6685_3','6685_2']:
    if site!='6685_0':app.selectbox(key='smoke_site').set_value(site)
    app.run()
    assert not app.exception,(site,[e.message for e in app.exception])
    assert any('No area is imputed' in e.value for e in app.caption)
    if site in ('6685_3','6685_2'):
        assert any('Partner denotes actin' in e.value for e in app.caption)
    print('PASS: exact MSA contact scope and all contact tabs for '+site,flush=True)
print('Legacy MSA contact smoke passed.')
standalone=AppTest.from_string('''
from msa_analysis import _msa_contact_analysis
_msa_contact_analysis(lambda title: 'cofilin' in title.lower(), 'cofilin_smoke', tabs=['C','D','E'])
''',default_timeout=180).run()
assert not standalone.exception,[e.message for e in standalone.exception]
assert any('No area is imputed' in e.value for e in standalone.caption)
print('PASS: standalone Cofilin contact analysis without site maps.')
