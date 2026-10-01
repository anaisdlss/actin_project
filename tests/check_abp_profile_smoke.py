"""Offline Streamlit smoke for bundled ABP profiles; run from a repo root."""
import json
import os
import sys
from pathlib import Path

from streamlit.testing.v1 import AppTest

ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)
sys.path.insert(0, str(ROOT / "script"))
from proteocast_results import result_file


APP = '''
import streamlit as st
from unittest.mock import patch
from proteocast_abp import _render_abp_proteocast
name = st.selectbox('Smoke ABP', ['Cofilin-1', 'Filamin-A'], key='smoke_abp')
# Keep the smoke offline; test an interior domain position, not an API response.
domains = [dict(name='Smoke domain', db='test', spans=[(2, 100)])]
with patch('proteocast_view.fetch_domains', return_value=domains):
    _render_abp_proteocast(name, [])
'''


def run(app, label):
    app.run()
    assert not app.exception, (label, [e.message for e in app.exception])
    print(f"PASS: {label}", flush=True)


app = AppTest.from_string(APP, default_timeout=180)
run(app, "ABP profile view starts")
if result_file(ROOT / "data/proteocast/abp", "Cofilin_1") is None:
    assert any("result unavailable" in element.value.lower() for element in list(app.info) + list(app.error))
    assert not app.get("plotly_chart")
    print("PASS: missing results retain the explicit unavailable message.")
else:
    for name, length in (("Cofilin-1", 166), ("Filamin-A", 2647)):
        app.selectbox(key="smoke_abp").set_value(name)
        run(app, f"{name} heatmap and continuous profile")
        charts = app.get("plotly_chart")
        assert len(charts) == 2, (name, len(charts))
        chart = json.loads(charts[1].proto.spec)
        assert len(chart["data"][0]["x"]) == length
        assert chart["data"][0]["x"] == list(range(1, length + 1))
        assert chart["layout"]["xaxis3"]["range"] == [0.5, length + 0.5]
        assert chart["layout"]["hovermode"] == "x unified"
        assert any(value is not None for value in chart["data"][1]["y"]), name
        assert any("Download ABP sensitivity" in button.label for button in app.get("download_button"))
    app.selectbox(key="smoke_abp").set_value("Cofilin-1")
    run(app, "return to Cofilin-1")
    app.toggle(key="pc_surface_Cofilin_1").set_value(True)
    run(app, "Cofilin-1 exposed-surface filter")
    print("All offline ABP profile smoke checks passed.")
