"""Exercise actual Streamlit selection callbacks and later manual changes."""
import json
import sys
import unittest
from pathlib import Path

from streamlit.testing.v1 import AppTest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "script"))

APP = '''
import streamlit as st
import pandas as pd
from unittest.mock import patch
import residue_explorer as view
pp = {
    "pos": pd.DataFrame({"canon": [3, 50, 51, 380],
                         "actin_aa": ["M", "G", "Q", "F"],
                         "n_abp": [0, 8, 4, 1]}),
    "res_abp": pd.DataFrame({"canon": [3, 50, 51, 380],
                             "abp": ["Example"] * 4,
                             "asa_max": [1., 8., 4., 2.]}),
}
with patch.object(view, "render_residue_fiche", lambda pp, p: st.text(f"detail:{p}")), \
     patch.object(view, "_render_actin_overview_3d", lambda pp, p, mode: st.text(f"surface:{p}")), \
     patch.object(view, "_render_actin_3d", lambda pp, p: st.text(f"surface:{p}")):
    RENDER(pp)
st.checkbox("Unrelated control", key="other")
'''


def select_chart_point(app, x):
    # AppTest does not expose a Plotly interaction wrapper. Send the same JSON
    # widget state as the frontend, so Streamlit invokes the real callback.
    chart_id = app.get("plotly_chart")[0].proto.id
    states = app._tree.get_widget_states()
    target = next((w for w in states.widgets if w.id == chart_id), None)
    if target is None:
        target = states.widgets.add(id=chart_id)
    target.string_value = json.dumps({"selection": {
        "points": [{"x": x, "curve_number": 0, "point_index": 1}],
        "point_indices": [1], "box": [], "lasso": [],
    }})
    app._run(states)
    assert not app.exception, [e.message for e in app.exception]


class ResidueSelectionTest(unittest.TestCase):
    def exercise(self, renderer, selector, clicked_position):
        app = AppTest.from_string(APP.replace("RENDER", renderer)).run()
        self.assertFalse(app.exception)
        select_chart_point(app, clicked_position)
        self.assertEqual(app.selectbox(key=selector).value, 50)  # G46
        app.selectbox(key=selector).set_value(51).run()
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox(key=selector).value, 51)
        self.assertIn("detail:51", [e.value for e in app.text])
        self.assertIn("surface:51", [e.value for e in app.text])
        app.checkbox(key="other").check().run()
        self.assertEqual(app.selectbox(key=selector).value, 51)
        # A subsequent explicit chart selection must still work.
        select_chart_point(app, 375 if isinstance(clicked_position, int) else "375")
        self.assertEqual(app.selectbox(key=selector).value, 380)
        app.selectbox(key=selector).set_value(3).run()
        self.assertEqual(app.selectbox(key=selector).value, 3)

    def test_overview_manual_choice_after_chart_click(self):
        self.exercise("view.render_actin_overview", "actin_ov_selbox", 46)

    def test_heatmap_manual_choice_after_chart_click(self):
        self.exercise("view.render_residue_tab", "explo_res_selbox", "46")


if __name__ == "__main__":
    unittest.main()
