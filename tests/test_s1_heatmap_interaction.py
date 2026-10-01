"""Display normalization, visible/exported scale, and fresh-event selection."""
import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from streamlit.testing.v1 import AppTest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "script"))
import s1_heatmaps as view


def sample_data(homo=True, hetero=True):
    return ([3, 6, 7], ["6685_1"] if homo else [],
            np.array([[10., 20., np.nan]]) if homo else np.empty((0, 3)),
            ["6685_0"] if hetero else [],
            np.array([[40., 80., 0.]]) if hetero else np.empty((0, 3)))


def chart_event(app, point):
    """Dispatch the frontend's Plotly JSON widget state to the real callback."""
    states = app._tree.get_widget_states()
    chart_id = app.get("plotly_chart")[0].proto.id
    target = next((w for w in states.widgets if w.id == chart_id), None)
    if target is None:
        target = states.widgets.add(id=chart_id)
    target.string_value = json.dumps({"selection": {
        "points": [point], "point_indices": [1], "box": [], "lasso": []}})
    app._run(states)
    assert not app.exception, [e.message for e in app.exception]


class S1ScaleTest(unittest.TestCase):
    def render(self, relative=False, homo=True, hetero=True):
        with patch.object(view.st, "markdown") as markdown, \
             patch.object(view.st, "caption"), \
             patch.object(view.st, "plotly_chart") as chart:
            view._render_s1_global_plotly(sample_data(homo, hetero), relative)
            return chart.call_args.args[0], markdown.call_args.args[0]

    def test_absolute_scale_fixed_percent_and_values_unchanged(self):
        fig, legend = self.render()
        self.assertIn("linear-gradient", legend)
        self.assertIn("rgb(255,255,204) 0%", legend)
        self.assertIn("rgb(128,0,38) 100%", legend)
        self.assertIn("<span>100%</span>", legend)
        for trace, expected in zip(fig.data[:2], ([10., 20., np.nan], [40., 80., 0.])):
            self.assertEqual((trace.zmin, trace.zmax), (0, 100))
            np.testing.assert_allclose(trace.z[0], expected, equal_nan=True)
            np.testing.assert_allclose(trace.customdata[0], expected, equal_nan=True)
            self.assertEqual(list(trace.x), [1, 2, 3])
        np.testing.assert_array_equal(fig.data[2].y, [2, 2, 0])

    def test_relative_scale_preserves_absolute_hover_values(self):
        fig, legend = self.render(relative=True)
        self.assertIn("<span>1</span>", legend)
        self.assertNotIn("<span>100%</span>", legend)
        for trace, expected in zip(fig.data[:2], ([.5, 1., np.nan], [.5, 1., 0.])):
            self.assertEqual((trace.zmin, trace.zmax), (0, 1))
            np.testing.assert_allclose(trace.z[0], expected, equal_nan=True)
        np.testing.assert_allclose(fig.data[1].customdata[0], [40., 80., 0.])
        self.assertIn("customdata", fig.data[1].hovertemplate)

    def test_one_native_export_scale_on_first_heatmap_in_all_cases(self):
        for homo, hetero in ((True, True), (True, False), (False, True)):
            with self.subTest(homo=homo, hetero=hetero):
                fig, legend = self.render(homo=homo, hetero=hetero)
                visible = [trace for trace in fig.data if trace.type == "heatmap" and trace.showscale]
                self.assertEqual(len(visible), 1)
                self.assertIs(visible[0], fig.data[0])
                bar = visible[0].colorbar
                self.assertEqual((bar.y, bar.yanchor, bar.lenmode, bar.len), (1, "top", "pixels", 180))
                self.assertIn("linear-gradient", legend)


GLOBAL_APP = '''
import numpy as np
import streamlit as st
import s1_heatmaps as view
data = ([3, 50, 51, 380], ["6685_1"], np.array([[1., 2., 3., 4.]]),
        ["6685_0"], np.array([[4., 3., 2., 1.]]))
view._render_s1_global_plotly(data, False, valid_clusters={"6685_1", "6685_0"})
st.selectbox("Cluster", ["6685_1", "6685_0"], key="sel_s1")
st.checkbox("Unrelated control", key="other")
'''

PATCH_APP = '''
import numpy as np
import pandas as pd
import streamlit as st
import s1_heatmaps as view
positions = [3, 50, 51, 380]
view._render_s1_position_detail((positions, pd.DataFrame(columns=["canon"]),
                                pd.DataFrame(columns=["canon"])), "6685_0")
values = np.array([1., 2., 3., 4.])
view._render_s1_patch_plotly((positions, values, [("c70", 1, "Partner", values)]), "6685_0")
st.checkbox("Unrelated control", key="other")
'''


class S1SelectionTest(unittest.TestCase):
    def test_global_manual_choice_survives_stale_chart_selection(self):
        app = AppTest.from_string(GLOBAL_APP).run()
        self.assertFalse(app.exception)
        chart_event(app, {"x": 46, "y": "6685_0", "curve_number": 1, "point_index": 1})
        self.assertEqual(app.selectbox(key="sel_s1").value, "6685_0")
        app.selectbox(key="sel_s1").set_value("6685_1").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox(key="sel_s1").value, "6685_1")
        app.checkbox(key="other").check().run()
        self.assertEqual(app.selectbox(key="sel_s1").value, "6685_1")
        chart_event(app, {"x": 47, "y": "6685_0", "curve_number": 1, "point_index": 2})
        self.assertEqual(app.selectbox(key="sel_s1").value, "6685_0")
        chart_event(app, {"x": 47, "y": 2, "curve_number": 2, "point_index": 2})
        self.assertEqual(app.selectbox(key="sel_s1").value, "6685_0")

    def test_patch_manual_choice_survives_stale_chart_selection(self):
        app = AppTest.from_string(PATCH_APP).run()
        self.assertFalse(app.exception)
        chart_event(app, {"x": 46, "y": "6685_0", "curve_number": 0, "point_index": 1})
        self.assertEqual(app.selectbox(key="s1posdet_6685_0").value, 50)
        app.selectbox(key="s1posdet_6685_0").set_value(51).run()
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox(key="s1posdet_6685_0").value, 51)
        app.checkbox(key="other").check().run()
        self.assertEqual(app.selectbox(key="s1posdet_6685_0").value, 51)
        chart_event(app, {"x": 375, "y": "6685_0", "curve_number": 0, "point_index": 374})
        self.assertEqual(app.selectbox(key="s1posdet_6685_0").value, 380)
        # A valid reference position absent from this detail's options is ignored.
        chart_event(app, {"x": 2, "y": "6685_0", "curve_number": 0, "point_index": 1})
        self.assertEqual(app.selectbox(key="s1posdet_6685_0").value, 380)

    def test_invalid_position_events_do_not_create_pending_selection(self):
        for x in (None, "unknown", 1.5, float("nan"), float("inf"), 0, 376):
            with self.subTest(x=x), patch.object(view.st, "session_state", {
                    "s1prof_6685_0": {"selection": {"points": [{"x": x}]}}}) as state:
                view._on_s1_patch_selection("6685_0")
                self.assertNotIn("_s1_click_6685_0", state)


if __name__ == "__main__":
    unittest.main()
