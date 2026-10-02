"""Shared display conventions without changing scientific values."""
import re
import numpy as np
import pandas as pd


def ordered_interactions(frame):
    """Keep identifiers intact and order numeric IDs numerically, with missing IDs last."""
    if 'interaction_id' not in frame:
        return frame
    result = frame.copy()
    values = pd.to_numeric(result.interaction_id, errors='coerce')
    if (result.interaction_id.isna() | values.notna()).all() and values.dropna().mod(1).eq(0).all():
        result['interaction_id'] = values.astype('Int64')
    return result.iloc[np.argsort(values.astype(float).fillna(np.inf).to_numpy(), kind='stable')].reset_index(drop=True)


def heatmap_cell_colors(trace, layout):
    from plotly.colors import get_colorscale, convert_colors_to_same_type
    z = np.asarray(trace.z, dtype=float)
    if z.ndim != 2 or not np.isfinite(z).any():
        return None
    axis = layout[trace.coloraxis] if trace.coloraxis else None
    scale = (axis.colorscale if axis is not None else trace.colorscale) or get_colorscale('Viridis')
    low = axis.cmin if axis is not None else trace.zmin
    high = axis.cmax if axis is not None else trace.zmax
    low = float(np.nanmin(z)) if low is None else low
    high = float(np.nanmax(z)) if high is None else high
    stops = np.array([item[0] for item in scale], dtype=float)
    rgb = np.array(convert_colors_to_same_type([item[1] for item in scale], colortype='tuple')[0])
    t = np.clip((z-low)/(high-low), 0, 1) if high > low else np.full(z.shape, .5)
    if (axis.reversescale if axis is not None else trace.reversescale):
        t = 1-t
    channels = np.stack([np.interp(t, stops, rgb[:,i]) for i in range(3)], axis=-1)
    channels = np.nan_to_num(channels * 255).round().astype(int)
    return np.array([[f'#{r:02x}{g:02x}{b:02x}' if np.isfinite(z[i,j]) else '#dddddd'
                      for j,(r,g,b) in enumerate(row)] for i,row in enumerate(channels)])


def cell_color_hover(figure):
    """Give heatmap tooltips the actual cell colour, rather than a miniature palette."""
    for trace in figure.data:
        if trace.type != 'heatmap' or trace.hoverinfo in ('skip', 'none'):
            continue
        if trace.hovertemplate and trace.hovertemplate.startswith('%{hovertext} '):
            continue
        colors = heatmap_cell_colors(trace, figure.layout)
        if colors is None:
            continue
        existing = trace.hovertext
        trace.hovertext = [[f'<span style="color:{color};font-size:20px">■</span>'
                            + ('<br>' + str(existing if isinstance(existing, str) else existing[i][j]) if existing is not None else '')
                            for j,color in enumerate(row)] for i,row in enumerate(colors)]
        trace.hovertemplate = '%{hovertext} ' + (trace.hovertemplate or ('%{text}<extra></extra>' if trace.text is not None else '%{y}<br>%{x}: %{z}<extra></extra>'))
    return figure


def plotly_chart(figure, **kwargs):
    import streamlit as st
    config = dict(kwargs.pop('config', {}) or {})
    config.update(displayModeBar=True, displaylogo=False, responsive=True, scrollZoom=True)
    return st.plotly_chart(cell_color_hover(figure), config=config, **kwargs)


def viewer_html(viewer):
    """Refit a 3D view when its own container changes width, including hidden views."""
    html = viewer._make_html()
    match = re.search(r'(viewer_\w+)\.render\(\);', html)
    if match:
        name = match.group(1)
        element = name.replace('viewer_', '3dmolviewer_')
        extra = f'''{name}.render();
        var lastWidth = -1;
        var resizeView = function() {{
          var element = document.getElementById('{element}');
          if (!element || element.clientWidth < 1 || element.clientWidth === lastWidth) return;
          lastWidth = element.clientWidth;
          {name}.resize(); {name}.zoomTo(); {name}.zoom(0.85); {name}.render();
        }};
        new ResizeObserver(resizeView).observe(document.getElementById('{element}'));
        resizeView();'''
        html = html[:match.start()] + extra + html[match.end():]
    return '<style>html,body{margin:0;width:100%;overflow:hidden;}div[id^="3dmolviewer_"]{max-width:100%;}</style>' + html


def color_legend(items):
    import html
    import streamlit as st
    st.markdown(' &nbsp; '.join(f'<span style="color:{color}">■</span> {html.escape(label)}'
                              for label,color in items), unsafe_allow_html=True)


def gradient_legend(label, colors, low, high, meaning):
    import html
    import streamlit as st
    st.markdown(f'<div style="max-width:420px;margin:8px auto"><b>{html.escape(label)}</b>'
                f'<div style="height:12px;background:linear-gradient(to right,{",".join(colors)});"></div>'
                f'<div style="display:flex;justify-content:space-between"><span>{low:.2f}</span><span>{high:.2f}</span></div></div>',
                unsafe_allow_html=True)
    st.caption(meaning)
