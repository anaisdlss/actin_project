"""Consistent position cursors for residue profiles and heatmaps."""


def position_hover(fig, axes=None, unified=True):
    """Share hover across vertically stacked position axes, keeping data intact.

    Matching axes created by make_subplots are not the same axis for Plotly hover.
    Reuse the bottom axis for traces in each column, preserving its title/range.
    Pass axes=['x'] for a figure with unrelated distributions in another column.
    Heatmaps keep their per-cell values; unified=False also suits sparse marks
    or domain intervals whose nearest endpoints are not simultaneous measures.
    """
    refs = list(dict.fromkeys(
        getattr(trace, 'xaxis', None) or 'x' for trace in fig.data
        if hasattr(trace, 'xaxis')))
    if axes is not None:
        refs = [ref for ref in refs if ref in axes]
    groups = {}
    for ref in refs:
        axis = fig.layout['xaxis' + ref[1:]]
        groups.setdefault(tuple(axis.domain or (0, 1)), []).append(ref)
    for group in groups.values():
        def bottom(ref):
            anchor = fig.layout['xaxis' + ref[1:]].anchor or 'y'
            yaxis = fig.layout['yaxis' + anchor[1:]]
            return (yaxis.domain or (0, 1))[0]
        shared = min(group, key=bottom)
        for trace in fig.data:
            if hasattr(trace, 'xaxis') and (trace.xaxis or 'x') in group:
                trace.xaxis = shared
                yref = getattr(trace, 'yaxis', None) or 'y'
                fig.layout['yaxis' + yref[1:]].anchor = shared
        fig.layout['xaxis' + shared[1:]].update(
            matches=None, showspikes=True, spikemode='across', spikesnap='data',
            spikedash='dot', spikecolor='#666666', spikethickness=1)
    fig.update_layout(hovermode='x unified' if unified else 'closest',
                      hoversubplots='axis', spikedistance=-1,
                      hoverlabel=dict(namelength=-1))
    return fig
