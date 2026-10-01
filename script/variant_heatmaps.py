"""Separate binary variant presence from quantitative substitution counts."""
import numpy as np
import pandas as pd
import plotly.graph_objects as go

AMINO_ACIDS = list('ACDEFGHIKLMNPQRSTVWY')
ABSENT_COLOR = '#F2F3F5'
PRESENT_COLOR = '#0072B2'


def presence_matrix(variants):
    """One selected gene/category: a substitution cell is present or absent."""
    matrix=pd.crosstab(variants.aa_alt,variants.position).reindex(
        index=AMINO_ACIDS,columns=range(1,376),fill_value=0)
    return matrix.gt(0).astype(int)


def presence_figure(matrix,gene,category):
    states=np.where(matrix.to_numpy().astype(bool),'Recorded in selected category','No record in selected category')
    figure=go.Figure(go.Heatmap(
        z=matrix.values,x=matrix.columns,y=matrix.index,zmin=0,zmax=1,
        colorscale=[[0,ABSENT_COLOR],[.4999,ABSENT_COLOR],[.5,PRESENT_COLOR],[1,PRESENT_COLOR]],
        showscale=False,customdata=states,
        hovertemplate='P60709 %{x} → %{y}<br>%{customdata}<extra></extra>'))
    for label,color in [('Not recorded (0)',ABSENT_COLOR),('Recorded (1)',PRESENT_COLOR)]:
        figure.add_trace(go.Scatter(x=[None],y=[None],mode='markers',name=label,
            marker=dict(symbol='square',size=12,color=color,line=dict(color='#7A7A7A',width=1)),
            hoverinfo='skip',showlegend=True))
    figure.update_layout(title=f'{gene} · {category}: individual substitutions',
        xaxis_title='Aligned P60709 position',yaxis_title='Alternate amino acid',height=480,
        margin=dict(l=60,r=20,t=70,b=100),
        legend=dict(orientation='h',x=0,y=-.23,xanchor='left',yanchor='top'))
    return figure


def substitution_count_trace(counts):
    """A continuous colour scale encodes actual integer counts, never severity."""
    maximum=max(1,int(counts.to_numpy().max()))
    ticks=list(range(maximum+1)) if maximum<=10 else sorted(set(np.linspace(0,maximum,6).round().astype(int)))
    return go.Heatmap(z=counts.values,x=counts.columns,y=counts.index,
        colorscale=[[0,ABSENT_COLOR],[.20,'#C6DBEF'],[.50,'#6BAED6'],[1,'#08519C']],
        zmin=0,zmax=maximum,
        colorbar=dict(title=dict(text='Substitutions<br>per position',side='right'),
                      tickmode='array',tickvals=ticks,tickformat='d',
                      len=.46,y=1,yanchor='top',x=1.02,thickness=16),
        hovertemplate='%{y}<br>P60709 %{x}<br>Distinct substitutions: %{z:d}<extra></extra>')
