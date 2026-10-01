"""Sensitivity of the connected-actin screen, without changing the dataset."""
import networkx as nx
import pandas as pd


def normalise_title(value):
    return '' if pd.isna(value) else ' '.join(str(value).strip().lower().split())


def compare_thresholds(entries, summary, thresholds=(3, 5)):
    allowed = ['homo', 'hetero']
    summary = summary[summary['Interface type'].isin(allowed)]
    actin_titles = set(summary.loc[pd.to_numeric(summary['Expect value'], errors='coerce').eq(0),
                                  'Result protein'].dropna().map(normalise_title))
    entries = entries[entries['Interface type'].isin(allowed)].copy()
    entries['pdb_id'] = entries['pdb_id'].astype(str).str.lower()
    rows = []
    for pdb, group in entries.groupby('pdb_id'):
        if pdb == '4b1z':  # same permanent exclusion as current pipeline
            continue
        graph = nx.Graph()
        titles = {}
        for _, row in group.iterrows():
            a,b = row['Interactor 1'],row['Interactor 2']
            graph.add_edge(a,b)
            titles[a],titles[b] = row['Interactor 1 title'],row['Interactor 2 title']
        actin = {node for node,title in titles.items() if normalise_title(title) in actin_titles}
        largest = max((len(c) for c in nx.connected_components(graph.subgraph(actin))),default=0)
        # Directly contacting proteins only; no non-actin/non-actin complex subunits.
        partners = sorted({str(titles[n]) for a in actin for n in graph.neighbors(a)
                           if n not in actin and pd.notna(titles[n])})
        rows.append(dict(pdb=pdb,largest_actin_component=largest,partner_names=partners))
    results=[]
    baseline = max(thresholds)
    base_rows=[r for r in rows if r['largest_actin_component']>=baseline]
    base_names={n for r in base_rows for n in r['partner_names']}
    for threshold in thresholds:
        selected=[r for r in rows if r['largest_actin_component']>=threshold]
        names={n for r in selected for n in r['partner_names']}
        results.append(dict(threshold=threshold,pdb_count=len(selected),partner_name_count=len(names),
                            additional_pdbs=[r['pdb'] for r in selected if r['largest_actin_component']<baseline],
                            additional_partner_names=sorted(names-base_names)))
    return results,rows


def render_threshold_comparison(read_csv):
    from pathlib import Path
    import streamlit as st
    entries = Path('data/raw/pdb_entry_results.csv')
    summary = Path('data/raw/ppi3d_actin_summary.csv')
    if not entries.exists() or not summary.exists():
        return
    with st.expander('Threshold sensitivity: 3 versus 5 connected actins'):
        st.caption('Exploratory comparison of the raw snapshot using the existing actin-name '
                   'annotation and connected-component rule. PDB 4b1z remains excluded. '
                   'Counts precede the domain/fragment and cluster-55649 filters. This does '
                   'not change the dataset used by the other sections. Partner counts are '
                   'distinct source names in direct contact with actin, not curated ABP identities.')
        results, rows = compare_thresholds(read_csv(str(entries)),read_csv(str(summary)))
        table = pd.DataFrame([{'Minimum connected actins':r['threshold'], 'PDBs':r['pdb_count'],
                               'Direct partner names':r['partner_name_count'],
                               'Additional PDBs vs 5':len(r['additional_pdbs']),
                               'Additional partner names vs 5':len(r['additional_partner_names'])}
                              for r in results])
        st.dataframe(table,hide_index=True,width='stretch')
        st.markdown('**Additional partner names at threshold 3**')
        st.dataframe(pd.DataFrame({'Source protein name':results[0]['additional_partner_names']}),
                     hide_index=True,width='stretch')
        details = pd.DataFrame(rows)
        details['partner_names'] = details['partner_names'].map(lambda x:' ; '.join(x))
        st.download_button('Download threshold comparison by PDB (CSV)',
                           details.to_csv(index=False).encode(),file_name='actin_threshold_comparison.csv',
                           mime='text/csv',key='threshold_comparison_download')
