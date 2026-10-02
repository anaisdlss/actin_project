"""Local calculation and source provenance controls in Documentation."""
from collections import deque
from pathlib import Path
import json
import subprocess
import sys
import pandas as pd
import streamlit as st


def render():
    root = Path(__file__).resolve().parents[1]
    report = root/'reports/local_rebuild'
    st.subheader('Rebuild scientific results from local sources')
    st.caption('Recalculate RSA, source-checked scientific tables, interface geometry, filament proximity, '
               'S1 figures and local FoldDisco controls. No remote job is submitted. '
               'Unchanged calculations are skipped after checking their inputs, code and outputs.')
    if (report/'last_run.json').exists():
        last = json.loads((report/'last_run.json').read_text())
        if last['state'] == 'failed':
            st.warning('Last rebuild stopped: '+last.get('error','See calculation log.')+
                       ' Run it again to resume; completed calculations are checked before reuse.')
        elif last['state'] == 'running':
            st.info('A rebuild was started. If it was interrupted, run it again to resume.')
    if (root/'data/.slim_deploy').exists():
        st.caption('Calculations run in the full local research project. This build displays its saved results.')
    elif st.button('Rebuild local scientific results', key='rebuild_local_science'):
        lines = deque(maxlen=15)
        with st.status('Checking and rebuilding local calculations…', expanded=True) as banner:
            output = st.empty()
            process = subprocess.Popen([sys.executable, '-u', 'tools/rebuild_local.py'], cwd=root,
                                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            for line in process.stdout:
                lines.append(line.rstrip())
                output.code('\n'.join(lines), language=None)
            if process.wait():
                banner.update(label='Rebuild stopped — see the reported error', state='error')
            else:
                st.cache_data.clear()
                banner.update(label='Registered local calculations rebuilt and verified', state='complete')
    with st.expander('Data origins and calculation receipts'):
        st.write('A CSV is a storage format. Its origin can be an external database, an external model, '
                 'a local calculation, or an undocumented historical export. A checksum identifies a file; '
                 'it does not establish the scientific validity or original source of its contents.')
        st.caption('The old conservation_vs_asa_per_position.csv has no recovered RSA recipe and is no longer '
                   'used for current RSA values. Variant snapshots retain their documented import origin; '
                   'their original database release dates remain unknown. ProteoCast scores remain external model outputs.')
        path = report/'csv_inventory.csv'
        if path.exists():
            st.caption('Inventory from the last local rebuild; run it again to check newly modified files.')
            data = pd.read_csv(path)
            unresolved = data.status.eq('unresolved historical provenance')
            st.write(f'{len(data)} CSV files inventoried; {int(unresolved.sum())} still need their historical provenance resolved.')
            st.dataframe(data[['file','status','producer','note']], hide_index=True, width='stretch')
            st.download_button('Download complete provenance inventory', path.read_bytes(),
                               file_name=path.name, mime='text/csv', key='provenance_inventory')
        else:
            st.info('Run the local rebuild to generate the provenance inventory.')
        st.caption('Detailed successful calculation receipts and failure logs are saved under reports/local_rebuild. '
                   'Imported source snapshots are kept intact. An unresolved historical file is not silently certified.')
