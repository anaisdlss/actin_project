"""Observed residue interaction types on a sequence-mapped actin reference."""
from io import StringIO
from pathlib import Path

import pandas as pd
import streamlit as st

CATEGORIES = {
    'No observed contact': '#D9D9D9',
    'Homo — actin only': '#E69F00',
    'Hetero — ABP only': '#0072B2',
    'Mixed — both': '#884EA0',
}


def contact_categories(records):
    """Use positive measured buried ASA, both actin sides; absence is observational."""
    retained = records[pd.to_numeric(records.asa, errors='coerce').gt(0)]
    homo = set(retained.loc[retained.kind.eq('homo'), 'position'].dropna().astype(int))
    hetero = set(retained.loc[retained.kind.eq('abp'), 'position'].dropna().astype(int))
    return {p: ('Mixed — both' if p in homo & hetero else
                'Homo — actin only' if p in homo else
                'Hetero — ABP only' if p in hetero else 'No observed contact')
            for p in range(1, 376)}


@st.cache_data(show_spinner=False)
def mapped_reference_atoms(text, reference):
    """Map atom serials by unambiguous sequence alignment, never by PDB number alone."""
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import protein_letters_3to1_extended
    from scientific_analysis import isoform_map
    structure = PDBParser(QUIET=True).get_structure('actin', StringIO(text))
    atoms = {}
    for chain in next(iter(structure)):
        residues = [r for r in chain if 'CA' in r and r.resname in protein_letters_3to1_extended]
        if len(residues) < 300:
            continue
        sequence = ''.join(protein_letters_3to1_extended[r.resname] for r in residues)
        mapping = isoform_map(sequence, reference)
        same = sum(sequence[i-1] == reference[p-1] for i, p in mapping.items())
        if len(mapping) < 300 or same / len(mapping) < .9:
            continue
        for i, p in mapping.items():
            atoms.setdefault(p, []).extend(a.serial_number for a in residues[i-1])
    return atoms


def atom_colors(text, mapped, categories):
    """Color one continuous surface; unmapped atoms retain the explicit gray key."""
    colors = {}
    for line in text.splitlines():
        if line.startswith(('ATOM  ', 'HETATM')):
            colors[int(line[6:11])] = CATEGORIES['No observed contact']
    for position, serials in mapped.items():
        color = CATEGORIES[categories.get(position, 'No observed contact')]
        colors.update({serial: color for serial in serials})
    return colors


def render_contact_surface(text, selected):
    if text is None:
        st.info('The actin reference structure is unavailable.')
        return
    import py3Dmol
    from footprint_comparison import FILES, load_footprints
    from scientific_analysis import sequence
    if not all(p.exists() for p in FILES):
        st.info('Interaction source tables are unavailable.')
        return
    categories = contact_categories(load_footprints(tuple(p.stat().st_mtime_ns for p in FILES)))
    mapped = mapped_reference_atoms(text, sequence(Path('data/P60709_ref.fasta')))
    if not mapped:
        st.warning('The reference structure cannot be mapped reliably to P60709.')
        return
    viewer = py3Dmol.view(width='100%', height=460)
    viewer.addModel(text, 'pdb')
    viewer.setStyle({}, {'cartoon': {'color': '#D9D9D9'}})
    viewer.addSurface(py3Dmol.SES, {'opacity': 1, 'colorscheme': {
        'prop': 'serial', 'map': atom_colors(text, mapped, categories)}})
    if selected in mapped:
        viewer.addStyle({'serial': mapped[selected]}, {'stick': {'color': '#111111', 'radius': .3}})
    viewer.zoomTo(); viewer.zoom(.82); viewer.setBackgroundColor('white')
    st.components.v1.html(viewer._make_html(), height=470, scrolling=False)
    st.markdown(' &nbsp; '.join(f'<span style="color:{color}">■</span> {label}' for label, color in CATEGORIES.items()), unsafe_allow_html=True)
    st.caption(f"Selected position: {selected} · {categories.get(selected, 'unmapped')}. Black sticks mark the selected residue. "
               'Categories combine observed positive buried-ASA contacts across the dataset; gray does not mean that binding is impossible. Unmapped atoms are also gray.')
