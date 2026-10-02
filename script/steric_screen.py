"""Rigid-placement proximity screening; not a binding/competition prediction."""
import numpy as np
from scipy.spatial import cKDTree
from Bio.SVDSuperimposer import SVDSuperimposer


def rigid_proximity(anchor, reference, partner_atoms, partner_residues, neighbor_atoms, minimum=300):
    positions=sorted(set(anchor)&set(reference))
    if len(positions)<minimum:raise ValueError('Fewer than 300 unambiguously aligned actin C-alpha positions.')
    if not len(partner_atoms) or not len(neighbor_atoms):raise ValueError('Missing protein heavy atoms.')
    fit=SVDSuperimposer();fit.set(np.array([reference[p] for p in positions]),np.array([anchor[p] for p in positions]));fit.run()
    rotation,translation=fit.get_rotran()
    transformed=np.asarray(partner_atoms)@rotation+translation
    distances=cKDTree(np.asarray(neighbor_atoms)).query(transformed)[0]
    result=dict(anchor_CA_count=len(positions),anchor_CA_RMSD_A=float(fit.get_rms()),
                minimum_neighbor_distance_A=float(distances.min()),ABP_heavy_atoms=len(partner_atoms),
                ABP_resolved_residues=len(set(partner_residues)))
    for cutoff in [1.5,2.,2.5]:
        mask=distances<cutoff;label=str(cutoff).replace('.','p')
        result[f'ABP_atoms_below_{label}_A']=int(mask.sum())
        result[f'ABP_residues_below_{label}_A']=len({p for p,m in zip(partner_residues,mask) if m})
    return result


def render_steric_screen():
    from pathlib import Path
    import pandas as pd
    import streamlit as st
    root=Path('reports/scientific_audit');path=root/'rigid_filament_proximity.csv'
    if not path.exists():return
    with st.expander('3D screening: ABP proximity to neighbouring actins'):
        st.caption('Each observed ABP–actin pair is placed on a central actin of 3J8A (tropomyosin filament) or 5YU8 (cofilactin) by fitting actin C-alpha atoms. Distances are measured between ABP heavy atoms and the other actins; the fitted actin and reference ABPs are excluded.')
        frame=pd.read_csv(path)
        abp=st.selectbox('ABP for 3D proximity screening',sorted(frame.ABP.unique(),key=str.casefold),key='steric_abp')
        selected=frame[frame.ABP.eq(abp)]
        cols=['site','reference_PDB','reference_anchor','state','anchor_CA_RMSD_A','minimum_neighbor_distance_A',
              'ABP_residues_below_1p5_A','ABP_residues_below_2p0_A','ABP_residues_below_2p5_A','source_PDB','source_actin','source_ABP']
        st.dataframe(selected[[c for c in cols if c in selected]],hide_index=True,width='stretch')
        st.caption('Counts are geometric warning signals at three explicit distance cutoffs, not validated clash scores. No flexibility, side-chain optimization, atom-specific radii or complete infinite filament is modeled. Missing atoms can hide overlaps. A low count does not prove compatibility; a high count needs structural inspection. These distances are independent of the ASA threshold above.')
        st.download_button('Download all rigid-placement measurements',frame.to_csv(index=False).encode(),
                           file_name='rigid_filament_proximity.csv',key='steric_csv')
        st.download_button('Download proximity-screen methods',(root/'rigid_filament_proximity_manifest.json').read_bytes(),
                           file_name='rigid_filament_proximity_methods.json',key='steric_manifest')
