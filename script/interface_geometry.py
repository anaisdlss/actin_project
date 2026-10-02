"""One-subunit superposition with independent neighbour displacement measurement."""
import numpy as np
from Bio.SVDSuperimposer import SVDSuperimposer


def compare_pair(anchor, neighbor, ref_anchor, ref_neighbor, minimum=300):
    """Fit anchor only. Try both reference orientations; never refit the neighbour."""
    best=None
    for flipped, (a, b) in enumerate(((ref_anchor, ref_neighbor), (ref_neighbor, ref_anchor))):
        ap=sorted(set(anchor)&set(a));bp=sorted(set(neighbor)&set(b))
        if len(ap)<minimum or len(bp)<minimum:continue
        fit=SVDSuperimposer();fit.set(np.array([a[p] for p in ap]),np.array([anchor[p] for p in ap]));fit.run()
        rotation,translation=fit.get_rotran()
        expected=np.array([b[p] for p in bp]);transformed=np.array([neighbor[p] for p in bp])@rotation+translation
        rms=float(np.sqrt(np.mean(np.sum((transformed-expected)**2,axis=1))))
        result=dict(anchor_CA_RMSD_A=float(fit.get_rms()),neighbor_CA_RMSD_A=rms,
                    n_anchor_CA=len(ap),n_neighbor_CA=len(bp),reference_reversed=bool(flipped))
        if best is None or rms<best['neighbor_CA_RMSD_A']:best=result
    return best
