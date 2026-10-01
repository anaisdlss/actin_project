"""Read raw PPI3D snapshots and compare ≥3/≥5 connected-actin screens."""
import sys
import json
import hashlib
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'script'))
from threshold_comparison import compare_thresholds
if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    paths=[ROOT/'data/raw/pdb_entry_results.csv', ROOT/'data/raw/ppi3d_actin_summary.csv']
    results,rows=compare_thresholds(*(pd.read_csv(p,sep=';',low_memory=False) for p in paths))
    report=dict(source_repository=ROOT.name,
                source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                method='Same actin title annotation (Expect value = 0) and graph connectivity as current pipeline; PDB 4b1z excluded. Counts precede the domain/fragment and cluster-55649 filters. Partners are distinct source names directly contacting an annotated actin, not curated ABP identities.',
                thresholds=results,structures=rows)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(results,ensure_ascii=False,indent=2))
