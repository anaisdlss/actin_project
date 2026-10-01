"""Regenerate S1 figures from the same exact-chain calculations used by the app.
Run from a repository root. Preserves overwritten files in a dated ZIP under reports/.
"""
import argparse,datetime,hashlib,json,os,sys,zipfile
from pathlib import Path
import streamlit  # before script/ on sys.path (streamlit.py shadows package)
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'script'))
from s1_heatmaps import _S1_GLOBAL_FILES,_build_s1_global_heatmap,_build_s1_patch_detail
import numbering


def run(scope='all'):
    os.chdir(ROOT)
    vis=ROOT/'data/visualisations';vis.mkdir(parents=True,exist_ok=True)
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    backups=ROOT/'reports/figure_backups';backups.mkdir(parents=True,exist_ok=True)
    old=[]
    if scope in ('all','global'):old += list(vis.glob('actin_s1_*.png'))
    if scope in ('all','clusters'):
        old += list((vis/'actin_s1_clusters').glob('*.png'))+list((vis/'actin_s1_clusters_by_c70').glob('*.png'))
        old += [ROOT/'data/filtered/actin_s1_canon_area_by_cluster.csv']
    with zipfile.ZipFile(backups/f'{stamp}.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in old:
            if p.exists():z.write(p,p.relative_to(ROOT))
    mtimes=tuple(Path(p).stat().st_mtime_ns for p in _S1_GLOBAL_FILES)
    pos,hl,hm,el,em=_build_s1_global_heatmap(mtimes)
    labels=[f'Homo: {x}' for x in hl]+[f'ABP: {x}' for x in el]
    matrix=np.vstack([hm,em]);outputs=[]
    def draw(mat,names,title,path,relative=False,maximum=None,value_label=None):
        fig,ax=plt.subplots(figsize=(18,max(2.3,.22*len(names)+1.6)))
        im=ax.imshow(mat,aspect='auto',cmap='Blues',vmin=0,vmax=maximum if maximum is not None else 1 if relative else 100,interpolation='none')
        ticks=[0]+list(range(24,375,25))
        ax.set_xticks(ticks);ax.set_xticklabels([str(i+1) for i in ticks],fontsize=8)
        ax.set_yticks(range(len(names)));ax.set_yticklabels(names,fontsize=7)
        ax.set_xlabel(numbering.AXIS_TITLE);ax.set_title(title,fontsize=11)
        fig.colorbar(im,ax=ax,label=value_label or ('Fraction of row maximum' if relative else 'Mean buried ASA (%)'),pad=.015)
        fig.tight_layout();path.parent.mkdir(parents=True,exist_ok=True)
        fig.savefig(path,dpi=120);plt.close(fig);outputs.append(path)
    if scope in ('all','global'):
        draw(matrix,labels,'S1 interface profiles — equal weight per C70',vis/'actin_s1_heatmap_absolute.png')
        maximum=matrix.max(axis=1,keepdims=True)
        relative=np.divide(matrix,maximum,out=np.zeros_like(matrix),where=maximum>0)
        draw(relative,labels,'S1 profiles — relative to each row maximum',vis/'actin_s1_all_equitable_heatmap.png',True)
        reference=[hl.index(p) for p in ['6685_1','6685_2','6685_3','6685_4'] if p in hl]
        for name,mat in [('homo',hm),('hetero',em)]:
            draw(np.any(mat>0,axis=0).astype(float).reshape(1,-1),[name],f'S1 {name}: observed positive buried ASA',vis/f'actin_s1_{name}_used_heatmap.png',maximum=1,value_label='Observed (0/1)')
        if reference:
            draw(hm[reference].mean(axis=0).reshape(1,-1),['6685_1–4'],'Reference sites 6685_1–4: equal mean of C70-weighted S1 profiles',vis/'actin_s1_homo_top4_mean_heatmap.png')
            draw(np.any(hm[reference]>0,axis=0).astype(float).reshape(1,-1),['6685_1–4'],'Reference sites 6685_1–4: observed positive buried ASA',vis/'actin_s1_top4homo_used_heatmap.png',maximum=1,value_label='Observed (0/1)')
        counts=(em>=.01).sum(axis=0).reshape(1,-1)
        draw(counts,['ABP sites'],'Number of S1 ABP sites with C70-weighted mean buried ASA ≥ 0.01%',vis/'actin_s1_hetero_count_heatmap.png',maximum=max(1,counts.max()),value_label='Number of binding sites')
        for name,mat in [('homo',hm),('abp',em)]:
            p=vis/f'actin_s1_{name}_p60709.csv';pd.DataFrame(mat,index=hl if name=='homo' else el,columns=range(1,376)).to_csv(p);outputs.append(p)
    if scope in ('all','clusters'):
        profiles=[];patches=sorted(set(hl+el),key=lambda s:int(s.split('_')[-1]))
        for i,patch in enumerate(patches):
            detail=_build_s1_patch_detail(patch,mtimes)
            if detail is None:raise ValueError(f'Missing profile {patch}')
            positions,profile,rows=detail
            if positions!=pos:raise ValueError('Different profile axes')
            profiles.append(profile)
            draw(profile.reshape(1,-1),[patch],f'{patch} — S1 profile, equal weight per C70',vis/'actin_s1_clusters'/f'{patch}.png')
            draw(np.array([r[3] for r in rows]),[f'{r[0]} (n={r[1]})' for r in rows],f'{patch} — C70 decomposition (all C70, including singletons)',vis/'actin_s1_clusters_by_c70'/f'{patch}.png')
            if i%25==0:print(f'{i+1}/{len(patches)} sites',flush=True)
        p=ROOT/'data/filtered/actin_s1_canon_area_by_cluster.csv'
        pd.DataFrame(profiles,index=pd.Index(patches,name='patch'),columns=positions).to_csv(p);outputs.append(p)
        # Remove obsolete S1 images only after backing them up and completing replacements.
        for p in old:
            if p.suffix=='.png' and p not in outputs:p.unlink()
    report={'generated_utc':stamp,'scope':scope,'reference':'P60709 1–375','S1_scope':'subunit_1 only, as in the app; both-side union is a separate footprint analysis',
            'aggregation':'Within C70: sum of per-interaction residue max ASA / total interactions; then equal-weight mean across C70. Absent interface residues contribute zero, not measured solvent accessibility.',
            'backup':str((backups/f'{stamp}.zip').relative_to(ROOT)),
            'sources':{p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in _S1_GLOBAL_FILES[1:]},
            'code':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'script/s1_heatmaps.py',ROOT/'script/numbering.py']},
            'outputs':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs}}
    (ROOT/'reports'/f's1_figures_{scope}.json').write_text(json.dumps(report,indent=2))
    print(f'Regenerated {len(outputs)} files',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--scope',choices=['all','global','clusters'],default='all')
    run(parser.parse_args().scope)
