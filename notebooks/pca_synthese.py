#!/usr/bin/env python
# coding: utf-8

import os as _os
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # cwd = racine projet (robuste, peu importe d'où on lance)


# # PCA de synthèse — montrer ce que disent les analyses
# 
# Condensé visuel des déterminants de la **conservation** d'une position d'actine.
# Unité = une position canonique (résidu MAFFT).
# 
# **Verdict honnête des tendances (testé sur les données) :**
# 
# | affirmation | verdict | chiffre |
# |---|---|---|
# | enfoui (rsa bas) = plus conservé | ✅ **fort** | rho = **−0.45** (p=2e-20) |
# | buried ASA% moyen ↑ = mieux conservé | ✅ faible | rho = +0.17, +0.07 à rsa égal |
# | + d'ABP différents ≠ + conservé | ✅ confirmé | brut −0.22 → **+0.13 à rsa égal** (ns) |
# | hydrophobe = plus conservé | ❌ **non** | rho = −0.02 (ns) |
# | aromatique = plus conservé | ⚠️ tendance ns | 4.54 vs 4.18, p=0.10 (n=31) |
# 
# → **Le déterminant dominant est l'enfouissement.** Le reste (nb de partenaires, hydrophobicité) n'ajoute rien une fois l'exposition contrôlée.

# In[81]:


import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import svd
from scipy.stats import spearmanr, mannwhitneyu
from pathlib import Path
_R = Path('.') if (Path('.')/'data').exists() else Path('.')

cv  = pd.read_csv(_R/'data/proteocast/conservation_vs_asa_per_position.csv').set_index('canon')
df3 = pd.read_csv(_R/'data/filtered/details/3.interface_residues.csv')
di  = pd.read_csv(_R/'data/filtered/details/1.interactions.csv')
da  = pd.read_csv(_R/'data/filtered/filtered_all_data.csv', low_memory=False)
pp  = pd.read_csv(_R/'data/filtered/proteins_per_pdb.csv'); ac = set(pp[pp.is_actin]['chain'].str.lower())

df3['canon'] = pd.to_numeric(df3['residue_number_canon_mafft'], errors='coerce')
df3['pct']   = pd.to_numeric(df3['buried_ASA_percent'].astype(str).str.replace('%','',regex=False), errors='coerce')
homo = set(di[di['chain_B_id'].str.lower().isin(ac)]['interaction_id'])
d = df3.dropna(subset=['canon']).copy(); d['canon'] = d.canon.astype(int)
d['aa'] = d['residue_name'].astype(str).str.upper().str.strip()

# acide amine majoritaire par position -> hydrophobicite (Kyte-Doolittle) + aromaticite
KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,
      'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
aa_can = d[d.aa.str.len()==1].groupby('canon')['aa'].agg(lambda s: s.mode().iloc[0])
cv['aa']             = aa_can.reindex(cv.index)
cv['hydrophobicite'] = cv['aa'].map(KD)
cv['aromatique']     = cv['aa'].isin(list('FWY')).astype(float)

# descripteurs d'interface cote actine (hetero)
het = d[~d.interaction_id.isin(homo)]
cv['asa_mean'] = het.groupby('canon')['pct'].mean().reindex(cv.index)     # buried ASA% moyen
mm = di.merge(da[['subunit_1','subunit_2','s2_sequence_cluster_40','s2_sequence_cluster_70']].drop_duplicates(),
              left_on=['chain_A_id','chain_B_id'], right_on=['subunit_1','subunit_2'], how='left')
cv['n_fam'] = het.assign(f=het.interaction_id.map(mm.set_index('interaction_id')['s2_sequence_cluster_40'])).groupby('canon')['f'].nunique().reindex(cv.index).fillna(0)
cv['n_abp'] = het.assign(f=het.interaction_id.map(mm.set_index('interaction_id')['s2_sequence_cluster_70'])).groupby('canon')['f'].nunique().reindex(cv.index).fillna(0)

P = cv
print(f'{len(P)} positions | descripteurs prets')
P[['conservation','rsa','asa_mean','combined_asa','homo_asa','hetero_asa','n_fam','n_abp','hydrophobicite','aromatique']].describe().round(2)


# # A — PCA des positions d'actine
# 
# Un point = **une position canonique** de l'actine. Message central : la conservation suit l'**enfouissement**, pas le nombre/type de partenaires.

# ## Figure 1 — Cercle des corrélations (le message en une image)
# 
# PCA sur les variables, conservation incluse. Flèches dans le **même sens** = corrélées ;
# **opposées** = anti-corrélées ; **perpendiculaires** = indépendantes.
# 🔴 conservation · 🔵 exposition (rsa, nb familles, enfoui ABP) · 🟢 enfouissement.

# In[82]:


vlab = {'conservation':'conservation','rsa':'exposition (rsa)','combined_asa':'enfouissement total',
        'homo_asa':'enfoui filament','hetero_asa':'enfoui ABP','n_fam':'nb familles ABP'}
grp  = {'conservation':'cons','rsa':'expo','n_fam':'expo','hetero_asa':'expo','combined_asa':'bury','homo_asa':'bury'}
gcol = {'cons':'#c0392b','expo':'#2980b9','bury':'#27ae60'}
feat = list(vlab); Q = P.dropna(subset=feat).copy()
X = Q[feat].values.astype(float); Xz = (X-X.mean(0))/(X.std(0)+1e-9)
U,S,Vt = svd(Xz, full_matrices=False); ev = (S**2)/np.sum(S**2); load = Vt.T*np.sqrt((S**2)/len(Q))

fig,ax = plt.subplots(figsize=(7.5,7.5)); th=np.linspace(0,2*np.pi,200)
ax.plot(np.cos(th),np.sin(th),color='#ccc',lw=1); ax.axhline(0,color='#eee',lw=.8); ax.axvline(0,color='#eee',lw=.8)
dy={'rsa':-0.10,'n_fam':+0.06,'homo_asa':+0.07,'combined_asa':-0.07}
for i,f in enumerate(feat):
    c=gcol[grp[f]]
    ax.arrow(0,0,load[i,0],load[i,1],head_width=.03,color=c,lw=2.2,length_includes_head=True)
    ax.text(load[i,0]*1.14,load[i,1]*1.14+dy.get(f,0),vlab[f],color=c,fontsize=10.5,ha='center',va='center',fontweight='bold')
ax.set_xlim(-1.3,1.3); ax.set_ylim(-1.3,1.3); ax.set_aspect('equal')
ax.set_xlabel(f'PC1 ({ev[0]*100:.0f}%)'); ax.set_ylabel(f'PC2 ({ev[1]*100:.0f}%)')
ax.set_title("Cercle des correlations\nsensibilité mutationnelle (rouge) s'oppose a l'exposition & au nb de partenaires (bleu),\ns'aligne sur l'enfouissement filament (vert)",fontsize=10.5)
plt.tight_layout(); plt.show()


# ## Figure 2 — PCA des positions, conservation en *overlay*
# 
# PCA sur descripteurs **structuraux seuls** (conservation exclue), points colorés ensuite.
# Le gradient **conservation** est l'**inverse** des gradients `rsa`/`n_fam` → la conservation suit l'axe d'exposition.

# In[83]:


# === PCA des positions d'actine — conservation INCLUSE comme descripteur ; couleur = role d'interface ===
# Non circulaire : la couleur (categorie d'interface) n'est pas la conservation.
sfeat = ['conservation','exposition','enfouissement','dist_nucleotide','n_abp']  # DETERMINANTS : exposition (intrinseque) / enfoui par le filament / enfoui par les ABP + promiscuite n_abp (cluster seq 70%)
def _catpos(r):
    if not r['at_interface']: return 'non-interface'
    if r['at_homo'] and r['at_hetero']: return 'actine-actine + ABP'
    if r['at_hetero']: return 'actine-ABP'
    if r['at_homo']: return 'actine-actine (filament)'
    return 'autre'
P['_cat']=P.apply(_catpos,axis=1)
P['exposition']=P['rsa']   # libelle clair pour rsa (accessibilite au solvant = exposition)
P['enfouissement']=P['combined_asa']   # enfouissement total (filament + ABP)
# === NOUVEAU descripteur : distance au nucleotide (ADP/ATP/ANP + Mg) = proximite du site fonctionnel ===
if 'dist_nucleotide' not in P.columns:
    import glob as _gl, os as _osn
    from pathlib import Path as _Pn
    _Rn=_Pn('..') if (_Pn('..')/'data').exists() else _Pn('.')
    _NUC={'ADP','ATP','ANP','MG'}
    def _pnuc(_f):
        aA={}; het=[]
        for ln in open(_f):
            if ln.startswith('ATOM') and ln[21]=='A':
                rs=ln[22:26].strip()
                if rs.lstrip('-').isdigit():
                    try: aA.setdefault(int(rs),[]).append((float(ln[30:38]),float(ln[38:46]),float(ln[46:54])))
                    except: pass
            elif ln.startswith('HETATM') and ln[17:20].strip() in _NUC:
                try: het.append((float(ln[30:38]),float(ln[38:46]),float(ln[46:54])))
                except: pass
        return aA,het
    _f3n=pd.read_csv(_Rn/'data/filtered/details/3.interface_residues.csv')
    _f3n['canon']=pd.to_numeric(_f3n['residue_number_canon_mafft'],errors='coerce')
    _f3n['rn']=pd.to_numeric(_f3n['residue_number_structure'],errors='coerce')
    _f3n=_f3n.dropna(subset=['canon','rn'])
    _din=pd.read_csv(_Rn/'data/filtered/details/1.interactions.csv').set_index('interaction_id')
    _rowsn=[]; _ndone=0
    for _f in _gl.glob(str(_Rn/'data/filtered/details/structures_files/pairwise/*.pdb')):
        _t=open(_f).read()
        if not any(n in _t for n in ('ADP','ATP','ANP')): continue
        _iid=int(_osn.path.basename(_f).split('_')[0])
        if _iid not in _din.index: continue
        _cA=_din.loc[_iid,'chain_A_id']
        _sub=_f3n[(_f3n.interaction_id==_iid)&(_f3n['chain']==_cA)]
        if len(_sub)<6: continue
        _of=(_sub['canon']-_sub['rn']).round()
        if _of.std()>0.5: continue                 # offset non constant (indels) -> on saute
        _off=int(_of.median()); _aA,_het=_pnuc(_f)
        if not _het or not _aA: continue
        _H=np.array(_het)
        for _rn,_xy in _aA.items():
            _d=np.sqrt(((np.array(_xy)[:,None,:]-_H[None,:,:])**2).sum(-1)).min()
            _rowsn.append((_rn+_off,_d))
        _ndone+=1
        if _ndone>=120: break
    _dnser=pd.DataFrame(_rowsn,columns=['canon','dist']).groupby('canon')['dist'].median()
    P['dist_nucleotide']=_dnser.reindex(P.index)
    P['dist_nucleotide']=P['dist_nucleotide'].fillna(P['dist_nucleotide'].median())
    print(f"dist_nucleotide ajoute : {_ndone} structures, {_dnser.size} positions canon")
R = P.dropna(subset=sfeat).copy()
Xs=R[sfeat].values.astype(float); Xsz=(Xs-Xs.mean(0))/(Xs.std(0)+1e-9)
from numpy.linalg import svd as _svd
Us,Ss,Vts=_svd(Xsz-Xsz.mean(0),full_matrices=False); pcs=Us*Ss; evs=(Ss**2)/np.sum(Ss**2)
_cor=np.array([[np.corrcoef(Xsz[:,j],pcs[:,k])[0,1] for k in range(2)] for j in range(len(sfeat))])

from matplotlib.lines import Line2D
_catcol={'actine-actine (filament)':'#D55E00','actine-ABP':'#56B4E9','actine-actine + ABP':'#009E73','non-interface':'#999999'}
_order=['actine-actine (filament)','actine-ABP','actine-actine + ABP','non-interface']
fig,ax=plt.subplots(figsize=(10,8))
for g in _order:
    m=(R['_cat']==g).values
    if mc:=m.any():
        ax.scatter(pcs[m,0],pcs[m,1],c=_catcol[g],s=34,edgecolor='k',lw=.25,alpha=.8,label=f'{g} (n={int(m.sum())})',zorder=2)
_f=0.8*np.abs(pcs).max()/np.abs(_cor).max()
for j,nm in enumerate(sfeat):
    _isc=(nm=='conservation')
    ax.arrow(0,0,_cor[j,0]*_f,_cor[j,1]*_f,color='#555',head_width=.10,lw=1.1,length_includes_head=True,alpha=.8,zorder=3)
    ax.annotate({'conservation':'sensibilité mutationnelle','enfouissement moyen (%ASA)':'engagement interaction (%ASA)'}.get(nm,nm),(_cor[j,0]*_f*1.07,_cor[j,1]*_f*1.07),fontsize=(8.5 if _isc else 7.5),
                color='#333',fontweight=('bold' if _isc else 'normal'),zorder=6)
ax.axhline(0,color='#eee',lw=.6); ax.axvline(0,color='#eee',lw=.6)
ax.set_xlabel(f'PC1 ({evs[0]*100:.0f}%)'); ax.set_ylabel(f'PC2 ({evs[1]*100:.0f}%)')
ax.legend(loc='best',fontsize=8)
ax.set_title("PCA des positions d'actine — determinants de la sensibilite mutationnelle ; couleur = role d'interface\n"
             f"{len(R)} positions ; PC1+PC2 = {(evs[0]+evs[1])*100:.0f}%")
plt.tight_layout(); plt.show()
print(f"{len(R)} positions | categories: {R['_cat'].value_counts().to_dict()}")


# ### Variante — positions, COEUR (>25% ASA)\n\nasa_mean et n_fam recalcules sur les seuls contacts ABP enfouis > 25% ASA.

# In[84]:


# === PCA des positions — COEUR (>25% ASA) : asa_mean et n_fam recalcules sur contacts enfouis ===
from pathlib import Path as _Pp
_Rp=_Pp('..') if (_Pp('..')/'data').exists() else _Pp('.')
_f3p=pd.read_csv(_Rp/'data/filtered/details/3.interface_residues.csv')
_dip=pd.read_csv(_Rp/'data/filtered/details/1.interactions.csv')
_dap=pd.read_csv(_Rp/'data/filtered/filtered_all_data.csv',low_memory=False)
_ppp=pd.read_csv(_Rp/'data/filtered/proteins_per_pdb.csv'); _acp=set(_ppp[_ppp.is_actin]['chain'].str.lower())
_f3p['pct']=pd.to_numeric(_f3p['buried_ASA_percent'].astype(str).str.replace('%','',regex=False),errors='coerce')
_f3p['canon']=pd.to_numeric(_f3p['residue_number_canon_mafft'],errors='coerce')
_homop=set(_dip[_dip['chain_B_id'].str.lower().isin(_acp)]['interaction_id'])
_dp=_f3p.dropna(subset=['canon','pct']).copy(); _dp['canon']=_dp.canon.astype(int)
_hetp=_dp[(~_dp.interaction_id.isin(_homop))&(_dp['pct']>25)]      # COEUR : contacts hetero enfouis >25%
_mp=_dip.merge(_dap[['subunit_1','subunit_2','s2_sequence_cluster_70']].drop_duplicates(),
               left_on=['chain_A_id','chain_B_id'],right_on=['subunit_1','subunit_2'],how='left')
Pc=P.copy()
Pc['asa_mean']=_hetp.groupby('canon')['pct'].mean().reindex(Pc.index)
Pc['n_abp']=_hetp.assign(f=_hetp.interaction_id.map(_mp.set_index('interaction_id')['s2_sequence_cluster_70'])).groupby('canon')['f'].nunique().reindex(Pc.index).fillna(0)
Pc['exposition']=Pc['rsa']   # libelle clair pour rsa
Pc['enfouissement moyen (%ASA)']=_dp[_dp['pct']>25].groupby('canon')['pct'].mean().reindex(Pc.index).fillna(0)   # COEUR : profondeur moyenne sur contacts >25% (filament + ABP) ; non enfouies = 0

sfeat=['conservation','exposition','enfouissement moyen (%ASA)','dist_nucleotide','n_abp']  # DETERMINANTS (idem cellule 6) ; n_abp = cluster seq 70% sur contacts COEUR
R=Pc.dropna(subset=sfeat).copy()
Xs=R[sfeat].values.astype(float); Xsz=(Xs-Xs.mean(0))/(Xs.std(0)+1e-9)
from numpy.linalg import svd as _svd
Us,Ss,Vts=_svd(Xsz-Xsz.mean(0),full_matrices=False); pcs=Us*Ss; evs=(Ss**2)/np.sum(Ss**2)
_cor=np.array([[np.corrcoef(Xsz[:,j],pcs[:,k])[0,1] for k in range(2)] for j in range(len(sfeat))])
from matplotlib.lines import Line2D
_catcol={'actine-actine (filament)':'#D55E00','actine-ABP':'#56B4E9','actine-actine + ABP':'#009E73','non-interface':'#999999'}
_order=['actine-actine (filament)','actine-ABP','actine-actine + ABP','non-interface']
fig,ax=plt.subplots(figsize=(10,8))
for g in _order:
    m=(R['_cat']==g).values
    if m.any():
        ax.scatter(pcs[m,0],pcs[m,1],c=_catcol[g],s=34,edgecolor='k',lw=.25,alpha=.8,label=f'{g} (n={int(m.sum())})',zorder=2)
_f=0.8*np.abs(pcs).max()/np.abs(_cor).max()
for j,nm in enumerate(sfeat):
    _isc=(nm=='conservation')
    ax.arrow(0,0,_cor[j,0]*_f,_cor[j,1]*_f,color='#555',head_width=.10,lw=1.1,length_includes_head=True,alpha=.8,zorder=3)
    ax.annotate({'conservation':'sensibilité mutationnelle','enfouissement moyen (%ASA)':'engagement interaction (%ASA)'}.get(nm,nm),(_cor[j,0]*_f*1.07,_cor[j,1]*_f*1.07),fontsize=(8.5 if _isc else 7.5),color='#333',fontweight=('bold' if _isc else 'normal'),zorder=6)
ax.axhline(0,color='#eee',lw=.6); ax.axvline(0,color='#eee',lw=.6)
ax.set_xlabel(f'PC1 ({evs[0]*100:.0f}%)'); ax.set_ylabel(f'PC2 ({evs[1]*100:.0f}%)')
ax.legend(loc='best',fontsize=8)
ax.set_title("PCA des positions d'actine — COEUR >25% ASA · determinants de la sensibilite mutationnelle\n"
             f"exposition / enfouissement total / distance nucleotide / promiscuite (n_abp) ; sensibilité mutationnelle = fleche ; {len(R)} positions ; PC1+PC2 = {(evs[0]+evs[1])*100:.0f}%")
plt.tight_layout(); plt.show()
print(f"{len(R)} positions (coeur) | categories: {R['_cat'].value_counts().to_dict()}")


# ### (Exploratoire) déterminants + variables intrinsèques — pour voir
# 
# Copie de la PCA des déterminants, **plus** les variables intrinsèques du résidu (hydrophobicité, aromaticité, charge, polarité) ajoutées comme descripteurs. Points recoloriés (filament rouge, ABP jaune, both vert, non-interface noir), labels et points agrandis.

# In[85]:


# === PCA des positions — COEUR >25% ASA + variables INTRINSEQUES (hydro/aromatique/charge/polaire) — POUR VOIR ===
from pathlib import Path as _Pp
_Rp=_Pp('..') if (_Pp('..')/'data').exists() else _Pp('.')
_f3p=pd.read_csv(_Rp/'data/filtered/details/3.interface_residues.csv')
_dip=pd.read_csv(_Rp/'data/filtered/details/1.interactions.csv')
_dap=pd.read_csv(_Rp/'data/filtered/filtered_all_data.csv',low_memory=False)
_ppp=pd.read_csv(_Rp/'data/filtered/proteins_per_pdb.csv'); _acp=set(_ppp[_ppp.is_actin]['chain'].str.lower())
_f3p['pct']=pd.to_numeric(_f3p['buried_ASA_percent'].astype(str).str.replace('%','',regex=False),errors='coerce')
_f3p['canon']=pd.to_numeric(_f3p['residue_number_canon_mafft'],errors='coerce')
_homop=set(_dip[_dip['chain_B_id'].str.lower().isin(_acp)]['interaction_id'])
_dp=_f3p.dropna(subset=['canon','pct']).copy(); _dp['canon']=_dp.canon.astype(int)
_hetp=_dp[(~_dp.interaction_id.isin(_homop))&(_dp['pct']>25)]      # COEUR : contacts hetero enfouis >25%
_mp=_dip.merge(_dap[['subunit_1','subunit_2','s2_sequence_cluster_70']].drop_duplicates(),
               left_on=['chain_A_id','chain_B_id'],right_on=['subunit_1','subunit_2'],how='left')
Pc=P.copy()
Pc['asa_mean']=_hetp.groupby('canon')['pct'].mean().reindex(Pc.index)
Pc['n_abp']=_hetp.assign(f=_hetp.interaction_id.map(_mp.set_index('interaction_id')['s2_sequence_cluster_70'])).groupby('canon')['f'].nunique().reindex(Pc.index).fillna(0)
Pc['exposition']=Pc['rsa']   # libelle clair pour rsa
Pc['enfouissement moyen (%ASA)']=_dp[_dp['pct']>25].groupby('canon')['pct'].mean().reindex(Pc.index).fillna(0)   # COEUR : profondeur moyenne sur contacts >25% (filament + ABP) ; non enfouies = 0


# --- descripteurs INTRINSEQUES (acide amine de la position via P60709) ---
_P60709b=('MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDD'
          'MEKIWHHTFYNELRVAPEEHPVLLTEAPLNPKANREKMTQIMFETFNTPAMYVAIQAVLSLYASGRTTGIVMDSGDGVTHTV'
          'PIYEGYALPHAILRLDLAGRDLTDYLMKILTERGYSFTTTAEREIVRDIKEKLCYVALDFEQEMATAASSSSLEKSYELPDG'
          'QVITIGNERFRCPEALFQPSFLGMESCGIHETTFNSIMKCDVDIRKDLYANTVLSGGTTMYPGIADRMQKEITALAPSTMKI'
          'KIIAPPERKYSVWIGGSILASLSTFQQMWISKQEYDESGPSIVHRKCF')
_KDb2={'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
_AROb=set('FWY'); _POLb=set('STNQHC'); _POSb=set('KR'); _NEGb=set('DE')
def _aab(r):
    try:
        i=int(r); return _P60709b[i-1] if 1<=i<=len(_P60709b) else None
    except Exception: return None
Pc['_aa']=Pc['Residue'].map(_aab) if 'Residue' in Pc.columns else None
Pc['hydrophobicite']=Pc['_aa'].map(lambda a:_KDb2.get(a,np.nan))
Pc['aromaticite']=Pc['_aa'].map(lambda a:1.0 if a in _AROb else 0.0)
Pc['charge']=Pc['_aa'].map(lambda a:1.0 if a in _POSb else (-1.0 if a in _NEGb else 0.0))
Pc['polarite']=Pc['_aa'].map(lambda a:1.0 if a in _POLb else 0.0)

sfeat=['conservation','exposition','enfouissement moyen (%ASA)','dist_nucleotide','n_abp','hydrophobicite','aromaticite','charge','polarite']  # + variables intrinseques (pour voir)
R=Pc.dropna(subset=sfeat).copy()
Xs=R[sfeat].values.astype(float); Xsz=(Xs-Xs.mean(0))/(Xs.std(0)+1e-9)
from numpy.linalg import svd as _svd
Us,Ss,Vts=_svd(Xsz-Xsz.mean(0),full_matrices=False); pcs=Us*Ss; evs=(Ss**2)/np.sum(Ss**2)
_cor=np.array([[np.corrcoef(Xsz[:,j],pcs[:,k])[0,1] for k in range(2)] for j in range(len(sfeat))])
from matplotlib.lines import Line2D
_catcol={'actine-actine (filament)':'#e74c3c','actine-ABP':'#f1c40f','actine-actine + ABP':'#009E73','non-interface':'#000000'}
_order=['actine-actine (filament)','actine-ABP','actine-actine + ABP','non-interface']
fig,ax=plt.subplots(figsize=(12,9.5))
for g in _order:
    m=(R['_cat']==g).values
    if m.any():
        ax.scatter(pcs[m,0],pcs[m,1],c=_catcol[g],s=75,edgecolor='k',lw=.35,alpha=.85,label=f'{g} (n={int(m.sum())})',zorder=5)
_f=0.8*np.abs(pcs).max()/np.abs(_cor).max()
for j,nm in enumerate(sfeat):
    _isc=(nm=='conservation')
    ax.arrow(0,0,_cor[j,0]*_f,_cor[j,1]*_f,color='#555',head_width=.10,lw=1.1,length_includes_head=True,alpha=.8,zorder=3)
    ax.annotate({'conservation':'sensibilité mutationnelle','enfouissement moyen (%ASA)':'engagement interaction (%ASA)'}.get(nm,nm),(_cor[j,0]*_f*1.07,_cor[j,1]*_f*1.07),fontsize=(15 if _isc else 13),color='#333',fontweight=('bold' if _isc else 'normal'),zorder=6)
ax.axhline(0,color='#eee',lw=.6); ax.axvline(0,color='#eee',lw=.6)
ax.set_xlabel(f'PC1 ({evs[0]*100:.0f}%)'); ax.set_ylabel(f'PC2 ({evs[1]*100:.0f}%)')
ax.legend(loc='best',fontsize=8)
ax.set_title("PCA des positions d'actine — COEUR >25% ASA + variables intrinseques (pour voir)\n"
             f"exposition / enfouissement total / distance nucleotide / promiscuite (n_abp) ; sensibilité mutationnelle = fleche ; {len(R)} positions ; PC1+PC2 = {(evs[0]+evs[1])*100:.0f}%")
plt.tight_layout(); plt.show()
print(f"{len(R)} positions (coeur) | categories: {R['_cat'].value_counts().to_dict()}")


# ### Variante — positions décrites par leurs propriétés PHYSICO-CHIMIQUES intrinsèques
# 
# Puisqu'un individu est **une position** (= un acide aminé de l'actine), on la décrit par ses propriétés **intrinsèques** : hydrophobicité, charge, aromaticité, polarité du résidu — définies pour **toutes** les positions (interface ou non) — plus les descripteurs **structuraux** (exposition, enfouissement moyen %ASA, distance au nucléotide, n_abp) et la **conservation** (incluse). Pas de chimie de contacts ici (relationnelle) : elle est réservée à la PCA des clusters d'interface.

# In[86]:


# === PCA des positions — descripteurs INTRINSEQUES du residu + structuraux + conservation ===
# Individu = une position d'actine. Descripteurs PROPRES a la position (definis partout) :
#   intrinseques : hydrophobicite, charge, aromaticite, polarite (de l'acide amine)
#   structuraux  : exposition (rsa), enfouissement moyen (%ASA), distance nucleotide, n_abp
#   sensibilité mutationnelle INCLUSE (fleche doree). Toutes les positions (interface + non-interface).
from pathlib import Path as _Pchem
_Rc = _Pchem('..') if (_Pchem('..')/'data').exists() else _Pchem('.')
_acset = set(pd.read_csv(_Rc/'data/filtered/proteins_per_pdb.csv').query('is_actin')['chain'].str.lower())

# --- acide amine de CHAQUE position (sequence actine de reference P60709) ---
_P60709 = ('MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDD'
           'MEKIWHHTFYNELRVAPEEHPVLLTEAPLNPKANREKMTQIMFETFNTPAMYVAIQAVLSLYASGRTTGIVMDSGDGVTHTV'
           'PIYEGYALPHAILRLDLAGRDLTDYLMKILTERGYSFTTTAEREIVRDIKEKLCYVALDFEQEMATAASSSSLEKSYELPDG'
           'QVITIGNERFRCPEALFQPSFLGMESCGIHETTFNSIMKCDVDIRKDLYANTVLSGGTTMYPGIADRMQKEITALAPSTMKI'
           'KIIAPPERKYSVWIGGSILASLSTFQQMWISKQEYDESGPSIVHRKCF')
_KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
_ARO = set('FWY'); _POL = set('STNQHC'); _POSc = set('KR'); _NEGc = set('DE')
_Pc = P.copy()
def _aa_of(r):
    try:
        i=int(r)
        return _P60709[i-1] if 1<=i<=len(_P60709) else None
    except Exception:
        return None
_Pc['_aa'] = _Pc['Residue'].map(_aa_of) if 'Residue' in _Pc.columns else None
# descripteurs intrinseques
_Pc['hydrophobicite'] = _Pc['_aa'].map(lambda a:_KD.get(a, np.nan))
_Pc['charge']         = _Pc['_aa'].map(lambda a: 1.0 if a in _POSc else (-1.0 if a in _NEGc else 0.0))
_Pc['aromaticite']    = _Pc['_aa'].map(lambda a: 1.0 if a in _ARO else 0.0)
_Pc['polarite']       = _Pc['_aa'].map(lambda a: 1.0 if a in _POL else 0.0)

# --- descripteurs structuraux (coeur >25% ASA) ---
_f3p = pd.read_csv(_Rc/'data/filtered/details/3.interface_residues.csv')
_f3p['pct']   = pd.to_numeric(_f3p['buried_ASA_percent'].astype(str).str.replace('%','',regex=False), errors='coerce')
_f3p['canon'] = pd.to_numeric(_f3p['residue_number_canon_mafft'], errors='coerce')
_dip = pd.read_csv(_Rc/'data/filtered/details/1.interactions.csv')
_dap = pd.read_csv(_Rc/'data/filtered/filtered_all_data.csv', low_memory=False)
_homop = set(_dip[_dip['chain_B_id'].str.lower().isin(_acset)]['interaction_id'])
_dp = _f3p.dropna(subset=['canon','pct']).copy(); _dp['canon'] = _dp.canon.astype(int)
_enf = _dp[_dp['pct']>25].groupby('canon')['pct'].mean()
_mp = _dip.merge(_dap[['subunit_1','subunit_2','s2_sequence_cluster_70']].drop_duplicates(),
                 left_on=['chain_A_id','chain_B_id'], right_on=['subunit_1','subunit_2'], how='left')
_hetp = _dp[(~_dp.interaction_id.isin(_homop))&(_dp['pct']>25)]
_nabp = _hetp.assign(f=_hetp.interaction_id.map(_mp.set_index('interaction_id')['s2_sequence_cluster_70'])).groupby('canon')['f'].nunique()

def _catp(r):
    if not r.get('at_interface', False): return 'non-interface'
    if r.get('at_homo', False) and r.get('at_hetero', False): return 'actine-actine + ABP'
    if r.get('at_hetero', False): return 'actine-ABP'
    if r.get('at_homo', False): return 'actine-actine (filament)'
    return 'autre'
_Pc['_cat'] = _Pc.apply(_catp, axis=1)
_Pc['exposition']                = _Pc['rsa']
_Pc['enfouissement moyen (%ASA)'] = _enf.reindex(_Pc.index).fillna(0.0)
_Pc['dist_nucleotide']           = _Pc['dist_nucleotide'] if 'dist_nucleotide' in _Pc.columns else np.nan
_Pc['n_abp']                     = _nabp.reindex(_Pc.index).fillna(0)
# COHERENCE Figure 2 : non-interface -> enfouissement & n_abp = 0 (suit l'etiquette de P)
_m_ni = (_Pc['_cat']=='non-interface')
_Pc.loc[_m_ni, ['enfouissement moyen (%ASA)','n_abp']] = 0.0

# classe chimique du residu (1 seule par acide amine) -> FORME du point (pas un descripteur)
def _chemclass(a):
    if a in _ARO: return 'aromatique'
    if a in _POSc: return 'charge +'
    if a in _NEGc: return 'charge -'
    if a in _POL: return 'polaire'
    return 'hydrophobe'
_Pc['_chem'] = _Pc['_aa'].map(_chemclass)
_feats_all = ['exposition','enfouissement moyen (%ASA)','dist_nucleotide','n_abp']   # descripteurs STRUCTURAUX seuls
Q = _Pc.dropna(subset=_feats_all+['conservation','_aa']).copy()
print(f'{len(Q)} positions | {len(_feats_all)} descripteurs structuraux + conservation ; chimie = forme du point')
print(f'roles : {Q["_cat"].value_counts().to_dict()}')

_catcol = {'actine-actine (filament)':'#D55E00','actine-ABP':'#f1c40f','actine-actine + ABP':'#009E73','non-interface':'#888888'}
_order = ['actine-actine (filament)','actine-ABP','actine-actine + ABP','non-interface']
_chemmark = {'aromatique':'^','charge +':'P','charge -':'v','polaire':'o','hydrophobe':'s'}
def _draw_pca(feats, ax, titre):
    X=Q[feats].values.astype(float); Xz=(X-X.mean(0))/(X.std(0)+1e-9)
    from numpy.linalg import svd as _svd
    U,S,Vt=_svd(Xz-Xz.mean(0),full_matrices=False); pcs=U*S; evs=(S**2)/np.sum(S**2)
    cor=np.array([[np.corrcoef(Xz[:,j],pcs[:,k])[0,1] for k in range(2)] for j in range(len(feats))])
    _catv=Q['_cat'].values; _chemv=Q['_chem'].values
    for g in _order:
        for cm,mk in _chemmark.items():
            m=(_catv==g)&(_chemv==cm)
            if m.any(): ax.scatter(pcs[m,0],pcs[m,1],c=_catcol[g],marker=mk,s=55,edgecolor='k',lw=.25,alpha=.8,zorder=2)
    _f=0.8*np.abs(pcs).max()/(np.abs(cor).max()+1e-9)
    for j,nm in enumerate(feats):
        _isc=(nm=='conservation')
        ax.arrow(0,0,cor[j,0]*_f,cor[j,1]*_f,color='#555',head_width=.10,lw=(1.6 if _isc else 1.1),length_includes_head=True,alpha=.85,zorder=4)
        ax.annotate({'conservation':'sensibilité mutationnelle','enfouissement moyen (%ASA)':'engagement interaction (%ASA)'}.get(nm,nm),(cor[j,0]*_f*1.07,cor[j,1]*_f*1.07),fontsize=(9 if _isc else 8),color='#333',fontweight=('bold' if _isc else 'normal'),zorder=6)
    ax.axhline(0,color='#eee',lw=.6); ax.axvline(0,color='#eee',lw=.6)
    ax.set_xlabel(f'PC1 ({evs[0]*100:.0f}%)'); ax.set_ylabel(f'PC2 ({evs[1]*100:.0f}%)')
    from matplotlib.lines import Line2D as _L2
    _lc=[_L2([0],[0],marker='o',ls='',mfc=_catcol[g],mec='k',ms=8,label=f'{g} (n={int((_catv==g).sum())})') for g in _order if (_catv==g).any()]
    _lm=[_L2([0],[0],marker=mk,ls='',mfc='#bbb',mec='k',ms=8,label=cm) for cm,mk in _chemmark.items()]
    _lg1=ax.legend(handles=_lc,loc='upper right',fontsize=7,title='couleur = role'); ax.add_artist(_lg1)
    ax.legend(handles=_lm,loc='lower right',fontsize=7,title='forme = chimie du residu')
    ax.set_title(titre,fontsize=10)

fig,ax=plt.subplots(figsize=(11,9))
_draw_pca(_feats_all+['conservation'], ax,
          "PCA des positions d'actine — descripteurs structuraux + sensibilité mutationnelle (fleche en gras)\ncouleur = role · forme = classe chimique du residu ; toutes positions")
plt.tight_layout(); plt.show()


# ## Figure 3 — Tableau de bord des déterminants (brut vs à exposition égale)
# 
# Pour chaque descripteur : corrélation de Spearman avec la conservation, **brute** et **à rsa égal**
# (résidualisation = on enlève l'effet de l'exposition). C'est la preuve directe, descripteur par descripteur.
# Une barre qui **change de signe** entre brut et « à rsa égal » = effet porté par l'exposition (artefact).

# In[87]:


def _resid(a,b):
    b=np.c_[np.ones(len(b)),b]; return a-b@np.linalg.lstsq(b,a,rcond=None)[0]
def rho_partial(col):
    s=P.dropna(subset=[col,'conservation','rsa'])
    return spearmanr(_resid(s[col].values,s['rsa'].values),_resid(s['conservation'].values,s['rsa'].values))[0]
def rho_brut(col):
    s=P.dropna(subset=[col,'conservation']); r,p=spearmanr(s[col],s['conservation']); return r,p

descr=[('rsa','exposition (rsa)'),('asa_mean','buried ASA% moyen'),('combined_asa','enfouissement total'),
       ('homo_asa','enfoui filament'),('n_fam','nb familles ABP'),('n_abp','nb ABP distincts'),
       ('hydrophobicite','hydrophobicite'),('aromatique','aromatique')]
rows=[]
for col,lab in descr:
    rb,pb=rho_brut(col); rp=(np.nan if col=='rsa' else rho_partial(col))
    rows.append((lab,rb,pb,rp))
T=pd.DataFrame(rows,columns=['lab','brut','p','partial']).iloc[::-1]
fig,ax=plt.subplots(figsize=(9,6)); y=np.arange(len(T)); h=0.38
cb=['#c0392b' if p<0.05 else '#bbbbbb' for p in T['p']]
ax.barh(y+h/2,T['brut'],height=h,color=cb,label='brut')
ax.barh(y-h/2,T['partial'],height=h,color='#2980b9',alpha=.85,label='à rsa égal')
for i,(_,r) in enumerate(T.iterrows()):
    s='***' if r['p']<1e-3 else '**' if r['p']<1e-2 else '*' if r['p']<.05 else 'ns'
    ax.text(r['brut']+(0.01 if r['brut']>=0 else -0.01),i+h/2,s,va='center',ha='left' if r['brut']>=0 else 'right',fontsize=8)
ax.set_yticks(y); ax.set_yticklabels(T['lab']); ax.axvline(0,color='k',lw=.8)
ax.set_xlabel('Spearman rho vs conservation'); ax.legend(loc='lower right',fontsize=9)
ax.set_title("Déterminants de la sensibilité mutationnelle : brut (rouge=signif) vs à exposition égale (bleu)\nseul l'enfouissement (rsa) reste fort ; nb de partenaires s'inverse ; hydrophobicité ≈ 0")
plt.tight_layout(); plt.show()
print(T.round(3).to_string(index=False))


# # B — PCA des interfaces
# 
# Un point = **un cluster d'interface C70**. Montre que les interfaces ont des **profils physico-chimiques différents**.

# ## Figure 4 — PCA des interfaces : les clusters diffèrent en physico-chimie (claim 6)
# 
# Ici un point = un **cluster d'interface C70** (≥ 3 interactions), décrit par la physico-chimie
# de ses résidus d'actine (hydrophobicité, charges +/−, aromaticité, enfouissement). Couleur =
# conservation moyenne du cluster, taille = nb d'interactions. Si les points **s'étalent** (pas un
# seul nuage compact), c'est que les interfaces n'ont **pas le même profil physico-chimique**.

# In[88]:


# === PCA des interfaces C70 : profils physico-chimiques cote actine ===
_cons = P['conservation']
_m = di.merge(da[['subunit_1','subunit_2','s2_actine','cluster_data_70']].drop_duplicates(),
              left_on=['chain_A_id','chain_B_id'], right_on=['subunit_1','subunit_2'], how='left')
_m['s2_actine'] = _m['s2_actine'].fillna(False).astype(bool)
_het = _m[~_m.s2_actine].dropna(subset=['cluster_data_70'])
_f3 = df3.copy()
_f3['aa']  = _f3['residue_name'].astype(str).str.upper().str.strip()
_f3['c70'] = _f3.interaction_id.map(_het.set_index('interaction_id')['cluster_data_70'])
_f3['act'] = _f3.interaction_id.map(_het.set_index('interaction_id')['chain_A_id'].str.lower())
_a = _f3[(_f3.chain.str.lower()==_f3.act) & _f3.c70.notna() & (_f3.aa.str.len()==1)].copy()
_a['kd']=_a.aa.map(KD); _a['pos']=_a.aa.isin(list('KR')); _a['neg']=_a.aa.isin(list('DE'))
_a['aro']=_a.aa.isin(list('FWY')); _a['cons']=_a.canon.map(_cons)
_g = _a.groupby('c70')
desc = pd.DataFrame({'hydrophobicite':_g['kd'].mean(),'charge +':_g['pos'].mean(),'charge -':_g['neg'].mean(),
                     'aromatique':_g['aro'].mean(),'enfoui (ASA%)':_g['pct'].mean(),
                     'conservation':_g['cons'].mean(),'n_intx':_g['interaction_id'].nunique()})
desc = desc[desc.n_intx>=3]
ifeat=['hydrophobicite','charge +','charge -','aromatique','enfoui (ASA%)']
Xi=desc[ifeat].values; Xiz=(Xi-Xi.mean(0))/(Xi.std(0)+1e-9)
Ui,Si,Vti=np.linalg.svd(Xiz,full_matrices=False); pci=Ui*Si; evi=(Si**2)/np.sum(Si**2)

fig,ax=plt.subplots(figsize=(9,7.5))
sc=ax.scatter(pci[:,0],pci[:,1],c=desc['conservation'],cmap='viridis',
              s=20+desc['n_intx']*2.5,edgecolor='k',lw=.3,alpha=.85)
fig.colorbar(sc,ax=ax,label='conservation moyenne du cluster')
_off={'aromatique':(0,-0.28),'enfoui (ASA%)':(0,0.22)}   # decalage anti-chevauchement
for i,f in enumerate(ifeat):
    ax.arrow(0,0,Vti[0,i]*3.2,Vti[1,i]*3.2,color='#c0392b',head_width=.13,lw=1.6,length_includes_head=True,alpha=.8)
    ox,oy=_off.get(f,(0,0))
    ax.text(Vti[0,i]*3.2*1.12+ox,Vti[1,i]*3.2*1.12+oy,f,color='#c0392b',fontsize=9.5,ha='center',fontweight='bold')
ax.axhline(0,color='#eee',lw=.8); ax.axvline(0,color='#eee',lw=.8)
ax.set_xlabel(f'PC1 ({evi[0]*100:.0f}%)'); ax.set_ylabel(f'PC2 ({evi[1]*100:.0f}%)')
ax.set_title(f"PCA des interfaces (cluster C70, n={len(desc)}) — profils physico-chimiques\n"
             "les clusters s'etalent : hydrophobes vs charges, enfouis vs superficiels (taille = nb interactions)",fontsize=10)
plt.tight_layout(); plt.show()
print(f"{len(desc)} clusters | etendue hydrophobicite {desc['hydrophobicite'].min():.1f}..{desc['hydrophobicite'].max():.1f} "
      f"| charge+ {desc['charge +'].min():.2f}..{desc['charge +'].max():.2f} | enfoui {desc['enfoui (ASA%)'].min():.0f}..{desc['enfoui (ASA%)'].max():.0f}%")


# ### Vue détaillée — PCA des interfaces C70 (actine S1 + partenaire S2, homo vs hétéro)

# In[89]:


# === PCA des interfaces au niveau cluster C70 (homo + hetero) ===
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path as _P
from scipy.spatial.distance import pdist, squareform
_R = _P('..') if _P('../data').exists() else _P('.')
_df3=pd.read_csv(_R/'data/filtered/details/3.interface_residues.csv')
_df3['pct']=pd.to_numeric(_df3['buried_ASA_percent'].astype(str).str.replace('%','',regex=False),errors='coerce')
_df3['a2']=pd.to_numeric(_df3['buried_ASA_Å²'].astype(str).str.replace('<0.1','0.05',regex=False),errors='coerce')
_df3['canon']=pd.to_numeric(_df3['residue_number_canon_mafft'],errors='coerce')
_df3=_df3[_df3.canon.notna()&_df3.pct.notna()].copy()
_di=pd.read_csv(_R/'data/filtered/details/1.interactions.csv')
_da=pd.read_csv(_R/'data/filtered/filtered_all_data.csv',low_memory=False)
_m=_di.merge(_da[['subunit_1','subunit_2','subunit_2_title','s2_actine','cluster_data_70']].drop_duplicates(),
             left_on=['chain_A_id','chain_B_id'],right_on=['subunit_1','subunit_2'],how='left').drop_duplicates('interaction_id')
_m=_m.dropna(subset=['cluster_data_70'])
_m['homo']=_m['s2_actine'].fillna(False).astype(bool)
_m['abp_title']=_m['subunit_2_title'].fillna('?').astype(str).str.replace(r'\s*\(.*?\)','',regex=True).str.strip().str[:22]
_clu2abp=_m[~_m['homo']].groupby('cluster_data_70')['abp_title'].agg(lambda x:x.value_counts().index[0])
_m['interface_area']=pd.to_numeric(_m['interface_area'].astype(str).str.replace('Å²','',regex=False).str.strip(),errors='coerce')
for _c in ['num_hbonds','num_salt_bridges','num_contacts']: _m[_c]=pd.to_numeric(_m[_c],errors='coerce')
_imap=_m.set_index('interaction_id'); _chA=_imap['chain_A_id'].str.lower(); _clu=_imap['cluster_data_70']
_d=_df3[_df3.interaction_id.isin(set(_m.interaction_id))].copy()
_d['_a']=_d.interaction_id.map(_chA); _d=_d[_d.chain.str.lower()==_d._a].copy()   # cote actine = chaine A (subunit_1)
_d['clu']=_d.interaction_id.map(_clu); _d['canon']=_d.canon.astype(int)
_cdf=pd.read_csv(_R/'data/proteocast/conservation_vs_asa_per_position.csv'); _cons=dict(zip(_cdf['canon'].astype(int),_cdf['conservation']))
_clsmap=dict(zip(_cdf['canon'].astype(int),_cdf['residue_class']))   # canon -> sensitive/tolerant
_KD={'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
def _md(s): mm=s.mode(); return mm.iloc[0] if len(mm) else None
_po=_d.groupby(['clu','canon']).agg(w=('pct','mean'),res=('residue_name',_md)).reset_index()
_po['kd']=_po['res'].map(lambda r:_KD.get(str(r),0.0)); _po['cons']=_po['canon'].map(lambda c:_cons.get(int(c),np.nan))
_po['ispos']=_po['res'].isin({'R','K'}).astype(float); _po['isneg']=_po['res'].isin({'D','E'}).astype(float); _po['isaro']=_po['res'].isin({'F','Y','W'}).astype(float)
def _wmi(g,c):
    w=g['w'].values; x=g[c].values; mk=~np.isnan(x)
    return (w[mk]*x[mk]).sum()/w[mk].sum() if w[mk].sum()>0 else np.nan
Fi=_po.groupby('clu').apply(lambda g:pd.Series({'conservation_actine':_wmi(g,'cons'),'hydrophobicite_actine':_wmi(g,'kd'),
    'charge_pos_actine':_wmi(g,'ispos'),'charge_neg_actine':_wmi(g,'isneg'),'aromaticite_actine':_wmi(g,'isaro'),'asa_moyen_actine':_wmi(g,'w')}))
Fi['n_residus']=_d.groupby('interaction_id').size().groupby(_clu).mean()
_ag=_m.groupby('cluster_data_70').agg(hb=('num_hbonds','sum'),sb=('num_salt_bridges','sum'),
    ct=('num_contacts','sum'),n=('interaction_id','nunique'),homo=('homo','mean'))
Fi['taille_interface_actine']=_d.groupby('interaction_id')['a2'].sum().groupby(_clu).mean()
Fi=Fi.join(_ag[['homo']])
_cc=pd.read_csv(_R/'data/filtered/details/4.inter-residue_contacts.csv'); _CHG={'R':1,'K':1,'D':-1,'E':-1}
_cc=_cc[_cc.interaction_id.isin(set(_m.interaction_id))].copy(); _cc['_clu']=_cc.interaction_id.map(_clu)
_cc['charge_compl']=-(_cc['residue_A_name'].map(lambda r:_CHG.get(str(r),0))*_cc['residue_B_name'].map(lambda r:_CHG.get(str(r),0)))
Fi['charge_compl']=_cc.groupby('_clu')['charge_compl'].mean()
Fi['pct_hbond']=100*_ag['hb']/_ag['ct']; Fi['pct_saltbridge']=100*_ag['sb']/_ag['ct']
# --- descripteurs cote PARTENAIRE (chaine B = S2 ; actine pour homo, ABP pour hetero) ---
_chB=_imap['chain_B_id'].str.lower()
_db=_df3[_df3.interaction_id.isin(set(_m.interaction_id))].copy()
_db['_b']=_db.interaction_id.map(_chB); _db=_db[_db.chain.str.lower()==_db._b].copy()
_db['clu']=_db.interaction_id.map(_clu)
_db['kd']=_db['residue_name'].map(lambda r:_KD.get(str(r),0.0))
_db['ispos']=_db['residue_name'].isin({'R','K'}).astype(float)
_db['isneg']=_db['residue_name'].isin({'D','E'}).astype(float)
_db['isaro']=_db['residue_name'].isin({'F','Y','W'}).astype(float)
def _wmb(g,col):
    w=g['pct'].values; x=g[col].values; mk=~np.isnan(x)
    return (w[mk]*x[mk]).sum()/w[mk].sum() if w[mk].sum()>0 else np.nan
Fb=_db.groupby('clu').apply(lambda g:pd.Series({'hydrophobicite_partenaire':_wmb(g,'kd'),
    'charge_pos_partenaire':_wmb(g,'ispos'),'charge_neg_partenaire':_wmb(g,'isneg'),'aromaticite_partenaire':_wmb(g,'isaro'),'asa_moyen_partenaire':g['pct'].mean()}))
Fi=Fi.join(Fb)
_featI=['conservation_actine','hydrophobicite_actine','charge_pos_actine','charge_neg_actine','aromaticite_actine','taille_interface_actine','asa_moyen_actine','pct_hbond','pct_saltbridge','hydrophobicite_partenaire','charge_pos_partenaire','charge_neg_partenaire','aromaticite_partenaire']  # conservation = DESCRIPTEUR (fleche) ; couleur = homo/hetero -> non circulaire
Fi=Fi.dropna(subset=_featI); _homo=Fi['homo']>0.5
print(f"{len(Fi)} clusters d'interface C70 ({int(_homo.sum())} homo / {int((~_homo).sum())} hetero)")

Xi=Fi[_featI].values.astype(float); Xzi=(Xi-Xi.mean(0))/(Xi.std(0)+1e-9)
Xci=Xzi-Xzi.mean(0); Ui,Si,Vti=np.linalg.svd(Xci,full_matrices=False); pcsi=Ui[:,:2]*Si[:2]; evi=(Si**2)/(Si**2).sum()*100
print(f"Variance PC1+PC2 = {evi[0]:.0f}% + {evi[1]:.0f}% = {evi[0]+evi[1]:.0f}%")
_Di=squareform(pdist(Xzi)); _lab=_homo.values.astype(int)
def _sili(D,lab):
    n=len(lab);s=np.zeros(n)
    for i in range(n):
        sm=(lab==lab[i]).copy();sm[i]=False
        a=D[i,sm].mean() if sm.any() else 0; b=min((D[i,lab==c].mean() for c in set(lab) if c!=lab[i]),default=0)
        s[i]=0 if max(a,b)==0 else (b-a)/max(a,b)
    return s.mean()
print(f"Silhouette separation homo vs hetero = {_sili(_Di,_lab):.3f}  (~0 -> interfaces homo non distinctes des hetero)")

_corI=np.array([[np.corrcoef(Xzi[:,j],pcsi[:,k])[0,1] for k in range(2)] for j in range(len(_featI))])
from matplotlib.lines import Line2D
_npc=_ag['n'].reindex(Fi.index).values.astype(float)            # nb d'interactions
_sz=30+(_npc/np.nanmax(_npc))*430                                # taille proportionnelle
_he=~_homo.values; _ho=_homo.values
_chet='#56B4E9'; _chomo='#D55E00'                                # couleur = type (Okabe-Ito)

fig,ax=plt.subplots(figsize=(11,8))
ax.scatter(pcsi[_he,0],pcsi[_he,1],c=_chet,s=_sz[_he],marker='o',edgecolor='k',lw=.3,alpha=.55,zorder=2)   # hetero
ax.scatter(pcsi[_ho,0],pcsi[_ho,1],c=_chomo,s=_sz[_ho],marker='o',edgecolor='k',lw=.5,alpha=.85,zorder=4)  # homo (meme marqueur)
_fI=0.8*np.abs(pcsi).max()/np.abs(_corI).max()
for j,nm in enumerate(_featI):
    _isc = (nm=='conservation_actine')                          # fleches uniformes, gris fonce
    ax.arrow(0,0,_corI[j,0]*_fI,_corI[j,1]*_fI,color='#555',
             head_width=.13,length_includes_head=True,lw=1.1,alpha=.8,zorder=3)
    ax.annotate({'conservation':'sensibilité mutationnelle','enfouissement moyen (%ASA)':'engagement interaction (%ASA)'}.get(nm,nm),(_corI[j,0]*_fI*1.06,_corI[j,1]*_fI*1.06),fontsize=(8 if _isc else 7),
                color='#333',fontweight=('bold' if _isc else 'normal'),zorder=7)
ax.legend(handles=[Line2D([0],[0],marker='o',ls='',mfc=_chet,mec='k',ms=9,label=f'hetero (ABP) n={int(_he.sum())}'),
                   Line2D([0],[0],marker='o',ls='',mfc=_chomo,mec='k',ms=9,label=f'homo (filament) n={int(_ho.sum())}')],
          loc='best',fontsize=9,title='couleur = type')
ax.axhline(0,color='#eee',lw=.6); ax.axvline(0,color='#eee',lw=.6)
ax.set_xlabel(f'PC1 ({evi[0]:.0f}%)'); ax.set_ylabel(f'PC2 ({evi[1]:.0f}%)')
ax.set_title("PCA des interfaces C70 — couleur = homo/hetero, taille = nb d'interactions\n"
             "conservation = DESCRIPTEUR (fleche, label en gras) -> pointe vers asa_moyen/taille = enfoui/grand ; homo et hetero se melangent")

# --- noms selectifs : 2 plus grosses interfaces filament + 1 label par ABP excentre ---
_dist2=np.hypot(pcsi[:,0],pcsi[:,1]); _idx2=list(Fi.index)
for i in np.where(_ho)[0][np.argsort(-_npc[np.where(_ho)[0]])][:2]:
    ax.annotate(f'{_idx2[i]}',(pcsi[i,0],pcsi[i,1]),fontsize=8,fontweight='bold',color='#5a1a00',
                xytext=(7,7),textcoords='offset points',zorder=9,
                bbox=dict(boxstyle='round,pad=.15',fc='white',ec='none',alpha=.8))
# (noms de proteines ABP retires pour l'instant)
plt.tight_layout(); plt.show()


# ### Variante — interface C70, COEUR (résidus > 25 % ASA)
# 
# Même PCA mais descripteurs calculés sur les seuls résidus enfouis > 25 % ASA (côtés actine et partenaire).

# In[90]:


# === PCA des interfaces au niveau cluster C70 (homo + hetero) ===
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path as _P
from scipy.spatial.distance import pdist, squareform
_R = _P('..') if _P('../data').exists() else _P('.')
_df3=pd.read_csv(_R/'data/filtered/details/3.interface_residues.csv')
_df3['pct']=pd.to_numeric(_df3['buried_ASA_percent'].astype(str).str.replace('%','',regex=False),errors='coerce')
_df3['a2']=pd.to_numeric(_df3['buried_ASA_Å²'].astype(str).str.replace('<0.1','0.05',regex=False),errors='coerce')
_df3['canon']=pd.to_numeric(_df3['residue_number_canon_mafft'],errors='coerce')
_df3=_df3[_df3.canon.notna()&_df3.pct.notna()].copy()
_di=pd.read_csv(_R/'data/filtered/details/1.interactions.csv')
_da=pd.read_csv(_R/'data/filtered/filtered_all_data.csv',low_memory=False)
_m=_di.merge(_da[['subunit_1','subunit_2','subunit_2_title','s2_actine','cluster_data_70']].drop_duplicates(),
             left_on=['chain_A_id','chain_B_id'],right_on=['subunit_1','subunit_2'],how='left').drop_duplicates('interaction_id')
_m=_m.dropna(subset=['cluster_data_70'])
_m['homo']=_m['s2_actine'].fillna(False).astype(bool)
_m['abp_title']=_m['subunit_2_title'].fillna('?').astype(str).str.replace(r'\s*\(.*?\)','',regex=True).str.strip().str[:22]
_clu2abp=_m[~_m['homo']].groupby('cluster_data_70')['abp_title'].agg(lambda x:x.value_counts().index[0])
_m['interface_area']=pd.to_numeric(_m['interface_area'].astype(str).str.replace('Å²','',regex=False).str.strip(),errors='coerce')
for _c in ['num_hbonds','num_salt_bridges','num_contacts']: _m[_c]=pd.to_numeric(_m[_c],errors='coerce')
_imap=_m.set_index('interaction_id'); _chA=_imap['chain_A_id'].str.lower(); _clu=_imap['cluster_data_70']
_d=_df3[_df3.interaction_id.isin(set(_m.interaction_id))].copy()
_d['_a']=_d.interaction_id.map(_chA); _d=_d[(_d.chain.str.lower()==_d._a)&(_d['pct']>25)].copy()   # cote actine COEUR >25% ASA
_d['clu']=_d.interaction_id.map(_clu); _d['canon']=_d.canon.astype(int)
_cdf=pd.read_csv(_R/'data/proteocast/conservation_vs_asa_per_position.csv'); _cons=dict(zip(_cdf['canon'].astype(int),_cdf['conservation']))
_clsmap=dict(zip(_cdf['canon'].astype(int),_cdf['residue_class']))   # canon -> sensitive/tolerant
_KD={'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
def _md(s): mm=s.mode(); return mm.iloc[0] if len(mm) else None
_po=_d.groupby(['clu','canon']).agg(w=('pct','mean'),res=('residue_name',_md)).reset_index()
_po['kd']=_po['res'].map(lambda r:_KD.get(str(r),0.0)); _po['cons']=_po['canon'].map(lambda c:_cons.get(int(c),np.nan))
_po['ispos']=_po['res'].isin({'R','K'}).astype(float); _po['isneg']=_po['res'].isin({'D','E'}).astype(float); _po['isaro']=_po['res'].isin({'F','Y','W'}).astype(float)
def _wmi(g,c):
    w=g['w'].values; x=g[c].values; mk=~np.isnan(x)
    return (w[mk]*x[mk]).sum()/w[mk].sum() if w[mk].sum()>0 else np.nan
Fi=_po.groupby('clu').apply(lambda g:pd.Series({'conservation_actine':_wmi(g,'cons'),'hydrophobicite_actine':_wmi(g,'kd'),
    'charge_pos_actine':_wmi(g,'ispos'),'charge_neg_actine':_wmi(g,'isneg'),'aromaticite_actine':_wmi(g,'isaro'),'asa_moyen_actine':_wmi(g,'w')}))
Fi['n_residus']=_d.groupby('interaction_id').size().groupby(_clu).mean()
_ag=_m.groupby('cluster_data_70').agg(hb=('num_hbonds','sum'),sb=('num_salt_bridges','sum'),
    ct=('num_contacts','sum'),n=('interaction_id','nunique'),homo=('homo','mean'))
Fi['taille_interface_actine']=_d.groupby('interaction_id')['a2'].sum().groupby(_clu).mean()
Fi=Fi.join(_ag[['homo']])
_cc=pd.read_csv(_R/'data/filtered/details/4.inter-residue_contacts.csv'); _CHG={'R':1,'K':1,'D':-1,'E':-1}
_cc=_cc[_cc.interaction_id.isin(set(_m.interaction_id))].copy(); _cc['_clu']=_cc.interaction_id.map(_clu)
_cc['charge_compl']=-(_cc['residue_A_name'].map(lambda r:_CHG.get(str(r),0))*_cc['residue_B_name'].map(lambda r:_CHG.get(str(r),0)))
Fi['charge_compl']=_cc.groupby('_clu')['charge_compl'].mean()
Fi['pct_hbond']=100*_ag['hb']/_ag['ct']; Fi['pct_saltbridge']=100*_ag['sb']/_ag['ct']
# --- descripteurs cote PARTENAIRE (chaine B = S2 ; actine pour homo, ABP pour hetero) ---
_chB=_imap['chain_B_id'].str.lower()
_db=_df3[_df3.interaction_id.isin(set(_m.interaction_id))].copy()
_db['_b']=_db.interaction_id.map(_chB); _db=_db[(_db.chain.str.lower()==_db._b)&(_db['pct']>25)].copy()  # partenaire COEUR >25%
_db['clu']=_db.interaction_id.map(_clu)
_db['kd']=_db['residue_name'].map(lambda r:_KD.get(str(r),0.0))
_db['ispos']=_db['residue_name'].isin({'R','K'}).astype(float)
_db['isneg']=_db['residue_name'].isin({'D','E'}).astype(float)
_db['isaro']=_db['residue_name'].isin({'F','Y','W'}).astype(float)
def _wmb(g,col):
    w=g['pct'].values; x=g[col].values; mk=~np.isnan(x)
    return (w[mk]*x[mk]).sum()/w[mk].sum() if w[mk].sum()>0 else np.nan
Fb=_db.groupby('clu').apply(lambda g:pd.Series({'hydrophobicite_partenaire':_wmb(g,'kd'),
    'charge_pos_partenaire':_wmb(g,'ispos'),'charge_neg_partenaire':_wmb(g,'isneg'),'aromaticite_partenaire':_wmb(g,'isaro'),'asa_moyen_partenaire':g['pct'].mean()}))
Fi=Fi.join(Fb)
_featI=['conservation_actine','charge_compl','pct_hbond','pct_saltbridge','hydrophobicite_partenaire','charge_pos_partenaire','charge_neg_partenaire','aromaticite_partenaire']  # ANGLE COMPLEMENTARITE/PARTENAIRE : on retire la physico-chimie cote actine (redondante avec le PCA des zones) ; charge_compl = complementarite electrostatique sur contacts reels
Fi=Fi.dropna(subset=_featI); _homo=Fi['homo']>0.5
print(f"{len(Fi)} clusters d'interface C70 ({int(_homo.sum())} homo / {int((~_homo).sum())} hetero)")

Xi=Fi[_featI].values.astype(float); Xzi=(Xi-Xi.mean(0))/(Xi.std(0)+1e-9)
Xci=Xzi-Xzi.mean(0); Ui,Si,Vti=np.linalg.svd(Xci,full_matrices=False); pcsi=Ui[:,:2]*Si[:2]; evi=(Si**2)/(Si**2).sum()*100
print(f"Variance PC1+PC2 = {evi[0]:.0f}% + {evi[1]:.0f}% = {evi[0]+evi[1]:.0f}%")
_Di=squareform(pdist(Xzi)); _lab=_homo.values.astype(int)
def _sili(D,lab):
    n=len(lab);s=np.zeros(n)
    for i in range(n):
        sm=(lab==lab[i]).copy();sm[i]=False
        a=D[i,sm].mean() if sm.any() else 0; b=min((D[i,lab==c].mean() for c in set(lab) if c!=lab[i]),default=0)
        s[i]=0 if max(a,b)==0 else (b-a)/max(a,b)
    return s.mean()
print(f"Silhouette separation homo vs hetero = {_sili(_Di,_lab):.3f}  (~0 -> interfaces homo non distinctes des hetero)")

_corI=np.array([[np.corrcoef(Xzi[:,j],pcsi[:,k])[0,1] for k in range(2)] for j in range(len(_featI))])
from matplotlib.lines import Line2D
_npc=_ag['n'].reindex(Fi.index).values.astype(float)            # nb d'interactions
_sz=30+(_npc/np.nanmax(_npc))*430                                # taille proportionnelle
_he=~_homo.values; _ho=_homo.values
_chet='#56B4E9'; _chomo='#D55E00'                                # couleur = type (Okabe-Ito)

fig,ax=plt.subplots(figsize=(11,8))
ax.scatter(pcsi[_he,0],pcsi[_he,1],c=_chet,s=_sz[_he],marker='o',edgecolor='k',lw=.3,alpha=.55,zorder=2)   # hetero
ax.scatter(pcsi[_ho,0],pcsi[_ho,1],c=_chomo,s=_sz[_ho],marker='o',edgecolor='k',lw=.5,alpha=.85,zorder=4)  # homo (meme marqueur)
_fI=0.8*np.abs(pcsi).max()/np.abs(_corI).max()
for j,nm in enumerate(_featI):
    _isc = (nm=='conservation_actine')                          # fleches uniformes, gris fonce
    ax.arrow(0,0,_corI[j,0]*_fI,_corI[j,1]*_fI,color='#555',
             head_width=.13,length_includes_head=True,lw=1.1,alpha=.8,zorder=3)
    ax.annotate({'conservation':'sensibilité mutationnelle','enfouissement moyen (%ASA)':'engagement interaction (%ASA)'}.get(nm,nm),(_corI[j,0]*_fI*1.06,_corI[j,1]*_fI*1.06),fontsize=(8 if _isc else 7),
                color='#333',fontweight=('bold' if _isc else 'normal'),zorder=7)
ax.legend(handles=[Line2D([0],[0],marker='o',ls='',mfc=_chet,mec='k',ms=9,label=f'hetero (ABP) n={int(_he.sum())}'),
                   Line2D([0],[0],marker='o',ls='',mfc=_chomo,mec='k',ms=9,label=f'homo (filament) n={int(_ho.sum())}')],
          loc='best',fontsize=9,title='couleur = type')
ax.axhline(0,color='#eee',lw=.6); ax.axvline(0,color='#eee',lw=.6)
ax.set_xlabel(f'PC1 ({evi[0]:.0f}%)'); ax.set_ylabel(f'PC2 ({evi[1]:.0f}%)')
ax.set_title("PCA des interfaces C70 (COEUR >25% ASA) — ANGLE COMPLEMENTARITE / PARTENAIRE (complementaire du PCA des zones)\n"
             "descripteurs = chimie du contact (H-bond, pont salin, complementarite de charge) + physico-chimie du PARTENAIRE ; sensibilité mutationnelle = fleche (gras)")

# --- noms selectifs : 2 plus grosses interfaces filament + 1 label par ABP excentre ---
_dist2=np.hypot(pcsi[:,0],pcsi[:,1]); _idx2=list(Fi.index)
for i in np.where(_ho)[0][np.argsort(-_npc[np.where(_ho)[0]])][:2]:
    ax.annotate(f'{_idx2[i]}',(pcsi[i,0],pcsi[i,1]),fontsize=8,fontweight='bold',color='#5a1a00',
                xytext=(7,7),textcoords='offset points',zorder=9,
                bbox=dict(boxstyle='round,pad=.15',fc='white',ec='none',alpha=.8))
# (noms de proteines ABP retires pour l'instant)
plt.tight_layout(); plt.show()


# ### PCA cluster binding-site — zones de liaison de l'actine
# 
# Un point = un **cluster de binding-site** (`s1_binding_site_cluster_data_70`, ex. 6685_0) = une **zone de l'actine** pouvant être liée par plusieurs ABP. Complète la PCA C70 (qui était par interface).

# In[91]:


# === PCA des clusters de BINDING-SITE de l'actine — TOUS les clusters ===
# Un point = une ZONE de liaison (s1_binding_site_cluster_data_70). Conservation = DESCRIPTEUR (fleche) ;
# couleur = homo (filament) / hetero (ABP) / mixte (les deux). Taille = nb d'interactions.
from pathlib import Path as _Pb
_Rb=_Pb('..') if (_Pb('..')/'data').exists() else _Pb('.')
_KDb={'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
_HYD=set('AVLIMFWC')   # residus hydrophobes
_dab=pd.read_csv(_Rb/'data/filtered/filtered_all_data.csv',low_memory=False)
_dib=pd.read_csv(_Rb/'data/filtered/details/1.interactions.csv')
_f3b=pd.read_csv(_Rb/'data/filtered/details/3.interface_residues.csv')
_cvb=pd.read_csv(_Rb/'data/proteocast/conservation_vs_asa_per_position.csv').set_index('canon')['conservation']
_mb=_dib.merge(_dab[['subunit_1','subunit_2','s2_actine','s1_binding_site_cluster_data_70','s2_sequence_cluster_70']].drop_duplicates(),
               left_on=['chain_A_id','chain_B_id'],right_on=['subunit_1','subunit_2'],how='left')
_mb['s2_actine']=_mb['s2_actine'].fillna(False).astype(bool)
_mb=_mb.dropna(subset=['s1_binding_site_cluster_data_70'])
# classe : homo (que filament) / hetero (que ABP) / mixte (les deux)
_fh=_mb.groupby('s1_binding_site_cluster_data_70')['s2_actine'].mean()
_zcls=_fh.map(lambda f:'homo' if f==1 else ('hetero' if f==0 else 'mixte'))
# nb d'ABP distincts (cluster sequence 70%) par zone, cote hetero
_nabp=_mb[~_mb.s2_actine].groupby('s1_binding_site_cluster_data_70')['s2_sequence_cluster_70'].nunique()
_iid2bs=_mb.set_index('interaction_id')['s1_binding_site_cluster_data_70']
_iid2act=_mb.set_index('interaction_id')['chain_A_id'].str.lower()
_f3b['pct']=pd.to_numeric(_f3b['buried_ASA_percent'].astype(str).str.replace('%','',regex=False),errors='coerce')
_f3b['canon']=pd.to_numeric(_f3b['residue_number_canon_mafft'],errors='coerce')
_f3b['aa']=_f3b['residue_name'].astype(str).str.upper().str.strip()
_f3b['bs']=_f3b.interaction_id.map(_iid2bs); _f3b['act']=_f3b.interaction_id.map(_iid2act)
_aB=_f3b[(_f3b.chain.str.lower()==_f3b.act)&_f3b.bs.notna()&(_f3b.aa.str.len()==1)].copy()
_aB['kd']=_aB.aa.map(_KDb); _aB['hyd']=_aB.aa.isin(_HYD); _aB['pos']=_aB.aa.isin(list('KR')); _aB['neg']=_aB.aa.isin(list('DE'))
_aB['aro']=_aB.aa.isin(list('FWY')); _aB['cons']=_aB.canon.map(_cvb)
_gB=_aB.groupby('bs')
_descB=pd.DataFrame({'conservation':_gB['cons'].mean(),'hydrophobe':_gB['hyd'].mean(),'charge +':_gB['pos'].mean(),
   'charge -':_gB['neg'].mean(),'aromatique':_gB['aro'].mean(),'enfoui (ASA%)':_gB['pct'].mean(),
   'n_intx':_gB['interaction_id'].nunique()})
_descB['n_abp']=_nabp.reindex(_descB.index).fillna(0)
_descB['classe']=_zcls.reindex(_descB.index)
_bf=['conservation','hydrophobe','charge +','charge -','aromatique','enfoui (ASA%)','n_abp']
_descB=_descB.dropna(subset=_bf)
_XB=_descB[_bf].values; _XBz=(_XB-_XB.mean(0))/(_XB.std(0)+1e-9)
_UB,_SB,_VtB=np.linalg.svd(_XBz-_XBz.mean(0),full_matrices=False); _pcB=_UB[:,:2]*_SB[:2]; _evB=(_SB**2)/np.sum(_SB**2)
_corB=np.array([[np.corrcoef(_XBz[:,j],_pcB[:,k])[0,1] for k in range(2)] for j in range(len(_bf))])

from matplotlib.lines import Line2D
_cc={'hetero':'#56B4E9','homo':'#D55E00','mixte':'#009E73'}
_cols=_descB['classe'].map(_cc).values
_sz=25+(_descB['n_intx'].values/_descB['n_intx'].max())*430
fig,ax=plt.subplots(figsize=(9.5,7.5))
ax.scatter(_pcB[:,0],_pcB[:,1],c=_cols,s=_sz,edgecolor='k',lw=.35,alpha=.8)
_fb=0.8*np.abs(_pcB).max()/np.abs(_corB).max()
for j,nm in enumerate(_bf):
    _isc=(nm=='conservation')
    ax.arrow(0,0,_corB[j,0]*_fb,_corB[j,1]*_fb,color='#555',head_width=.07,lw=1.1,length_includes_head=True,alpha=.8,zorder=3)
    ax.annotate({'conservation':'sensibilité mutationnelle','enfouissement moyen (%ASA)':'engagement interaction (%ASA)'}.get(nm,nm),(_corB[j,0]*_fb*1.08,_corB[j,1]*_fb*1.08),fontsize=(8.5 if _isc else 7.5),
                color='#333',fontweight=('bold' if _isc else 'normal'),zorder=6)
# noms des 4 clusters principaux (interfaces filament : 6685_1..4)
_idxB=list(_descB.index)
for b in ['6685_1','6685_2','6685_3','6685_4']:
    if b in _idxB:
        i=_idxB.index(b)
        ax.annotate(b,(_pcB[i,0],_pcB[i,1]),fontsize=8,fontweight='bold',color='#7a2e00',
                    xytext=(6,5),textcoords='offset points',zorder=7,
                    bbox=dict(boxstyle='round,pad=.13',fc='white',ec='none',alpha=.8))
ax.legend(handles=[Line2D([0],[0],marker='o',ls='',mfc=_cc[g],mec='k',ms=9,
                          label=f'{g} (n={int((_descB.classe==g).sum())})') for g in ['hetero','homo','mixte'] if (_descB.classe==g).any()],
          loc='best',fontsize=9,title='couleur = type de zone')
ax.axhline(0,color='#eee',lw=.6); ax.axvline(0,color='#eee',lw=.6)
ax.set_xlabel(f'PC1 ({_evB[0]*100:.0f}%)'); ax.set_ylabel(f'PC2 ({_evB[1]*100:.0f}%)')
ax.set_title(f"PCA des clusters de BINDING-SITE de l'actine — TOUS ({len(_descB)} zones)\n"
             "sensibilité mutationnelle = descripteur ; couleur = homo/hetero/mixte ; taille = nb d'interactions")
plt.tight_layout(); plt.show()
print(f"{len(_descB)} zones | classes: {_descB['classe'].value_counts().to_dict()} | var PC1+PC2={(_evB[0]+_evB[1])*100:.0f}%")


# ### Variante — binding-site, clusters contenant une **myosine** en évidence
# 
# Même PCA que ci-dessus (mêmes coordonnées, mêmes flèches), mais les clusters de binding-site dont au moins un partenaire est une **myosine** (hors tropomyosine) sont colorés en **rose**. *(À exécuter après la cellule PCA binding-site ci-dessus : réutilise `_descB`, `_pcB`, `_evB`, `_corB`, `_bf`, `_sz`, `_mb`, `_dab`.)*

# In[92]:


# === PCA binding-site — MEME PCA, coloree selon la presence d'une MYOSINE ===
from matplotlib.lines import Line2D as _L2
_t2b = _dab[['subunit_2','subunit_2_title']].drop_duplicates()
_mbm = _mb.merge(_t2b, on='subunit_2', how='left')
_is_myo = lambda t: ('myosin' in str(t).lower()) and ('tropomyosin' not in str(t).lower())
_mbm['is_myo'] = _mbm['subunit_2_title'].map(_is_myo)
_has_myo = _mbm.groupby('s1_binding_site_cluster_data_70')['is_myo'].any()
_myo = _has_myo.reindex(_descB.index).fillna(False).values
# nom(s) de myosine par cluster (pour les labels)
_myo_names = (_mbm[_mbm['is_myo']].groupby('s1_binding_site_cluster_data_70')['subunit_2_title']
              .apply(lambda x: ' / '.join(sorted(set(x)))).to_dict())

fig,ax=plt.subplots(figsize=(10.5,8))
ax.scatter(_pcB[~_myo,0],_pcB[~_myo,1],c='#cfd8dc',s=_sz[~_myo],edgecolor='k',lw=.3,alpha=.55)
ax.scatter(_pcB[_myo,0],_pcB[_myo,1],c='#d81b60',s=_sz[_myo],edgecolor='k',lw=.6,alpha=.95,zorder=4)
_fb=0.8*np.abs(_pcB).max()/np.abs(_corB).max()
for j,nm in enumerate(_bf):
    ax.arrow(0,0,_corB[j,0]*_fb,_corB[j,1]*_fb,color='#888',head_width=.07,lw=1,length_includes_head=True,alpha=.7,zorder=3)
    ax.annotate({'conservation':'sensibilité mutationnelle','enfouissement moyen (%ASA)':'engagement interaction (%ASA)'}.get(nm,nm),(_corB[j,0]*_fb*1.08,_corB[j,1]*_fb*1.08),fontsize=7.5,color='#333',zorder=6)
_idxB=list(_descB.index)
_txt=[]
for i2,b in enumerate(_idxB):
    if _myo[i2]:
        _lbl=f"{b}\n{_myo_names.get(b,'')}"
        _txt.append(ax.annotate(_lbl,(_pcB[i2,0],_pcB[i2,1]),fontsize=7,color='#880e4f',fontweight='bold',
                    xytext=(5,4),textcoords='offset points',zorder=7,
                    bbox=dict(boxstyle='round,pad=.12',fc='white',ec='#d81b60',lw=.4,alpha=.85)))
try:
    from adjustText import adjust_text; adjust_text(_txt,ax=ax,only_move={'text':'xy'})
except Exception: pass
ax.legend(handles=[_L2([0],[0],marker='o',ls='',mfc='#d81b60',mec='k',ms=10,label=f'avec myosine (n={int(_myo.sum())})'),
                   _L2([0],[0],marker='o',ls='',mfc='#cfd8dc',mec='k',ms=10,label=f'autre ABP (n={int((~_myo).sum())})')],
          loc='best',fontsize=10,title='couleur')
ax.axhline(0,color='#eee',lw=.6); ax.axvline(0,color='#eee',lw=.6)
ax.set_xlabel(f'PC1 ({_evB[0]*100:.0f}%)'); ax.set_ylabel(f'PC2 ({_evB[1]*100:.0f}%)')
ax.set_title(f"PCA binding-site — clusters AVEC myosine ({int(_myo.sum())}) ; label = cluster + proteines ; taille = nb d'interactions")
plt.tight_layout(); plt.show()
for b in _idxB:
    if b in _myo_names: print(f'{b}: {_myo_names[b]}')


# ### Variante — binding-site, COEUR seulement (résidus > 25 % ASA)
# 
# Même PCA mais les descripteurs de chaque zone ne sont calculés que sur les **résidus réellement enfouis (> 25 % ASA)** — on retire les contacts superficiels (effleurements).

# In[93]:


# === PCA des clusters de BINDING-SITE de l'actine — COEUR (residus >25% ASA) ===
# Un point = une ZONE de liaison (s1_binding_site_cluster_data_70). Conservation = DESCRIPTEUR (fleche) ;
# couleur = homo (filament) / hetero (ABP) / mixte (les deux). Taille = nb d'interactions.
from pathlib import Path as _Pb
_Rb=_Pb('..') if (_Pb('..')/'data').exists() else _Pb('.')
_KDb={'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
_HYD=set('AVLIMFWC')   # residus hydrophobes
_dabcr=pd.read_csv(_Rb/'data/filtered/filtered_all_data.csv',low_memory=False)
_dibcr=pd.read_csv(_Rb/'data/filtered/details/1.interactions.csv')
_f3bcr=pd.read_csv(_Rb/'data/filtered/details/3.interface_residues.csv')
_cvbcr=pd.read_csv(_Rb/'data/proteocast/conservation_vs_asa_per_position.csv').set_index('canon')['conservation']
_mbcr=_dibcr.merge(_dabcr[['subunit_1','subunit_2','s2_actine','s1_binding_site_cluster_data_70','s2_sequence_cluster_70']].drop_duplicates(),
               left_on=['chain_A_id','chain_B_id'],right_on=['subunit_1','subunit_2'],how='left')
_mbcr['s2_actine']=_mbcr['s2_actine'].fillna(False).astype(bool)
_mbcr=_mbcr.dropna(subset=['s1_binding_site_cluster_data_70'])
# classe : homo (que filament) / hetero (que ABP) / mixte (les deux)
_fhcr=_mbcr.groupby('s1_binding_site_cluster_data_70')['s2_actine'].mean()
_zclscr=_fhcr.map(lambda f:'homo' if f==1 else ('hetero' if f==0 else 'mixte'))
# nb d'ABP distincts (cluster sequence 70%) par zone, cote hetero
_nabpcr=_mbcr[~_mbcr.s2_actine].groupby('s1_binding_site_cluster_data_70')['s2_sequence_cluster_70'].nunique()
_iid2bscr=_mbcr.set_index('interaction_id')['s1_binding_site_cluster_data_70']
_iid2actcr=_mbcr.set_index('interaction_id')['chain_A_id'].str.lower()
_f3bcr['pct']=pd.to_numeric(_f3bcr['buried_ASA_percent'].astype(str).str.replace('%','',regex=False),errors='coerce')
_f3bcr['canon']=pd.to_numeric(_f3bcr['residue_number_canon_mafft'],errors='coerce')
_f3bcr['aa']=_f3bcr['residue_name'].astype(str).str.upper().str.strip()
_f3bcr['bs']=_f3bcr.interaction_id.map(_iid2bscr); _f3bcr['act']=_f3bcr.interaction_id.map(_iid2actcr)
_aBcr=_f3bcr[(_f3bcr.chain.str.lower()==_f3bcr.act)&_f3bcr.bs.notna()&(_f3bcr.aa.str.len()==1)&(_f3bcr['pct']>25)].copy()  # COEUR : residus enfouis >25% ASA seulement
_aBcr['kd']=_aBcr.aa.map(_KDb); _aBcr['hyd']=_aBcr.aa.isin(_HYD); _aBcr['pos']=_aBcr.aa.isin(list('KR')); _aBcr['neg']=_aBcr.aa.isin(list('DE'))
_aBcr['aro']=_aBcr.aa.isin(list('FWY')); _aBcr['cons']=_aBcr.canon.map(_cvbcr)
_posBcr=_aBcr.sort_values('pct').drop_duplicates(['bs','canon'],keep='last')   # 1 ligne/position = l'occurrence la + enfouie (ASA max) -> moyenne PAR POSITION, pas par interaction
_gBcr=_posBcr.groupby('bs')
_descBcr=pd.DataFrame({'conservation':_gBcr['cons'].mean(),'%_hydrophobe':_gBcr['hyd'].mean(),'%_charge +':_gBcr['pos'].mean(),
   '%_charge -':_gBcr['neg'].mean(),'%_aromatique':_gBcr['aro'].mean(),'enfouissement moyen (%ASA)':_gBcr['pct'].mean(),
   'n_intx':_aBcr.groupby('bs')['interaction_id'].nunique()})   # taille = nb d'interactions, sur TOUTES les occurrences (pas dedupe)
_descBcr['n_abp']=_nabpcr.reindex(_descBcr.index).fillna(0)
_descBcr['n_res']=_gBcr.size().reindex(_descBcr.index)   # taille empreinte = nb de positions coeur distinctes (apres dedup par position)
_descBcr['classe']=_zclscr.reindex(_descBcr.index)
_bfcr=['conservation','%_hydrophobe','%_charge +','%_charge -','%_aromatique','enfouissement moyen (%ASA)','n_res','n_abp']
_descBcr=_descBcr.dropna(subset=_bfcr)
_XBcr=_descBcr[_bfcr].values; _XBcrz=(_XBcr-_XBcr.mean(0))/(_XBcr.std(0)+1e-9)
_UBcr,_SBcr,_VtBcr=np.linalg.svd(_XBcrz-_XBcrz.mean(0),full_matrices=False); _pcBcr=_UBcr[:,:2]*_SBcr[:2]; _evBcr=(_SBcr**2)/np.sum(_SBcr**2)
_corBcr=np.array([[np.corrcoef(_XBcrz[:,j],_pcBcr[:,k])[0,1] for k in range(2)] for j in range(len(_bfcr))])

from matplotlib.lines import Line2D
_cccr={'hetero':'#56B4E9','homo':'#D55E00','mixte':'#009E73'}
_colscr=_descBcr['classe'].map(_cccr).values
_szcr=110+(_descBcr['n_intx'].values/_descBcr['n_intx'].max())*950
fig,ax=plt.subplots(figsize=(13,10.5))
ax.scatter(_pcBcr[:,0],_pcBcr[:,1],c=_colscr,s=_szcr,edgecolor='k',lw=.5,alpha=.9,zorder=5)
_fbcr=0.8*np.abs(_pcBcr).max()/np.abs(_corBcr).max()
for j,nm in enumerate(_bfcr):
    _isc=(nm=='conservation')
    ax.arrow(0,0,_corBcr[j,0]*_fbcr,_corBcr[j,1]*_fbcr,color='#444',head_width=.14,lw=2.3,length_includes_head=True,alpha=.85,zorder=3)
    ax.annotate({'conservation':'sensibilité mutationnelle','enfouissement moyen (%ASA)':'engagement interaction (%ASA)'}.get(nm,nm),(_corBcr[j,0]*_fbcr*1.10,_corBcr[j,1]*_fbcr*1.10),fontsize=(15 if _isc else 13),
                color='#333',fontweight=('bold' if _isc else 'normal'),zorder=6)
# === noms des clusters : code + ABP(s) presents (familles courtes), sans fond blanc ===
import re as _re2, textwrap as _tw
def _shortabp(n):                                  # nom de famille court (cas speciaux + 1er mot)
    n=str(n).strip()
    if 'complex subunit 5-like' in n: return 'ARPC5L'
    m=_re2.search(r'complex subunit (\d+[A-Za-z]?)',n)
    if m: return 'ARPC'+m.group(1)
    if _re2.match(r'Actin-related protein 3',n): return 'Arp3'
    if _re2.match(r'Actin-related protein 2',n): return 'Arp2'
    if 'formin' in n.lower(): return _re2.sub(r'-\d+$','',n)        # Inverted formin-2 -> Inverted formin
    if 'tropomyosin' in n.lower(): return 'Tropomyosin'
    if 'unconventional' in n.lower() and 'myosin' in n.lower(): return 'Myosin A'   # non-conv. -> distincte
    if 'myosin' in n.lower(): return 'Myosin'                       # toutes les autres myosines
    if 'capping' in n.lower(): return 'Capping'                     # F-actin-capping... -> Capping
    if 'filamin' in n.lower(): return 'Filamin'                     # Filamin-A -> Filamin
    if 'diaphanous' in n.lower(): return 'Diaphanous-'+(_re2.search(r'homolog (\d+)',n).group(1) if _re2.search(r'homolog (\d+)',n) else '1')
    if 'catenin' in n.lower(): return 'Catenin'                     # Catenin alpha-1 = Alpha-catenin-like
    if _re2.match(r'Epididymis secretory protein Li',n): return 'ESP-Li '+(_re2.search(r'(\d+)$',n).group(1) if _re2.search(r'(\d+)$',n) else '')
    w=_re2.split(r'[ ,]',n)[0]
    return _re2.sub(r'-\d+$','',w)
_nmsrc=_dabcr[~_dabcr['s2_actine'].fillna(False)].copy()
_nmsrc['t']=_nmsrc['subunit_2_title'].fillna('').str.replace(r'\s*\(.*','',regex=True).str.strip().map(_shortabp)
_z2fams=_nmsrc[_nmsrc['t']!=''].groupby('s1_binding_site_cluster_data_70')['t'].agg(lambda x:[f for f,_ in x.value_counts().items()])
def _zlabel(bs):                                   # homo -> code seul ; sinon code + familles d'ABP (1/ligne, max 5)
    if (bs not in _descBcr.index) or (_descBcr.loc[bs,'classe']=='homo'): return bs
    _f=list(_z2fams.get(bs,[]))
    if 'Myosin' in _f and 'Myosin A' in _f: _f=[x for x in _f if x!='Myosin A']   # evite Myosin + Myosin A
    if not _f: return bs
    return bs+'\n'+'\n'.join(_f[:5])+(f'\n+{len(_f)-5}' if len(_f)>5 else '')
_idxBcr=list(_descBcr.index)
_rad=np.hypot(_pcBcr[:,0],_pcBcr[:,1])
_RLAB=np.percentile(_rad,88)   # <<< seuil de rayon : baisse -> plus de labels, monte -> moins
_top7=set(_descBcr['n_intx'].sort_values(ascending=False).head(7).index)
_main=['6685_1','6685_2','6685_3','6685_4']
_txt=[]
for b in _main:                                        # 4 clusters principaux (gras)
    if b in _idxBcr:
        i=_idxBcr.index(b)
        _txt.append(ax.text(_pcBcr[i,0],_pcBcr[i,1],_zlabel(b),fontsize=11,fontweight='bold',color='#7a2e00',zorder=7))
for k2,bs in enumerate(_idxBcr):                       # peripherie (rayon) + 7 plus liees
    if bs in _main: continue
    if (_rad[k2]>=_RLAB) or (bs in _top7):
        _txt.append(ax.text(_pcBcr[k2,0],_pcBcr[k2,1],_zlabel(bs),fontsize=8,color='#222',zorder=7))
from adjustText import adjust_text                      # anti-chevauchement des labels de clusters uniquement
adjust_text(_txt,ax=ax,arrowprops=dict(arrowstyle='-',lw=0.5,color='#999'),
            expand_text=(1.05,1.25),expand_points=(1.02,1.05),force_text=(0.4,0.6),force_points=(0.1,0.2),lim=400)
_lgc=ax.legend(handles=[Line2D([0],[0],marker='o',ls='',mfc=_cccr[g],mec='k',ms=12,
                        label=f'{g} (n={int((_descBcr.classe==g).sum())})') for g in ['hetero','homo','mixte'] if (_descBcr.classe==g).any()],
               loc='upper left',fontsize=12,title='couleur = type de zone',title_fontsize=12)
ax.add_artist(_lgc)   # garder la legende couleur quand on ajoute celle de la taille
# legende TAILLE = nb d'interactions (memes valeurs de reference que le mapping _szcr)
_vmax=_descBcr['n_intx'].max()
_refsz=sorted(set([1,int(_vmax*0.25),int(_vmax*0.5),int(_vmax)]))
_hsz=[Line2D([0],[0],marker='o',ls='',mfc='#bbbbbb',mec='k',ms=np.sqrt(60+(v/_vmax)*750)*0.55,label=str(v)) for v in _refsz]
ax.legend(handles=_hsz,loc='upper left',bbox_to_anchor=(0.0,0.83),fontsize=9,title="taille = nb d'interactions",title_fontsize=10,
          labelspacing=1.3,borderpad=0.7,handletextpad=1.0,framealpha=.9)
ax.axhline(0,color='#eee',lw=.6); ax.axvline(0,color='#eee',lw=.6)
ax.set_xlabel(f'PC1 ({_evBcr[0]*100:.0f}%)',fontsize=13); ax.set_ylabel(f'PC2 ({_evBcr[1]*100:.0f}%)',fontsize=13); ax.tick_params(labelsize=12)
ax.set_title(f"PCA des clusters de BINDING-SITE de l'actine — COEUR >25% ASA ({len(_descBcr)} zones)\n"
             "sensibilité mutationnelle = descripteur ; couleur = homo/hetero/mixte ; taille = nb d'interactions",fontsize=14)
plt.tight_layout(); plt.show()
print(f"{len(_descBcr)} zones | classes: {_descBcr['classe'].value_counts().to_dict()} | var PC1+PC2={(_evBcr[0]+_evBcr[1])*100:.0f}%")


# ### Variante — PCA des clusters d'interface + descripteurs relationnels (liaison H, pont salin)
# 
# Mêmes individus (zones de liaison, cœur >25% ASA) et descripteurs intrinsèques, **plus** deux descripteurs **relationnels** qui n'ont de sens qu'au niveau d'une interface : le **% de contacts en liaison H** et le **% en pont salin**. (Pas de paire aromatique/polaire.)

# In[94]:


# === PCA des clusters de BINDING-SITE — COEUR >25% ASA + descripteurs RELATIONNELS (liaison H, pont salin) ===
# Un point = une ZONE de liaison (s1_binding_site_cluster_data_70). Conservation = DESCRIPTEUR (fleche) ;
# couleur = homo (filament) / hetero (ABP) / mixte (les deux). Taille = nb d'interactions.
from pathlib import Path as _Pb
_Rb=_Pb('..') if (_Pb('..')/'data').exists() else _Pb('.')
_KDb={'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
_HYD=set('AVLIMFWC')   # residus hydrophobes
_dabcr=pd.read_csv(_Rb/'data/filtered/filtered_all_data.csv',low_memory=False)
_dibcr=pd.read_csv(_Rb/'data/filtered/details/1.interactions.csv')
_f3bcr=pd.read_csv(_Rb/'data/filtered/details/3.interface_residues.csv')
_cvbcr=pd.read_csv(_Rb/'data/proteocast/conservation_vs_asa_per_position.csv').set_index('canon')['conservation']
_mbcr=_dibcr.merge(_dabcr[['subunit_1','subunit_2','s2_actine','s1_binding_site_cluster_data_70','s2_sequence_cluster_70']].drop_duplicates(),
               left_on=['chain_A_id','chain_B_id'],right_on=['subunit_1','subunit_2'],how='left')
_mbcr['s2_actine']=_mbcr['s2_actine'].fillna(False).astype(bool)
_mbcr=_mbcr.dropna(subset=['s1_binding_site_cluster_data_70'])
# classe : homo (que filament) / hetero (que ABP) / mixte (les deux)
_fhcr=_mbcr.groupby('s1_binding_site_cluster_data_70')['s2_actine'].mean()
_zclscr=_fhcr.map(lambda f:'homo' if f==1 else ('hetero' if f==0 else 'mixte'))
# nb d'ABP distincts (cluster sequence 70%) par zone, cote hetero
_nabpcr=_mbcr[~_mbcr.s2_actine].groupby('s1_binding_site_cluster_data_70')['s2_sequence_cluster_70'].nunique()
_iid2bscr=_mbcr.set_index('interaction_id')['s1_binding_site_cluster_data_70']
_iid2actcr=_mbcr.set_index('interaction_id')['chain_A_id'].str.lower()
_f3bcr['pct']=pd.to_numeric(_f3bcr['buried_ASA_percent'].astype(str).str.replace('%','',regex=False),errors='coerce')
_f3bcr['canon']=pd.to_numeric(_f3bcr['residue_number_canon_mafft'],errors='coerce')
_f3bcr['aa']=_f3bcr['residue_name'].astype(str).str.upper().str.strip()
_f3bcr['bs']=_f3bcr.interaction_id.map(_iid2bscr); _f3bcr['act']=_f3bcr.interaction_id.map(_iid2actcr)
_aBcr=_f3bcr[(_f3bcr.chain.str.lower()==_f3bcr.act)&_f3bcr.bs.notna()&(_f3bcr.aa.str.len()==1)&(_f3bcr['pct']>25)].copy()  # COEUR : residus enfouis >25% ASA seulement
_aBcr['kd']=_aBcr.aa.map(_KDb); _aBcr['hyd']=_aBcr.aa.isin(_HYD); _aBcr['pos']=_aBcr.aa.isin(list('KR')); _aBcr['neg']=_aBcr.aa.isin(list('DE'))
_aBcr['aro']=_aBcr.aa.isin(list('FWY')); _aBcr['cons']=_aBcr.canon.map(_cvbcr)
_posBcr=_aBcr.sort_values('pct').drop_duplicates(['bs','canon'],keep='last')   # 1 ligne/position = l'occurrence la + enfouie (ASA max) -> moyenne PAR POSITION, pas par interaction
_gBcr=_posBcr.groupby('bs')
_descBcr=pd.DataFrame({'conservation':_gBcr['cons'].mean(),'%_hydrophobe':_gBcr['hyd'].mean(),'%_charge +':_gBcr['pos'].mean(),
   '%_charge -':_gBcr['neg'].mean(),'%_aromatique':_gBcr['aro'].mean(),'enfouissement moyen (%ASA)':_gBcr['pct'].mean(),
   'n_intx':_aBcr.groupby('bs')['interaction_id'].nunique()})   # taille = nb d'interactions, sur TOUTES les occurrences (pas dedupe)
_descBcr['n_abp']=_nabpcr.reindex(_descBcr.index).fillna(0)
_descBcr['n_res']=_gBcr.size().reindex(_descBcr.index)   # taille empreinte = nb de positions coeur distinctes (apres dedup par position)
_descBcr['classe']=_zclscr.reindex(_descBcr.index)
# --- contacts COEUR par binding-site : % liaison H et % pont salin (descripteurs RELATIONNELS) ---
_ccB=pd.read_csv(_Rb/'data/filtered/details/4.inter-residue_contacts.csv',low_memory=False)
def _nB(s): return pd.to_numeric(s.astype(str).str.replace('%','',regex=False).str.replace('Å²','',regex=False).str.strip(),errors='coerce')
_ccB['asa_A']=_nB(_ccB['asa_pct_A']); _ccB['asa_B']=_nB(_ccB['asa_pct_B'])
_ccB['is_hb']=_ccB['contact_type'].astype(str).str.contains('H-bond').astype(int)
_ccB['is_sb']=_ccB['contact_type'].astype(str).str.contains('Salt').astype(int)
_ccB['bs']=_ccB.interaction_id.map(_iid2bscr); _ccB['act']=_ccB.interaction_id.map(_iid2actcr)
_ccB['act_asa']=np.where(_ccB['chain_A_id'].str.lower()==_ccB['act'],_ccB['asa_A'],
                  np.where(_ccB['chain_B_id'].str.lower()==_ccB['act'],_ccB['asa_B'],np.nan))
_ccBc=_ccB[_ccB['bs'].notna()&(_ccB['act_asa']>25)]   # COEUR : residu actine enfoui >25% ASA
_hbsb=_ccBc.groupby('bs').agg(hb=('is_hb','mean'),sb=('is_sb','mean'))*100
_descBcr['% liaison H']=_hbsb['hb'].reindex(_descBcr.index).fillna(0)
_descBcr['% pont salin']=_hbsb['sb'].reindex(_descBcr.index).fillna(0)
_bfcr=['conservation','%_hydrophobe','%_charge +','%_charge -','%_aromatique','% liaison H','% pont salin','enfouissement moyen (%ASA)','n_res','n_abp']
_descBcr=_descBcr.dropna(subset=_bfcr)
_XBcr=_descBcr[_bfcr].values; _XBcrz=(_XBcr-_XBcr.mean(0))/(_XBcr.std(0)+1e-9)
_UBcr,_SBcr,_VtBcr=np.linalg.svd(_XBcrz-_XBcrz.mean(0),full_matrices=False); _pcBcr=_UBcr[:,:2]*_SBcr[:2]; _evBcr=(_SBcr**2)/np.sum(_SBcr**2)
_corBcr=np.array([[np.corrcoef(_XBcrz[:,j],_pcBcr[:,k])[0,1] for k in range(2)] for j in range(len(_bfcr))])

from matplotlib.lines import Line2D
_cccr={'hetero':'#009E73','homo':'#D55E00','mixte':'#F0E442'}
_colscr=_descBcr['classe'].map(_cccr).values
_szcr=110+(_descBcr['n_intx'].values/_descBcr['n_intx'].max())*950
fig,ax=plt.subplots(figsize=(13,10.5))
ax.scatter(_pcBcr[:,0],_pcBcr[:,1],c=_colscr,s=_szcr,edgecolor='k',lw=.5,alpha=.9,zorder=5)
_fbcr=0.8*np.abs(_pcBcr).max()/np.abs(_corBcr).max()
for j,nm in enumerate(_bfcr):
    _isc=(nm=='conservation')
    ax.arrow(0,0,_corBcr[j,0]*_fbcr,_corBcr[j,1]*_fbcr,color='#444',head_width=.14,lw=2.3,length_includes_head=True,alpha=.85,zorder=3)
    ax.annotate({'conservation':'sensibilité mutationnelle','enfouissement moyen (%ASA)':'engagement interaction (%ASA)'}.get(nm,nm),(_corBcr[j,0]*_fbcr*1.10,_corBcr[j,1]*_fbcr*1.10),fontsize=(15 if _isc else 13),
                color='#333',fontweight='bold',zorder=8,bbox=dict(boxstyle='round,pad=0.12',fc='white',ec='none',alpha=0.75))
# === noms des clusters : code + ABP(s) presents (familles courtes), sans fond blanc ===
import re as _re2, textwrap as _tw
def _shortabp(n):                                  # nom de famille court (cas speciaux + 1er mot)
    n=str(n).strip()
    if 'complex subunit 5-like' in n: return 'ARPC5L'
    m=_re2.search(r'complex subunit (\d+[A-Za-z]?)',n)
    if m: return 'ARPC'+m.group(1)
    if _re2.match(r'Actin-related protein 3',n): return 'Arp3'
    if _re2.match(r'Actin-related protein 2',n): return 'Arp2'
    if 'formin' in n.lower(): return _re2.sub(r'-\d+$','',n)        # Inverted formin-2 -> Inverted formin
    if 'tropomyosin' in n.lower(): return 'Tropomyosin'
    if 'unconventional' in n.lower() and 'myosin' in n.lower(): return 'Myosin A'   # non-conv. -> distincte
    if 'myosin' in n.lower(): return 'Myosin'                       # toutes les autres myosines
    if 'capping' in n.lower(): return 'Capping'                     # F-actin-capping... -> Capping
    if 'filamin' in n.lower(): return 'Filamin'                     # Filamin-A -> Filamin
    if 'diaphanous' in n.lower(): return 'Diaphanous-'+(_re2.search(r'homolog (\d+)',n).group(1) if _re2.search(r'homolog (\d+)',n) else '1')
    if 'catenin' in n.lower(): return 'Catenin'                     # Catenin alpha-1 = Alpha-catenin-like
    if _re2.match(r'Epididymis secretory protein Li',n): return 'ESP-Li '+(_re2.search(r'(\d+)$',n).group(1) if _re2.search(r'(\d+)$',n) else '')
    w=_re2.split(r'[ ,]',n)[0]
    return _re2.sub(r'-\d+$','',w)
_nmsrc=_dabcr[~_dabcr['s2_actine'].fillna(False)].copy()
_nmsrc['t']=_nmsrc['subunit_2_title'].fillna('').str.replace(r'\s*\(.*','',regex=True).str.strip().map(_shortabp)
_z2fams=_nmsrc[_nmsrc['t']!=''].groupby('s1_binding_site_cluster_data_70')['t'].agg(lambda x:[f for f,_ in x.value_counts().items()])
def _zlabel(bs):                                   # homo -> code seul ; sinon code + familles d'ABP (1/ligne, max 5)
    if (bs not in _descBcr.index) or (_descBcr.loc[bs,'classe']=='homo'): return bs
    _f=list(_z2fams.get(bs,[]))
    if 'Myosin' in _f and 'Myosin A' in _f: _f=[x for x in _f if x!='Myosin A']   # evite Myosin + Myosin A
    if not _f: return bs
    return bs+'\n'+'\n'.join(_f[:5])+(f'\n+{len(_f)-5}' if len(_f)>5 else '')
_idxBcr=list(_descBcr.index)
_rad=np.hypot(_pcBcr[:,0],_pcBcr[:,1])
_RLAB=np.percentile(_rad,88)   # <<< seuil de rayon : baisse -> plus de labels, monte -> moins
_top7=set(_descBcr['n_intx'].sort_values(ascending=False).head(7).index)
_main=['6685_1','6685_2','6685_3','6685_4']
_txt=[]
for b in _main:                                        # 4 clusters principaux (gras)
    if b in _idxBcr:
        i=_idxBcr.index(b)
        _txt.append(ax.text(_pcBcr[i,0],_pcBcr[i,1],_zlabel(b),fontsize=11,fontweight='bold',color='#7a2e00',zorder=7))
for k2,bs in enumerate(_idxBcr):                       # peripherie (rayon) + 7 plus liees
    if bs in _main: continue
    if (_rad[k2]>=_RLAB) or (bs in _top7):
        _txt.append(ax.text(_pcBcr[k2,0],_pcBcr[k2,1],_zlabel(bs),fontsize=9,color='#222',zorder=7))
from adjustText import adjust_text                      # anti-chevauchement des labels de clusters uniquement
adjust_text(_txt,ax=ax,arrowprops=dict(arrowstyle='-',lw=0.5,color='#999'),
            expand_text=(1.05,1.25),expand_points=(1.02,1.05),force_text=(0.4,0.6),force_points=(0.1,0.2),lim=400)
_hcol=[Line2D([0],[0],marker='o',ls='',mfc=_cccr[g],mec='k',ms=12,
       label=f'{g} (n={int((_descBcr.classe==g).sum())})') for g in ['hetero','homo','mixte'] if (_descBcr.classe==g).any()]
_lgc=ax.legend(handles=_hcol,loc='upper left',fontsize=12,title='couleur = type de zone',title_fontsize=12)
ax.add_artist(_lgc)   # garder la legende couleur quand on ajoute celle de la taille
# legende TAILLE = nb d'interactions (memes valeurs de reference que le mapping _szcr)
_vmax=_descBcr['n_intx'].max()
_refsz=sorted(set([1,int(_vmax*0.25),int(_vmax*0.5),int(_vmax)]))
_hsz=[Line2D([0],[0],marker='o',ls='',mfc='#bbbbbb',mec='k',ms=np.sqrt(60+(v/_vmax)*750)*0.55,label=str(v)) for v in _refsz]
ax.legend(handles=_hsz,loc='upper left',bbox_to_anchor=(0.0,0.83),fontsize=9,title="taille = nb d'interactions",title_fontsize=10,
          labelspacing=1.3,borderpad=0.7,handletextpad=1.0,framealpha=.9)
ax.axhline(0,color='#eee',lw=.6); ax.axvline(0,color='#eee',lw=.6)
ax.set_xlabel(f'PC1 ({_evBcr[0]*100:.0f}%)',fontsize=13); ax.set_ylabel(f'PC2 ({_evBcr[1]*100:.0f}%)',fontsize=13); ax.tick_params(labelsize=12)
ax.set_title(f"PCA des clusters de BINDING-SITE — COEUR >25% ASA + liaison H / pont salin ({len(_descBcr)} zones)\n"
             "sensibilité mutationnelle = descripteur ; couleur = homo/hetero/mixte ; taille = nb d'interactions",fontsize=14)
plt.tight_layout(); plt.show()
print(f"{len(_descBcr)} zones | classes: {_descBcr['classe'].value_counts().to_dict()} | var PC1+PC2={(_evBcr[0]+_evBcr[1])*100:.0f}%")


# ### Variante — binding-site COEUR, clusters à **myosine entourés**
# 
# Même PCA COEUR que ci-dessus (couleurs homo/hétéro/mixte conservées), mais les clusters dont un partenaire est une **myosine** sont **entourés d'un cercle rose** et annotés avec leur **cluster + nom(s) de myosine**. *(À exécuter après la cellule COEUR ci-dessus : réutilise `_pcBcr`, `_descBcr`, `_szcr`, `_corBcr`, `_bfcr`, `_evBcr`, `_idxBcr`, `_mbcr`, `_dabcr`, `_cccr`.)*

# In[95]:


# === PCA binding-site COEUR — clusters AVEC myosine ENTOURES + noms ===
from matplotlib.lines import Line2D
_t2c=_dabcr[['subunit_2','subunit_2_title']].drop_duplicates()
_mc=_mbcr.merge(_t2c,on='subunit_2',how='left')
_ismyo=lambda t:('myosin' in str(t).lower()) and ('tropomyosin' not in str(t).lower())
_mc['is_myo']=_mc['subunit_2_title'].map(_ismyo)
_myoname=(_mc[_mc.is_myo].groupby('s1_binding_site_cluster_data_70')['subunit_2_title']
          .apply(lambda x:' / '.join(sorted(set(x)))).to_dict())
_myoC=np.array([b in _myoname for b in _idxBcr])

fig,ax=plt.subplots(figsize=(13,10.5))
_colsc=_descBcr['classe'].map(_cccr).values
ax.scatter(_pcBcr[:,0],_pcBcr[:,1],c=_colsc,s=_szcr,edgecolor='k',lw=.5,alpha=.55,zorder=2)
# anneaux roses autour des clusters myosine
ax.scatter(_pcBcr[_myoC,0],_pcBcr[_myoC,1],s=_szcr[_myoC]*2.4,facecolor='none',edgecolor='#d81b60',lw=2.6,zorder=5)
_fb2=0.8*np.abs(_pcBcr).max()/np.abs(_corBcr).max()
for j,nm in enumerate(_bfcr):
    ax.arrow(0,0,_corBcr[j,0]*_fb2,_corBcr[j,1]*_fb2,color='#888',head_width=.12,lw=1.6,length_includes_head=True,alpha=.7,zorder=3)
    ax.annotate({'conservation':'sensibilité mutationnelle','enfouissement moyen (%ASA)':'engagement interaction (%ASA)'}.get(nm,nm),(_corBcr[j,0]*_fb2*1.1,_corBcr[j,1]*_fb2*1.1),fontsize=12,color='#333',zorder=6)
_txtc=[]
for k,b in enumerate(_idxBcr):
    if _myoC[k]:
        _txtc.append(ax.text(_pcBcr[k,0],_pcBcr[k,1],f"{b}\n{_myoname[b]}",fontsize=8,fontweight='bold',color='#880e4f',zorder=8))
try:
    from adjustText import adjust_text
    adjust_text(_txtc,ax=ax,arrowprops=dict(arrowstyle='-',lw=.6,color='#d81b60'),force_text=(.4,.6),lim=400)
except Exception: pass
ax.legend(handles=[Line2D([0],[0],marker='o',ls='',mfc=_cccr[g],mec='k',ms=11,label=f'{g} (n={int((_descBcr.classe==g).sum())})')
                   for g in ['hetero','homo','mixte'] if (_descBcr.classe==g).any()]
          +[Line2D([0],[0],marker='o',ls='',mfc='none',mec='#d81b60',mew=2.4,ms=15,label=f'avec myosine (n={int(_myoC.sum())})')],
          loc='upper left',fontsize=11,title='couleur = type de zone · cercle = myosine',title_fontsize=11)
ax.axhline(0,color='#eee',lw=.6); ax.axvline(0,color='#eee',lw=.6)
ax.set_xlabel(f'PC1 ({_evBcr[0]*100:.0f}%)',fontsize=13); ax.set_ylabel(f'PC2 ({_evBcr[1]*100:.0f}%)',fontsize=13)
ax.set_title(f"PCA binding-site COEUR — clusters avec MYOSINE entoures ({int(_myoC.sum())})",fontsize=14)
plt.tight_layout(); plt.show()
for b in _idxBcr:
    if b in _myoname: print(f'{b}: {_myoname[b]}')


# # C — Chimie des contacts : résidus d'actine conservés vs tolérants
# 
# Un point = **un contact** résidu-résidu. On compare la **fraction de chaque type d'appariement** (charges complémentaires, pont salin, liaison H, packing hydrophobe, stacking aromatique, paire polaire) entre les contacts portés par des résidus d'actine **sensibles** (conservés) vs **tolérants**. Test khi². Complète les PCA zones/positions sans redondance.
# 

# In[96]:


# === PCA des contacts residu-residu — couleur = classe (sensitive/tolerant) — AUTONOME ===
# === PCA des contacts residu-residu - descripteurs par cote A/B ===
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path as _P
_R=_P('..') if _P('../data').exists() else _P('.')
_c=pd.read_csv(_R/'data/filtered/details/4.inter-residue_contacts.csv',low_memory=False)
def _num(s): return pd.to_numeric(s.astype(str).str.replace('Å²','',regex=False).str.replace('%','',regex=False).str.strip(),errors='coerce')
_c['asa_A']=_num(_c['asa_pct_A']); _c['asa_B']=_num(_c['asa_pct_B'])   # %ASA (0-100)
_KD={'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
_CHG={'R':1,'K':1,'D':-1,'E':-1}; _ARO={'F','Y','W'}; _POL={'S','T','N','Q','H','C'}
_a=_c['residue_A_name'].astype(str); _b=_c['residue_B_name'].astype(str)
_c['hydrophobie_A']=_a.map(lambda r:_KD.get(r,np.nan)); _c['hydrophobie_B']=_b.map(lambda r:_KD.get(r,np.nan))
_c['polaire_A']=_a.isin(_POL).astype(int); _c['polaire_B']=_b.isin(_POL).astype(int)
_qA=_a.map(lambda r:_CHG.get(r,0)); _qB=_b.map(lambda r:_CHG.get(r,0)); _c['compl_charge']=-(_qA*_qB)
_c['aromatique_A']=_a.isin(_ARO).astype(int); _c['aromatique_B']=_b.isin(_ARO).astype(int)
_c['charge_A']=_qA.values; _c['charge_B']=_qB.values             # charges nettes separees
# conservation cote actine
_pp=pd.read_csv(_R/'data/filtered/proteins_per_pdb.csv'); _ac=set(_pp[_pp.is_actin]['chain'].str.lower())
_cdf=pd.read_csv(_R/'data/proteocast/conservation_vs_asa_per_position.csv'); _cons=dict(zip(_cdf['canon'].astype(int),_cdf['conservation']))
_cnA=pd.to_numeric(_c['residue_A_canon_mafft'],errors='coerce'); _cnB=pd.to_numeric(_c['residue_B_canon_mafft'],errors='coerce')
_isA=_c['chain_A_id'].str.lower().isin(_ac).values; _isB=_c['chain_B_id'].str.lower().isin(_ac).values
_c['iface']=np.where(_isA&_isB,'actine-actine',np.where(_isA^_isB,'actine-ABP','autre'))  # type de contact homo/hetero
def _cv(cn,ia): return _cons.get(int(cn),np.nan) if (ia and pd.notna(cn)) else np.nan
_csA=np.array([_cv(cn,ia) for cn,ia in zip(_cnA,_isA)]); _csB=np.array([_cv(cn,ia) for cn,ia in zip(_cnB,_isB)])
_c['conservation_actine']=np.nanmean(np.vstack([_csA,_csB]),axis=0)
def _typ(t):
    if pd.isna(t): return 'vdW'
    if 'Salt' in t and 'H-bond' in t: return 'H-bond+sel'
    if 'Salt' in t: return 'pont salin'
    if 'H-bond' in t: return 'H-bond'
    return 'vdW'
_c['type']=_c['contact_type'].map(_typ)
_feat=['asa_A','asa_B','hydrophobie_A','hydrophobie_B','polaire_A','polaire_B','charge_A','charge_B','compl_charge','aromatique_A','aromatique_B','conservation_actine']
_c=_c.dropna(subset=[f for f in _feat if f!='conservation_actine']).reset_index(drop=True)  # retire contacts a residu inconnu (NaN) -> PCA stable

# === Variante : PCA des contacts SANS conservation, couleur = classe ProteoCast (actine) ===
# Cellule AUTONOME (le bloc data ci-dessous reconstruit _c, _feat, _cdf, _ac)
_clsmap=dict(zip(_cdf['canon'].astype(int),_cdf['residue_class']))   # canon -> sensitive / tolerant
_cnA=pd.to_numeric(_c['residue_A_canon_mafft'],errors='coerce'); _cnB=pd.to_numeric(_c['residue_B_canon_mafft'],errors='coerce')
_isA=_c['chain_A_id'].str.lower().isin(_ac).values; _isB=_c['chain_B_id'].str.lower().isin(_ac).values
_clA=_cnA.map(lambda x:_clsmap.get(int(x)) if pd.notna(x) else None)
_clB=_cnB.map(lambda x:_clsmap.get(int(x)) if pd.notna(x) else None)
_clsv=np.array(['inconnu']*len(_c),dtype=object)
_useB=_isB & _clB.notna().values; _clsv[_useB]=_clB[_useB].astype(str).values
_useA=_isA & _clA.notna().values; _clsv[_useA]=_clA[_useA].astype(str).values   # cote actine A prioritaire
import collections as _co; print('classes des contacts :', dict(_co.Counter(_clsv)))

# PCA sans le descripteur conservation_actine
# descripteurs de PAIRING (propres au contact, NON redondants avec PCA zones/positions) :
_c['match_hydrophobe']=((_c['hydrophobie_A']>0)&(_c['hydrophobie_B']>0)).astype(int)
_c['pi_aromatique']=(_c['aromatique_A'].astype(bool)&_c['aromatique_B'].astype(bool)).astype(int)
_c['match_polaire']=(_c['polaire_A'].astype(bool)&_c['polaire_B'].astype(bool)).astype(int)
_c['is_hbond']=_c['type'].isin(['H-bond','H-bond+sel']).astype(int)
_c['is_saltbridge']=_c['type'].isin(['pont salin','H-bond+sel']).astype(int)
_c['enfoui_contact']=_c[['asa_A','asa_B']].mean(axis=1)
_c['compl']=(_c['compl_charge']>0).astype(int)
# --- n_prot : nb de familles (seq cluster 70, ACTINE comprise) en contact avec la position actine ---
_din=pd.read_csv(_R/'data/filtered/details/1.interactions.csv'); _dan=pd.read_csv(_R/'data/filtered/filtered_all_data.csv',low_memory=False)
_fam=_din.merge(_dan[['subunit_1','subunit_2','s2_sequence_cluster_70']].drop_duplicates(),left_on=['chain_A_id','chain_B_id'],right_on=['subunit_1','subunit_2'],how='left').drop_duplicates('interaction_id')
_i2f=dict(zip(_fam.interaction_id,_fam.s2_sequence_cluster_70))
_cnA2=pd.to_numeric(_c['residue_A_canon_mafft'],errors='coerce'); _cnB2=pd.to_numeric(_c['residue_B_canon_mafft'],errors='coerce')
_isBn=_c['chain_B_id'].str.lower().isin(_ac).values; _fvn=_c.interaction_id.map(_i2f)
from collections import defaultdict as _ddn
_Sn=_ddn(set)
for _cn,_ib,_cnb,_ff in zip(_cnA2,_isBn,_cnB2,_fvn):
    if pd.notna(_cn): _Sn[int(_cn)].add('actine' if _ib else _ff)
    if _ib and pd.notna(_cnb): _Sn[int(_cnb)].add('actine')
_npr={k:len([x for x in v if pd.notna(x)]) for k,v in _Sn.items()}
_c['n_prot']=_cnA2.map(lambda x:_npr.get(int(x),np.nan) if pd.notna(x) else np.nan)
# === QUANTITATIF : type d'appariement (%) + nb de familles en contact (n_prot) ===
from scipy.stats import chi2_contingency, mannwhitneyu
_msS=_clsv=='sensitive'; _msT=_clsv=='tolerant'; _ns=int(_msS.sum()); _nt=int(_msT.sum())
_c['compl']=(_c['compl_charge']>0).astype(int)
_types=[('compl','charges complémentaires'),('is_saltbridge','pont salin'),('is_hbond','liaison H'),
        ('match_hydrophobe','contact hydrophobe'),('pi_aromatique','paire aromatique'),('match_polaire','paire polaire')]
_lab=[l for _,l in _types]
_fs=[100*_c.loc[_msS,k].mean() for k,_ in _types]; _ft=[100*_c.loc[_msT,k].mean() for k,_ in _types]
_pv=[chi2_contingency([[int(_c.loc[_msS,k].sum()),_ns-int(_c.loc[_msS,k].sum())],[int(_c.loc[_msT,k].sum()),_nt-int(_c.loc[_msT,k].sum())]])[1] for k,_ in _types]
_x=np.arange(len(_types)); _w=0.38
fig,(ax,ax2)=plt.subplots(1,2,figsize=(13,6),gridspec_kw={'width_ratios':[4.2,1]})
_b1=ax.bar(_x-_w/2,_fs,_w,label=f'sensitive (n={_ns})',color='#D55E00'); _b2=ax.bar(_x+_w/2,_ft,_w,label=f'tolerant (n={_nt})',color='#0072B2')
for _bars in (_b1,_b2):
    for _r in _bars: ax.annotate(f'{_r.get_height():.1f}',(_r.get_x()+_r.get_width()/2,_r.get_height()),ha='center',va='bottom',fontsize=8)
for _i,_pp2 in enumerate(_pv):
    _s='***' if _pp2<1e-3 else '**' if _pp2<1e-2 else '*' if _pp2<.05 else 'ns'
    ax.annotate(_s,(_x[_i],max(_fs[_i],_ft[_i])+0.9),ha='center',fontsize=10,color='#333',fontweight='bold')
ax.set_xticks(_x); ax.set_xticklabels(_lab,rotation=18,ha='right',fontsize=9); ax.set_ylabel('% des contacts (cote actine)',fontsize=11)
ax.set_ylim(0,max(max(_fs),max(_ft))*1.18); ax.legend(fontsize=10); ax.grid(axis='y',alpha=.3)
_posdf=pd.DataFrame({'canon':_cnA2.values,'n_prot':_c['n_prot'].values,'cls':_clsv}).dropna(subset=['canon','n_prot']).drop_duplicates('canon')   # 1 valeur PAR POSITION (evite la pseudo-replication)
_es=_posdf.loc[_posdf.cls=='sensitive','n_prot'].mean(); _et=_posdf.loc[_posdf.cls=='tolerant','n_prot'].mean()
_pn=mannwhitneyu(_posdf.loc[_posdf.cls=='sensitive','n_prot'],_posdf.loc[_posdf.cls=='tolerant','n_prot']).pvalue
_bn=ax2.bar([0,1],[_es,_et],color=['#D55E00','#0072B2'])
for _r in _bn: ax2.annotate(f'{_r.get_height():.1f}',(_r.get_x()+_r.get_width()/2,_r.get_height()),ha='center',va='bottom',fontsize=9)
ax2.set_xticks([0,1]); ax2.set_xticklabels(['sens.','tol.'],fontsize=9); ax2.set_ylabel('nb familles en contact (moy./position)',fontsize=10)
ax2.set_ylim(0,max(_es,_et)*1.28); ax2.set_title('n_prot',fontsize=10); ax2.grid(axis='y',alpha=.3)
ax2.annotate('ns' if _pn>=.05 else f'p={_pn:.0e}',(0.5,max(_es,_et)*1.13),ha='center',fontsize=9,fontweight='bold',color='#333')
fig.suptitle("Type d'appariement + nb familles en contact (n_prot) — SENSIBLES vs TOLERANTS — TOUS contacts (brut)",fontsize=11); plt.tight_layout(); plt.show()
print(f"enfoui_contact moyen : sensitive {_c.loc[_msS,'enfoui_contact'].mean():.1f}%  vs  tolerant {_c.loc[_msT,'enfoui_contact'].mean():.1f}%")
print(f"n_prot moyen : sensitive {_es:.1f}  vs  tolerant {_et:.1f}  (Mann-Whitney p={_pn:.1e}) -> {'ns (la promiscuite ne separe pas)' if _pn>=.05 else 'significatif'}")


# ### Variante — COEUR (les 2 résidus > 25 % ASA)
# 
# Mêmes barres, mais en ne gardant que les contacts où **les deux** résidus sont enfouis > 25 % ASA (vrais contacts de cœur). Version COEUR → **résultats** ; version brute ci-dessus → annexe/robustesse.
# 

# In[97]:


# === PCA des contacts residu-residu — couleur = classe (sensitive/tolerant) — AUTONOME ===
# === PCA des contacts residu-residu - descripteurs par cote A/B ===
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path as _P
_R=_P('..') if _P('../data').exists() else _P('.')
_c=pd.read_csv(_R/'data/filtered/details/4.inter-residue_contacts.csv',low_memory=False)
def _num(s): return pd.to_numeric(s.astype(str).str.replace('Å²','',regex=False).str.replace('%','',regex=False).str.strip(),errors='coerce')
_c['asa_A']=_num(_c['asa_pct_A']); _c['asa_B']=_num(_c['asa_pct_B'])   # %ASA (0-100)
_KD={'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
_CHG={'R':1,'K':1,'D':-1,'E':-1}; _ARO={'F','Y','W'}; _POL={'S','T','N','Q','H','C'}
_a=_c['residue_A_name'].astype(str); _b=_c['residue_B_name'].astype(str)
_c['hydrophobie_A']=_a.map(lambda r:_KD.get(r,np.nan)); _c['hydrophobie_B']=_b.map(lambda r:_KD.get(r,np.nan))
_c['polaire_A']=_a.isin(_POL).astype(int); _c['polaire_B']=_b.isin(_POL).astype(int)
_qA=_a.map(lambda r:_CHG.get(r,0)); _qB=_b.map(lambda r:_CHG.get(r,0)); _c['compl_charge']=-(_qA*_qB)
_c['aromatique_A']=_a.isin(_ARO).astype(int); _c['aromatique_B']=_b.isin(_ARO).astype(int)
_c['charge_A']=_qA.values; _c['charge_B']=_qB.values             # charges nettes separees
# conservation cote actine
_pp=pd.read_csv(_R/'data/filtered/proteins_per_pdb.csv'); _ac=set(_pp[_pp.is_actin]['chain'].str.lower())
_cdf=pd.read_csv(_R/'data/proteocast/conservation_vs_asa_per_position.csv'); _cons=dict(zip(_cdf['canon'].astype(int),_cdf['conservation']))
_cnA=pd.to_numeric(_c['residue_A_canon_mafft'],errors='coerce'); _cnB=pd.to_numeric(_c['residue_B_canon_mafft'],errors='coerce')
_isA=_c['chain_A_id'].str.lower().isin(_ac).values; _isB=_c['chain_B_id'].str.lower().isin(_ac).values
_c['iface']=np.where(_isA&_isB,'actine-actine',np.where(_isA^_isB,'actine-ABP','autre'))  # type de contact homo/hetero
def _cv(cn,ia): return _cons.get(int(cn),np.nan) if (ia and pd.notna(cn)) else np.nan
_csA=np.array([_cv(cn,ia) for cn,ia in zip(_cnA,_isA)]); _csB=np.array([_cv(cn,ia) for cn,ia in zip(_cnB,_isB)])
_c['conservation_actine']=np.nanmean(np.vstack([_csA,_csB]),axis=0)
def _typ(t):
    if pd.isna(t): return 'vdW'
    if 'Salt' in t and 'H-bond' in t: return 'H-bond+sel'
    if 'Salt' in t: return 'pont salin'
    if 'H-bond' in t: return 'H-bond'
    return 'vdW'
_c['type']=_c['contact_type'].map(_typ)
_feat=['asa_A','asa_B','hydrophobie_A','hydrophobie_B','polaire_A','polaire_B','charge_A','charge_B','compl_charge','aromatique_A','aromatique_B','conservation_actine']
_c=_c.dropna(subset=[f for f in _feat if f!='conservation_actine']).reset_index(drop=True)
_c=_c[(_c['asa_A']>25)&(_c['asa_B']>25)].reset_index(drop=True)  # COEUR : les 2 residus enfouis >25% ASA  # retire contacts a residu inconnu (NaN) -> PCA stable

# === Variante : PCA des contacts SANS conservation, couleur = classe ProteoCast (actine) ===
# Cellule AUTONOME (le bloc data ci-dessous reconstruit _c, _feat, _cdf, _ac)
_clsmap=dict(zip(_cdf['canon'].astype(int),_cdf['residue_class']))   # canon -> sensitive / tolerant
_cnA=pd.to_numeric(_c['residue_A_canon_mafft'],errors='coerce'); _cnB=pd.to_numeric(_c['residue_B_canon_mafft'],errors='coerce')
_isA=_c['chain_A_id'].str.lower().isin(_ac).values; _isB=_c['chain_B_id'].str.lower().isin(_ac).values
_clA=_cnA.map(lambda x:_clsmap.get(int(x)) if pd.notna(x) else None)
_clB=_cnB.map(lambda x:_clsmap.get(int(x)) if pd.notna(x) else None)
_clsv=np.array(['inconnu']*len(_c),dtype=object)
_useB=_isB & _clB.notna().values; _clsv[_useB]=_clB[_useB].astype(str).values
_useA=_isA & _clA.notna().values; _clsv[_useA]=_clA[_useA].astype(str).values   # cote actine A prioritaire
import collections as _co; print('classes des contacts :', dict(_co.Counter(_clsv)))

# PCA sans le descripteur conservation_actine
# descripteurs de PAIRING (propres au contact, NON redondants avec PCA zones/positions) :
_c['match_hydrophobe']=((_c['hydrophobie_A']>0)&(_c['hydrophobie_B']>0)).astype(int)             # packing hydrophobe
_c['pi_aromatique']=(_c['aromatique_A'].astype(bool)&_c['aromatique_B'].astype(bool)).astype(int) # stacking aromatique
_c['match_polaire']=(_c['polaire_A'].astype(bool)&_c['polaire_B'].astype(bool)).astype(int)       # paire polaire-polaire
_c['is_hbond']=_c['type'].isin(['H-bond','H-bond+sel']).astype(int)
_c['is_saltbridge']=_c['type'].isin(['pont salin','H-bond+sel']).astype(int)
_c['enfoui_contact']=_c[['asa_A','asa_B']].mean(axis=1)                                           # profondeur du contact
_featc=['compl_charge','match_hydrophobe','pi_aromatique','match_polaire','is_hbond','is_saltbridge','enfoui_contact']  # PAIRING : nature de l'appariement A-B (pas la physico-chimie d'un seul residu)
# --- n_prot : nb de familles (seq cluster 70, ACTINE comprise) en contact avec la position actine ---
_din=pd.read_csv(_R/'data/filtered/details/1.interactions.csv'); _dan=pd.read_csv(_R/'data/filtered/filtered_all_data.csv',low_memory=False)
_fam=_din.merge(_dan[['subunit_1','subunit_2','s2_sequence_cluster_70']].drop_duplicates(),left_on=['chain_A_id','chain_B_id'],right_on=['subunit_1','subunit_2'],how='left').drop_duplicates('interaction_id')
_i2f=dict(zip(_fam.interaction_id,_fam.s2_sequence_cluster_70))
_cnA2=pd.to_numeric(_c['residue_A_canon_mafft'],errors='coerce'); _cnB2=pd.to_numeric(_c['residue_B_canon_mafft'],errors='coerce')
_isBn=_c['chain_B_id'].str.lower().isin(_ac).values; _fvn=_c.interaction_id.map(_i2f)
from collections import defaultdict as _ddn
_Sn=_ddn(set)
for _cn,_ib,_cnb,_ff in zip(_cnA2,_isBn,_cnB2,_fvn):
    if pd.notna(_cn): _Sn[int(_cn)].add('actine' if _ib else _ff)
    if _ib and pd.notna(_cnb): _Sn[int(_cnb)].add('actine')
_npr={k:len([x for x in v if pd.notna(x)]) for k,v in _Sn.items()}
_c['n_prot']=_cnA2.map(lambda x:_npr.get(int(x),np.nan) if pd.notna(x) else np.nan)
# === QUANTITATIF : type d'appariement (%) + nb de familles en contact (n_prot) ===
from scipy.stats import chi2_contingency, mannwhitneyu
_msS=_clsv=='sensitive'; _msT=_clsv=='tolerant'; _ns=int(_msS.sum()); _nt=int(_msT.sum())
_c['compl']=(_c['compl_charge']>0).astype(int)
_types=[('compl','charges complémentaires'),('is_saltbridge','pont salin'),('is_hbond','liaison H'),
        ('match_hydrophobe','contact hydrophobe'),('pi_aromatique','paire aromatique'),('match_polaire','paire polaire')]
_lab=[l for _,l in _types]
_fs=[100*_c.loc[_msS,k].mean() for k,_ in _types]; _ft=[100*_c.loc[_msT,k].mean() for k,_ in _types]
_pv=[chi2_contingency([[int(_c.loc[_msS,k].sum()),_ns-int(_c.loc[_msS,k].sum())],[int(_c.loc[_msT,k].sum()),_nt-int(_c.loc[_msT,k].sum())]])[1] for k,_ in _types]
_x=np.arange(len(_types)); _w=0.38
fig,(ax,ax2)=plt.subplots(1,2,figsize=(13,6),gridspec_kw={'width_ratios':[4.2,1]})
_b1=ax.bar(_x-_w/2,_fs,_w,label=f'sensitive (n={_ns})',color='#D55E00'); _b2=ax.bar(_x+_w/2,_ft,_w,label=f'tolerant (n={_nt})',color='#0072B2')
for _bars in (_b1,_b2):
    for _r in _bars: ax.annotate(f'{_r.get_height():.1f}',(_r.get_x()+_r.get_width()/2,_r.get_height()),ha='center',va='bottom',fontsize=8)
for _i,_pp2 in enumerate(_pv):
    _s='***' if _pp2<1e-3 else '**' if _pp2<1e-2 else '*' if _pp2<.05 else 'ns'
    ax.annotate(_s,(_x[_i],max(_fs[_i],_ft[_i])+0.9),ha='center',fontsize=10,color='#333',fontweight='bold')
ax.set_xticks(_x); ax.set_xticklabels(_lab,rotation=18,ha='right',fontsize=9); ax.set_ylabel('% des contacts (cote actine)',fontsize=11)
ax.set_ylim(0,max(max(_fs),max(_ft))*1.18); ax.legend(fontsize=10); ax.grid(axis='y',alpha=.3)
_posdf=pd.DataFrame({'canon':_cnA2.values,'n_prot':_c['n_prot'].values,'cls':_clsv}).dropna(subset=['canon','n_prot']).drop_duplicates('canon')   # 1 valeur PAR POSITION (evite la pseudo-replication)
_es=_posdf.loc[_posdf.cls=='sensitive','n_prot'].mean(); _et=_posdf.loc[_posdf.cls=='tolerant','n_prot'].mean()
_pn=mannwhitneyu(_posdf.loc[_posdf.cls=='sensitive','n_prot'],_posdf.loc[_posdf.cls=='tolerant','n_prot']).pvalue
_bn=ax2.bar([0,1],[_es,_et],color=['#D55E00','#0072B2'])
for _r in _bn: ax2.annotate(f'{_r.get_height():.1f}',(_r.get_x()+_r.get_width()/2,_r.get_height()),ha='center',va='bottom',fontsize=9)
ax2.set_xticks([0,1]); ax2.set_xticklabels(['sens.','tol.'],fontsize=9); ax2.set_ylabel('nb familles en contact (moy./position)',fontsize=10)
ax2.set_ylim(0,max(_es,_et)*1.28); ax2.set_title('n_prot',fontsize=10); ax2.grid(axis='y',alpha=.3)
ax2.annotate('ns' if _pn>=.05 else f'p={_pn:.0e}',(0.5,max(_es,_et)*1.13),ha='center',fontsize=9,fontweight='bold',color='#333')
fig.suptitle("Type d'appariement + nb familles en contact (n_prot) — SENSIBLES vs TOLERANTS — COEUR >25% ASA",fontsize=11); plt.tight_layout(); plt.show()
print(f"enfoui_contact moyen : sensitive {_c.loc[_msS,'enfoui_contact'].mean():.1f}%  vs  tolerant {_c.loc[_msT,'enfoui_contact'].mean():.1f}%")
print(f"n_prot moyen : sensitive {_es:.1f}  vs  tolerant {_et:.1f}  (Mann-Whitney p={_pn:.1e}) -> {'ns (la promiscuite ne separe pas)' if _pn>=.05 else 'significatif'}")


# In[98]:


# === PCA des contacts COEUR >25% ASA + descripteur n_prot (familles en contact, actine comprise) ===
import warnings; warnings.filterwarnings('ignore')
from pathlib import Path as _P
_R=_P('..') if _P('../data').exists() else _P('.')
_c=pd.read_csv(_R/'data/filtered/details/4.inter-residue_contacts.csv',low_memory=False)
_da=pd.read_csv(_R/'data/filtered/filtered_all_data.csv',low_memory=False)
_di=pd.read_csv(_R/'data/filtered/details/1.interactions.csv')
_pp=pd.read_csv(_R/'data/filtered/proteins_per_pdb.csv'); _ac=set(_pp[_pp.is_actin]['chain'].str.lower())
_cdf=pd.read_csv(_R/'data/proteocast/conservation_vs_asa_per_position.csv'); _clsmap=dict(zip(_cdf['canon'].astype(int),_cdf['residue_class']))
def _num(s): return pd.to_numeric(s.astype(str).str.replace('%','',regex=False).str.replace('Å²','',regex=False).str.strip(),errors='coerce')
_c['asa_A']=_num(_c['asa_pct_A']); _c['asa_B']=_num(_c['asa_pct_B'])
_KD={'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
_CHG={'R':1,'K':1,'D':-1,'E':-1}; _ARO={'F','Y','W'}; _POL={'S','T','N','Q','H','C'}
_a=_c['residue_A_name'].astype(str); _b=_c['residue_B_name'].astype(str)
_c['hydrophobie_A']=_a.map(lambda r:_KD.get(r,np.nan)); _c['hydrophobie_B']=_b.map(lambda r:_KD.get(r,np.nan))
_c['polaire_A']=_a.isin(_POL).astype(int); _c['polaire_B']=_b.isin(_POL).astype(int)
_qA=_a.map(lambda r:_CHG.get(r,0)); _qB=_b.map(lambda r:_CHG.get(r,0)); _c['compl_charge']=(-(_qA*_qB)).values
_c['aromatique_A']=_a.isin(_ARO).astype(int); _c['aromatique_B']=_b.isin(_ARO).astype(int)
_c['charge_A']=_qA.values; _c['charge_B']=_qB.values
_c['cnA']=pd.to_numeric(_c['residue_A_canon_mafft'],errors='coerce'); _c['cnB']=pd.to_numeric(_c['residue_B_canon_mafft'],errors='coerce')
_c=_c.dropna(subset=['hydrophobie_A','hydrophobie_B']).reset_index(drop=True)
_c=_c[(_c.asa_A>25)&(_c.asa_B>25)].reset_index(drop=True)   # COEUR : les 2 residus enfouis >25%
# --- n_prot : nb de familles (sequence cluster 70, ACTINE comprise) en contact avec la position actine ---
_fam=_di.merge(_da[['subunit_1','subunit_2','s2_sequence_cluster_70']].drop_duplicates(),
               left_on=['chain_A_id','chain_B_id'],right_on=['subunit_1','subunit_2'],how='left').drop_duplicates('interaction_id')
_i2f=dict(zip(_fam.interaction_id,_fam.s2_sequence_cluster_70)); _fv=_c.interaction_id.map(_i2f)
_isB=_c['chain_B_id'].str.lower().isin(_ac).values    # actine = chain_A ; isB=True -> contact homo (actine-actine)
from collections import defaultdict as _dd
_S=_dd(set)
for _cn,_ib,_cnb,_f in zip(_c.cnA,_isB,_c.cnB,_fv):
    if pd.notna(_cn): _S[int(_cn)].add('actine' if _ib else _f)
    if _ib and pd.notna(_cnb): _S[int(_cnb)].add('actine')
_nprot={k:len([x for x in v if pd.notna(x)]) for k,v in _S.items()}
_c['n_prot']=_c.cnA.map(lambda x:_nprot.get(int(x),np.nan) if pd.notna(x) else np.nan)
# --- PCA ---
_feat=['asa_A','asa_B','hydrophobie_A','hydrophobie_B','polaire_A','polaire_B','charge_A','charge_B','compl_charge','aromatique_A','aromatique_B','n_prot']
_c=_c.dropna(subset=_feat).reset_index(drop=True)
_clsv=_c.cnA.map(lambda x:_clsmap.get(int(x)) if pd.notna(x) else None).fillna('inconnu').values
_X=_c[_feat].values.astype(float); _Xz=(_X-_X.mean(0))/(_X.std(0)+1e-9); _Xc=_Xz-_Xz.mean(0)
_cov=np.cov(_Xc.T); _val,_vec=np.linalg.eigh(_cov); _o=_val.argsort()[::-1]; _val,_vec=_val[_o],_vec[:,_o]
_pcs=_Xc@_vec[:,:2]; _ev=_val/_val.sum()*100
_cor=np.array([[np.corrcoef(_Xz[:,j],_pcs[:,k])[0,1] for k in range(2)] for j in range(len(_feat))])
_ccol={'sensitive':'#D55E00','tolerant':'#0072B2','inconnu':'#dddddd'}
np.random.seed(0); fig,ax=plt.subplots(figsize=(13,10))
for g in ['inconnu','tolerant','sensitive']:
    sub=np.where(_clsv==g)[0]
    if len(sub)==0: continue
    if len(sub)>18000: sub=np.random.choice(sub,18000,replace=False)
    ax.scatter(_pcs[sub,0],_pcs[sub,1],s=7,alpha=.25 if g=='inconnu' else .4,color=_ccol[g],label=f'{g} (n={int((_clsv==g).sum())})',edgecolor='none')
_fk=0.6*np.abs(_pcs).max()/np.abs(_cor).max()
for j,nm in enumerate(_feat):
    ax.annotate('',xy=(_cor[j,0]*_fk,_cor[j,1]*_fk),xytext=(0,0),arrowprops=dict(arrowstyle='-|>',color='k',lw=1.3,mutation_scale=12),zorder=5)
    ax.annotate({'conservation':'sensibilité mutationnelle','enfouissement moyen (%ASA)':'engagement interaction (%ASA)'}.get(nm,nm),(_cor[j,0]*_fk*1.08,_cor[j,1]*_fk*1.08),fontsize=10,fontweight='bold',zorder=6)
ax.set_xlabel(f'PC1 ({_ev[0]:.0f}%)',fontsize=12); ax.set_ylabel(f'PC2 ({_ev[1]:.0f}%)',fontsize=12)
ax.legend(markerscale=3,fontsize=10,loc='best',title='classe ProteoCast (actine)')
ax.set_title('PCA des contacts COEUR >25% ASA + n_prot (nb familles seq70 en contact, actine comprise)\ncouleur = classe (sensitive/tolerant)',fontsize=11)
plt.tight_layout(); plt.show()
print(f"Variance PC1+PC2 = {_ev[0]:.0f}% + {_ev[1]:.0f}% = {_ev[0]+_ev[1]:.0f}% ({len(_feat)} descripteurs dont n_prot)")
print(f"n_prot moyen (par contact) : sensible {_c['n_prot'][_clsv=='sensitive'].mean():.1f} vs tolerant {_c['n_prot'][_clsv=='tolerant'].mean():.1f}")


# In[99]:


# === Version BARRES de la PCA des contacts COEUR : effet de CHAQUE descripteur (sensible vs tolerant) ===
# A executer juste apres la cellule PCA des contacts COEUR (reutilise _c, _feat, _clsv).
from scipy.stats import mannwhitneyu
_mS=_clsv=='sensitive'; _mT=_clsv=='tolerant'
_rows=[]
for _ff in _feat:
    _av=_c.loc[_mS,_ff].astype(float); _bv=_c.loc[_mT,_ff].astype(float)
    _d=(_av.mean()-_bv.mean())/(np.sqrt((_av.var()+_bv.var())/2)+1e-9)   # Cohen's d (taille d'effet)
    _pp=mannwhitneyu(_av,_bv).pvalue
    _rows.append((_ff,_av.mean(),_bv.mean(),_d,_pp))
_bd=pd.DataFrame(_rows,columns=['descripteur','moy_sensible','moy_tolerant','effet_d','p']).sort_values('effet_d')
fig,ax=plt.subplots(figsize=(9,7))
_col=['#D55E00' if x>0 else '#0072B2' for x in _bd['effet_d']]
ax.barh(_bd['descripteur'],_bd['effet_d'],color=_col,edgecolor='k',lw=.3)
_mx=np.abs(_bd['effet_d']).max()
for _y,(_dd2,_pp2) in enumerate(zip(_bd['effet_d'],_bd['p'])):
    _s='***' if _pp2<1e-3 else '**' if _pp2<1e-2 else '*' if _pp2<.05 else 'ns'
    ax.annotate(_s,(_dd2+(0.012*_mx*8 if _dd2>=0 else -0.012*_mx*8),_y),va='center',ha='left' if _dd2>=0 else 'right',fontsize=9,fontweight='bold')
ax.axvline(0,color='k',lw=.8); ax.set_xlim(-_mx*1.3,_mx*1.3)
ax.set_xlabel("taille d'effet (Cohen's d)  —  >0 : plus chez SENSIBLE (orange)   |   <0 : plus chez TOLERANT (bleu)",fontsize=9)
ax.set_title("Contacts COEUR >25% ASA : pouvoir discriminant de chaque descripteur (SENSIBLE vs TOLERANT)\n"
             "barres = ecart standardise ; * p<.05  ** p<.01  *** p<.001 (Mann-Whitney)",fontsize=10)
ax.grid(axis='x',alpha=.3); plt.tight_layout(); plt.show()
print("(n enorme -> p toujours petit ; lire la TAILLE D'EFFET. A-side & n_prot = proprietes de position -> pseudo-replication)")
print(_bd.round(3).to_string(index=False))


# # D — PCA des ABP
# 
# Un point = **un ABP**. Deux vues : (1) **empreinte** — ABP regroupés par leur *site* de liaison sur l'actine ; (2) **biophysique** — ABP décrits par la physico-chimie de leur interface.

# In[100]:


# Heatmap ABP x positions actine (= streamlit) + dendrogramme (pour Word)
import seaborn as sns, matplotlib as mpl
from scipy.cluster.hierarchy import leaves_list
from pathlib import Path as _P
_R = _P('..') if _P('../data').exists() else _P('.')
df3 = pd.read_csv(_R/'data/filtered/details/3.interface_residues.csv')
df3['buried_ASA_percent'] = pd.to_numeric(df3['buried_ASA_percent'].astype(str).str.replace('%','',regex=False), errors='coerce')
df3['canon'] = pd.to_numeric(df3['residue_number_canon_mafft'], errors='coerce')
df3 = df3[df3.canon.notna() & df3.buried_ASA_percent.notna()].copy()
di = pd.read_csv(_R/'data/filtered/details/1.interactions.csv')
da = pd.read_csv(_R/'data/filtered/filtered_all_data.csv', low_memory=False)
pp = pd.read_csv(_R/'data/filtered/proteins_per_pdb.csv')
actin_ch = set(pp[pp.is_actin]['chain'].str.lower())
homo = set(di[di['chain_B_id'].str.lower().isin(actin_ch)]['interaction_id'])
het = di[~di.interaction_id.isin(homo)].merge(
    da[['subunit_1','subunit_2','subunit_2_title','s2_actine','cluster_data_70']],
    left_on=['chain_A_id','chain_B_id'], right_on=['subunit_1','subunit_2'], how='left').drop_duplicates('interaction_id')
het = het[het['s2_actine'].fillna(False)==False].copy()
import re as _re
def _arp_name(n):
    if 'complex subunit 5-like' in n: return 'ARPC5L'
    mc = _re.search(r'complex subunit (\d+[A-Za-z]?)', n)
    if mc: return 'ARPC' + mc.group(1)
    if _re.match(r'Actin-related protein 3(\b|,|$)', n): return 'Arp3'
    if _re.match(r'Actin-related protein 2(\b|,|$)', n): return 'Arp2'
    return n
het['abp'] = het['subunit_2_title'].fillna('Unknown').str.replace(r'\s*\(.*?\)','',regex=True).str.strip().str[:50].map(_arp_name)
s1ch=het.set_index('interaction_id')['chain_A_id'].str.lower(); abpn=het.set_index('interaction_id')['abp']; c70n=het.set_index('interaction_id')['cluster_data_70']
d3 = df3[df3.interaction_id.isin(set(het.interaction_id))].copy()
d3['_s1c']=d3.interaction_id.map(s1ch); d3=d3[d3.chain.str.lower()==d3._s1c].copy()
d3['abp']=d3.interaction_id.map(abpn); d3['canon']=d3.canon.astype(int); d3['c70']=d3.interaction_id.map(c70n)
abp_freq=d3.groupby('abp')['interaction_id'].nunique().sort_values(ascending=False)
ac70n=d3.groupby(['abp','c70'])['interaction_id'].nunique()
ag=d3.groupby(['abp','c70','canon'])['buried_ASA_percent'].sum().reset_index(name='s')
ag['a70']=ag.apply(lambda r:r['s']/max(ac70n.get((r.abp,r.c70),1),1),axis=1)
nc70=d3.groupby('abp')['c70'].nunique()
eq=ag.groupby(['abp','canon'])['a70'].sum().reset_index(); eq['v']=eq['a70']/eq['abp'].map(nc70)
pivot=eq.pivot(index='abp',columns='canon',values='v').loc[abp_freq.index].fillna(0)   # noms seuls (pas de parentheses)

from scipy.cluster.hierarchy import linkage as _lkh
_Zbin=_lkh((pivot.values>0).astype(float), method='average', metric='jaccard')  # regroupement BINAIRE (presence/absence)
g=sns.clustermap(pivot, cmap='YlOrRd', figsize=(24,10.5), row_cluster=True, col_cluster=False, row_linkage=_Zbin,
                 dendrogram_ratio=(0.12,0.0), xticklabels=False, cbar_pos=None)  # couleurs=%ASA, lignes groupees en binaire
order=g.dendrogram_row.reordered_ind
g.ax_heatmap.yaxis.set_ticks_position('right'); g.ax_heatmap.yaxis.set_label_position('right')
g.ax_heatmap.set_yticks(np.arange(len(pivot))+0.5)
g.ax_heatmap.set_yticklabels([pivot.index[i] for i in order], fontsize=9, rotation=0)
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_position([0.04, 0.06, 0.55, 0.85])   # heatmap a gauche
g.ax_row_dendrogram.set_position([0.74, 0.06, 0.12, 0.85]); g.ax_row_dendrogram.invert_xaxis()  # X=0.74 -> dendro colle aux noms (baisser X = plus proche)
g.ax_heatmap.set_xlabel('Positions actine (canoniques)')
# legende gradient : barre horizontale propre en haut a gauche
cax=g.fig.add_axes([0.05,0.945,0.16,0.012])
sm=mpl.cm.ScalarMappable(cmap='YlOrRd', norm=mpl.colors.Normalize(pivot.values.min(), pivot.values.max()))
g.fig.colorbar(sm, cax=cax, orientation='horizontal', label='% ASA moyen equitable')
cax.tick_params(labelsize=7); cax.xaxis.set_label_position('top'); cax.xaxis.set_ticks_position('top')
plt.show()


# ### Choix du nombre de clusters K — methode du coude\n\nAvant de fixer K=3 pour lempreinte ABP, on le justifie : inertie intra-cluster (coude) et silhouette en fonction de K, sur la meme donnee (profil %ASA, distance correlation).

# In[101]:


# === Methode du coude (+ silhouette) : choix du nombre de clusters K pour l'empreinte ABP ===
# Meme donnee que la heatmap/PCoA : presence/absence (ASA>0), distance de Jaccard, linkage average.
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage as _lk, fcluster as _fc
_Bp=(pivot.values>0).astype(float)
_De=np.nan_to_num(squareform(pdist(_Bp,metric='jaccard')))
_ne=_De.shape[0]; _Je=np.eye(_ne)-np.ones((_ne,_ne))/_ne
_Ge=-0.5*_Je@(_De**2)@_Je; _ve,_ce=np.linalg.eigh(_Ge); _oe=_ve.argsort()[::-1]; _ve,_ce=_ve[_oe],_ce[:,_oe]
_pv=np.maximum(_ve,0); _embe=_ce[:,:min(8,_ne)]*np.sqrt(_pv[:min(8,_ne)])
_Ze=_lk(_Bp,method='average',metric='jaccard')
def _sile(D,lab):
    s=[]
    for i in range(len(lab)):
        sm=[D[i,j] for j in range(len(lab)) if lab[j]==lab[i] and j!=i]; a=np.mean(sm) if sm else 0
        ot=[np.mean([D[i,j] for j in range(len(lab)) if lab[j]==c]) for c in set(lab) if c!=lab[i]]; b=min(ot) if ot else 0
        s.append(0 if max(a,b)==0 else (b-a)/max(a,b))
    return np.mean(s)
from scipy.cluster.vq import kmeans2 as _kmE   # coude par K-MEANS (et non CAH)
_Ks=list(range(1,11)); _wcss=[]; _silh=[]
for K in _Ks:
    _bk=None
    for _sd in range(20):
        _cc,_ll=_kmE(_embe,K,minit='++',seed=_sd)
        if len(set(_ll))<K and K>1: continue
        _ss=((_embe-_cc[_ll])**2).sum()
        if _bk is None or _ss<_bk[0]: _bk=(_ss,_ll)
    _wcss.append(_bk[0]); _silh.append(_sile(_De,_bk[1]) if K>1 else np.nan)
_drops=[_wcss[i]-_wcss[i+1] for i in range(len(_wcss)-1)]
_Kopt=int(np.argmax(_drops))+2   # coude = K juste apres la plus grosse chute de WCSS
fig,ax1=plt.subplots(figsize=(7,4.5))
ax1.plot(_Ks,_wcss,'o-',color='#1f77b4'); _=0; ax1.set_xlabel('nombre de clusters K'); ax1.set_ylabel('inertie intra-cluster (WCSS)')
ax1.set_title(f"Methode du coude\nplus grosse chute de WCSS -> coude a K={_Kopt}"); plt.tight_layout(); plt.show()
print('WCSS:',[round(w,2) for w in _wcss]); print('silhouette:',[round(s,3) if not np.isnan(s) else None for s in _silh])
print(f'-> K={_Kopt} (coude). NB binaire : silhouette faible, clusters desequilibres -> regroupement binaire peu net.')


# In[102]:


try:
    from adjustText import adjust_text
except Exception:
    def adjust_text(*a, **k): pass
# helper importe depuis abp_analysis (rend la cellule autonome)
import re as _re
def _arp_name(n):
    if 'complex subunit 5-like' in n: return 'ARPC5L'
    mc = _re.search(r'complex subunit (\d+[A-Za-z]?)', n)
    if mc: return 'ARPC' + mc.group(1)
    if _re.match(r'Actin-related protein 3(\b|,|$)', n): return 'Arp3'
    if _re.match(r'Actin-related protein 2(\b|,|$)', n): return 'Arp2'
    return n

# k-means (k=2) + PCA ; couleur=cluster, forme=cote filament (+/-/partout)
from scipy.cluster.vq import kmeans2
from matplotlib.lines import Line2D
from pathlib import Path as _P
_R = _P('..') if _P('../data').exists() else _P('.')

# --- tendance cote filament par ABP (via actin_filament_positions.csv) ---
fil = pd.read_csv(_R/'data/filtered/actin_filament_positions.csv'); c2l = dict(zip(fil['chain'].str.lower(), fil['label']))
_da = pd.read_csv(_R/'data/filtered/filtered_all_data.csv', low_memory=False); _da['s2_actine']=_da['s2_actine'].fillna(False).astype(bool)
_di = pd.read_csv(_R/'data/filtered/details/1.interactions.csv').merge(
    _da[['subunit_1','subunit_2','subunit_2_title','s2_actine']].drop_duplicates(), left_on=['chain_A_id','chain_B_id'], right_on=['subunit_1','subunit_2'], how='left')
_di['s2_actine']=_di['s2_actine'].fillna(False).astype(bool); _di=_di[~_di['s2_actine']].copy()
_di['abp']=_di['subunit_2_title'].fillna('Unknown').str.replace(r'\s*\(.*?\)','',regex=True).str.strip().str[:50].map(_arp_name)
_di['fil']=_di['chain_A_id'].str.lower().map(c2l)
def _simp(l):
    if pd.isna(l): return None
    if l=='+': return '+'          # bout barbe strict
    if l=='-': return '-'          # bout pointe strict
    return 'side'                  # side, +2, +3, -2, -3 -> lateral
_di['grp']=_di['fil'].map(_simp)
# tendance par presence (ensemble des positions), identique a streamlit _pos_border_color :
# + (avec ou sans side) -> barbe ; - (avec ou sans side) -> pointe ; side seul ou les 2 bouts -> lateral
posset=_di.dropna(subset=['grp']).groupby('abp')['grp'].apply(set)
def _tend(s):
    hp='+' in s; hm='-' in s
    if hp and not hm: return '+'
    if hm and not hp: return '-'
    return 'partout'
tendance={a:_tend(posset.get(a,set())) for a in posset.index}
tend_arr=np.array([tendance.get(a,'partout') for a in pivot.index])

# --- BINAIRE : presence/absence (ASA>0) ; distance de Jaccard (= meme regroupement que la heatmap) ---
from scipy.spatial.distance import pdist, squareform
_Bp=(pivot.values>0).astype(float)
_Dj=np.nan_to_num(squareform(pdist(_Bp,metric='jaccard')))
_n=_Dj.shape[0]; _Jc=np.eye(_n)-np.ones((_n,_n))/_n
_G=-0.5*_Jc@(_Dj**2)@_Jc
_val,_vec=np.linalg.eigh(_G); _ord=_val.argsort()[::-1]; _val,_vec=_val[_ord],_vec[:,_ord]
_posv=np.maximum(_val,0); pcs=_vec[:,:2]*np.sqrt(_posv[:2]); ev=_posv/_posv.sum()*100
# clusters = K-MEANS sur l'embedding PCoA binaire (Jaccard), K issu du coude
from scipy.cluster.vq import kmeans2 as _km
K = 4   # <<< NOMBRE DE CLUSTERS (A) — modifie ce chiffre (coude conseille K=4)
_kd=min(6,_n); _emb=_vec[:,:_kd]*np.sqrt(_posv[:_kd])
_best=None
for _sd in range(40):
    _cc,_ll=_km(_emb,K,minit='++',seed=_sd)
    if len(set(_ll))<K: continue
    _ss=((_emb-_cc[_ll])**2).sum()
    if _best is None or _ss<_best[0]: _best=(_ss,_ll)
labels=_best[1]

# --- taille = conservation moyenne de l'interface actine (ponderee par enfouissement) ---
_cdf=pd.read_csv(_R/'data/proteocast/conservation_vs_asa_per_position.csv')
_cmap=dict(zip(_cdf['canon'].astype(int), _cdf['conservation']))
_cvec=np.array([_cmap.get(int(c),np.nan) for c in pivot.columns]); _msk=~np.isnan(_cvec)
_W=pivot.values.astype(float)
_num=(_W[:,_msk]*_cvec[_msk]).sum(1); _den=_W[:,_msk].sum(1)
cons_abp=np.where(_den>0,_num/np.where(_den>0,_den,1),np.nan)   # conservation moyenne ponderee de l'interface
cmin=np.nanmin(cons_abp); cmax=np.nanmax(cons_abp)
_t=(cons_abp-cmin)/(cmax-cmin); sizes=20+620*np.power(np.clip(_t,0,1),2.5); sizes=np.where(np.isnan(sizes),20,sizes)

cols=['#1f77b4','#d62728','#2ca02c','#ff7f0e']; shp={'+':'^','-':'v','partout':'o'}
fig,ax=plt.subplots(figsize=(15,11))
for cl in range(K):
    for t,mk in shp.items():
        m=(labels==cl)&(tend_arr==t)
        if m.any(): ax.scatter(pcs[m,0],pcs[m,1],marker=mk,color=cols[cl],s=sizes[m],edgecolor='k',lw=0.3)  # taille = conservation
import textwrap as _tw
_wrapn=lambda s: '\n'.join(_tw.wrap(str(s), 18))  # nom complet, replie a 18 car.
texts=[ax.text(pcs[i,0],pcs[i,1],_wrapn(nm),fontsize=6.5,ha='center') for i,nm in enumerate(pivot.index)]
adjust_text(texts,ax=ax,arrowprops=dict(arrowstyle='-',lw=0.4,color='grey'),
            expand_text=(1.06,1.12),expand_points=(1.06,1.12),
            force_text=(0.3,0.5),force_points=(0.12,0.2),lim=200)
sh=[Line2D([0],[0],marker=m,color='grey',ls='',ms=10,label=l) for l,m in [('barbé (+)','^'),('pointé (−)','v'),('partout','o')]]
co=[Line2D([0],[0],marker='s',color=cols[c],ls='',ms=10,label=f'cluster {c+1}') for c in range(K)]
ax.set_xlabel(f'PCoA1 ({ev[0]:.0f}%)'); ax.set_ylabel(f'PCoA2 ({ev[1]:.0f}%)')
ax.set_title('Empreinte des ABP — BINAIRE (presence/absence, Jaccard) + clustering K-MEANS\ncouleur=cluster · forme=côté filament · taille=conservation')
plt.tight_layout(); plt.show()
for cl in range(K):
    print(f'Cluster {cl+1} (n={int((labels==cl).sum())}):', ', '.join(pivot.index[labels==cl])); print()
_cs=pd.Series(cons_abp,index=pivot.index).dropna().sort_values(ascending=False)
print('Conservation interface — top 5 :'); print(_cs.head().round(3).to_string())
print('\nConservation interface — bottom 5 :'); print(_cs.tail().round(3).to_string())


# ### Comparaison — Option B (continu, corrélation)\n\nMême ABP mais regroupés par **intensité** denfouissement (%ASA continu) au lieu de présence/absence. À comparer avec la PCoA binaire ci-dessus.

# In[108]:


# === COMPARAISON — Option B : empreinte CONTINUE (%ASA) + heatmap correlation + k-means ===
# Reutilise pivot, tend_arr, sizes, shp, cols, adjust_text des cellules precedentes.
import seaborn as _sns
from scipy.spatial.distance import pdist as _pd2, squareform as _sq2
from scipy.cluster.hierarchy import linkage as _lk2
from scipy.cluster.vq import kmeans2 as _km2
import matplotlib as _mpl
from matplotlib.lines import Line2D

# --- (1) Heatmap CONTINUE, lignes groupees par correlation (= distance de l'option B) ---
_ZBh=_lk2(pivot.values, method='average', metric='correlation')
gB=_sns.clustermap(pivot, cmap='YlOrRd', figsize=(18,9), row_cluster=True, col_cluster=False,
                   row_linkage=_ZBh, dendrogram_ratio=(0.10,0.0), xticklabels=False, cbar_pos=None)
gB.ax_heatmap.set_yticks(np.arange(len(pivot))+0.5)
gB.ax_heatmap.set_yticklabels([pivot.index[i] for i in gB.dendrogram_row.reordered_ind], fontsize=12, rotation=0)
gB.ax_heatmap.yaxis.set_ticks_position('right'); gB.ax_heatmap.yaxis.set_label_position('right')  # noms a DROITE
gB.ax_heatmap.set_ylabel('')   # enlever le label 'abp'
_hp=gB.ax_heatmap.get_position()
gB.ax_heatmap.set_position([_hp.x0, _hp.y0, 0.58, _hp.height])              # retrecit la heatmap -> place pour les noms
gB.ax_row_dendrogram.set_position([0.92, _hp.y0, 0.2, _hp.height]); gB.ax_row_dendrogram.invert_xaxis()  # x=0.78 (position) · 0.16 = LONGUEUR du dendrogramme
gB.ax_heatmap.set_xlabel('Positions actine (canoniques)', fontsize=13)
gB.fig.suptitle("Option B — Heatmap %ASA, lignes groupees par CORRELATION (intensite d'enfouissement)",y=1.01,fontsize=12)
# --- barre de couleur (gradient % ASA) ---
_smB=_mpl.cm.ScalarMappable(cmap='YlOrRd', norm=_mpl.colors.Normalize(pivot.values.min(), pivot.values.max()))
_caxB=gB.fig.add_axes([0.6,0.00,0.26,0.028])   # [x, y, LARGEUR, HAUTEUR] de la barre
gB.fig.colorbar(_smB, cax=_caxB, orientation='horizontal', label='% ASA enfoui (moyen)')
_caxB.tick_params(labelsize=11); _caxB.xaxis.label.set_size(12)
plt.show()

# --- (2) PCoA correlation + clustering K-MEANS (sur l'embedding) ---
_Db=np.nan_to_num(_sq2(_pd2(pivot.values,metric='correlation')))
_nb=_Db.shape[0]; _Jb=np.eye(_nb)-np.ones((_nb,_nb))/_nb
_Gb=-0.5*_Jb@(_Db**2)@_Jb; _vb,_cb=np.linalg.eigh(_Gb); _ob=_vb.argsort()[::-1]; _vb,_cb=_vb[_ob],_cb[:,_ob]
_pvb=np.maximum(_vb,0); pcsB=_cb[:,:2]*np.sqrt(_pvb[:2]); evB=_pvb/_pvb.sum()*100
# --- coude propre a B (distance correlation) : choisit KB ---
from scipy.cluster.vq import kmeans2 as _km2b
_kdB=min(8,_nb); _embB_all=_cb[:,:_kdB]*np.sqrt(_pvb[:_kdB])
def _silb(D,lab):
    s=[]
    for i in range(len(lab)):
        sm=[D[i,j] for j in range(len(lab)) if lab[j]==lab[i] and j!=i]; a=np.mean(sm) if sm else 0
        ot=[np.mean([D[i,j] for j in range(len(lab)) if lab[j]==cc]) for cc in set(lab) if cc!=lab[i]]; b=min(ot) if ot else 0
        s.append(0 if max(a,b)==0 else (b-a)/max(a,b))
    return np.mean(s)
_KsB=list(range(1,11)); _wB=[]; _sB=[]
for _K in _KsB:
    _bb=None
    for _sd in range(20):
        _cc,_ll=_km2b(_embB_all,_K,minit='++',seed=_sd)
        if len(set(_ll))<_K and _K>1: continue
        _ss=((_embB_all-_cc[_ll])**2).sum()
        if _bb is None or _ss<_bb[0]: _bb=(_ss,_ll)
    _wB.append(_bb[0]); _sB.append(_silb(_Db,_bb[1]) if _K>1 else np.nan)
_dropB=[_wB[i]-_wB[i+1] for i in range(len(_wB)-1)]
KB = 3   # <<< NOMBRE DE CLUSTERS (B) — modifie ce chiffre (coude auto=2, K=3 a du sens bio)
figc,axe1=plt.subplots(figsize=(7,4))
axe1.plot(_KsB,_wB,'o-',color='#1f77b4'); # (pas de ligne K)
axe1.set_xlabel('K'); axe1.set_ylabel('inertie intra-cluster'); axe1.set_title('Option B — methode du coude (correlation)'); axe1.set_xticks(_KsB)
plt.tight_layout(); plt.show()
# KB defini par le coude ci-dessus
_kd=min(6,_nb); _embB=_cb[:,:_kd]*np.sqrt(_pvb[:_kd])
_best=None
for sd in range(40):
    cc,ll=_km2(_embB,KB,minit='++',seed=sd)
    if len(set(ll))<KB: continue
    ss=((_embB-cc[ll])**2).sum()
    if _best is None or ss<_best[0]: _best=(ss,ll)
labB=_best[1]
fig,ax=plt.subplots(figsize=(15,11))
for cl in range(KB):
    for t,mk in shp.items():
        m=(labB==cl)&(tend_arr==t)
        if m.any(): ax.scatter(pcsB[m,0],pcsB[m,1],marker=mk,color=cols[cl%len(cols)],s=sizes[m],edgecolor='k',lw=0.3)  # taille = conservation
import textwrap as _tw
_wrapn=lambda s: '\n'.join(_tw.wrap(str(s), 18))  # nom complet, replie a 18 car.
texts=[ax.text(pcsB[i,0],pcsB[i,1],_wrapn(nm),fontsize=9,ha='center') for i,nm in enumerate(pivot.index)]
adjust_text(texts,ax=ax,arrowprops=dict(arrowstyle='-',lw=0.4,color='grey'),
            expand_text=(1.06,1.12),expand_points=(1.06,1.12),
            force_text=(0.3,0.5),force_points=(0.12,0.2),lim=200)
ax.set_xlabel(f'PCoA1 ({evB[0]:.0f}%)'); ax.set_ylabel(f'PCoA2 ({evB[1]:.0f}%)')
from matplotlib.lines import Line2D as _L2
_hsh=[_L2([0],[0],marker=m,color='grey',ls='',ms=10,label=l) for l,m in [('barbé (+)','^'),('pointé (−)','v'),('partout','o')]]
_hcl=[_L2([0],[0],marker='s',color=cols[k%len(cols)],ls='',ms=10,label=f'cluster {k+1}') for k in range(KB)]
_hsz=[_L2([0],[0],marker='o',color='grey',ls='',ms=ss,label=lb) for ss,lb in [(5,'peu sensible'),(12,'très sensible')]]
_lgc=ax.legend(handles=_hcl,loc='upper left',fontsize=8,title='couleur = cluster',framealpha=.9)
ax.add_artist(_lgc)
_lgf=ax.legend(handles=_hsh,loc='upper right',fontsize=8,title='forme = côté filament',framealpha=.9)
ax.add_artist(_lgf)
ax.legend(handles=_hsz,loc='lower left',fontsize=8,title='taille = sensibilité mutationnelle',framealpha=.9)
ax.set_title(f"Option B — empreinte CONTINUE (%ASA, correlation) + clustering K-MEANS (K={KB})\n"
             "regroupe par INTENSITE d'enfouissement · forme=côté filament · taille=conservation")
plt.tight_layout(); plt.show()
for cl in range(KB):
    print(f'[B kmeans] Cluster {cl+1} (n={int((labB==cl).sum())}):', ', '.join(pivot.index[labB==cl])[:150])
print(f'\nB : PCoA1+2 = {evB[0]+evB[1]:.0f}% | k-means sur intensite %ASA. A = memes sites (binaire) ; B = meme intensite (continu).')


# ### La conservation differe-t-elle entre clusters dempreinte ? (justifie le facteur taille)\n\nLe marqueur = conservation moyenne de linterface ABP. Test : la conservation depend-elle du cluster (option B) ?

# In[104]:


# === La conservation differe-t-elle entre clusters d'empreinte (B) ? (justifie le facteur taille) ===
from scipy.stats import kruskal, mannwhitneyu
from itertools import combinations
_cv=np.asarray(cons_abp); _lb=np.asarray(labB)
_grps=[_cv[(_lb==k)&~np.isnan(_cv)] for k in range(KB)]
_H,_pk=kruskal(*_grps)
_lblnames={0:'cluster 1',1:'cluster 2',2:'cluster 3',3:'cluster 4'}
fig,ax=plt.subplots(figsize=(7,5))
bp=ax.boxplot(_grps,labels=[f'{_lblnames[k]}\n(n={len(_grps[k])})' for k in range(KB)],
              patch_artist=True,showmeans=True,meanprops=dict(marker='D',mfc='white',mec='k',ms=6))
for b,c in zip(bp['boxes'],cols): b.set_facecolor(c); b.set_alpha(.6)
ax.set_ylabel('conservation moyenne de l interface (ponderee enfouissement)')
ax.set_title(f"Conservation par cluster d'empreinte ABP (option B, K={KB})\nKruskal-Wallis : H={_H:.2f}, p={_pk:.3f}")
plt.tight_layout(); plt.show()
print(f'Kruskal-Wallis : H={_H:.2f}, p={_pk:.3f}'+('  -> SIGNIFICATIF' if _pk<0.05 else '  -> NS'))
for a,b in combinations(range(KB),2):
    print(f'  cluster{a+1} (med={np.median(_grps[a]):.2f}) vs cluster{b+1} (med={np.median(_grps[b]):.2f}) : MW p={mannwhitneyu(_grps[a],_grps[b])[1]:.3f}')
print("\nLecture : le facteur TAILLE (conservation) porte un vrai signal -> a garder.")


# In[105]:


# === PCA biophysique (ordination) : 11 descripteurs d'interface cote actine ===
import warnings; warnings.filterwarnings('ignore')
import re as _re
from pathlib import Path as _P
from matplotlib.lines import Line2D
from scipy.spatial.distance import pdist, squareform
_R = _P('..') if _P('../data').exists() else _P('.')

def _arpn(n):
    if 'complex subunit 5-like' in n: return 'ARPC5L'
    mc=_re.search(r'complex subunit (\d+[A-Za-z]?)',n)
    if mc: return 'ARPC'+mc.group(1)
    if _re.match(r'Actin-related protein 3(\b|,|$)',n): return 'Arp3'
    if _re.match(r'Actin-related protein 2(\b|,|$)',n): return 'Arp2'
    return n
_df3=pd.read_csv(_R/'data/filtered/details/3.interface_residues.csv')
_df3['pct']=pd.to_numeric(_df3['buried_ASA_percent'].astype(str).str.replace('%','',regex=False),errors='coerce')
_df3['a2']=pd.to_numeric(_df3['buried_ASA_Å²'].astype(str).str.replace('<0.1','0.05',regex=False),errors='coerce')
_df3['canon']=pd.to_numeric(_df3['residue_number_canon_mafft'],errors='coerce')
_df3=_df3[_df3.canon.notna()&_df3.pct.notna()].copy()
_di=pd.read_csv(_R/'data/filtered/details/1.interactions.csv')
_da=pd.read_csv(_R/'data/filtered/filtered_all_data.csv',low_memory=False)
_pp=pd.read_csv(_R/'data/filtered/proteins_per_pdb.csv'); _ac=set(_pp[_pp.is_actin]['chain'].str.lower())
_homo=set(_di[_di['chain_B_id'].str.lower().isin(_ac)]['interaction_id'])
_het=_di[~_di.interaction_id.isin(_homo)].merge(
    _da[['subunit_1','subunit_2','subunit_2_title','s2_actine','cluster_data_70','s1_binding_site_cluster_data_70']],
    left_on=['chain_A_id','chain_B_id'],right_on=['subunit_1','subunit_2'],how='left').drop_duplicates('interaction_id')
_het=_het[_het['s2_actine'].fillna(False)==False].copy()
_het['abp']=_het['subunit_2_title'].fillna('Unknown').str.replace(r'\s*\(.*?\)','',regex=True).str.strip().str[:50].map(_arpn)
# un "site" = cluster de site de liaison S1 (fallback = cluster_data_70 si manquant)
_het['site']=_het['s1_binding_site_cluster_data_70'].fillna('c70_'+_het['cluster_data_70'].astype(str))
_im=_het.set_index('interaction_id'); _s1=_im['chain_A_id'].str.lower(); _ab=_im['abp']; _st=_im['site']
_dd=_df3[_df3.interaction_id.isin(set(_het.interaction_id))].copy()
_dd['_c']=_dd.interaction_id.map(_s1); _dd=_dd[_dd.chain.str.lower()==_dd._c].copy()   # cote actine
_dd['abp']=_dd.interaction_id.map(_ab); _dd['site']=_dd.interaction_id.map(_st); _dd['canon']=_dd.canon.astype(int)
_cdf=pd.read_csv(_R/'data/proteocast/conservation_vs_asa_per_position.csv'); _cons=dict(zip(_cdf['canon'].astype(int),_cdf['conservation']))
_KD={'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
def _mode(s): m=s.mode(); return m.iloc[0] if len(m) else None
# --- etage 1 : par (abp, site, canon) -> enfouissement moyen de la position dans le site ---
_pos=_dd.groupby(['abp','site','canon']).agg(w=('pct','mean'),a2=('a2','mean'),res=('residue_name',_mode)).reset_index()
_pos['kd']=_pos['res'].map(lambda r:_KD.get(str(r),0.0)); _pos['cons']=_pos['canon'].map(lambda c:_cons.get(int(c),np.nan))
_pos['ispos']=_pos['res'].isin({'R','K'}).astype(float); _pos['isneg']=_pos['res'].isin({'D','E'}).astype(float); _pos['isaro']=_pos['res'].isin({'F','Y','W'}).astype(float)
def _wm(g,col):
    w=g['w'].values; x=g[col].values; m=~np.isnan(x)
    return (w[m]*x[m]).sum()/w[m].sum() if w[m].sum()>0 else np.nan
# --- etage 1bis : par (abp, site) moyennes ponderees enfouissement ; etage 2 : moyenne entre sites ---
_site=_pos.groupby(['abp','site']).apply(lambda g:pd.Series({
    'conservation_actine':_wm(g,'cons'),'hydrophobicite':_wm(g,'kd'),'charge_pos':_wm(g,'ispos'),
    'charge_neg':_wm(g,'isneg'),'aromaticite':_wm(g,'isaro'),'asa_moyen':_wm(g,'w'),'taille_interface':g['a2'].sum()})).reset_index()
F=_site.groupby('abp')[['conservation_actine','hydrophobicite','charge_pos','charge_neg','aromaticite','asa_moyen']].mean()
F['taille_interface']=_site.groupby('abp')['taille_interface'].max()   # la plus grosse interface de l'ABP (pas la moyenne/somme)
# ratios H-bond / pont salin : par site puis moyenne entre sites
_hb=_het.groupby(['abp','site'])[['num_contacts','num_hbonds','num_salt_bridges']].sum().reset_index()
_hb['pct_hbond']=100*_hb['num_hbonds']/_hb['num_contacts']; _hb['pct_saltbridge']=100*_hb['num_salt_bridges']/_hb['num_contacts']
F=F.join(_hb.groupby('abp')[['pct_hbond','pct_saltbridge']].mean())
# position ordinale (aussi utilisee en couleur)
_fil=pd.read_csv(_R/'data/filtered/actin_filament_positions.csv'); _c2l=dict(zip(_fil['chain'].str.lower(),_fil['label']))
_het['fil']=_het['chain_A_id'].str.lower().map(_c2l)
_het['grp']=_het['fil'].map(lambda l: l if l in ('+','-') else ('side' if pd.notna(l) else None))
_ps=_het.dropna(subset=['grp']).groupby('abp')['grp'].apply(set)
F['position']=pd.Series({a:(1.0 if ('+' in _ps.get(a,set()) and '-' not in _ps.get(a,set())) else (-1.0 if ('-' in _ps.get(a,set()) and '+' not in _ps.get(a,set())) else 0.0)) for a in F.index})
F['n_sites']=_het.groupby('abp')['cluster_data_70'].nunique()   # nb de clusters d'interface C70 de l'ABP
_abs1=_het.dropna(subset=['s1_binding_site_cluster_data_70']).groupby('abp')['s1_binding_site_cluster_data_70'].apply(set)
_s1a={}
for a,ss in _abs1.items():
    for s in ss: _s1a.setdefault(s,set()).add(a)
F['n_competiteurs']=pd.Series({a:len({x for s in _abs1.get(a,set()) for x in _s1a.get(s,set())}-{a}) for a in F.index})
F=F[['conservation_actine','hydrophobicite','charge_pos','charge_neg','aromaticite','taille_interface','pct_hbond','pct_saltbridge','position','n_sites','n_competiteurs']]
_miss=F[F.isna().any(axis=1)].index.tolist(); F=F.fillna(F.median(numeric_only=True))
print(f"{len(F)} ABP, {F.shape[1]} descripteurs (agreg. 2 etages : par site -> pondere enfouissement -> moyenne entre sites). Imputes: {_miss if _miss else 'aucun'}")

# z-score + PCA
X=F.values.astype(float); Xz=(X-X.mean(0))/(X.std(0)+1e-9)
Xc=Xz-Xz.mean(0); U,S,Vt=np.linalg.svd(Xc,full_matrices=False); pcs=U[:,:2]*S[:2]; ev=(S**2)/(S**2).sum()*100
print(f"Variance PC1+PC2 = {ev[0]:.0f}% + {ev[1]:.0f}% = {ev[0]+ev[1]:.0f}%")

# couleur = cluster K-MEANS de l'Option B (empreinte continue %ASA) -> memes 3 groupes/couleurs
_pcol={1.0:'#E69F00',-1.0:'#0072B2',0.0:'#999999'}; _plab={1.0:'barbé (+)',-1.0:'pointé (−)',0.0:'latéral'}
_abp2cl=dict(zip(pivot.index,labB)); _clset=set(pivot.index)   # cluster de l'Option B (cellule k-means)
_pc=np.array([cols[_abp2cl[a]%len(cols)] if a in _clset else '#cccccc' for a in F.index])

# === Biplot : individus (couleur = position) + variables (fleches) ===
cor=np.array([[np.corrcoef(Xz[:,j],pcs[:,k])[0,1] for k in range(2)] for j in range(F.shape[1])])
_f=0.72*np.abs(pcs).max()/np.abs(cor).max()
fig,axb=plt.subplots(figsize=(11,9))
axb.scatter(pcs[:,0],pcs[:,1],c=_pc,s=55,edgecolor='k',lw=.3,zorder=2)
# labels anti-chevauchement : ABP peripheriques, en reservant les emplacements des labels de fleches
_dist=np.hypot(pcs[:,0],pcs[:,1]); _thr=np.percentile(_dist,55)
_xr=pcs[:,0].max()-pcs[:,0].min(); _yr=pcs[:,1].max()-pcs[:,1].min()
_placed=[(cor[j,0]*_f*1.06,cor[j,1]*_f*1.06) for j in range(F.shape[1])]
def _free(x,y,mx=0.06,my=0.045):
    return all(not(abs(x-px)<mx*_xr and abs(y-py)<my*_yr) for px,py in _placed)
for i in np.argsort(-_dist):
    if _dist[i]<_thr: break
    x,y=pcs[i,0],pcs[i,1]
    if not _free(x,y): continue
    _placed.append((x,y)); sd=-1 if x>0 else 1
    axb.annotate(F.index[i][:16],(x,y),fontsize=6.5,zorder=4,ha='right' if sd<0 else 'left',
                 xytext=(6*sd,3),textcoords='offset points',bbox=dict(boxstyle='round,pad=0.1',fc='white',ec='none',alpha=0.6))
for j,nm in enumerate(F.columns):
    axb.arrow(0,0,cor[j,0]*_f,cor[j,1]*_f,color='#c0392b',head_width=.12,length_includes_head=True,lw=1.4,zorder=3)
    axb.annotate(nm,(cor[j,0]*_f*1.06,cor[j,1]*_f*1.06),fontsize=8,color='#7b241c',ha='center',zorder=3)
axb.axhline(0,color='grey',lw=.5); axb.axvline(0,color='grey',lw=.5)
axb.set_xlabel(f'PC1 ({ev[0]:.0f}%)'); axb.set_ylabel(f'PC2 ({ev[1]:.0f}%)')
axb.set_title('Biplot biophysique des ABP : les profils d\'interface se séparent en physico-chimie\ncouleur = cluster K-MEANS (Option B) · flèches = descripteurs')
axb.legend(handles=[Line2D([0],[0],marker='o',ls='',mfc=cols[k%len(cols)],mec='k',ms=9,label=f'cluster {k+1}') for k in range(KB)],loc='best',fontsize=8,title='couleur = cluster (Option B)')
plt.tight_layout(); plt.show()

# silhouette (confirme que c'est une ordination, pas des clusters nets)
from scipy.cluster.vq import kmeans2
_D=squareform(pdist(Xz))
def _sil(D,lab):
    n=len(lab); s=np.zeros(n)
    for i in range(n):
        sm=(lab==lab[i]).copy(); sm[i]=False
        a=D[i,sm].mean() if sm.any() else 0
        b=min((D[i,lab==c].mean() for c in set(lab) if c!=lab[i]),default=0)
        s[i]=0 if max(a,b)==0 else (b-a)/max(a,b)
    return s.mean()
_bb=None
for sd in range(60):
    c,l=kmeans2(Xz,2,minit='++',seed=sd)
    if len(set(l))<2: continue
    ss=((Xz-c[l])**2).sum()
    if _bb is None or ss<_bb[0]: _bb=(ss,l)
print(f"Silhouette k-means k=2 = {_sil(_D,_bb[1]):.3f}  (vs ~0.30 pour la PCA d'empreinte -> ici gradient, pas clusters nets)")


# ### Variante — ABP biophysique, COEUR (résidus > 25 % ASA)
# 
# Même PCA des ABP mais descripteurs d'interface calculés sur les seuls résidus d'actine enfouis > 25 % ASA.

# In[106]:


# === PCA biophysique (ordination) : 11 descripteurs d'interface cote actine ===
import warnings; warnings.filterwarnings('ignore')
import re as _re
from pathlib import Path as _P
from matplotlib.lines import Line2D
from scipy.spatial.distance import pdist, squareform
_R = _P('..') if _P('../data').exists() else _P('.')

def _arpn(n):
    if 'complex subunit 5-like' in n: return 'ARPC5L'
    mc=_re.search(r'complex subunit (\d+[A-Za-z]?)',n)
    if mc: return 'ARPC'+mc.group(1)
    if _re.match(r'Actin-related protein 3(\b|,|$)',n): return 'Arp3'
    if _re.match(r'Actin-related protein 2(\b|,|$)',n): return 'Arp2'
    return n
_df3=pd.read_csv(_R/'data/filtered/details/3.interface_residues.csv')
_df3['pct']=pd.to_numeric(_df3['buried_ASA_percent'].astype(str).str.replace('%','',regex=False),errors='coerce')
_df3['a2']=pd.to_numeric(_df3['buried_ASA_Å²'].astype(str).str.replace('<0.1','0.05',regex=False),errors='coerce')
_df3['canon']=pd.to_numeric(_df3['residue_number_canon_mafft'],errors='coerce')
_df3=_df3[_df3.canon.notna()&_df3.pct.notna()].copy()
_di=pd.read_csv(_R/'data/filtered/details/1.interactions.csv')
_da=pd.read_csv(_R/'data/filtered/filtered_all_data.csv',low_memory=False)
_pp=pd.read_csv(_R/'data/filtered/proteins_per_pdb.csv'); _ac=set(_pp[_pp.is_actin]['chain'].str.lower())
_homo=set(_di[_di['chain_B_id'].str.lower().isin(_ac)]['interaction_id'])
_het=_di[~_di.interaction_id.isin(_homo)].merge(
    _da[['subunit_1','subunit_2','subunit_2_title','s2_actine','cluster_data_70','s1_binding_site_cluster_data_70']],
    left_on=['chain_A_id','chain_B_id'],right_on=['subunit_1','subunit_2'],how='left').drop_duplicates('interaction_id')
_het=_het[_het['s2_actine'].fillna(False)==False].copy()
_het['abp']=_het['subunit_2_title'].fillna('Unknown').str.replace(r'\s*\(.*?\)','',regex=True).str.strip().str[:50].map(_arpn)
# un "site" = cluster de site de liaison S1 (fallback = cluster_data_70 si manquant)
_het['site']=_het['s1_binding_site_cluster_data_70'].fillna('c70_'+_het['cluster_data_70'].astype(str))
_im=_het.set_index('interaction_id'); _s1=_im['chain_A_id'].str.lower(); _ab=_im['abp']; _st=_im['site']
_dd=_df3[_df3.interaction_id.isin(set(_het.interaction_id))].copy()
_dd['_c']=_dd.interaction_id.map(_s1); _dd=_dd[(_dd.chain.str.lower()==_dd._c)&(_dd['pct']>25)].copy()   # cote actine COEUR >25% ASA
_dd['abp']=_dd.interaction_id.map(_ab); _dd['site']=_dd.interaction_id.map(_st); _dd['canon']=_dd.canon.astype(int)
_cdf=pd.read_csv(_R/'data/proteocast/conservation_vs_asa_per_position.csv'); _cons=dict(zip(_cdf['canon'].astype(int),_cdf['conservation']))
_KD={'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
def _mode(s): m=s.mode(); return m.iloc[0] if len(m) else None
# --- etage 1 : par (abp, site, canon) -> enfouissement moyen de la position dans le site ---
_pos=_dd.groupby(['abp','site','canon']).agg(w=('pct','mean'),a2=('a2','mean'),res=('residue_name',_mode)).reset_index()
_pos['kd']=_pos['res'].map(lambda r:_KD.get(str(r),0.0)); _pos['cons']=_pos['canon'].map(lambda c:_cons.get(int(c),np.nan))
_pos['ispos']=_pos['res'].isin({'R','K'}).astype(float); _pos['isneg']=_pos['res'].isin({'D','E'}).astype(float); _pos['isaro']=_pos['res'].isin({'F','Y','W'}).astype(float)
def _wm(g,col):
    w=g['w'].values; x=g[col].values; m=~np.isnan(x)
    return (w[m]*x[m]).sum()/w[m].sum() if w[m].sum()>0 else np.nan
# --- etage 1bis : par (abp, site) moyennes ponderees enfouissement ; etage 2 : moyenne entre sites ---
_site=_pos.groupby(['abp','site']).apply(lambda g:pd.Series({
    'conservation_actine':_wm(g,'cons'),'hydrophobicite':_wm(g,'kd'),'charge_pos':_wm(g,'ispos'),
    'charge_neg':_wm(g,'isneg'),'aromaticite':_wm(g,'isaro'),'asa_moyen':_wm(g,'w'),'taille_interface':g['a2'].sum()})).reset_index()
F=_site.groupby('abp')[['conservation_actine','hydrophobicite','charge_pos','charge_neg','aromaticite','asa_moyen']].mean()
F['taille_interface']=_site.groupby('abp')['taille_interface'].max()   # la plus grosse interface de l'ABP (pas la moyenne/somme)
# ratios H-bond / pont salin : par site puis moyenne entre sites
_hb=_het.groupby(['abp','site'])[['num_contacts','num_hbonds','num_salt_bridges']].sum().reset_index()
_hb['pct_hbond']=100*_hb['num_hbonds']/_hb['num_contacts']; _hb['pct_saltbridge']=100*_hb['num_salt_bridges']/_hb['num_contacts']
F=F.join(_hb.groupby('abp')[['pct_hbond','pct_saltbridge']].mean())
# position ordinale (aussi utilisee en couleur)
_fil=pd.read_csv(_R/'data/filtered/actin_filament_positions.csv'); _c2l=dict(zip(_fil['chain'].str.lower(),_fil['label']))
_het['fil']=_het['chain_A_id'].str.lower().map(_c2l)
_het['grp']=_het['fil'].map(lambda l: l if l in ('+','-') else ('side' if pd.notna(l) else None))
_ps=_het.dropna(subset=['grp']).groupby('abp')['grp'].apply(set)
F['position']=pd.Series({a:(1.0 if ('+' in _ps.get(a,set()) and '-' not in _ps.get(a,set())) else (-1.0 if ('-' in _ps.get(a,set()) and '+' not in _ps.get(a,set())) else 0.0)) for a in F.index})
F['n_sites']=_het.groupby('abp')['cluster_data_70'].nunique()   # nb de clusters d'interface C70 de l'ABP
_abs1=_het.dropna(subset=['s1_binding_site_cluster_data_70']).groupby('abp')['s1_binding_site_cluster_data_70'].apply(set)
_s1a={}
for a,ss in _abs1.items():
    for s in ss: _s1a.setdefault(s,set()).add(a)
F['n_competiteurs']=pd.Series({a:len({x for s in _abs1.get(a,set()) for x in _s1a.get(s,set())}-{a}) for a in F.index})
F=F[['conservation_actine','hydrophobicite','charge_pos','charge_neg','aromaticite','taille_interface','pct_hbond','pct_saltbridge','position','n_sites','n_competiteurs']]
_miss=F[F.isna().any(axis=1)].index.tolist(); F=F.fillna(F.median(numeric_only=True))
print(f"{len(F)} ABP, {F.shape[1]} descripteurs (agreg. 2 etages : par site -> pondere enfouissement -> moyenne entre sites). Imputes: {_miss if _miss else 'aucun'}")

# z-score + PCA
X=F.values.astype(float); Xz=(X-X.mean(0))/(X.std(0)+1e-9)
Xc=Xz-Xz.mean(0); U,S,Vt=np.linalg.svd(Xc,full_matrices=False); pcs=U[:,:2]*S[:2]; ev=(S**2)/(S**2).sum()*100
print(f"Variance PC1+PC2 = {ev[0]:.0f}% + {ev[1]:.0f}% = {ev[0]+ev[1]:.0f}%")

# couleur = cluster K-MEANS de l'Option B (empreinte continue %ASA) -> memes 3 groupes/couleurs
_pcol={1.0:'#E69F00',-1.0:'#0072B2',0.0:'#999999'}; _plab={1.0:'barbé (+)',-1.0:'pointé (−)',0.0:'latéral'}
_abp2cl=dict(zip(pivot.index,labB)); _clset=set(pivot.index)   # cluster de l'Option B (cellule k-means)
_pc=np.array([cols[_abp2cl[a]%len(cols)] if a in _clset else '#cccccc' for a in F.index])

# === Biplot : individus (couleur = position) + variables (fleches) ===
cor=np.array([[np.corrcoef(Xz[:,j],pcs[:,k])[0,1] for k in range(2)] for j in range(F.shape[1])])
_f=0.72*np.abs(pcs).max()/np.abs(cor).max()
fig,axb=plt.subplots(figsize=(11,9))
axb.scatter(pcs[:,0],pcs[:,1],c=_pc,s=55,edgecolor='k',lw=.3,zorder=2)
# labels anti-chevauchement : ABP peripheriques, en reservant les emplacements des labels de fleches
_dist=np.hypot(pcs[:,0],pcs[:,1]); _thr=np.percentile(_dist,55)
_xr=pcs[:,0].max()-pcs[:,0].min(); _yr=pcs[:,1].max()-pcs[:,1].min()
_placed=[(cor[j,0]*_f*1.06,cor[j,1]*_f*1.06) for j in range(F.shape[1])]
def _free(x,y,mx=0.06,my=0.045):
    return all(not(abs(x-px)<mx*_xr and abs(y-py)<my*_yr) for px,py in _placed)
for i in np.argsort(-_dist):
    if _dist[i]<_thr: break
    x,y=pcs[i,0],pcs[i,1]
    if not _free(x,y): continue
    _placed.append((x,y)); sd=-1 if x>0 else 1
    axb.annotate(F.index[i][:16],(x,y),fontsize=6.5,zorder=4,ha='right' if sd<0 else 'left',
                 xytext=(6*sd,3),textcoords='offset points',bbox=dict(boxstyle='round,pad=0.1',fc='white',ec='none',alpha=0.6))
for j,nm in enumerate(F.columns):
    axb.arrow(0,0,cor[j,0]*_f,cor[j,1]*_f,color='#c0392b',head_width=.12,length_includes_head=True,lw=1.4,zorder=3)
    axb.annotate(nm,(cor[j,0]*_f*1.06,cor[j,1]*_f*1.06),fontsize=8,color='#7b241c',ha='center',zorder=3)
axb.axhline(0,color='grey',lw=.5); axb.axvline(0,color='grey',lw=.5)
axb.set_xlabel(f'PC1 ({ev[0]:.0f}%)'); axb.set_ylabel(f'PC2 ({ev[1]:.0f}%)')
axb.set_title('Biplot biophysique des ABP (COEUR >25% ASA) : les profils d\'interface se séparent en physico-chimie\ncouleur = cluster K-MEANS (Option B) · flèches = descripteurs')
axb.legend(handles=[Line2D([0],[0],marker='o',ls='',mfc=cols[k%len(cols)],mec='k',ms=9,label=f'cluster {k+1}') for k in range(KB)],loc='best',fontsize=8,title='couleur = cluster (Option B)')
plt.tight_layout(); plt.show()

# silhouette (confirme que c'est une ordination, pas des clusters nets)
from scipy.cluster.vq import kmeans2
_D=squareform(pdist(Xz))
def _sil(D,lab):
    n=len(lab); s=np.zeros(n)
    for i in range(n):
        sm=(lab==lab[i]).copy(); sm[i]=False
        a=D[i,sm].mean() if sm.any() else 0
        b=min((D[i,lab==c].mean() for c in set(lab) if c!=lab[i]),default=0)
        s[i]=0 if max(a,b)==0 else (b-a)/max(a,b)
    return s.mean()
_bb=None
for sd in range(60):
    c,l=kmeans2(Xz,2,minit='++',seed=sd)
    if len(set(l))<2: continue
    ss=((Xz-c[l])**2).sum()
    if _bb is None or ss<_bb[0]: _bb=(ss,l)
print(f"Silhouette k-means k=2 = {_sil(_D,_bb[1]):.3f}  (vs ~0.30 pour la PCA d'empreinte -> ici gradient, pas clusters nets)")


# ## Comment lire — correspondance figure ⇄ test
# 
# - **Fig 1** : flèche `conservation` opposée à `rsa` & `nb familles` (qui sont collées) → enfouissement gouverne, nb de partenaires va avec l'exposition.
# - **Fig 2** : le gradient de couleur conservation est l'inverse de celui de rsa/n_fam.
# - **Fig 3** : `rsa` est la seule barre forte et stable ; `nb familles`/`nb ABP` **changent de signe** à rsa égal (artefact) ; `hydrophobicité` ≈ 0.
# 
# **Phrase défendable :** la conservation d'une position d'actine est gouvernée par son **enfouissement** ;
# ni le nombre de partenaires, ni l'hydrophobicité n'ajoutent d'effet propre une fois l'exposition contrôlée.
