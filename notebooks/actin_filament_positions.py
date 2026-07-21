#!/usr/bin/env python
# coding: utf-8

import os as _os
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # cwd = racine projet (robuste, peu importe d'où on lance)


# # Caractérisation des positions des actines dans les filaments
# 
# Pour chaque actine de chaque PDB (filaments de 4+ actines uniquement), on attribue un label :
# - **(+)** : extrémité barbée — profil {6685_1, 6685_3}, manque 6685_2 et 6685_4
# - **(+2)** : adjacent à (+), 3 interactions actin–actin (moins de voisins)
# - **(+3)** : adjacent à (+), 4 interactions actin–actin (plus de voisins)
# - **(-)** : extrémité pointée — profil {6685_2, 6685_4}, manque 6685_1 et 6685_3
# - **(-2)** : adjacent à (-), 3 interactions actin–actin (moins de voisins)
# - **(-3)** : adjacent à (-), 4 interactions actin–actin (plus de voisins)
# - **(side)** : tout le reste (intérieur, adjacent aux deux extrémités, composant ≤ 3)
# 
# Filaments de ≤ 3 actines exclus.
# 
# **Validation** : profiline et formines (INF2) se trouvent côté (+) barbé dans 9azp.

# In[1]:


import pandas as pd
from collections import defaultdict

df = pd.read_csv('data/filtered/filtered_all_data.csv', low_memory=False)
pp = pd.read_csv('data/filtered/proteins_per_pdb.csv')

print('Interactions totales :', len(df))
print('Protéines par PDB :', len(pp))


# In[2]:


# Interactions actin-actin uniquement (AUCUNE exclusion : chaque actine est
# assignée à son/ses méga-cluster(s), même structures atypiques)
aa = df[df['s1_actine'].fillna(False) & df['s2_actine'].fillna(False)].copy()
aa = aa[['pdb_id', 'subunit_1', 'subunit_2', 's1_supercluster', 's2_supercluster']].drop_duplicates()
print('Interactions actin-actin :', len(aa), '| PDB :', aa['pdb_id'].nunique())


# In[3]:


# Pour chaque actine : superclusters occupés et nombre d'interactions actin-actin
# chain → {occupied_superclusters, n_actin, neighbors}

# ── Repliement des superclusters variants sur leur canonique ─────────────────
# Pour l'assignation des positions, seules 4 interfaces existent :
#   6685_1, 6685_2, 6685_3, 6685_4
# Les superclusters minoritaires (6685_23, 6685_109, 6685_274) sont des
# variantes géométriques des mêmes contacts physiques. On l'identifie via
# l'appariement des deux faces de chaque contact actine-actine :
#   6685_109 <-> 6685_2  (comme 6685_1)  → 6685_109 ≡ 6685_1
#   6685_274 <-> 6685_2  (comme 6685_1)  → 6685_274 ≡ 6685_1
#   6685_23  <-> 6685_4  (comme 6685_3)  → 6685_23  ≡ 6685_3
# Sans ce repliement, une actine dont un contact tombe dans une variante paraît
# « manquer » une interface canonique et bascule du mauvais côté (+/-).
CANON_SC = {'6685_109': '6685_1', '6685_274': '6685_1', '6685_23': '6685_3'}

def _canon(sc_str):
    return {CANON_SC.get(x, x) for x in str(sc_str).split(';')}

# Structure : par (pdb_id, chain)
occupied_sc  = defaultdict(set)   # (pdb, chain) → set de superclusters
n_actin_int  = defaultdict(int)   # (pdb, chain) → nb d'interactions actin-actin
neighbors    = defaultdict(set)   # (pdb, chain) → set de chains voisines actin

for _, row in aa.iterrows():
    pdb = row['pdb_id']
    s1, s2 = row['subunit_1'], row['subunit_2']
    sc1, sc2 = row['s1_supercluster'], row['s2_supercluster']

    # s1 : son binding site = sc1
    if pd.notna(sc1):
        occupied_sc[(pdb, s1)].update(_canon(sc1))
    n_actin_int[(pdb, s1)] += 1
    neighbors[(pdb, s1)].add(s2)

    # s2 : son binding site = sc2
    if pd.notna(sc2):
        occupied_sc[(pdb, s2)].update(_canon(sc2))
    n_actin_int[(pdb, s2)] += 1
    neighbors[(pdb, s2)].add(s1)

print('Actines avec au moins 1 interaction actin-actin :', len(n_actin_int))


# In[4]:


# Actines de référence (TOUTES, aucune exclusion de PDB)
actin_chains = pp[pp['is_actin']].copy()
actin_chains['key'] = list(zip(actin_chains['pdb_id'], actin_chains['chain']))
print('Actines uniques :', len(actin_chains))

from collections import deque

def _connected_components(nodes, adj):
    seen = set()
    components = []
    for n in nodes:
        if n in seen:
            continue
        comp = set()
        queue = deque([n])
        while queue:
            cur = queue.popleft()
            if cur in seen:
                continue
            seen.add(cur)
            comp.add(cur)
            for nb in adj.get(cur, set()):
                if nb not in seen:
                    queue.append(nb)
        components.append(comp)
    return components

# Exclusion : composants ≤ 3 (filaments trop courts)
forced_side = set()
component_size = {}
for pdb, grp in actin_chains.groupby('pdb_id'):
    chains = list(grp['chain'])
    adj_pdb = {c: neighbors.get((pdb, c), set()) for c in chains}
    comps = _connected_components(chains, adj_pdb)
    for comp in comps:
        sz = len(comp)
        for c in comp:
            component_size[(pdb, c)] = sz
            if sz <= 3:
                forced_side.add((pdb, c))

print(f'Actines forcées side (composant ≤ 3) : {len(forced_side)}')
from collections import Counter
print('Distribution tailles composants:', Counter(component_size.values()))


# In[5]:


# ── Étape 1 : label (-) et (+) ───────────────────────────────────────────────
label = {}

for _, row in actin_chains.iterrows():
    key = (row['pdb_id'], row['chain'])
    osc = occupied_sc[key]

    if key in forced_side:
        label[key] = 'side'
        continue

    # (+) barbée : profil {6685_1, 6685_3} — manque 6685_2 et 6685_4
    is_p = ('6685_2' not in osc) and ('6685_4' not in osc)
    # (-) pointée : profil {6685_2, 6685_4} — manque 6685_1 et 6685_3
    is_m = ('6685_1' not in osc) and ('6685_3' not in osc)

    if is_m and is_p:
        label[key] = 'side'
    elif is_m:
        label[key] = '-'
    elif is_p:
        label[key] = '+'
    else:
        label[key] = None

n_m = sum(1 for v in label.values() if v == '-')
n_p = sum(1 for v in label.values() if v == '+')
n_side = sum(1 for v in label.values() if v == 'side')
n_unlabeled = sum(1 for v in label.values() if v is None)
print(f'(-): {n_m}  (+): {n_p}  side: {n_side}  non-labelisés: {n_unlabeled}')


# In[6]:


# ── Étape 2 : -2/-3/+2/+3/side sur les non-labelisés ────────────────────────
for key, lbl in label.items():
    if lbl is not None:
        continue
    pdb = key[0]
    nbrs = neighbors[key]
    n = n_actin_int[key]

    adj_m = any(label.get((pdb, nb)) == '-' for nb in nbrs)
    adj_p = any(label.get((pdb, nb)) == '+' for nb in nbrs)

    if adj_m and adj_p:
        label[key] = 'side'
    elif adj_m:
        if n == 3:
            label[key] = '-2'
        elif n == 4:
            label[key] = '-3'
        else:
            label[key] = 'side'
    elif adj_p:
        if n == 3:
            label[key] = '+2'
        elif n == 4:
            label[key] = '+3'
        else:
            label[key] = 'side'
    else:
        label[key] = 'side'

from collections import Counter
print('Distribution des labels :')
print(Counter(label.values()))


# In[7]:


# ── Sauvegarde CSV ────────────────────────────────────────────────────────────
rows = []
for _, row in actin_chains.iterrows():
    key = (row['pdb_id'], row['chain'])
    rows.append({
        'pdb_id': row['pdb_id'],
        'chain': row['chain'],
        'protein': row['protein'],
        'n_actin_interactions': n_actin_int[key],
        'occupied_superclusters': ','.join(sorted(occupied_sc[key])) if occupied_sc[key] else '',
        'component_size': component_size.get(key, 0),
        'label': label[key],
    })

out = pd.DataFrame(rows).sort_values(['pdb_id', 'chain']).reset_index(drop=True)
out.to_csv('data/filtered/actin_filament_positions.csv', index=False)
print('Sauvegardé : data/filtered/actin_filament_positions.csv')
print(f'{len(out)} lignes')
out.head(20)


# In[8]:


# ── Validation automatique de la convention +/- ──────────────────────────────
# Protéines de référence biologiquement connues :
#   barbée (+) : Profilin-1, F-actin-capping protein alpha/beta
#   pointée (-) : Tropomodulin-1, Tropomodulin-1.9
_BARBED_REF  = ['Profilin-1', 'F-actin-capping protein subunit alpha-1',
                'F-actin-capping protein subunit beta']
_POINTED_REF = ['Tropomodulin-1', 'Tropomodulin-1.9']

_label_map = {(r['pdb_id'], r['chain']): r['label'] for _, r in out.iterrows()}
_pp_abp = pp[~pp['is_actin']].copy()
_df_hetero = df[(df['s1_actine'].fillna(False) & ~df['s2_actine'].fillna(False)) |
                (~df['s1_actine'].fillna(False) & df['s2_actine'].fillna(False))].copy()

errors = []
for ref_list, expected, side_name in [
    (_BARBED_REF,  {'+', '+2', '+3'}, 'barbée (+)'),
    (_POINTED_REF, {'-', '-2', '-3'}, 'pointée (-)'),
]:
    for prot in ref_list:
        chains = _pp_abp[_pp_abp['protein'] == prot]['chain'].tolist()
        if not chains:
            continue
        for ch in chains:
            pdb = ch.split('_')[0]
            rows = _df_hetero[
                (_df_hetero['subunit_1'] == ch) | (_df_hetero['subunit_2'] == ch)
            ]
            for _, r in rows.iterrows():
                actin_ch = r['subunit_1'] if r['s2_actine'] else r['subunit_2']
                lbl = _label_map.get((pdb, actin_ch))
                if lbl and lbl not in expected and lbl != 'side':
                    errors.append(
                        f'  {prot} ({ch}) → actine {actin_ch} labellisée "{lbl}" '
                        f'mais attendu côté {side_name}'
                    )

if errors:
    print('❌ CONVENTION INVALIDE — inverser is_p / is_m dans cell-5 :')
    for e in errors[:10]:
        print(e)
    raise AssertionError('Convention +/- incorrecte')
else:
    print('✓ Convention validée automatiquement :')
    print(f'  Protéines barbées {_BARBED_REF} → toutes côté (+)')
    print(f'  Protéines pointées {_POINTED_REF} → toutes côté (-)')


# In[9]:


# Vérification visuelle sur un PDB exemple
pdb_ex = '7tpt'
ex = out[out['pdb_id'] == pdb_ex][['chain', 'n_actin_interactions', 'occupied_superclusters', 'label']]
print(f'Vérification sur {pdb_ex} :')
print(ex.to_string(index=False))

