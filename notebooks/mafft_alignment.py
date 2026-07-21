#!/usr/bin/env python
# coding: utf-8

import os as _os
_os.chdir(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # cwd = racine projet (robuste, peu importe d'où on lance)


# # Alignements multiples MAFFT par groupe de séquences
# 
# **Pipeline** :
# 1. Chargement de la table d'interactions exportée depuis `cluster_interaction_analysis.ipynb`
# 2. Création des fichiers FASTA :
#    - S1 (actine) : 1–2 groupes selon `s1_sequence_cluster_70`
#    - S2 (partenaires) : par patch × sous-groupe `s2_sequence_cluster_70`
# 3. Alignements MAFFT (`--auto`)
# 4. Analyse par position : composition en acides aminés + fréquence d'interface
# 5. Export JSON par patch → input de `pymol_bfactor.ipynb`
# 
# **Prérequis** : avoir exécuté `cluster_interaction_analysis.ipynb` (cellules d'export en fin de notebook)

# In[1]:


import pandas as pd
import numpy as np
import subprocess
import json
import re
from pathlib import Path
from collections import Counter

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ── Chemins ───────────────────────────────────────────────────────────────
DATA_DIR  = Path('data')
EXPORT_DIR = DATA_DIR / 'exports'
FASTA_DIR  = EXPORT_DIR / 'fasta'
ALN_DIR    = DATA_DIR / 'alignments'

for d in [FASTA_DIR, ALN_DIR]:
    d.mkdir(parents=True, exist_ok=True)

print('Chemins OK')


# ## 1. Chargement des données

# In[2]:


DATA_DIR     = Path('data')
DETAILS_PATH = DATA_DIR / 'filtered' / 'details'

df_all     = pd.read_csv(DATA_DIR / 'filtered/filtered_all_data.csv')
df_summary = pd.read_csv(DATA_DIR / 'filtered/filtered_summary.csv')  # contient déjà interaction_id

df_merged = df_all.merge(
    df_summary,
    left_on=['subunit_1', 'subunit_2'],
    right_on=['Protein ID', 'Interactor ID'],
    how='left'
)
print(f"df_all    : {len(df_all)}")
print(f"df_merged : {len(df_merged)}")
print(f"Groupes s1_sequence_cluster_70 : {df_merged['s1_sequence_cluster_70'].nunique()}")


# ## 2. Création des fichiers FASTA
# 
# Un fichier par groupe `s1_sequence_cluster_70` contenant toutes les séquences uniques (S1 + S2).

# In[3]:


from collections import defaultdict

# Mapping subunit → titre protéine
subunit_title = {}
for _, row in df_merged.iterrows():
    for subunit_col, title_col in [
        ('subunit_1', 'subunit_1_title'),
        ('subunit_2', 'subunit_2_title'),
    ]:
        subunit = row[subunit_col]
        title   = str(row.get(title_col, '')).strip()
        if title and title != 'nan':
            subunit_title[subunit] = title

# Aplatir : une entrée par séquence (S1 ou S2) avec son cluster
seq_entries = []  # (cluster_id, sequence, subunit)
for _, row in df_merged.iterrows():
    for seq_col, subunit_col, cluster_col in [
        ('s1_sequence', 'subunit_1', 's1_sequence_cluster_70'),
        ('s2_sequence', 'subunit_2', 's2_sequence_cluster_70'),
    ]:
        seq     = str(row.get(seq_col, '')).strip()
        cluster = row.get(cluster_col)
        subunit = row[subunit_col]
        if seq and seq != 'nan' and pd.notna(cluster):
            seq_entries.append((int(cluster), seq, subunit))

# Grouper par cluster → séquences uniques
cluster_data = defaultdict(lambda: defaultdict(list))  # cluster_id → seq → [subunits]
for cluster_id, seq, subunit in seq_entries:
    if subunit not in cluster_data[cluster_id][seq]:
        cluster_data[cluster_id][seq].append(subunit)

# Écrire un FASTA par cluster
n_files = 0
for cluster_id, seq_to_subunits in sorted(cluster_data.items()):
    records = [
        SeqRecord(
            Seq(seq),
            id='|'.join(subunits),
            description=' | '.join(dict.fromkeys(
                subunit_title.get(s, '') for s in subunits
                if subunit_title.get(s, '')
            )),
        )
        for seq, subunits in seq_to_subunits.items()
    ]
    if len(records) >= 2:
        path = FASTA_DIR / f"cluster_{cluster_id}.fasta"
        SeqIO.write(records, path, 'fasta')
        n_files += 1
        print(f"  Cluster {cluster_id} : {len(records)} séquences uniques")
    else:
        print(f"  Cluster {cluster_id} : {len(records)} séquence (ignoré — MAFFT nécessite ≥ 2)")

print(f'\nFASTA : {n_files} groupe(s) → {FASTA_DIR}/')


# ## 3. Alignements MAFFT
# 
# MAFFT est appelé avec `--auto` (choisit automatiquement l'algorithme selon la taille).

# In[4]:


import sys
import tempfile
import shutil

# ── Binaires MAFFT (bypass du wrapper bash qui se bloque sur ce système) ──
_MAFFT_PREFIX = Path(
    '/Users/user/Desktop/stage/actin_project/.pixi/envs/default/libexec/mafft'
)

# Vérification que les binaires sont disponibles
for _bin in ('disttbfast', 'f2cl', 'version'):
    if not (_MAFFT_PREFIX / _bin).exists():
        raise EnvironmentError(
            f'Binaire MAFFT manquant : {_MAFFT_PREFIX / _bin}\n'
            'Vérifier que l\'environnement pixi est activé : `pixi shell`\n'
            'Sinon : `pixi install` pour installer les dépendances.'
        )

_version = subprocess.check_output(
    [str(_MAFFT_PREFIX / 'version')], stderr=subprocess.STDOUT
).decode().strip()
print(f'MAFFT binaries : v{_version}')


def run_mafft(fasta_path: Path, aln_path: Path) -> bool:
    """Aligne les séquences du fichier FASTA et écrit le résultat dans aln_path.

    Utilise les binaires MAFFT directement (disttbfast + f2cl) pour contourner
    le wrapper bash qui se bloque sur ce système.
    Algorithme : FFT-NS-2 (progression 2 cycles, rapide et robuste).
    """
    tmpdir = Path(tempfile.mkdtemp())
    try:
        # Préparer infile (même traitement que le wrapper bash)
        infile = tmpdir / 'infile'
        infile.write_text(fasta_path.read_text().replace('\r', '\n') + '\n')

        # Étape 1 : disttbfast — calcul des distances + alignement progressif
        # -E 2  : 2 cycles (FFT-NS-2)
        # -C 1  : 1 thread (évite -C 0-0 qui cause le blocage)
        # -F    : FFT activé
        # -J d  : algorithme NJ pour l'arbre guide
        r1 = subprocess.run(
            [str(_MAFFT_PREFIX / 'disttbfast'),
             '-q', '0', '-E', '2', '-V', '-1.53',
             '-s', '0.0', '-W', '6', '-C', '1', '-F', '-J', 'd'],
            stdin=open(infile), stdout=open(tmpdir / 'pre', 'w'),
            stderr=subprocess.DEVNULL, cwd=tmpdir, timeout=120
        )
        if r1.returncode != 0 or (tmpdir / 'pre').stat().st_size == 0:
            return False

        # Étape 2 : f2cl — convertit le format interne → FASTA
        r2 = subprocess.run(
            [str(_MAFFT_PREFIX / 'f2cl'), '-n', '-1', '-f', '-l', '60'],
            stdin=open(tmpdir / 'pre'), stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL, cwd=tmpdir, timeout=30
        )
        if r2.returncode != 0 or not r2.stdout.strip():
            return False

        aln_path.write_bytes(r2.stdout)
        return True

    except subprocess.TimeoutExpired:
        return False
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


fastas = sorted(FASTA_DIR.glob('*.fasta'))
ok = 0
for fasta in fastas:
    aln = ALN_DIR / fasta.name.replace('.fasta', '.aln')
    print(f'  → {fasta.stem} ...', end=' ', flush=True)
    if run_mafft(fasta, aln):
        ok += 1
        print('✓')
    else:
        print('✗')

print(f'\nAlignements : {ok}/{len(fastas)} → {ALN_DIR}/')


# ## 4. Mapping positions MAFFT → tables 3 et 4
# 
# Pour chaque séquence alignée, on construit le mapping :
# `position dans la séquence originale → colonne dans l'alignement`
# 
# Ce mapping est appliqué aux tables 3 (résidus d'interface) et 4 (contacts inter-résidus)
# pour ajouter une colonne `canon_mafft` : position canonique commune à toutes les interactions du cluster.

# In[6]:


def build_seq_to_col_mapping(aligned_seq: str) -> dict:
    """
    Pour une séquence alignée (avec gaps '-'), retourne le mapping
    position_séquence (1-based) → colonne_alignement (1-based).
    """
    mapping = {}
    seq_pos = 0
    for col, aa in enumerate(aligned_seq, start=1):
        if aa != '-':
            seq_pos += 1
            mapping[seq_pos] = col
    return mapping


# ── Mappings depuis les alignements MAFFT ────────────────────────────────
subunit_to_mapping = {}
for aln_path in sorted(ALN_DIR.glob('*.aln')):
    try:
        alignment = AlignIO.read(str(aln_path), 'fasta')
    except Exception:
        continue
    for record in alignment:
        mapping = build_seq_to_col_mapping(str(record.seq))
        for subunit in record.id.split('|'):
            subunit_to_mapping[subunit] = mapping

# ── Mapping identité pour les clusters à 1 séquence (pas d'alignement) ──
# canon_mafft = position dans la séquence (pas de gaps)
seq_by_subunit = {}
for _, row in df_merged.iterrows():
    seq_by_subunit[row['subunit_1']] = str(row['s1_sequence'])
    seq_by_subunit[row['subunit_2']] = str(row['s2_sequence'])

n_identity = 0
for subunit, seq in seq_by_subunit.items():
    if subunit not in subunit_to_mapping and seq and seq.upper() != 'NAN':
        subunit_to_mapping[subunit] = {i: i for i in range(1, len(seq) + 1)}
        n_identity += 1

print(f'Mappings depuis alignements : {len(subunit_to_mapping) - n_identity} subunits')
print(f'Mappings identité (1 séquence) : {n_identity} subunits')

# ── Charger tables 3 et 4 ────────────────────────────────────────────────
df_res3 = pd.read_csv(DATA_DIR / 'filtered/details/3.interface_residues.csv')
df_res4 = pd.read_csv(DATA_DIR / 'filtered/details/4.inter-residue_contacts.csv')
df_int1 = pd.read_csv(DATA_DIR / 'filtered/details/1.interactions.csv')

# chain_A (= S1) par interaction
chain_a_map = dict(zip(df_int1['interaction_id'], df_int1['chain_A_id']))

# ── Table 4 : residue_A_sequence → canon_mafft ───────────────────────────
def map_t4(row):
    subunit = chain_a_map.get(int(row['interaction_id']))
    mapping = subunit_to_mapping.get(subunit, {})
    seq_pos = row['residue_A_sequence']
    if pd.isna(seq_pos):
        return pd.NA
    return mapping.get(int(seq_pos), pd.NA)

df_res4['residue_A_canon_mafft'] = df_res4.apply(map_t4, axis=1)
n4 = df_res4['residue_A_canon_mafft'].notna().sum()
print(f'Table 4 : {n4}/{len(df_res4)} contacts avec residue_A_canon_mafft')

# ── Table 3 : residue_number_sequence → canon_mafft (chain S1 uniquement) ─
def map_t3(row):
    iid = int(row['interaction_id'])
    subunit = chain_a_map.get(iid)
    if row['chain'] != subunit:
        return pd.NA  # chaîne S2 → pas de mapping S1
    mapping = subunit_to_mapping.get(subunit, {})
    seq_pos = row['residue_number_sequence']
    if pd.isna(seq_pos):
        return pd.NA
    return mapping.get(int(seq_pos), pd.NA)

df_res3['residue_number_canon_mafft'] = df_res3.apply(map_t3, axis=1)
n3 = df_res3['residue_number_canon_mafft'].notna().sum()
print(f'Table 3 : {n3}/{len(df_res3)} résidus S1 avec residue_number_canon_mafft')

# ── Sauvegarder ───────────────────────────────────────────────────────────
df_res4.to_csv(DATA_DIR / 'filtered/details/4.inter-residue_contacts.csv', index=False)
df_res3.to_csv(DATA_DIR / 'filtered/details/3.interface_residues.csv', index=False)
print('Tables 3 et 4 mises à jour avec canon_mafft')


# ## 4. Analyse des alignements
# 
# Pour chaque colonne de l'alignement :
# - **Composition** : fréquence de chaque acide aminé (les « différentes versions » du nœud)
# - **Fréquence d'interface** : proportion des interactions du patch où cette position canonique est en contact

# In[7]:


def alignment_composition(aln_path: Path) -> dict:
    """
    Pour chaque colonne de l'alignement, retourne la composition en acides aminés.
    Returns: {col_index: {'variants': {aa: freq}, 'n_seqs': int, 'gap_fraction': float}}
    """
    try:
        alignment = AlignIO.read(str(aln_path), 'fasta')
    except Exception:
        return {}

    n_seqs = len(alignment)
    composition = {}
    for col in range(alignment.get_alignment_length()):
        residues = [str(alignment[i].seq[col]).upper() for i in range(n_seqs)]
        non_gap = [r for r in residues if r not in ('-', 'X', '*')]
        if not non_gap:
            continue
        counts = Counter(non_gap)
        total = len(non_gap)
        composition[col] = {
            'variants': {aa: round(c / total, 3) for aa, c in counts.most_common()},
            'n_seqs': n_seqs,
            'gap_fraction': round(1 - total / n_seqs, 3)
        }
    return composition


compositions = {}
for aln in sorted(ALN_DIR.glob('*.aln')):
    cluster_id = aln.stem.replace('cluster_', '')
    comp = alignment_composition(aln)
    compositions[cluster_id] = comp
    print(f'  Cluster {cluster_id} : {len(comp)} colonnes alignées')

print(f'\n{len(compositions)} alignement(s) analysé(s)')

