# Onglet « Comparaison structurale » — explication pour le tuteur

## 1. La question scientifique

Sur l'actine, plusieurs protéines partenaires (**ABP** = Actin-Binding Proteins)
se fixent au **même site de surface**. Question :

> Les ABP qui touchent le **même site** de l'actine partagent-ils le **même
> repliement / domaine** (elles sont apparentées), ou est-ce une **convergence**
> de familles **différentes** vers la même zone ?

## 2. D'où viennent les données

- **PPI3D** (bioinformatics.lt) : base d'interfaces protéine–protéine **dérivée du
  PDB** (structures 3D expérimentales). On en extrait toutes les interactions
  impliquant l'actine → `all_data.csv`, puis filtrées → `filtered_all_data.csv`.
- **Fichiers de structure 3D** : les assemblages biologiques (`.cif`/`.pdb`)
  téléchargés depuis PPI3D pour chaque complexe (`structures_files/assembly/`).
- **Résidus d'interface** : PPI3D fournit, par interaction, les résidus en contact
  et l'aire enfouie (`3.interface_residues.csv`, `4.inter-residue_contacts.csv`).

## 3. Le pipeline informatique (étape 9 « Analyses structurales ABP »)

Pour chaque **site de liaison** (cluster de surface d'actine), on prend les ABP qui
le touchent et on les compare structuralement (`script/abp_site_domain/`) :

| Étape | Outil | Ce que ça produit |
|---|---|---|
| Isoler la chaîne ABP de chaque complexe | **Biopython** | 1 fichier `.pdb` par ABP |
| Repliement (protéine entière) | **Foldseek** (easy-search, TM-align) | TM-score + %identité toutes paires |
| Repliement de la **zone d'interface** | **tmtools** (TM-align en Python) | TM-score de la seule région de contact |
| Domaines protéiques | **API InterPro** (via UniProt/SIFTS) | domaines Pfam de chaque ABP |
| Structure secondaire du contact | **PyMOL** (`dss`) | hélice / brin / boucle |
| Empreinte sur l'actine | résidus canoniques **MAFFT** + **Jaccard** | résidus d'actine touchés + recouvrement |

## 4. Signification des variables

### Table 1 — « Comparaison structurale des ABP de ce site » (une ligne par ABP)

- **ABP** : nom de la protéine partenaire de l'actine.
- **Famille** : regroupement automatique des ABP. Deux ABP sont dans la même
  famille s'ils partagent **≥ 1 domaine Pfam** OU le **même UniProt** (composantes
  connexes d'un graphe). Ex. `F04 · Vinculin family`.
- **Interface (struct. 2ndaire)** : structure secondaire **dominante** des résidus
  qui touchent l'actine, avec les proportions `H` (hélice) / `β` (brin) / `L`
  (boucle) — calculée par PyMOL `dss`. Dit **comment** l'ABP agrippe l'actine.
- **Domaines Pfam** : domaines protéiques annotés (InterPro/Pfam).
- **UniProt** : accession UniProt de la protéine.

### Table 2 — « Comparaison 2 à 2 des ABP du site » (une ligne par paire)

- **ABP A / ABP B** : la paire comparée.
- **Même famille** : partagent-ils un domaine Pfam / le même UniProt ? (Oui / Non).
- **TM interface (±6)** : **TM-score** (TM-align) de la **zone d'interface** élargie
  de 6 résidus de contexte de chaque côté. Mesure la similarité de repliement de la
  **seule région qui touche l'actine**.
- **TM structure entière** : **TM-score** sur la **protéine complète**.
  → Rappel : **TM-score ≥ 0,5 = même repliement** ; proche de 0 = repliements sans
  rapport (le TM-score va de 0 à 1, indépendant de la taille).
- **%id séquence** : **identité de séquence** sur la région alignée (Foldseek).
  `—` = trop peu de similarité pour aligner.
- **Empreinte actine (% communs)** : fraction de **résidus d'actine communs** aux
  deux empreintes (rapportée à la **plus petite** des deux empreintes). Répond à :
  « touchent-elles **les mêmes résidus** d'actine, ou seulement une région voisine ? »

## 5. Le résultat clé (l'histoire à raconter)

- **Même famille** (ex. Cofilin-1 vs Cofilin-2) : TM entière ≈ 0,97, %id ≈ 79 %,
  empreinte ≈ 94 % → ce sont vraiment les mêmes protéines/repliement (attendu).
- **Familles DIFFÉRENTES** au même site (ex. Vinculin vs Cofilin vs Formin) :
  TM faible (0,2–0,3 → **repliements sans rapport**), %id faible/`—`, **MAIS**
  empreinte actine élevée (70–94 % de résidus communs).

**Conclusion** : des protéines **structuralement non apparentées** touchent **les
mêmes résidus** de l'actine → la convergence est **positionnelle**, pas structurale.
C'est la **surface exposée de l'actine** (accessibilité au solvant) qui dicte le
site de fixation, pas une parenté de repliement des ABP.

## 6. Détails techniques (entrées / TM-align / InterPro)

### Entrée
Fichiers de **structure 3D** de PPI3D : l'**assemblage biologique** (`.cif`/`.pdb`)
de chaque complexe (`structures_files/assembly/`), avec les coordonnées atomiques.

### Quelle partie de la protéine on compare (Biopython découpe le PDB)
- **Protéine entière** (`02_extract_chains`) : on isole la **chaîne ABP** (bonne
  lettre de chaîne, uniquement les acides aminés — pas eau/ligands) → `.pdb`
  mono-chaîne.
- **Zone d'interface** (`05_extract_interface`) : on garde **seulement les résidus
  ABP qui touchent l'actine**, listés par PPI3D dans `3.interface_residues.csv`
  (`residue_number_structure`), éventuellement élargis d'un **contexte PAD** de
  ±0/6/12/20/35 résidus (`PAD=0` = contacts stricts ; « ±6 » = 6 résidus de contexte).
  → la partie à comparer est définie par **les contacts fournis par PPI3D**.

### Mécanique du TM-align (`08_interface_tm_sweep`, tmtools)
1. `load_ca` lit les **coordonnées 3D des Cα** + la **séquence** de chaque `.pdb`.
2. `tm_align(coords1, coords2, seq1, seq2)` **superpose** rigidement les deux
   structures pour **maximiser le TM-score** (0→1, indépendant de la taille).
3. On garde `max(tm_norm_chain1, tm_norm_chain2)`.
→ Input du TM-align = **coordonnées Cα + séquence** ; comparaison purement
géométrique de repliement.

### PDB → UniProt → InterPro (`03_interpro`) — ce n'est PAS la même chose
```
PDB + chaîne  ──SIFTS (API PDBe)──►  accession UniProt  ──API InterPro──►  domaines Pfam
```
- **SIFTS** (EBI/PDBe) : mappe la chaîne PDB vers son **accession UniProt** (identité
  de la protéine) → colonne **UniProt**.
- **InterPro** (EBI) : à partir de cet UniProt, renvoie les **domaines** Pfam/InterPro
  → colonne **Domaines Pfam**.
- **UniProt = l'identité** (la clé) ; **InterPro/Pfam = les domaines** annotés pour
  cette identité. La colonne **Famille** relie ensuite les ABP partageant ≥1 domaine
  Pfam ou le même UniProt.
