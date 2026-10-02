# User guide — Actin–ABP interaction analysis

This guide explains **what the app shows, where the data comes from, and what
each number means**. Short captions live in the app itself; this document is the
reference for the methods and the vocabulary.

---

## 1. What this app is

A structural analysis of **actin** and its **actin-binding proteins (ABPs)**,
built from the retained snapshot of actin co-structures retrieved from
[PPI3D](https://bioinformatics.lt/ppi3d) (actin = UniProt **P60709**).

For each 3D structure the app looks at **who touches actin, where, and how** —
which residues form the interface, how buried they are, how the binding sites
group together, where partner footprints overlap, which partners occur in the same structures, and the model-derived mutational sensitivity of contacted positions.

---

## 2. Where the data comes from (pipeline)

The local project separates two operations under **Documentation → Data management**:

- **Rebuild local scientific results** computes documented RSA, the ABP motif catalogue, variant/sensitivity/interface summaries, interface geometry, filament proximity, S1 figures and local FoldDisco controls from the installed sources. The procedure checks file contents, code and software versions, skips unchanged results, records each successful calculation and resumes after failure. It does not submit external jobs.
- **Run / update** refreshes PPI3D through the existing nine-stage acquisition pipeline, then runs the local scientific rebuild. If the remote service cannot be checked, the update stops with an error and retains installed data. A network failure is not reported as a successful refresh.

**Data origins and calculation receipts** lists every CSV, its checksum and the provenance that could actually be established. A recovered producer script is distinct from a verified historical run. Imported variant snapshots have an identified import source but unknown original release dates. ProteoCast scores and database assertions remain external inputs. Undocumented historical files are listed explicitly; they are not certified by their presence.

Current residue RSA is recomputed from experimental **7PDZ chain I**, with three explicit contexts: isolated chain in the same conformation, six-actin fragment, and fragment plus capping proteins. It is **not a mean over all actins or all PDBs**. The isolated-chain context defines the RSA ≥ 0.2 baseline in mutational-sensitivity comparisons. Missing/modified/incomplete residues retain missing RSA. The old `conservation_vs_asa_per_position.csv` is retained for historical comparison, but current readers and the migrated analysis scripts no longer use it as a measurement source.

Source summary:

1. **PPI3D** — all 3D interactions involving actin (summary + structures +
   inter-residue contacts). We keep only assemblies with **≥ 5 connected actin
   subunits** (i.e. real filament/oligomer contexts, not isolated chains).
2. **Interface residues & contacts** — computed by PPI3D: for every pair of
   chains in contact, the residues at the interface, their **buried surface
   area**, and the residue-residue **contact area**.
3. **MAFFT** — multiple sequence alignments per sequence cluster, to map every
   structure's numbering onto a common **canonical numbering**.
4. **Clustering** — interactions are grouped structurally (see *C70 clusters*).
5. **B-factors** — interface metrics written into PDB files for the 3D viewers.
6–7. **Visualisations** — heatmaps and PyMOL scripts per binding site.
8–9. **ABP analyses** — competition/cooperation, Foldseek/InterPro domains,
   TM-align, interface footprints, physicochemistry.

> The **shared / Cloud** version ships this data pre-computed: the pipeline and
> the ProteoCast computation are disabled there (they only run on the full local
> project).

### How long it takes (local only)

- **Regenerating the full dataset** (*Run / update*) takes **≈ 1 hour** on a
  fresh clone — it downloads every actin co-structure from PPI3D and runs
  MAFFT + the structural analyses. It is **resumable** (steps already done are
  skipped), so you can stop and come back.
- **ProteoCast** (per-ABP mutational landscape) is a separate, opt-in
  computation run on [proteocast.ijm.fr](https://proteocast.ijm.fr), **one job
  per ABP**. Computing **all** ABPs can take **several hours** (some large
  proteins take ~20 min each). It is resumable too. Missing results and failed
  jobs have separate diagnostics; a failure alone does not show that a protein
  cannot be computed.

---

## 3. Key terms (glossary)

| Term | Meaning |
|---|---|
| **ABP** | Actin-binding protein — any non-actin partner in contact with actin. |
| **S1 / S2** | The two sides of an interaction. By convention **S1 = actin**, S2 = the partner (or the other actin, for actin–actin contacts). |
| **Interface residue** | A residue that loses solvent-accessible surface upon binding (i.e. it touches the partner). |
| **Buried %ASA** | For one residue: the fraction of its surface **buried** at the interface. **0 % = fully exposed, 100 % = fully buried.** The more buried, the more central to the contact. In heatmaps: pale = little buried, dark red = strongly buried. |
| **Buried contact area (Å²)** | The physical contact surface between **two** specific residues (one on each side). |
| **Binding site (S1) cluster** | Groups of actin surface patches used by partners (labels like `6685_2`). Two partners on the same binding-site cluster use the **same region of actin**. |
| **C70 cluster** (`cluster_data_70`) | Interface cluster from PPI3D (for example `0_7797_0`). The 70 refers to the upstream sequence-identity clustering level, not 70% identical contacts. Interfaces are then grouped by contact-area similarity. Membership does not establish identical geometry. |
| **Actin residue number (UniProt P60709)** | Actin residues are numbered as in human beta-actin (UniProt P60709, 1–375), whatever the PDB's own numbering. Internally, the app aligns every actin chain (MAFFT) and converts the alignment column to the P60709 residue; the few alignment columns absent from P60709 (N-terminal insertion of alpha-actins, one internal insertion) are shown as "ins.". Partner (ABP) positions are still alignment columns of their own sequence cluster. |
| **ProteoCast sensitivity** | Negative mean of the 20 supplied substitution scores at a position, including the unchanged amino acid. A model-derived measure of mutational constraint, not sequence identity or a clinical classification. |
| **Competition network** | Overlap of observed ABP footprints. This suggests possible competition but does not prove a steric clash or competitive binding. |
| **Cooperation network** | Co-presence of ABPs in a retained PDB entry. This describes structural co-occurrence, not demonstrated cooperative binding. |
| **Footprint** | The set of actin residues a given ABP contacts. |
| **Representative pair** | For a cluster, the single most frequent structure shown in the 3D viewer / sequences (so you look at one clear example, not an average). |

---

## 4. Navigation and scientific pages

The sidebar follows the eleven sections of the scientific brief. Selecting a
section opens its own page. Documentation displays the guide directly. Where
a section contains several analyses, the **View** control selects which one
is rendered. Main tables, selected networks and structural evidence are visible
without an extra disclosure control. Residue, protein, cluster
and explicitly keyed analysis choices are retained during navigation in the
same browser session. Reloading the browser starts a new Streamlit session.

| Page | Available views and purpose |
|---|---|
| **Documentation** | Guide displayed directly; **Data management** groups cache reload and, in the full project, data updates and ProteoCast job diagnostics. |
| **Summary tables** | **Structures** (retained PDB explorer), **Source tables**, **Residue numbering**, **Dataset checks**. |
| **Actin use at a residue level** | **Residue explorer** with contacts and a 3D surface; **Binding-site heatmap** across homo/hetero sites. |
| **Actin-actin interfaces** | **Binding sites**, **Compare footprints**, **Structural evidence**. Mixed sites are also included. |
| **ABP-actin interfaces** | **By protein** (clusters and structures), **Overview and heatmap**, **Binding sites**. |
| **Comparative analyses of binding sites** | **Binding-site clusters**, **Interaction clusters**, **ABP networks**, **ABP pairs and sequences**. |
| **Actin mutational sensitivity** | **Overview**, **By ABP footprint**, **Solvent accessibility**. |
| **Human actin variants** | Gene/category selection, substitution maps, footprint summaries and exploratory associations. |
| **ABP mutational sensitivity** | **Profiles**, **3D structure**, **Sequence alignments**. Missing results are identified explicitly. |
| **Physico-chemical properties of the interface** | **Actin surface chemistry**, **By binding site** (contact chemistry, alignments and structural comparisons). |
| **Homolog search** | Historical results and tracked searches for single or combined ABP motifs; motif coordinates and provenance are retained. New submissions run in the full project. |

Selecting a binding site in **Actin–actin interfaces** or **ABP–actin interfaces** displays its contact details, heatmap, network and 3D view on that same page. Clicking a
binding-site heatmap cell opens that cluster with the clicked residue selected
when that position has contact details in the cluster.
In the global binding-site network, a site node selects its details and an ABP
node opens the protein page. Global cluster networks and catalogues are grouped
under **All binding sites: network and table** or **All interaction clusters:
network and table**, above the selected-cluster results.

The PDB menu follows the current retained interaction dataset automatically.
Missing titles are recovered from local metadata when possible. Missing titles
or structure files alone do not remove an otherwise retained structure.
The structure explorer combines chain/interaction networks, interface surfaces,
sequences and cluster assignments. Orange identifies actin and blue identifies
ABPs in the assembly view; the selected-pair view has its own labelled colours.

Position plots use a vertical hover guide. Aligned numerical tracks share a
tooltip to read their values at the same residue. Heatmaps retain the hovered
cell's value. Downloads and saved ProteoCast figures are grouped below the
profile; calculation controls are kept in Documentation. Loading a page does
not start the pipeline or submit a new ProteoCast job.

---

## 5. What each readout indicates (quick reference)

Because the interface is kept clean, the meaning of every element is listed here.

**Counts and headers**
- *"N actin residues · M partner proteins · n=X interactions"* — how much data
  the current view is built on: X = number of structural interactions pooled,
  N/M = distinct residues/partners involved.
- *"2,152 rows · 23 columns"* — size of the underlying table.
- *"N interactions"* next to a cluster — how many solved interactions fall in it
  (bigger = more frequently observed geometry).

**Colours**
- **Buried %ASA scale** (pale yellow → dark red): how buried a residue is at the
  interface. Pale = barely touching, dark red = deeply buried = central to the
  contact. Used in the residue networks and the interface sequences.
- **3D surfaces**: yellow = the selected chain, blue = its partner, grey = the
  rest of the assembly. In binding-site 3D, the actin surface is shaded by the
  buried %ASA of the contacted patch.
- **ABP network node colour** = ABP family; **node size** = share of its binding
  sites that are contested (competition view).
- **Purple gradient** (specificity views) = how many partners contact a given
  actin position (bright = shared by many, pale = specific to one).

**Individual readouts**
- **Bipartite / radial residue networks** — each residue is a node coloured by
  buried %ASA; an edge is a residue-residue contact. Hover a node for its
  canonical position, %ASA and interaction count.
- **Interface sequences coloured by %ASA** — the linear version of the same
  information: each interface residue is shaded by how buried it is.
- **Mutational sensitivity profiles** — model-derived sensitivity along actin, with the selected ABP footprint and structural accessibility shown separately. Comparisons and p-values are exploratory; residues are not independent observations.
- **Residue sensitivity vs surface mean** — located in **Actin mutational sensitivity**. The selected residue is retained when opening **Mutational sensitivity of this residue** from the contact explorer. The legacy surface mean uses RSA ≥ 0.2; that source's monomer/filament provenance remains unresolved.
- **Interaction-type surface** — orange: positive buried-ASA observations with actin only; blue: ABPs only; purple: both; gray: no observed positive contact or unmapped atoms. Both sides of actin–actin interactions are included. Black sticks locate the selected residue. Absence of an observed contact does not establish absence of binding.
- **ABP landscape** — scores cover the whole supplied query, with observed contacts and annotated domains on aligned tracks. Missing scores remain missing.
- **Competition / Cooperation networks** — observed footprint overlap / co-presence, respectively. Neither establishes a functional mechanism by itself.
- **FoldDisco results** — coverage, raw score and RMSD describe a structural-motif match. Historical normalized scores use the stated source-labelled or best-hit denominator, can exceed one, and are not probabilities. A low RMSD on few residues is inconclusive.
- **Combined FoldDisco motifs** — choose additional sites observed on one ABP chain. Coordinates from different structures are never joined. Edit the PDB residue list explicitly when necessary. The [public server](https://github.com/soedinglab/MMseqs2-App/blob/master/frontend/FoldDiscoSearch.vue) limits motifs to 32 residues; larger motifs remain exportable for a local search. Searches never launch just by changing a selector.
- **Tracked searches** — new results retain the exact submitted PDB, residue list, database names, ticket and timestamps. Complete with no returned alignments, failed, pending and unsupported are distinct states. Saved historical results remain separate. Returned hits may be capped by the service and are not proof of homology or actin binding.
- **Interface geometry audit** — all retained homo pairs with sufficient local coordinates are compared to 3J8A (with tropomyosin) and 5YU8 (with cofilin). One subunit is fitted; the second is measured under that transform. Results are summarized once per PDB and remain descriptive until biological interpretation is reviewed.

---

## 6. Notes & caveats

- **ProteoCast results may be missing** because of service failures, input
  limitations or incomplete jobs. Consult the recorded diagnostics before
  retrying; the cause cannot be inferred from a missing score file alone.
- Numbers are **structure-derived** (from deposited PDB co-structures): they
  describe the interfaces that have actually been solved, not every possible
  interaction.
- PDB **4b1z** is deliberately excluded from all analyses.


### Display and interpretation update — 2 October 2026

- Interaction IDs in source tables sort numerically. Data identifiers and joins remain unchanged.
- ProteoCast residue summaries are labelled **mutational sensitivity**, defined as minus the mean of the 20 supplied substitution scores. They are neither sequence-identity percentages nor clinical severity grades. Original source filenames are retained for compatibility.
- Heatmap tooltips include the hovered cell colour. Binary maps show recorded/not-recorded states; quantitative ASA maps retain numerical percentages. Plot zoom controls stay visible.
- Three-dimensional panels fit their container when resized. Binding-site views can reconstruct observed actin–ABP pairs from the local assemblies and contact tables, even without pre-generated cluster PDBs. Partner superpositions require at least 100 matching actin Cα pairs, 80% matching resolved residues and RMSD ≤ 5 Å; they do not demonstrate co-binding.
- ABP 3D sensitivity requires complete ProteoCast score grids and exact agreement between query and structure numbering/identity. Otherwise the original AlphaFold model is explicitly labelled **pLDDT confidence**. Arbitrary PDB B-factors are never assumed to be sensitivity.
- Contact-region zoom stays visible, disabled with an explanation when contacts cannot be mapped. Its bounds are observed contact positions with a margin, not a predicted domain boundary.
- FoldDisco query positions refer to the selected ABP's PDB chain. The resolved sequence and 3D motif identify the selected residues. Distant sequence positions may be neighbours in 3D. The public server limit remains 32 residues; the app never truncates motifs automatically.

These display changes do not resolve the remaining scientific validations: biological interpretation of major/minor interfaces, provenance of the legacy RSA reference, unavailable ProteoCast results, and positive/negative control validation of FoldDisco candidates. Sequence conservation, model sensitivity, structural similarity and clinical annotation remain distinct concepts.


## Accessing FoldDisco results

In **Homolog search**, select an ABP and its actin site. **Saved results** shows a tracked result when available, with a downloadable CSV, links to PDB/AlphaFold records and the exact query provenance. Editing the preparation below does not overwrite results. **Historical results** contains the original exports; exact historical query provenance is incomplete. An empty inventory does not mean a negative search.

**Local validation: known ABPs and observed actin contacts** displays the completed searches against representative ABP chains already in this dataset. This includes motifs longer than the public server's 32-residue limit. It is a finite local control panel, not a search of the full PDB or AlphaFold databases. The contact columns test whether matched target residues fall on the target's observed actin interface at the same site or at any recorded site. Partial self-matches are flagged. These checks do not establish evolutionary homology or negative-control specificity.

**Prepare or edit a new motif** shows the resolved sequence, selected residues and 3D motif. The research version can submit a new motif of at most 32 residues; the public version displays saved results. Export both the motif text and its matching prepared PDB: a numeric source chain is renamed A in this export because the FoldDisco query parser only accepts letter chain prefixes. Coordinates and residue numbers are unchanged. Invalid old queries are retained for audit and are not negative results.

In **Actin-actin interfaces → Compare footprints**, the **3D screening** panel reports distances from rigidly placed ABPs to neighbouring actins in two finite reference filaments. Proximity counts are warning signals, not validated predictions of steric competition. Their cutoffs, fit quality, source structures and limitations accompany the table and downloads.
