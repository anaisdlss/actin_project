# User guide — Actin–ABP interaction analysis

This guide explains **what the app shows, where the data comes from, and what
each number means**. Short captions live in the app itself; this document is the
reference for the methods and the vocabulary.

---

## 1. What this app is

A structural analysis of **actin** and its **actin-binding proteins (ABPs)**,
built from every 3D co-structure of actin available in
[PPI3D](https://bioinformatics.lt/ppi3d) (actin = UniProt **P60709**).

For each 3D structure the app looks at **who touches actin, where, and how** —
which residues form the interface, how buried they are, how the binding sites
group together, whether partners compete or cooperate, and how conserved the
contacted positions are.

---

## 2. Where the data comes from (pipeline)

The full dataset is regenerated locally by a **9-step pipeline** (button
*Run / update*). Summary of the sources:

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
| **C70 cluster** (`cluster_data_70`) | Structural clustering of whole **interactions** at a 70 % threshold — interactions in the same C70 cluster have the **same 3D interface geometry** (labels like `0_7797_0`). |
| **Actin residue number (UniProt P60709)** | Actin residues are numbered as in human beta-actin (UniProt P60709, 1–375), whatever the PDB's own numbering. Internally, the app aligns every actin chain (MAFFT) and converts the alignment column to the P60709 residue; the few alignment columns absent from P60709 (N-terminal insertion of alpha-actins, one internal insertion) are shown as "ins.". Partner (ABP) positions are still alignment columns of their own sequence cluster. |
| **Conservation** | Evolutionary conservation of an actin position (ProteoCast/GEMME): higher = more conserved = less tolerant to mutation. |
| **Competition** | Two ABPs **compete** if their footprints on actin overlap (their C70 clusters cover the same region beyond a chosen %). |
| **Cooperation** | Two ABPs **cooperate** if they are **co-present in the same PDB** (they coexist on actin at the same time). |
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
| **Actin binding sites conservation** | **Overview**, **By ABP footprint**, **Solvent accessibility**. |
| **Human actin variants** | Gene/category selection, substitution maps, footprint summaries and exploratory associations. |
| **ABP conservation** | **Profiles**, **3D structure**, **Sequence alignments**. Missing results are identified explicitly. |
| **Physico-chemical properties of the interface** | **Actin surface chemistry**, **By binding site** (contact chemistry, alignments and structural comparisons). |
| **Homolog search** | Existing FoldDisco candidates, query motifs and annotations for the selected ABP. |

**Explore selected cluster** opens the matching comparative view. Clicking a
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
- **Conservation plot (ProteoCast on actin)** — grey line = conservation along
  the whole actin sequence; red dots = the positions the ABP contacts. Tells you
  whether a partner binds **conserved** (functionally important) or **variable**
  regions of actin.
- **Footprint vs surface** (ProteoCast panel) — *higher* / *lower*: is the
  actin footprint of this ABP more or less conserved than the rest of the actin
  surface? A Mann-Whitney p-value quantifies it.
- **Residue conservation … vs mean surface** — for one actin position: its
  conservation and how far it sits above/below the average surface residue.
- **ProteoCast mutational landscape** — per position of the ABP, how sensitive
  it is to mutation (dark = deleterious/constrained). The green track marks the
  positions that touch actin, so you see if the binding interface is under
  constraint.
- **Competition / Cooperation networks** — an edge means two ABPs compete
  (overlapping footprints) or cooperate (co-present in a PDB).
- **FoldDisco "same motif" reading** — whether two ABPs (or an ABP vs the PDB)
  share the same 3D interface geometry: coverage (% shared), normalised score
  (quality 0-1), RMSD (fit). A low RMSD on few residues is *inconclusive*.

---

## 6. Notes & caveats

- **ProteoCast results may be missing** because of service failures, input
  limitations or incomplete jobs. Consult the recorded diagnostics before
  retrying; the cause cannot be inferred from a missing score file alone.
- Numbers are **structure-derived** (from deposited PDB co-structures): they
  describe the interfaces that have actually been solved, not every possible
  interaction.
- PDB **4b1z** is deliberately excluded from all analyses.
