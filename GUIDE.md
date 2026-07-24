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

**ProteoCast** (per-ABP mutational landscape) is a separate, opt-in computation
run on [proteocast.ijm.fr](https://proteocast.ijm.fr) (one job per ABP).

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
| **Canonical position** | Residue numbering in the shared MSA reference — lets you compare the **same** position across different structures/species (independent of each PDB's own numbering). |
| **Conservation** | Evolutionary conservation of an actin position (ProteoCast/GEMME): higher = more conserved = less tolerant to mutation. |
| **Competition** | Two ABPs **compete** if their footprints on actin overlap (their C70 clusters cover the same region beyond a chosen %). |
| **Cooperation** | Two ABPs **cooperate** if they are **co-present in the same PDB** (they coexist on actin at the same time). |
| **Footprint** | The set of actin residues a given ABP contacts. |
| **Representative pair** | For a cluster, the single most frequent structure shown in the 3D viewer / sequences (so you look at one clear example, not an average). |

---

## 4. The app, section by section

### Data download *(local only)*
Runs / updates the 9-step pipeline. Resumable. Hidden on the shared version.

### Filtered data & Valid PDB structures (explorer)
Pick a PDB. You get:
- **Interaction network** — the chains of that structure and who touches whom
  (orange = actin, green = ABP). Click an edge (a pair) or a node (a chain).
- **3D visualisation — interface contacts** — the structure in surface;
  the selected chain is **yellow**, its partner **blue**, the rest grey.
- **Sequences — interface residues** — the two sequences with their interface
  residues highlighted (yellow = S1, blue = S2). Hover a highlighted residue to
  read its **canonical position**.
- **Interaction clusters per pair** — for the selected pair, the C70 cluster and
  the S1/S2 binding-site clusters it belongs to.

### Interaction clusters
Actin's binding sites, grouped. Two selectors:
- **Patch S1 binding site** → the *Interactive network — actin residues ↔
  partners* (each actin residue coloured by buried %ASA) + an **Interface 3D**
  of that binding site.
- **Patch cluster_data_70** (C70) → a **bipartite network** (actin residues ↔
  partner residues, coloured by buried %ASA) + the **3D interface of the
  representative pair** + the **two interface sequences coloured by buried
  %ASA** (pale → dark red; hover = position + %ASA).

### ABP
- **Overview** — table of all ABPs.
- **ABP networks** — **Competition** (footprints overlap) or **Cooperation**
  (co-present in a PDB). Node colour = ABP family, node size = proportion of
  binding sites in conflict.
- **Per-ABP detail** — for one ABP: its interaction clusters, actin–actin homo
  interactions on the same PDBs, Foldseek interface-motif discovery, and its
  **ProteoCast** panel.
- **ProteoCast panel** — the ABP's **mutational landscape** (per-position
  substitution sensitivity / conservation) and the **actin-side conservation**
  of the residues it contacts (is the footprint more or less conserved than the
  rest of the actin surface?), plus a 3D structure coloured by mutational
  sensitivity.

### MSA — Interface proteins
Per family (myosins, tropomyosins…), the aligned sequences with interface
columns highlighted, actin canonical positions, and a *vs C. elegans* identity
comparison.

### Interactive explorer
Paste a query actin sequence (which positions vary → which ABPs are involved),
or compare an ABP pair (footprint overlap), on the per-residue passport table.

---

## 5. Reading the main visualisations

- **Bipartite / radial residue networks** — each residue is a node coloured by
  **buried %ASA** (pale yellow = low, dark red = high). Edge = a contact. Hover a
  node for its canonical position, %ASA and interaction count.
- **Interface sequences coloured by %ASA** — the linear view of the same
  information: each interface residue is shaded by how buried it is at the
  contact.
- **3D surfaces** — the actual structure; colours mark the selected chain / its
  partner / the contacted region.
- **Conservation plot (ProteoCast on actin)** — grey line = conservation along
  the whole actin sequence; red dots = the positions this ABP contacts. Lets you
  see if a partner binds conserved or variable regions.

---

## 6. Notes & caveats

- **A few ProteoCast entries can't be computed** by proteocast.ijm.fr (very large
  proteins, fusion constructs, or too weak an MSA) — this is expected, not an
  app bug.
- Numbers are **structure-derived** (from deposited PDB co-structures): they
  describe the interfaces that have actually been solved, not every possible
  interaction.
- PDB **4b1z** is deliberately excluded from all analyses.
