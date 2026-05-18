# Actin–ABP Interaction Analysis Pipeline

Automated structural analysis of actin–actin binding protein (ABP) interactions extracted from the [PPI3D](https://bioinformatics.lt/ppi3d) database (UniProt entry **P60709**, Actin cytoplasmic 1).

The pipeline retrieves all available 3D co-structures, filters assemblies containing at least four connected actin subunits, clusters interface residues at 70 % similarity (C70), computes buried surface area (B-factor encoding), and produces interactive visualizations via a Streamlit web interface and PyMOL session scripts.

---

## Requirements

| Tool | Purpose |
|------|---------|
| [pixi](https://pixi.sh) | Environment and dependency management (installed below) |
| [PyMOL](https://pymol.org) | 3D structure visualization |
| Internet access | Data retrieval from PPI3D and RCSB PDB |

Python and all Python dependencies are managed automatically by pixi — no manual installation required.

---

## Installation

### 1 — Install pixi

**macOS / Linux**

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

**Windows** (PowerShell)

```powershell
iwr -useb https://pixi.sh/install.ps1 | iex
```

Close and reopen the terminal before proceeding.

### 2 — Clone the repository

```bash
git clone https://github.com/anaisdlss/actin_project.git
cd actin_project
```

> Windows: use Git Bash or PowerShell with [Git for Windows](https://git-scm.com/download/win).

### 3 — Install the environment

```bash
pixi install
```

This downloads Python and all required libraries. Allow a few minutes on first run.

---

## Running the Pipeline

**macOS**

```bash
caffeinate -i pixi run python -m script.data_extract.pipeline_data
```

**Linux / Windows**

```bash
pixi run python -m script.data_extract.pipeline_data
```

> `caffeinate` (macOS only) prevents the system from sleeping during execution.  
> Estimated duration: **30–60 minutes** depending on network speed.

The pipeline is incremental: steps whose outputs are already up to date are skipped automatically on subsequent runs.

| Step | Description |
|------|-------------|
| 1/13 | Download PPI3D interaction summary |
| 2/13 | Download PDB entry metadata |
| 3/13 | Download full cluster table |
| 4/13 | Filter structures (≥ 4 connected actin subunits) |
| 5/13 | Download per-residue interface details |
| 6/13 | MAFFT sequence alignment by cluster |
| 7/13 | C70 interaction cluster analysis |
| 8/13 | Compute C70 interface B-factors |
| 9/13 | Generate global PyMOL script (all clusters) |
| 10/13 | Generate per-binding-site PyMOL scripts |
| 11/13 | Detailed C70 interface analysis |
| 12/13 | S1 binding site heatmap |
| 13/13 | Compute S1 cluster B-factors |

---

## Web Interface

```bash
pixi run streamlit run script/streamlit.py
```

Opens at `http://localhost:8501`. The interface provides interactive bipartite networks, cluster statistics, sequence alignments, and 3D B-factor viewers.

---

## PyMOL Visualizations

Scripts are generated under:

```
data/filtered/details/structures_files/bfactor_c70_interface/
```

**Global view** — all C70 cluster representatives superimposed on the reference actin (8iah chain L):

```
view_full_surface.pml
```

**Per-binding-site view** — one script per S1 actin binding site cluster (red nodes in the Streamlit network):

```
by_s1_cluster/<cluster_name>.pml
```

To execute a script, use `File > Run Script…` in PyMOL, or type in the PyMOL console:

```
@/path/to/script.pml
```

**Color scheme:** green gradient = ABP partner (heterologous interaction) · pink gradient = actin partner (homologous interaction) · color intensity encodes the percentage of buried accessible surface area at the interface.
