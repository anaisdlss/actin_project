# Actin–ABP interaction analysis pipeline

Automated structural analysis of the interactions between actin and its
actin-binding proteins (ABPs), extracted from the
[PPI3D](https://bioinformatics.lt/ppi3d) database
(UniProt entry **P60709**, Actin cytoplasmic 1).

The pipeline retrieves every available 3D co-structure, keeps the assemblies
that contain at least four connected actin subunits, clusters interface
residues at 70 % similarity (C70), computes the buried surface area at each
interface (encoded as a B-factor), and produces interactive visualisations
through a Streamlit web interface and PyMOL session scripts.

---

## Installation

> **Supported systems:** macOS (Intel and Apple Silicon) and Linux.
> Windows is not supported because one dependency (MAFFT) is unavailable on
> that platform.
> Make sure **pixi** is installed on your machine (see
> [pixi.sh](https://pixi.sh)).

### 0 — Open a terminal

Every command below is run in a terminal:

- **macOS:** `Cmd + Space` → type `Terminal` → Enter
- **Linux:** `Ctrl + Alt + T`

### 1 — Clone the repository

```bash
cd ~/Desktop
git clone https://github.com/anaisdlss/actin_project.git
cd actin_project
```

### 2 — Install the environment

```bash
pixi install
```

The first run can take a few minutes.

### 3 — Launch the web interface

```bash
pixi run streamlit run script/streamlit.py
```

A page opens automatically in your browser (otherwise, go to
`http://localhost:8501`).

---

## Generating the data

The repository ships **code only** — the `data/` folder is not included
(it is git-ignored) because everything can be regenerated locally.

Once the interface is open, go to the **Data download** section and click
**Run / update**. The pipeline runs entirely from the interface — no extra
command is needed. It takes **30 to 60 minutes** depending on your connection.
Steps that are already done are skipped automatically on later runs.

---

## PyMOL visualisations

The scripts are generated automatically by the pipeline in:

```
data/filtered/details/structures_files/bfactor_c70_interface/
```

**Global view** — every C70 cluster representative superimposed on the
reference actin (8iah chain L):

```
view_full_surface.pml
```

**Per binding-site view** — one script per S1 binding-site cluster (red nodes
in the Streamlit network):

```
by_s1_cluster/<cluster_name>.pml
```

To run a script, use `File > Run Script…` in PyMOL, or type in the PyMOL
console:

```
@/path/to/the/script.pml
```

**Colour code:** green gradient = ABP partner (heterologous interaction) ·
pink gradient = actin partner (homologous interaction) · intensity encodes
the percentage of buried accessible surface at the interface.

---

## Installing PyMOL

PyMOL is only required for the 3D visualisations. Download it from
[pymol.org](https://pymol.org) and follow the installation instructions for
your operating system.
