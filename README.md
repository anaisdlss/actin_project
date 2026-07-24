# Actin–ABP interaction analysis

Structural analysis of actin and its actin-binding proteins (ABPs), built from
every 3D co-structure in [PPI3D](https://bioinformatics.lt/ppi3d) (actin,
UniProt **P60709**). The pipeline keeps assemblies with ≥ 5 connected actin
subunits, clusters interface residues, computes buried surface at each
interface, and shows everything in an interactive Streamlit app.

> **Systems:** macOS (Intel / Apple Silicon) and Linux. Windows is not
> supported (the MAFFT dependency is unavailable there).

**[→ USER GUIDE (GUIDE.md)](GUIDE.md)** — what the app shows, where the data
comes from, and what every number / colour means (also shown inside the app,
top section *Documentation*).

All code lives in **`script/`**; `data/` is regenerated locally (see step 3).

## 1. Install

Requires **pixi** ([pixi.sh](https://pixi.sh)). Install it first (macOS / Linux):

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

Then restart your terminal and set up the project:

```bash
git clone https://github.com/anaisdlss/actin_project.git
cd actin_project
pixi install        # recreates the environment (a few minutes)
```

## 2. Launch the app

```bash
pixi run streamlit run script/streamlit.py
```

It opens in your browser (otherwise: `http://localhost:8501`).

## 3. Generate the data

The repo ships **code only** — `data/` is git-ignored and regenerated locally.

In the app, open **Data download** and click **Run / update**. It runs the whole
pipeline (9 steps, ~1 h on a fresh clone). It is **resumable**: already-computed
steps are skipped, so you can quit and come back — it continues where it stopped.

## ProteoCast (optional)

Computing the per-ABP mutational landscape (**ProteoCast**) is separate from
`Run / update` and takes **several hours** (one job per ABP on
[proteocast.ijm.fr](https://proteocast.ijm.fr)). In the **ABP** section →
*ABP ProteoCast*, click **Compute all missing ProteoCast**. It is resumable.

## PyMOL (optional)

Only needed for the 3D scripts the pipeline writes under
`data/filtered/details/structures_files/bfactor_c70_interface/`. Install from
[pymol.org](https://pymol.org), then `File > Run Script…` (or `@/path/to.pml`).

## Shareable version (Streamlit Cloud)

A slim, read-only build (data pre-bundled, pipeline disabled) can be deployed to
Streamlit Community Cloud: `python script/make_slim_deploy.py` creates a
self-contained `deploy/` folder to push to a separate repo.
