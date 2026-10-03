# Actin–ABP interaction analysis

## Reproducible local calculations

In the full research project, open **Documentation → Data management → Rebuild local scientific results**, or run `pixi run rebuild`. This offline workflow rebuilds fifteen registered scientific calculations and verifies source/code/output checksums before reusing a result. `pixi run data-status` checks freshness; `pixi run data-inventory` inventories CSV provenance. The public app displays saved results and does not launch these local tools.

Successful calculation receipts, logs and the CSV inventory are saved in `reports/local_rebuild/`. An interrupted or failed step is never accepted merely because an old output exists. Remote PPI3D updates remain a separate action; imported ProteoCast and variant source snapshots are not represented as locally generated measurements. Historical source gaps remain visible in the inventory.

The RSA displayed by the app now comes from the reproducible 7PDZ chain I calculation, with explicit isolated/fragment/capped-fragment contexts. It does not average unrelated structures or silently reuse the legacy CSV.


Structural analysis of actin and its actin-binding proteins (ABPs), built from
every 3D co-structure in [PPI3D](https://bioinformatics.lt/ppi3d) (actin,
UniProt **P60709**). The pipeline keeps assemblies with ≥ 5 connected actin
subunits, clusters interface residues, computes buried surface at each
interface, and shows everything in an interactive Streamlit app.

> **Systems:** macOS (Intel / Apple Silicon) and Linux. Windows is not
> supported (the MAFFT dependency is unavailable there).

**[→ USER GUIDE (GUIDE.md)](GUIDE.md)** — what the app shows, where the data
comes from, and what every number / colour means (also shown inside the app,
*Documentation*).

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

The sidebar opens one scientific page at a time. Where several analyses are
available, use **View** to choose between them. Selections are retained while moving between pages
in the same session; calculations and dataset updates are under Documentation.

## 3. Generate the data

The repo ships **code only** — `data/` is git-ignored and regenerated locally.

In the app, open **Documentation → Data management** and click **Run / update**. It runs the whole
pipeline (9 steps, ~1 h on a fresh clone). It is **resumable**: already-computed
steps are skipped, so you can quit and come back — it continues where it stopped.

## ProteoCast (optional)

Computing the per-ABP mutational landscape (**ProteoCast**) is separate from
`Run / update` and takes **several hours** (one job per ABP on
[proteocast.ijm.fr](https://proteocast.ijm.fr)). Open **Documentation → Data management**, then use **Compute missing ProteoCast**. Results remain under **ABP conservation**. It is resumable.

## PyMOL (optional)

Only needed for the 3D scripts the pipeline writes under
`data/filtered/details/structures_files/bfactor_c70_interface/`. Install from
[pymol.org](https://pymol.org), then `File > Run Script…` (or `@/path/to.pml`).

## Shareable version (Streamlit Cloud)

A slim, read-only build (data pre-bundled, pipeline disabled) can be deployed to
Streamlit Community Cloud: `python script/make_slim_deploy.py` creates a
self-contained `deploy/` folder to push to a separate repo.

## Scientific app reorganization

See [diagnostic and roadmap](SCIENTIFIC_APP_ROADMAP.md) for the September 2026 reorganization, scientific checks still required and the local/public synchronization strategy.

Regression checks (in an environment with the app dependencies installed): `python -m unittest discover -s tests -v`.

## Updating the public dataset

From the full local environment, rebuild first, then preview a checked sync:

```sh
python tools/rebuild_local.py
python tools/sync_cloud_dataset.py --source /path/to/actin_project --destination /path/to/actin-abp-app
```

Add `--apply` to copy the verified snapshot. This preserves public-only imported
results, omits local calculation work files and keeps original local data intact.
Run `python tools/check_dataset.py` and `python tools/check_app.py` in the public
checkout before publication. Use this checked synchronization for updates rather
than the historical `script/make_slim_deploy.py` export.

The public `master` branch is the publication branch. The local research project
continues on `codex/scientific-app-foundations`. A successful push alone does not
prove that the remote Streamlit runtime has rebuilt successfully.
