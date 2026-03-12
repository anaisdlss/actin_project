# Projet interaction actin-ABP

## 1. Installation du projet et environnement

Assurez vous que pixi et python >3.11 (anaconda) soient installés.

Clonez le repertoire, allez au dossier du projet, et activez l'environnement pixi:

```bash
git clone https://github.com/anaisdlss/actin_project.git
cd actin_project
pixi install
pixi shell
```
## 2. Generate data

Using the provided scripts, it is possible to retrieve data from the PPI3D database using the UniProt identifier P60709.
This generates the summary results as well as the detailed information for each interaction.
The clustered results and their corresponding details are also retrieved.
Finally, additional information about the PDB entries associated with these results is extracted.

For reproducibility purposes, the generated identifiers (e.g. job IDs) are stored and reused until a new update of the PPI3D database occurs.
When the database is updated, it is recommended to rerun the scripts in order to regenerate the results and obtain the most up-to-date data.

The scripts in the script/ directory allow retrieval of interaction data from the PPI3D database using the UniProt identifier P60709 (Actin).

The pipeline is composed of four main steps.

### Data structure

```
data/
├── ppi3d_actin_summary.csv
├── pdb_entry_results.csv
├── details/
│   ├── 1.interactions.csv
│   ├── 2.proteins.csv
│   ├── 3.interface_residues.csv
│   ├── 4.inter-residue_contacts.csv
│   ├── 6.meta_alignement.csv
│   ├── 7.alignment_sequences.csv
│   ├── 8.structures.csv
│   ├── structures_files/
│       ├── pairwise/
│       ├── assembly/
│       ├── pymol/
├── clusters/
│   ├── clusters_summary.csv
│   ├── pdb_entry_cluster.csv
│   └── details/
│       ├── ...

```

### 1. Retrieve summary results

```bash
python script/get_summary_results.py
```
This script queries the PPI3D database and generates the summary of all interactions involving actin.
Output:
```code
data/ppi3d_actin_summary.csv
data/metadata.json
```
The metadata.json file stores the job identifier returned by PPI3D to ensure reproducibility.

### 2. Retrieve interaction details
```bash
python script/get_interaction_details.py
```
This script downloads the detailed information for each interaction found in the summary.
Outputs are stored in:

```code
data/details/
```
Structure files (PDB / CIF) are stored in:
```code
data/details/structures_files/
```

### 3. Retrieve clusters
```bash
python script/get_cluster.py
```

This script extracts the clustered interaction results from PPI3D.
Outputs:
```code
data/clusters/clusters_summary.csv
data/clusters/details/
```

### 4. Retrieve PDB entry interactions
```bash
python script/get_pdb_entries.py
```
This script retrieves all interactions present in each PDB entry identified in the summary results.
The script will ask whether to use:
```code
cluster summary or default summary
```
Outputs:
```code
data/pdb_entry_results.csv
```
or
```code
data/clusters/pdb_entry_results.csv
```
depending on the selected mode.



