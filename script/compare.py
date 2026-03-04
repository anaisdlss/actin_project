import pandas as pd

# fichiers
script_file = "ppi3d_actin_results.csv"
site_file = "ppi3d_internet.tsv"

# ouvrir les csv
df_script = pd.read_csv(script_file)
df_site = pd.read_csv(site_file, sep="\t")
# récupérer les PDB ID
script_pdb = set(df_script["PDB ID"].astype(str).str.strip().str.lower())
site_pdb = set(df_site["PDB ID"].astype(str).str.strip().str.lower())

print("PDB uniques dans script :", len(script_pdb))
print("PDB uniques dans site :", len(site_pdb))
print()

# différences
missing_in_script = site_pdb - script_pdb
missing_in_site = script_pdb - site_pdb
common = script_pdb & site_pdb

print("PDB communs :", len(common))
print()

print("PDB présents sur internet mais PAS dans le script :")
print(sorted(missing_in_script))
print()

print("PDB présents dans le script mais PAS sur internet :")
print(sorted(missing_in_site))
