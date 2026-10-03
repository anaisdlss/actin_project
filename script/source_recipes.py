"""Recovered producer recipes, distinct from proof of a historical execution."""
from pathlib import Path

RECIPES = {
    'data/exports/abp_site_domain/abp_interpro.csv': 'script/abp_site_domain/03_interpro.py',
    'data/exports/abp_site_domain/abp_master.csv': 'script/abp_site_domain/04_crosstab.py',
    'data/exports/abp_site_domain/site_cluster_summary.csv': 'script/abp_site_domain/04_crosstab.py',
    'data/exports/abp_site_domain/interface_tm_sweep.csv': 'script/abp_site_domain/08_interface_tm_sweep.py',
    'data/exports/abp_site_domain/familles.csv': 'script/abp_site_domain/11_family_regroup.py',
    'data/exports/abp_site_domain/convergences_inter_familles.csv': 'script/abp_site_domain/11_family_regroup.py',
    'data/exports/abp_site_domain/actin_footprint_overlap.csv': 'script/abp_site_domain/13_actin_footprint_overlap.py',
    'data/exports/abp_site_domain/folddisco_interface.csv': 'script/abp_site_domain/13_folddisco_interface.py',
    'data/exports/abp_site_domain/folddisco_discovery.csv': 'script/abp_site_domain/14_folddisco_discovery.py',
    'data/exports/abp_site_domain/folddisco_names.csv': 'script/abp_site_domain/14b_resolve_names.py',
    'data/exports/abp_site_domain/interface_chemistry.csv': 'script/abp_site_domain/16_interface_chemistry.py',
    'data/exports/abp_site_domain/interface_secondary_structure.csv': 'script/abp_site_domain/15_interface_secondary_structure.py',
    'data/exports/abp_site_domain/actin_residue_determinants.csv': 'script/abp_site_domain/18_site_determinants.py',
    'data/exports/abp_site_domain/interface_domain_by_abp.csv': 'script/abp_site_domain/19_interface_domain_family.py',
    'data/exports/abp_site_domain/interface_family_refined.csv': 'script/abp_site_domain/19_interface_domain_family.py',
    'data/filtered/proteins_per_pdb.csv': 'script/graphe_filter.py',
    'data/filtered/filtered_all_data.csv': 'script/graphe_filter.py; script/binding_site_superclusters.py',
    'data/filtered/filtered_summary.csv': 'script/graphe_filter.py',
    'data/filtered/filtered_pdb_entry.csv': 'script/graphe_filter.py',
    'data/filtered/actin_filament_positions.csv': 'script/actin_filament_positions.py',
    'data/raw/all_data.csv': 'script/data_extract/get_cluster_table.py',
    'data/raw/pdb_entry_results.csv': 'script/data_extract/get_pdb_entries.py',
    'data/raw/ppi3d_actin_summary.csv': 'script/data_extract/get_summary_results.py',
    'data/exports/interactions_for_mafft.csv': 'script/cluster_interaction_analysis.py',
    'data/filtered/patches_infos_s1_binding_site.csv': 'script/cluster_interaction_analysis.py; script/interface_analysis_s1.py',
    'data/filtered/patches_infos_cluster_data_70.csv': 'script/cluster_interaction_analysis.py',
    'data/filtered/binding_site_superclusters.csv': 'script/binding_site_superclusters.py',
    'data/filtered/s1_cluster_reference.csv': 'script/interface_analysis_s1.py',
    'data/filtered/s1_superclusters.csv': 'script/data_analysis.py',
    'data/proteocast/abp_inputs/manifest.csv': 'script/proteocast_prep_manifest.py',
    **{f'data/filtered/details/{name}': 'script/data_extract/get_interaction_details.py; script/mafft_pipeline.py'
       for name in ['1.interactions.csv', '2.proteins.csv', '3.interface_residues.csv',
                    '4.inter-residue_contacts.csv', '5.ligands.csv', '6.meta_alignement.csv',
                    '7.alignment_sequences.csv', '8.structures.csv']},
}


def recovered_origin(root, relative):
    if Path(relative).name.startswith('filament_accessibility_7pdz_I'):
        return ('archived previous reference', 'tools/calculate_filament_accessibility.py (Git history)',
                'Bovine 7PDZ replaced by human 8DNH. Original method receipt retained; excluded from current RSA.')
    if relative in RECIPES:
        recipe = RECIPES[relative]
        if all((Path(root)/name).exists() for name in recipe.split('; ')):
            return ('producer recipe found; historical build unverified', recipe,
                    'The producing code is present. No original input/code/output receipt proves how this historical copy was built.')
    if relative.startswith('data/proteocast/abp/') and relative.endswith('/4.query_ProteoCast.csv'):
        return ('external model output; verify job provenance', 'ProteoCast',
                'Imported/downloaded ProteoCast mutation scores. These are not locally recomputed; query sequence validation remains required.')
    if relative.startswith('data/proteocast/abp/') and relative.endswith('/rsa_values.csv'):
        return ('external bundle; RSA recipe unresolved', 'ProteoCast bundle',
                'RSA distributed with the model result. Its calculation must not be confused with the locally documented human 8DNH RSA.')
    if relative.startswith('data/proteocast/abp/') and Path(relative).name in {
            '8.query_Segmentation.csv', '14.query_GEMME_pLDDT.csv', 'mapping_fasta_pdb.csv'}:
        return ('external bundle; original execution unverified', 'ProteoCast supplementary output',
                'File type documented by the ProteoCast bundle README_static.md. No original job/version inferred from the filename; retain as imported evidence.')
    return None
