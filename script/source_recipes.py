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
    'data/exports/abp_site_domain/actin_residue_determinants.csv': 'script/abp_site_domain/18_site_determinants.py',
    'data/exports/abp_site_domain/interface_domain_by_abp.csv': 'script/abp_site_domain/19_interface_domain_family.py',
    'data/exports/abp_site_domain/interface_family_refined.csv': 'script/abp_site_domain/19_interface_domain_family.py',
    'data/filtered/proteins_per_pdb.csv': 'script/graphe_filter.py',
    'data/filtered/filtered_all_data.csv': 'script/graphe_filter.py; script/binding_site_superclusters.py',
    'data/filtered/filtered_summary.csv': 'script/graphe_filter.py',
    'data/filtered/filtered_pdb_entry_results.csv': 'script/graphe_filter.py',
    'data/filtered/actin_filament_positions.csv': 'script/actin_filament_positions.py',
    **{f'data/filtered/details/{name}': 'script/data_extract/get_interaction_details.py; script/mafft_pipeline.py'
       for name in ['1.interactions.csv', '2.proteins.csv', '3.interface_residues.csv',
                    '4.inter-residue_contacts.csv', '5.ligands.csv', '6.meta_alignement.csv']},
}


def recovered_origin(root, relative):
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
                'RSA distributed with the model result. Its calculation must not be confused with the locally documented 7PDZ RSA.')
    return None
