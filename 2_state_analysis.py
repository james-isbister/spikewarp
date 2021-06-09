'''
2_state_analysis.py

This script can be run once the first script has completed for the two "analyses_to_process" options:
- "CustomDBSCANExtra"
- "Unclustered"



### Configuration

1. Update the variable "analysis_1_output" on line 45 of this script to the path originally specified for the variable "figures_root_directory" in "recipes/investigation_recipe.json". For example:
analysis_1_output = "/Users/james/IsbisterEtAl-SciRep2021-Figures/",

2. Update the variable "number_of_processors" on line 46 to the number of threads you would like to use (initially set to 8).



### Options

1. To produce figures illustrating example correlations between clusters line 47 may be uncommented.


### Run

The script can be run with the following command:
python 2_state_analysis.py

'''


import numpy as np
import itertools
import multiprocessing
import pickle
import matplotlib.pyplot as plt

import spikewarp as sw


######################################
# CONFIGURATION
######################################

analysis_1_output = [PATH_TO_BE_SET_BY_USER]
processes=8
plot_opt = sw.new_default_plot_opt()
# plot_opt['plot_cluster_pair_multi_stats'] = True



######################################
# SETUP
######################################
time_window_dir = "5.4-55.4/"
fig_dir = analysis_1_output + time_window_dir + 'StateAnalysis/'
pickle_root = analysis_1_output + time_window_dir + "PairwiseAnalysis/"
# clustering_types = ['CustomDBSCAN', 'CustomDBSCANExtra', 'CustomGaussianMixtureInfCrit', 'Unclustered']
clustering_types = ['CustomDBSCANExtra', 'Unclustered']
state_types = ['Original', 'Adjusted', "MeanFilledPCA"]
cluster_pair_types = ["triplet", "quadruplet", "all_clusters", "including_non_stationary"]



######################################
# STEP 1 - PROCESS DIFFERENT CLUSTERING TYPES
######################################

# 1.1 Prepare data holders and multiprocessing_pool
meta_dict_by_clustering_type = {}
exp_cond_units_holders_by_clustering_type = {}
exp_cond_pop_state_procs_by_clustering_type = {}
exp_cond_schs_for_conds_by_clustering_type = {}
dh_by_clustering_type = {}
results = []
multiprocessing_pool = multiprocessing.Pool(processes=processes)

# 1.2 Run analysis for each clustering type
for clustering_type in clustering_types:

	dh = sw.create_dir_holder(fig_dir, clustering_type)
	dh_by_clustering_type[clustering_type] = dh

	meta_dict_by_clustering_type[clustering_type] = sw.create_meta_dict(cluster_pair_types, state_types)
	exp_cond_units_holders_by_clustering_type[clustering_type] = {}
	exp_cond_pop_state_procs_by_clustering_type[clustering_type] = {}
	exp_cond_schs_for_conds_by_clustering_type[clustering_type] = {}

	clustering_pickle_file = pickle_root + clustering_type + "/single_cluster_holders_by_exp_cond_and_shuffle_option.pickle"
	with open(clustering_pickle_file, 'rb') as handle:
		single_cluster_holders_by_exp_cond_and_shuffle_option = pickle.load(handle)

	for exp_key, exp_dict in single_cluster_holders_by_exp_cond_and_shuffle_option.items():
		for cond_key, cond_dict in exp_dict.items():
			results.extend([multiprocessing_pool.apply_async(sw.ClusteringTypeIterator, args=(clustering_type, plot_opt, state_types, cond_dict['normal-1'], dh, "E" + exp_key + "_C" + cond_key, None, cluster_pair_types))])
	
# 1.3 Unpack parallel results
list_of_results = [r.get() for r in results]
multiprocessing_pool.close(); multiprocessing_pool.join()
sw.unpack_parralel_results(list_of_results, meta_dict_by_clustering_type, exp_cond_units_holders_by_clustering_type, exp_cond_pop_state_procs_by_clustering_type, exp_cond_schs_for_conds_by_clustering_type, clustering_types, state_types, dh_by_clustering_type)

# 1.4 Create misc meta plots
sw.normal_histo_plot([meta_dict_by_clustering_type['CustomDBSCANExtra']['including_non_stationary']['angle_conf_intervals'], meta_dict_by_clustering_type['Unclustered']['including_non_stationary']['angle_conf_intervals']], fig_dir + 'Confidence_interval_widths', bins=18, histo_range=[0.0, 90.0], x_axis_left_buffer=0.01, x_axis_label="Confidence interval width", y_axis_label="Frequency", custom_x_tick_locators=[90.0, 10.0], custom_y_tick_locators=[5, 5], alpha=0.78, labels=['CustomDBSCANExtra', 'Unclustered'])

including_non_stationary_clust_d = sw.check_intersect_and_union_for_cluster_pair_type("CustomDBSCANExtra", "Unclustered", meta_dict_by_clustering_type, 'including_non_stationary', 'unique_pair_strings', optional_figure_path=fig_dir + 'including_non_stationary_ClustVenn.pdf')
all_clust_d = sw.check_intersect_and_union_for_cluster_pair_type("CustomDBSCANExtra", "Unclustered", meta_dict_by_clustering_type, 'all_clusters', 'unique_pair_strings', optional_figure_path=fig_dir + 'AllClustVenn.pdf')
corr_quadruplet_d  = sw.check_intersect_and_union_for_cluster_pair_type("CustomDBSCANExtra", "Unclustered", meta_dict_by_clustering_type, 'quadruplet', 'quadruplet_keys_ENOUGH_QUAD_CONJ_CORR', optional_figure_path=fig_dir + 'ENOUGH_QUAD_CONJ_CORR_Venn.pdf')
quadruplet_d  = sw.check_intersect_and_union_for_cluster_pair_type("CustomDBSCANExtra", "Unclustered", meta_dict_by_clustering_type, 'quadruplet', 'quadruplet_keys_ENOUGH_QUAD_CONJ', optional_figure_path=fig_dir + 'ENOUGH_QUAD_CONJ.pdf')



######################################
# STEP 2 - CREATE CustomDBSCANExtra and UNCLUSTERED POOL AND RUN ANALYSIS
######################################

# 2.1 Create pool of clusters and unit holders by condition
exp_cond_unit_holders_conj_all, schs_for_conds_conj_all, exp_cond_key_conj_all, meta_dict_conj, dh_by_conj = sw.create_conj_pool(exp_cond_units_holders_by_clustering_type, exp_cond_pop_state_procs_by_clustering_type, exp_cond_schs_for_conds_by_clustering_type, "CustomDBSCANExtra", "Unclustered", fig_dir, cluster_pair_types, state_types)

# 2.2 Prepare data holders and multiprocessing_pool
multiprocessing_pool = multiprocessing.Pool(processes=processes)
results = []

# 2.3 Run analysis for pool by condition
for exp_cond_unit_holder, schs_for_cond, exp_cond_key in zip(exp_cond_unit_holders_conj_all, schs_for_conds_conj_all, exp_cond_key_conj_all):
	results.extend([multiprocessing_pool.apply_async(sw.ClusteringTypeIterator, args=('Conj', plot_opt, state_types, schs_for_cond, dh_by_conj['Conj'], exp_cond_key + "_Conj", exp_cond_unit_holder, cluster_pair_types))])

# 2.3 Unpack parallel results
list_of_results = [r.get() for r in results]
multiprocessing_pool.close(); multiprocessing_pool.join()
sw.unpack_parralel_results(list_of_results, meta_dict_conj, {}, {}, {}, ['Conj'], state_types, dh_by_conj)

# 2.4 Create misc meta plots
meta_dicts_for_comp = [meta_dict_conj['Conj'], meta_dict_by_clustering_type["CustomDBSCANExtra"], meta_dict_by_clustering_type["Unclustered"]]
labels = ['Conj', 'CustomDBSCANExtra', 'Unclustered']
sw.normal_histo_plot([meta_dict['quadruplet']['cluster_pairs_state_corr_pvals'] for meta_dict in meta_dicts_for_comp], fig_dir + 'conj_dbscan_and_unclustered_cluster_pairs_state_corr_pvals', bins=20, histo_range=[0.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="p value", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[5, 5], alpha=0.78, labels=labels)
sw.normal_histo_plot([meta_dict['quadruplet']['cluster_pairs_state_corr_rvals'] for meta_dict in meta_dicts_for_comp], fig_dir + 'conj_dbscan_and_unclustered_cluster_pairs_state_corr_rvals', bins=40, histo_range=[-1.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="r value", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[5, 5], alpha=0.78, labels=labels)
sw.normal_histo_plot([meta_dict['quadruplet']['cluster_pairs_state_corr_r2vals'] for meta_dict in meta_dicts_for_comp], fig_dir + 'conj_dbscan_and_unclustered_cluster_pairs_state_corr_r2vals', bins=20, histo_range=[0.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$r^2$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[5, 5], alpha=0.78)
sw.normal_histo_plot([meta_dict['quadruplet']['difference_prediction_from_state_pvalues'] for meta_dict in meta_dicts_for_comp], fig_dir + 'conj_dbscan_and_unclustered_cluster_pairs_NON45_difference_prediction_from_state_pvalues', bins=20, histo_range=[0.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="p value", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[5, 5], alpha=0.78)
sw.normal_histo_plot([meta_dict['quadruplet']['difference_prediction_from_state_rsquared_values'] for meta_dict in meta_dicts_for_comp], fig_dir + 'conj_dbscan_and_unclustered_cluster_pairs_NON45_difference_prediction_from_state_rsquared_values', bins=20, histo_range=[0.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$r^2$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[5, 5], alpha=0.78)

if hasattr(sw, "placeholder_8"):
	sw.placeholder_8(meta_dicts_for_comp, fig_dir, meta_dict_by_clustering_type)


