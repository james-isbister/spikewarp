import spikewarp as sw

def create_meta_dict(cluster_pair_types, state_types):
	
	"""
	Creation of meta dictionary for storage of analysis results
	"""

	meta_dict = {}	
	
	for cluster_pair_type in cluster_pair_types:
		meta_dict[cluster_pair_type] = {}
		meta_dict[cluster_pair_type] = {
			'cluster_pairs_state_corr_pvals': [],
			'cluster_pairs_state_corr_rvals': [],
			'cluster_pairs_state_corr_r2vals': [],
			'explained_variance_ratio_by_PCA_component': [],
			'difference_prediction_from_state_pvalues': [],
			'difference_prediction_from_state_rvalues': [],
			'difference_prediction_from_state_rsquared_values': [],
			'difference_prediction_from_state_FA_predicted_sds': [],
			'difference_prediction_from_state_conj_orig_sds': [],
			'difference_prediction_from_state_cluster_orig_sds': [],
			'difference_prediction_from_state_state_predic_error_sds': [],

			'sch_single_unit_prediction_from_state_pvalues': [],
			'sch_single_unit_prediction_from_state_rvalues': [],
			'sch_single_unit_prediction_from_state_rsquared_values': [],
			'sch_single_unit_prediction_from_state_FA_predicted_sds': [],
			'sch_single_unit_prediction_from_state_conj_orig_sds': [],
			'sch_single_unit_prediction_from_state_cluster_orig_sds': [],
			'sch_single_unit_prediction_from_state_state_predic_error_sds': [],

			'proportion_of_pairwise_conj_trials_by_single_unit_ENOUGH_QUAD_CONJ':[],
			'proportion_of_pairwise_conj_trials_by_single_unit_ENOUGH_QUAD_CONJ_CORR':[],
			'proportion_of_pairwise_conj_trials_by_single_unit_ALL':[],
			'unique_pair_strings':[],
			'quadruplet_keys_ENOUGH_QUAD_CONJ_CORR':[],
			'quadruplet_keys_ENOUGH_QUAD_CONJ':[],
			'draw_dicts':[],
			'pairwise_state_slopes_Original':[],
			'pairwise_state_offsets_Original':[],
			'pairwise_state_slopes_Adjusted':[],
			'pairwise_state_offsets_Adjusted':[],

			'angle_conf_intervals':[],
		}

	meta_dict['stimulation_frequencies_of_condition'] = []


	if hasattr(sw, "placeholder_8"):
		sw.placeholder_8(meta_dict, state_types)

	
	return meta_dict



class PPDirectoryHolder(object):
	def __init__(self):
		
		"""
		Class for creation and storage of directories
		"""

def create_dir_holder(fig_dir, clustering_type):

	"""
	Function for creation and intiailisation of directory holder
	"""

	dh = PPDirectoryHolder()

	dh.clustering_fig_dir = fig_dir + clustering_type + "/"

	dh.meta_stats_fig_dir = dh.clustering_fig_dir + "MetaStats/"
	dh.pca_explained_var_meta_stats = dh.meta_stats_fig_dir + "PCA/"
	dh.state_corr_meta_stats = dh.meta_stats_fig_dir + "StateCorr/"
	dh.diff_corr_meta_stats = dh.meta_stats_fig_dir + "DiffCorr/"
	dh.pairs_of_clusters_times_and_state_corrs_fig_dir = dh.clustering_fig_dir + "PairsOfClusters/"
	dh.pairs_of_clusters_state_comp_fig_dir = dh.pairs_of_clusters_times_and_state_corrs_fig_dir + "StateComparison/"
	dh.single_pairwise_st_dewarp_fig_dir = dh.pairs_of_clusters_times_and_state_corrs_fig_dir + "SinglePairSTDewarp/"
	dh.single_pairwise_diff_dewarp_fig_dir = dh.pairs_of_clusters_times_and_state_corrs_fig_dir + "SinglePairDiffDewarp/"

	sw.makedirs(dh.meta_stats_fig_dir)
	sw.makedirs(dh.pca_explained_var_meta_stats)
	sw.makedirs(dh.state_corr_meta_stats)
	sw.makedirs(dh.diff_corr_meta_stats)
	sw.makedirs(dh.pairs_of_clusters_times_and_state_corrs_fig_dir)
	sw.makedirs(dh.pairs_of_clusters_state_comp_fig_dir)
	sw.makedirs(dh.single_pairwise_st_dewarp_fig_dir)
	sw.makedirs(dh.single_pairwise_diff_dewarp_fig_dir)
	
	if hasattr(sw, "placeholder_7"):
		sw.placeholder_7(dh)
	
	return dh



"""
Default plotting options
"""

import copy

default_plot_opt = {'plot_dewarp_and_state_dependence': False,
			'plot_pairwise_dewarp': True,
			'plot_multi_neuron_timewarp_single_images': True,
			'plot_multi_neuron_timewarp_videos': False,
			'plot_pairwise_corr_mat':True,
			'plot_multi_quad_timewarp_videos':False, 
			'make_single_quad_vid':False,
			'plot_cluster_pair_multi_stats':True}

def new_default_plot_opt():

	"""
	Function returns a copy of the default plotting options
	"""

	return copy.copy(default_plot_opt)



