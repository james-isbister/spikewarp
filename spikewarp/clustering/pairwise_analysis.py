import numpy as np
import multiprocessing
import time

import spikewarp as sw

def run_pairwise_analysis(random_seed_counter,
						experiment_code,
						data_loader,
						analysis):

	"""
	This function performs the parallel clustering (over parameter range) for cell pairs where both spiked above atleast twice above threshold
	"""
  	
  	# Start timer
	single_experiement_pairwise_start = time.time()

	# Use custom condition indices if specified
	condition_indices_to_process = data_loader.condition_indices_to_process
	if ("quick_custom_conditions_to_process_for_experiment" in analysis):
		condition_indices_to_process = analysis["quick_custom_conditions_to_process_for_experiment"]

	# Prepare data holders and multiprocessing pool
	mcah_for_each_shuffle_type = {shuffle_option_string:sw.MetaClusterAnalysisHolder(shuffle_option_string) for shuffle_option_string in analysis['option_helper'].shuffle_option_strings}
	tw_files_by_stwcboo_and_shuffle_option = prepare_tw_files_by_stwcboo_and_shuffle_option_dict(analysis)
	single_cluster_holders_by_condition_and_shuffle_option = prepare_single_cluster_holders_by_condition_and_shuffle_option_for_parallel_unpack(condition_indices_to_process, analysis)
	pool = multiprocessing.Pool(processes=analysis['number_of_processors']); results = []

	for condition_index in condition_indices_to_process: 

		# Get SCSUOs for condition index
		list_of_queried_scsuo = [scsuo for scsuo in data_loader.list_of_all_scsuos if (scsuo.condition_index == condition_index)]
		number_of_units = len(list_of_queried_scsuo)
		if (number_of_units > 0):
			
			# Find cell pairs where both spiked above atleast twice above threshold
			number_of_trials_for_condition = list_of_queried_scsuo[0].number_of_trials_for_condition
			did_spike_atleast_once_for_each_gid_BY_TRIAL = np.zeros((number_of_trials_for_condition, number_of_units), dtype=bool)
			for unit_index in range(number_of_units):
				scsuo = list_of_queried_scsuo[unit_index]
				did_spike_atleast_once_for_each_gid_BY_TRIAL[:, unit_index] = scsuo.true_if_spiked_atleast_once_on_trial.flatten()	
			pairwise_did_both_spike_trial_count_matrix = np.dot(did_spike_atleast_once_for_each_gid_BY_TRIAL.T.astype(int), did_spike_atleast_once_for_each_gid_BY_TRIAL.astype(int))
			true_where_did_both_spike_count_above_threshold = np.logical_and(pairwise_did_both_spike_trial_count_matrix >= 20, pairwise_did_both_spike_trial_count_matrix >= 20)
			pairs_of_cell_indices_where_both_did_spike_above_threshold = np.asarray(np.nonzero(np.triu(true_where_did_both_spike_count_above_threshold, 1)))
			number_of_neuron_pairs_where_both_did_spike_count_above_threshold = pairs_of_cell_indices_where_both_did_spike_above_threshold.shape[1]

			# Draw plots and update meta data
			sw.draw_basic_similarity_matrix(pairwise_did_both_spike_trial_count_matrix, analysis['directory_holder'].pairwise_reliabilities_directory + experiment_code + "_" + str(condition_index) + "_pairwise_did_both_spike_trial_count_matrix")
			indices_above_diagonal = np.triu_indices(number_of_units, k=1)
			all_unique_pairwise_did_both_spike_trial_count = pairwise_did_both_spike_trial_count_matrix[indices_above_diagonal[0], indices_above_diagonal[1]]
			sw.normal_histo_plot([all_unique_pairwise_did_both_spike_trial_count], analysis['directory_holder'].pairwise_reliabilities_directory + experiment_code + "_" + str(condition_index) + "_did_both_spike_on_trial_count_histogram", bins=number_of_trials_for_condition, histo_range=[0.0, float(number_of_trials_for_condition + 1)])
			mcah_for_each_shuffle_type['normal-1'].extend_pairwise_reliabilities(all_unique_pairwise_did_both_spike_trial_count.astype(float)/list_of_queried_scsuo[0].number_of_trials_for_condition, all_unique_pairwise_did_both_spike_trial_count)
			mcah_for_each_shuffle_type['normal-1'].plot_pairwise_reliability_plots(analysis['directory_holder'].pairwise_reliabilities_directory, experiment_code)

			# Iterate through cell pairs where both spiked atleast twice above threshold
			for threshold_pair_index in range(number_of_neuron_pairs_where_both_did_spike_count_above_threshold):
				cell_0_index = pairs_of_cell_indices_where_both_did_spike_above_threshold[0, threshold_pair_index]
				cell_1_index = pairs_of_cell_indices_where_both_did_spike_above_threshold[1, threshold_pair_index]
				scsuo_0 = list_of_queried_scsuo[cell_0_index]
				scsuo_1 = list_of_queried_scsuo[cell_1_index]
				random_seed_counter = random_seed_counter + 1
				async_object = pool.apply_async(sw.ClusteringParamIterator, args=(experiment_code,
																				random_seed_counter,
																				scsuo_0, 
																				scsuo_1, 
																				condition_index, 
																				analysis))

				results.append(async_object)

	

	# Unpack parallel results
	result_objects = [r.get() for r in results]; pool.close(); pool.join()
	spike_pairs_single_clusterings_by_shuffle_option = [result_object.single_clusterings_by_shuffle_option for result_object in result_objects]
	unpack_single_cluster_holders(spike_pairs_single_clusterings_by_shuffle_option, mcah_for_each_shuffle_type, single_cluster_holders_by_condition_and_shuffle_option)
	add_appropriate_pairwise_response_dist_files_from_result_objects(analysis, tw_files_by_stwcboo_and_shuffle_option, result_objects)

	# End timer and return
	print("Experiment: " + str(experiment_code) + ", Pairwise Analysis End: " + str(time.time() - single_experiement_pairwise_start))
	return mcah_for_each_shuffle_type, random_seed_counter, tw_files_by_stwcboo_and_shuffle_option, single_cluster_holders_by_condition_and_shuffle_option 


def add_appropriate_pairwise_response_dist_files_from_result_objects(analysis, 
																	tw_files_by_stwcboo_and_shuffle_option,
																	result_objects):

	"""
	Function adds the appropriate pairwise response distribution files to list for later gif creation
	"""

	if (analysis['do_plot_super_time_warp_cluster_busters'] & analysis['do_plot_super_time_warp_cluster_buster_videos_ALL_PAIRS']):
		for result_object in result_objects:
			for shuffle_option in analysis['option_helper'].shuffle_options:
				
				shuffle_option_string = shuffle_option[0] + str(shuffle_option[1])

				for stwcboo_index, stwcboo in enumerate(analysis['super_time_warp_cluster_buster_option_objects']):
					tw_files_by_stwcboo_and_shuffle_option[shuffle_option_string][stwcboo_index].extend(result_object.tw_files_by_stwcboo_by_shuffle_option[shuffle_option_string][stwcboo_index])

					

def prepare_tw_files_by_stwcboo_and_shuffle_option_dict(analysis):
	
	tw_files_by_stwcboo_and_shuffle_option = {}
	for shuffle_option in analysis['option_helper'].shuffle_options:
		shuffle_option_string = shuffle_option[0] + str(shuffle_option[1])
		tw_files_by_stwcboo_and_shuffle_option[shuffle_option_string] = [[] for stwcboo in analysis['super_time_warp_cluster_buster_option_objects']]

	return tw_files_by_stwcboo_and_shuffle_option


def prepare_single_cluster_holders_by_condition_and_shuffle_option_for_parallel_unpack(condition_indices_to_process, analysis):

	single_cluster_holders_by_condition_and_shuffle_option = {}
	for condition_index in condition_indices_to_process: 
		single_cluster_holders_by_condition_and_shuffle_option[str(condition_index)] = {}
		for shuffle_option in analysis['option_helper'].shuffle_options:
			shuffle_option_string = shuffle_option[0] + str(shuffle_option[1])		
			single_cluster_holders_by_condition_and_shuffle_option[str(condition_index)][shuffle_option_string] = []

	return single_cluster_holders_by_condition_and_shuffle_option


def unpack_single_cluster_holders(spike_pairs_single_clusterings_by_shuffle_option, mcah_for_each_shuffle_type, single_cluster_holders_by_condition_and_shuffle_option):
	
	for spike_pair_single_clusterings_by_shuffle_option in spike_pairs_single_clusterings_by_shuffle_option:
		for shuffle_option_string, single_clusterings_for_shuffle_option in spike_pair_single_clusterings_by_shuffle_option.items():

			for single_clustering in single_clusterings_for_shuffle_option:
				mcah_for_each_shuffle_type[single_clustering.shuffle_option_string].extend_standard_cluster_arrays(single_clustering)

				for single_cluster_holder in single_clustering.secondary_single_cluster_holders:

					if (single_cluster_holder.final_stage_2_cluster):

						single_cluster_holders_by_condition_and_shuffle_option[str(single_clustering.condition_index)][shuffle_option_string].extend([single_cluster_holder])
