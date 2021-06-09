import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import warnings
import math

import spikewarp as sw

def process_schs(schs_for_cond, exp_cond_str):

	# Prepare single experimental condition stationary and made stationary clusters for analysis

	unique_cond_units = []
	unique_unit_info = []
	unique_pair_strings = []
	non_stationary_to_pop = []

	units_holder = None
	if hasattr(sw, "UnitsHolder"):
		units_holder = UnitsHolder(exp_cond_str)

	for sch_index, sch in enumerate(schs_for_cond):

		sch.cluster_type = ''

		if (len(sch.stage_2_cluster.analysis_dicts['OriginalTestsPassed']['PCA_predicited_state_for_cluster_trials']) > 0):

			sch.cluster_type = 'OriginalTestsPassed'
			sch.trials = sch.stage_2_cluster.real_trial_indices_of_cluster
			sch.spikes = sch.stage_2_cluster.cluster_spike_pairs
			sch.num_spike_pairs = sch.spikes.shape[0]

		
		elif (len(sch.stage_2_cluster.analysis_dicts['SelectivelyDifferencedBoxJenkinsTestsPassed']['PCA_predicited_state_for_cluster_trials']) > 0):

			sch.cluster_type = 'SelectivelyDifferencedBoxJenkinsTestsPassed'
			sch.trials = sch.stage_2_cluster.real_trial_indices_of_cluster[1:]
			sch.spikes = sch.stage_2_cluster.selective_differences_scaled_down + sch.stage_1_cluster.mean
			sch.num_spike_pairs = sch.spikes.shape[0]

		sch.unit_indices = [sch.scsuo_0_channel_unit_pair_index, sch.scsuo_1_channel_unit_pair_index]
		sch.unique_pair_str = exp_cond_str + "_" + str(np.min(sch.unit_indices)) + "_" + str(np.max(sch.unit_indices))

		if (sch.cluster_type != ''):

			sch.proportion_of_conjunctive_trials_in_cluster = sch.num_spike_pairs/sch.number_of_conjunctive_trials

			sch.predicted_states = np.asarray(sch.stage_2_cluster.analysis_dicts[sch.cluster_type]['PCA_predicited_state_for_cluster_trials']).flatten()				
			sch.distances_from_factor_line = sch.stage_2_cluster.analysis_dicts[sch.cluster_type]['distances_from_factor_line']
			sch.cluster_correlation_r_value = sch.stage_2_cluster.analysis_dicts[sch.cluster_type]['LR_rvalue']
			sch.cluster_correlation_p_value = sch.stage_2_cluster.analysis_dicts[sch.cluster_type]['LR_pvalue']
			sch.predicted_state_dependent_differences_sd = sch.stage_2_cluster.analysis_dicts[sch.cluster_type]['FA_1sd_estimated_difference_area_BS_mean']
			sch.Conjunctive_original_FS_diff_SD = np.std(sch.spike_pair_distribution_spikes[:, 1] - sch.spike_pair_distribution_spikes[:, 0])
			sch.is_PCA_BS_empirical_pvalue_different_from_45 = sch.stage_2_cluster.analysis_dicts[sch.cluster_type]['is_PCA_BS_empirical_pvalue_different_from_45']
			sch.FA_N_sds = [sch.stage_2_cluster.analysis_dicts[sch.cluster_type]['FA_N0_sd_BS_mean'], sch.stage_2_cluster.analysis_dicts[sch.cluster_type]['FA_N1_sd_BS_mean']]
			sch.Conj_original_N_sds = np.std(sch.spike_pair_distribution_spikes, axis=0)

			v = sch.predicted_states
			sch.normalised_predicted_states = v / np.std(v)
			sch.normalised_predicted_states_first_correction = sch.normalised_predicted_states
			
			unique_pair_strings.append(sch.unique_pair_str)
			sch.unit_pair_str = str(sch.scsuo_0_channel_unit_pair_index) + str(sch.scsuo_1_channel_unit_pair_index)

			if (units_holder != None):
				units_holder.placeholder_4(sch)		
				

		else:
			non_stationary_to_pop.append(sch_index)
		
	if (units_holder != None):
		units_holder.placeholder_5(schs_for_cond)

	
	return non_stationary_to_pop, unique_pair_strings, units_holder