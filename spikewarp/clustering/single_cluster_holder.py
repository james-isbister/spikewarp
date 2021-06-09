from scipy import linalg
from scipy import stats

from sklearn.metrics import r2_score # Temp used to avoid deprecation warning

import math
import statsmodels.api as sm
import numpy as np

import spikewarp as sw

"""
Cluster options
"""

class SingleClusterHolderOptions(object):

	def __init__(self):
	
		self.intersection_dist = 3.0
		self.inner_bound = 3.0
		self.outlier_bound = 4.0
		self.rings = [self.inner_bound, self.outlier_bound]

		self.minimum_elements_in_stage_1_cluster = 30


"""
Class for creating and managing Stage 1 and 2 clusters for single original cluster
"""

class SingleClusterHolder(object):

	def __init__(self, single_cluster_holder_options, original_cluster_indices, first_spike_pairs_for_original_cluster, spike_pair_distribution_spikes, number_of_conjunctive_trials, number_of_trials_for_condition, non_singleton_index, scsuo_0_channel_unit_pair_index, scsuo_1_channel_unit_pair_index, condition_index, real_trial_indices_on_which_both_neurons_spiked, stimulation_frequency_of_condition, neurophys_layer_strings_of_unit_pair, exc_inh_types_of_unit_pair, barrels_of_unit_pair):

		self.single_cluster_holder_options = single_cluster_holder_options
		
		self.original_cluster_indices = original_cluster_indices
		self.first_spike_pairs_for_original_cluster = first_spike_pairs_for_original_cluster
		self.spike_pair_distribution_spikes = spike_pair_distribution_spikes

		self.number_of_conjunctive_trials = number_of_conjunctive_trials
		self.number_of_trials_for_condition = number_of_trials_for_condition
		self.non_singleton_index = non_singleton_index

		self.scsuo_0_channel_unit_pair_index = scsuo_0_channel_unit_pair_index
		self.scsuo_1_channel_unit_pair_index = scsuo_1_channel_unit_pair_index

		self.condition_index = condition_index

		self.real_trial_indices_on_which_both_neurons_spiked = real_trial_indices_on_which_both_neurons_spiked

		self.stimulation_frequency_of_condition = stimulation_frequency_of_condition
		self.neurophys_layer_strings_of_unit_pair = neurophys_layer_strings_of_unit_pair
		self.exc_inh_types_of_unit_pair = exc_inh_types_of_unit_pair
		self.barrels_of_unit_pair = barrels_of_unit_pair

		self.interesting_to_plot_super_time_warp_cluster_buster = False

		self.p_value_threshold = 0.005

		self.final_stage_1_cluster = False
		self.final_stage_2_cluster = False

		self.stage_1_cluster = None
		self.stage_2_cluster = None



	def create_intersection_ellipse_and_stage_1_cluster(self, size_limit_of_new_cluster, with_custom_stage_1_cluster_expansion=True):

		mean_point_first_spike_times_for_original_cluster = [np.mean(self.first_spike_pairs_for_original_cluster[:, 0]), np.mean(self.first_spike_pairs_for_original_cluster[:, 1])]
		std_first_spike_neuron_0_for_cluster = np.std(self.first_spike_pairs_for_original_cluster[:, 0])
		std_first_spike_neuron_1_for_cluster = np.std(self.first_spike_pairs_for_original_cluster[:, 1])
		self.flat_intersection_ellipse = sw.create_ellipse(mean_point_first_spike_times_for_original_cluster, max(0.5, std_first_spike_neuron_0_for_cluster*self.single_cluster_holder_options.intersection_dist), max(0.5, std_first_spike_neuron_1_for_cluster*self.single_cluster_holder_options.intersection_dist), 0.0) 

		if ((std_first_spike_neuron_0_for_cluster > 0.0) & (std_first_spike_neuron_1_for_cluster > 0.0)):

			self.stage_1_cluster = sw.Cluster()
			if (with_custom_stage_1_cluster_expansion):
				self.stage_1_cluster.create_from_pairs(self.single_cluster_holder_options, 
											self.first_spike_pairs_for_original_cluster,
											self.spike_pair_distribution_spikes, 
											self.number_of_conjunctive_trials, 
											self.number_of_trials_for_condition,
											self.real_trial_indices_on_which_both_neurons_spiked)
			else:
				self.stage_1_cluster.create_nonflat_cluster_from_cluster_spikes(self.single_cluster_holder_options, 
																				self.first_spike_pairs_for_original_cluster, 
																				self.original_cluster_indices, 
																				self.number_of_conjunctive_trials, 
																				self.number_of_trials_for_condition,
																				self.real_trial_indices_on_which_both_neurons_spiked)

			
			if ((self.stage_1_cluster.number_of_samples_in_ellipse >= self.single_cluster_holder_options.minimum_elements_in_stage_1_cluster) & (self.stage_1_cluster.number_of_samples_in_ellipse < size_limit_of_new_cluster)):

				self.final_stage_1_cluster = True
				self.stage_1_cluster.create_stage_1_cluster_statistics_dictionary()

				if (self.stage_1_cluster.analysis_dicts['FlatClusterStats']['FlatCluster_LR_pvalue'] < self.p_value_threshold):
					self.final_stage_2_cluster = True
					self.interesting_to_plot_super_time_warp_cluster_buster = True
				


	def calculate_stage_2_cluster_with_time_series_and_factor_analysis(self, custom_stage_2_estimation=False):

		number_of_bootstrap_iterations = 10000

		self.calculate_bootstrap_pca_cluster(number_of_bootstrap_iterations, "Original", custom_stage_2_estimation=custom_stage_2_estimation) # Creates Stage 2 cluster!
		self.calculate_bootstrap_pca_cluster(number_of_bootstrap_iterations, "Differenced", spikes_to_use=self.stage_2_cluster.pairs_differenced)
		if (self.stage_2_cluster.use_selective_differences):
			self.calculate_bootstrap_pca_cluster(number_of_bootstrap_iterations, "SelectivelyDifferenced", spikes_to_use=self.stage_2_cluster.selective_differences_scaled_down)
			self.calculate_bootstrap_pca_cluster(number_of_bootstrap_iterations, "SelectivelyDifferencedBoxJenkins", spikes_to_use=self.stage_2_cluster.selective_differences_arima_residuals)

		sw.bootstrap_factor_analysis(self.stage_2_cluster, 200, "Original")	
		sw.bootstrap_factor_analysis(self.stage_2_cluster, 200, "SelectivelyDifferenced")	
		sw.bootstrap_factor_analysis(self.stage_2_cluster, 200, "SelectivelyDifferencedBoxJenkins")	



	def calculate_bootstrap_pca_cluster(self, number_of_bootstrap_iterations, iteration_key, custom_stage_2_estimation=False, spikes_to_use=None):

		if (iteration_key == "Original"):

			number_of_bootstrap_steps = 1
			if (custom_stage_2_estimation):
				number_of_bootstrap_steps = 3

			previous_step_cluster = self.stage_1_cluster
			for bootstrap_step in range(number_of_bootstrap_steps):

				bootstrap_pca_data_dict = sw.bootstrap_PCA_body(previous_step_cluster.cluster_spike_pairs, number_of_bootstrap_iterations, single_cluster_holder_options=self.single_cluster_holder_options)

				cluster = sw.Cluster()
				cluster.analysis_dicts = self.stage_1_cluster.analysis_dicts.copy()

				if (custom_stage_2_estimation):
					

					if (bootstrap_step < number_of_bootstrap_steps - 1):
						cluster.create_from_covariance_matrix(self.single_cluster_holder_options, 
							[bootstrap_pca_data_dict['BS_PCA_mean_of_mean0'], 
							bootstrap_pca_data_dict['BS_PCA_mean_of_mean1']], 
							bootstrap_pca_data_dict['BS_PCA_mean_covariance_matrix'], 
							self.spike_pair_distribution_spikes, 
							self.number_of_conjunctive_trials, 
							self.number_of_trials_for_condition, 
							self.real_trial_indices_on_which_both_neurons_spiked,
							make_differencing_calculations=False)

					else:
						cluster.create_from_final_spikes_and_covariance(self.single_cluster_holder_options, 
						previous_step_cluster.cluster_spike_pairs, 
						previous_step_cluster.cluster_indices, 
						[bootstrap_pca_data_dict['BS_PCA_mean_of_mean0'], 
						bootstrap_pca_data_dict['BS_PCA_mean_of_mean1']], 
						bootstrap_pca_data_dict['BS_PCA_mean_covariance_matrix'], 
						self.number_of_conjunctive_trials, 
						self.number_of_trials_for_condition, 
						self.real_trial_indices_on_which_both_neurons_spiked,
						make_differencing_calculations=True)



				else:
					cluster.create_from_final_spikes_and_covariance(self.single_cluster_holder_options, 
						self.stage_1_cluster.cluster_spike_pairs, 
						self.stage_1_cluster.cluster_indices, 
						[bootstrap_pca_data_dict['BS_PCA_mean_of_mean0'], 
						bootstrap_pca_data_dict['BS_PCA_mean_of_mean1']], 
						bootstrap_pca_data_dict['BS_PCA_mean_covariance_matrix'], 
						self.number_of_conjunctive_trials, 
						self.number_of_trials_for_condition, 
						self.real_trial_indices_on_which_both_neurons_spiked,
						make_differencing_calculations=(bootstrap_step == number_of_bootstrap_steps - 1))


				previous_step_cluster = cluster
				if (bootstrap_step == number_of_bootstrap_steps - 1):
					self.stage_2_cluster = cluster

		else:
			bootstrap_pca_data_dict = sw.bootstrap_PCA_body(spikes_to_use, number_of_bootstrap_iterations, single_cluster_holder_options=self.single_cluster_holder_options)


		suffixes = ["", "TestsPassed", "TestsPassedAndNormal", "TestsPassedActuallyDifferenced"]
		for suffix in suffixes:
			analysis_dict_key = iteration_key + suffix
			if ((analysis_dict_key in self.stage_2_cluster.analysis_dicts.keys()) and (not self.stage_2_cluster.analysis_dicts[analysis_dict_key]['is_empty'])):
				for key, value in bootstrap_pca_data_dict.items():
					self.stage_2_cluster.analysis_dicts[analysis_dict_key][key] = value

