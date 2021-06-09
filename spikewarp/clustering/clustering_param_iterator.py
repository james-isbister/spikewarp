import copy
import glob
import random
import math

import spikewarp as sw

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kde

from shapely.geometry import Point

from numpy.random import RandomState

"""
Class for iterating over clustering parameter range and creating surrogate angle control datasets
"""

class ClusteringParamIterator(object):

	def __init__(self, 
				experiment_code,
				random_seed_counter, 
				scsuo_0, 
				scsuo_1, 
				condition_index,
				analysis_recipe):

		self.experiment_code = experiment_code
		self.random_seed_counter = random_seed_counter
		self.scsuo_0 = scsuo_0
		self.scsuo_1 = scsuo_1
		self.condition_index = condition_index
		self.analysis_recipe = analysis_recipe
		self.all_value_combinations = copy.copy(self.analysis_recipe['option_helper'].all_value_combinations)
		self.channel_unit_suffix_0 = '_U' + str(self.scsuo_0.channel_unit_pair)
		self.channel_unit_suffix_1 = '_U' + str(self.scsuo_1.channel_unit_pair)
		self.experiment_condition_units_string = self.experiment_code + '_' + str(self.condition_index) + self.channel_unit_suffix_0 + self.channel_unit_suffix_1
		true_for_trials_both_spiking = np.logical_and(self.scsuo_0.true_if_spiked_atleast_once_on_trial, self.scsuo_1.true_if_spiked_atleast_once_on_trial)
		trial_indices_on_which_both_neurons_spiked = np.nonzero(true_for_trials_both_spiking)
		self.real_trial_indices_on_which_both_neurons_spiked = trial_indices_on_which_both_neurons_spiked[0]
		self.first_spike_pairs_original = np.stack((self.scsuo_0.first_spikes_for_condition_trials[trial_indices_on_which_both_neurons_spiked], self.scsuo_1.first_spikes_for_condition_trials[trial_indices_on_which_both_neurons_spiked]), axis=1)
		self.single_clusterings_by_shuffle_option = {}
		self.tw_files_by_stwcboo_by_shuffle_option = {}
		for shuffle_option in self.analysis_recipe['option_helper'].shuffle_options:
			shuffle_option_string = shuffle_option[0] + str(shuffle_option[1])
			self.single_clusterings_by_shuffle_option[shuffle_option_string] = []
			self.tw_files_by_stwcboo_by_shuffle_option[shuffle_option_string] = [[] for stwcboo in self.analysis_recipe['super_time_warp_cluster_buster_option_objects']]

		self.main_processing_loop()
		sw.make_pairwise_distribution_analysis_plots(self.analysis_recipe, 
													self.first_spike_pairs_original, 
													self.experiment_condition_units_string,
													self.single_clusterings_by_shuffle_option, 
													self.condition_index, 
													self.experiment_code, 
													self.tw_files_by_stwcboo_by_shuffle_option)

		

	



	def main_processing_loop(self):

		zip_objects = []
		for value_combinations_for_shuffle_type in self.analysis_recipe['option_helper'].all_value_combinations_as_shuffle_type_organised_list:
			number_of_value_combinations_for_shuffle_type = len(value_combinations_for_shuffle_type)

			spikes_to_use_for_each_shuffle_type_value_combination = [self.first_spike_pairs_original for combination_index in range(number_of_value_combinations_for_shuffle_type)]
			still_need_to_be_rotated_for_shuffle_type_each_combination = [True for combination_index in range(number_of_value_combinations_for_shuffle_type)]
			original_clustering_param_0s_for_each_shuffle_type_rotated_combination = [-1.0 for combination_index in range(number_of_value_combinations_for_shuffle_type)]
			size_limit_of_new_cluster_for_shuffle_type_each_combination = [100000 for combination_index in range(number_of_value_combinations_for_shuffle_type)]

			zip_object = list(zip(value_combinations_for_shuffle_type, spikes_to_use_for_each_shuffle_type_value_combination, still_need_to_be_rotated_for_shuffle_type_each_combination, original_clustering_param_0s_for_each_shuffle_type_rotated_combination, size_limit_of_new_cluster_for_shuffle_type_each_combination))
			zip_objects.append(zip_object)


		for zip_object in zip_objects:

			new_zip_object = []
			previous_single_clusterings_for_shuffle_type = []
			total_angled_cluster_count = 0

			number_of_value_combinations = len(zip_object)
			for value_comb_index, (value_combination, spikes_to_use, still_need_to_be_rotated, original_clustering_param_0s_for_rotated_combination, size_limit_of_new_cluster) in enumerate(zip_object):

				# 1. Do clustering for clustering parameter
				single_clustering = sw.SingleClustering([self.scsuo_0, self.scsuo_1], 
														spikes_to_use, 
														self.analysis_recipe['clustering_type'],
														self.analysis_recipe['custom_checks'],
														value_combination,
														self.experiment_condition_units_string,
														self.random_seed_counter,
														original_clustering_param_0s_for_rotated_combination,
														size_limit_of_new_cluster,
														self.condition_index,
														self.real_trial_indices_on_which_both_neurons_spiked)

				single_clustering.initial_statistics_and_clustering()
				

				# 2. Create Stage 1 clusters and optionally pop intersecting
				single_clustering.create_stage_1_clusters_and_optionally_pop_intersecting(pop_intersecting=self.analysis_recipe['custom_checks'])

				# 3. Check clustering not the same as previous param value
				if ((previous_single_clusterings_for_shuffle_type == []) or ((previous_single_clusterings_for_shuffle_type != []) and (np.array_equal(single_clustering.cluster_indices_for_non_outliers, previous_single_clusterings_for_shuffle_type[-1].cluster_indices_for_non_outliers) == False))):
					single_clustering.do_use_clusters_in_analysis = True

					# 4. prepare Stage 1 and Stage 2 data dictionaries, create Stage 2 cluster and do time series and factor analysis
					single_clustering.prepare_data_dictionaries_for_stage_1_and_2_and_calculate_stage_2_clusters_and_time_series_and_factor_analysis(total_angled_cluster_count, previous_single_clusterings_for_shuffle_type, custom_stage_2_estimation=self.analysis_recipe['custom_stage_2_estimation'], custom_stage_2_repeat_and_shared_pop=self.analysis_recipe['custom_stage_2_repeat_and_shared_pop'])
					total_angled_cluster_count += single_clustering.final_angled_cluster_count

					
				previous_single_clusterings_for_shuffle_type.append(single_clustering)

				
				if ((single_clustering.shuffle_option[0] not in ["sample_correlated_cluster_pca_ellipse_rotated_to_45", "45cluster_with_second_cluster", "rotate_correlatd_clusters_to_45_degrees", "shuffle_correlated_clusters", "flatten_and_sample_correlated_clusters"]) or (still_need_to_be_rotated == False)):

					if (single_clustering.do_use_clusters_in_analysis):
						self.single_clusterings_by_shuffle_option[single_clustering.shuffle_option_string].append(single_clustering)
					
				elif (single_clustering.do_use_clusters_in_analysis):
					new_zip_object = self.create_modified_control_distributions(single_clustering, spikes_to_use, new_zip_object)

				if (value_comb_index == number_of_value_combinations - 1):
					zip_objects.append(new_zip_object)

							 


	def create_45_cluster(self, single_cluster_holder):

		cov_matrix = np.copy(single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['BS_PCA_mean_covariance_matrix'])
		mean = [single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['BS_PCA_mean_of_mean0'], single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['BS_PCA_mean_of_mean1']]
		average_sd = (single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['BS_PCA_mean_component_0_sd'] + single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['BS_PCA_mean_component_1_sd']) / 2.0
		cov_matrix[0, 0] = average_sd**2; cov_matrix[1, 1] = average_sd**2; cov_matrix[0, 1] = abs(cov_matrix[0, 1]); cov_matrix[1, 0] = abs(cov_matrix[1, 0])

		return mean, cov_matrix



	def create_modified_control_distributions(self, single_clustering, spikes_to_use, new_zip_object):

		spikes_to_use_for_rotation = []

		if (single_clustering.shuffle_option[0] in ["sample_correlated_cluster_pca_ellipse_rotated_to_45", "flatten_and_sample_correlated_clusters", "45cluster_with_second_cluster"]):
			np.random.seed(self.random_seed_counter)
			prng = RandomState(self.random_seed_counter)

		size_limit_of_new_cluster = 10000

		for single_cluster_holder in single_clustering.secondary_single_cluster_holders:

			if (single_cluster_holder.final_stage_2_cluster):

				if (spikes_to_use_for_rotation == []):
					spikes_to_use_for_rotation = np.copy(spikes_to_use)

				if (single_clustering.shuffle_option[0] == "sample_correlated_cluster_pca_ellipse_rotated_to_45"):


					mean, cov_matrix = self.create_45_cluster(single_cluster_holder)

					flat_cov_matrix = cov_matrix.copy()
					flat_cov_matrix[0, 1] = 0.0; flat_cov_matrix[1, 0] = 0.0;
										
					indices_within_flat_delete_bound = sw.indices_within_x_mahalanobis_distances(mean, flat_cov_matrix, spikes_to_use, 6.0)
					indices_within_delete_bound = indices_within_flat_delete_bound
					
					number_of_indices_within_delete_bound = len(indices_within_delete_bound)

					if (sw.is_pos_def(cov_matrix)):
						sampled_from_cov_matrix = prng.multivariate_normal(mean, cov_matrix, [number_of_indices_within_delete_bound])
						spikes_to_use_for_rotation[indices_within_delete_bound, :] = sampled_from_cov_matrix
					else:
						spikes_to_use_for_rotation = []



				elif (single_clustering.shuffle_option[0] == "45cluster_with_second_cluster"):
					
					spikes_to_use_for_rotation = np.zeros((single_cluster_holder.stage_2_cluster.number_of_samples_in_ellipse * 2, 2))

					mean, cov_matrix = self.create_45_cluster(single_cluster_holder)
					if (sw.is_pos_def(cov_matrix)):
						sampled_from_cov_matrix = prng.multivariate_normal(mean, cov_matrix, [single_cluster_holder.stage_2_cluster.number_of_samples_in_ellipse])
						spikes_to_use_for_rotation[:single_cluster_holder.stage_2_cluster.number_of_samples_in_ellipse, :] = sampled_from_cov_matrix

						a0 = np.random.normal(single_cluster_holder.stage_2_cluster.analysis_dicts['FlatClusterStats']['FlatCluster_FS_Mean0'] + 10.0 * math.sqrt(cov_matrix[0,0]), single_cluster_holder.stage_2_cluster.analysis_dicts['FlatClusterStats']['FlatCluster_N0_FS_SD'], [single_cluster_holder.stage_2_cluster.number_of_samples_in_ellipse])
						a1 = np.random.normal(single_cluster_holder.stage_2_cluster.analysis_dicts['FlatClusterStats']['FlatCluster_FS_Mean1'] + 10.0 * math.tan(math.radians(15)) * math.sqrt(cov_matrix[0,0]), single_cluster_holder.stage_2_cluster.analysis_dicts['FlatClusterStats']['FlatCluster_N1_FS_SD'], [single_cluster_holder.stage_2_cluster.number_of_samples_in_ellipse])
						spikes_to_use_for_rotation[single_cluster_holder.stage_2_cluster.number_of_samples_in_ellipse:, :] = np.stack((a0, a1), axis=1)
					else:
						spikes_to_use_for_rotation = []	


				elif (single_clustering.shuffle_option[0] == "rotate_correlatd_clusters_to_45_degrees"):

					if (single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['is_PCA_BS_empirical_pvalue_different_from_45']):

						cluster_mean_point = np.mean(spikes_to_use_for_rotation[single_cluster_holder.stage_1_cluster.cluster_indices, :], axis=0)

						spikes_to_use_for_rotation[single_cluster_holder.stage_1_cluster.cluster_indices, :] = spikes_to_use_for_rotation[single_cluster_holder.stage_1_cluster.cluster_indices, :] - cluster_mean_point

						theta = np.radians(single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['BS_PCA_mean_angle_difference_from_45'])
						c, s = np.cos(theta), np.sin(theta)
						R = np.array(((c,-s), (s, c)))
						spikes_to_use_for_rotation[single_cluster_holder.stage_1_cluster.cluster_indices, :] = np.matmul(R, spikes_to_use_for_rotation[single_cluster_holder.stage_1_cluster.cluster_indices, :].T).T

						spikes_to_use_for_rotation[single_cluster_holder.stage_1_cluster.cluster_indices, :] = spikes_to_use_for_rotation[single_cluster_holder.stage_1_cluster.cluster_indices, :] + cluster_mean_point

						size_limit_of_new_cluster = single_cluster_holder.stage_1_cluster.number_of_samples_in_ellipse + 5

					else:
						spikes_to_use_for_rotation = []


				elif (single_clustering.shuffle_option[0] == "shuffle_correlated_clusters"):

					spikes_to_be_shuffled = spikes_to_use_for_rotation[single_cluster_holder.stage_1_cluster.cluster_indices, 1]
					random.Random(self.random_seed_counter).shuffle(spikes_to_be_shuffled)
					spikes_to_use_for_rotation[single_cluster_holder.stage_1_cluster.cluster_indices, 1] = spikes_to_be_shuffled
					size_limit_of_new_cluster = single_cluster_holder.stage_1_cluster.number_of_samples_in_ellipse + 5
					

				elif (single_clustering.shuffle_option[0] == "flatten_and_sample_correlated_clusters"):

					a0 = np.random.normal(single_cluster_holder.stage_2_cluster.analysis_dicts['FlatClusterStats']['FlatCluster_FS_Mean0'], single_cluster_holder.stage_2_cluster.analysis_dicts['FlatClusterStats']['FlatCluster_N0_FS_SD'], [single_cluster_holder.stage_1_cluster.number_of_samples_in_ellipse])
					a1 = np.random.normal(single_cluster_holder.stage_2_cluster.analysis_dicts['FlatClusterStats']['FlatCluster_FS_Mean1'], single_cluster_holder.stage_2_cluster.analysis_dicts['FlatClusterStats']['FlatCluster_N1_FS_SD'], [single_cluster_holder.stage_1_cluster.number_of_samples_in_ellipse])
					spikes_to_use_for_rotation[single_cluster_holder.stage_1_cluster.cluster_indices, :] = np.stack((a0, a1), axis=1)
					size_limit_of_new_cluster = single_cluster_holder.stage_1_cluster.number_of_samples_in_ellipse + 5

		if (spikes_to_use_for_rotation != []):

			temp_option_helper = sw.OptionHelper([single_clustering.shuffle_option], self.analysis_recipe['option_helper'].info_types, self.analysis_recipe['option_helper'].clustering_param_0s)

			self.all_value_combinations.extend(temp_option_helper.all_value_combinations)
			rotated_spikes_to_use_for_each_value_combination = [spikes_to_use_for_rotation for combination_index in range(temp_option_helper.total_number_of_combinations)]
			rotated_still_need_to_be_rotated_for_each_value_combination = [False for combination_index in range(temp_option_helper.total_number_of_combinations)]
			rotated_original_clustering_param_0s_for_each_rotated_combination = [single_clustering.clustering_param_0 for combination_index in range(temp_option_helper.total_number_of_combinations)]
			size_limit_of_new_cluster_for_each_rotated_combination = [size_limit_of_new_cluster for combination_index in range(temp_option_helper.total_number_of_combinations)]
			
			for temp_combination_index in range(temp_option_helper.total_number_of_combinations):
				new_zip_object.append((temp_option_helper.all_value_combinations[temp_combination_index], rotated_spikes_to_use_for_each_value_combination[temp_combination_index], rotated_still_need_to_be_rotated_for_each_value_combination[temp_combination_index], rotated_original_clustering_param_0s_for_each_rotated_combination[temp_combination_index], size_limit_of_new_cluster_for_each_rotated_combination[temp_combination_index]))


		return new_zip_object
		
