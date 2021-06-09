import numpy as np
import random
from scipy.stats import linregress

import warnings

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from sklearn.cluster import DBSCAN
from sklearn.mixture import GaussianMixture

import spikewarp as sw




"""
Class for managing single clustering for single clustering parameter combination and generating correlation control datasets
"""
class SingleClustering(object):

	def __init__(self, 
				scsuos, 
				first_spike_pairs_original, 
				clustering_type,
				custom_checks,
				value_combination,
				experiment_condition_units_string,
				random_seed_counter,
				original_clustering_param_0s_for_rotated_combination,
				size_limit_of_new_cluster,
				condition_index,
				real_trial_indices_on_which_both_neurons_spiked):


		self.first_spike_pairs_original = first_spike_pairs_original
		self.number_of_conjunctive_trials = self.first_spike_pairs_original.shape[0]
		self.clustering_type = clustering_type
		self.custom_checks = custom_checks
		self.value_combination = value_combination
		self.experiment_condition_units_string = experiment_condition_units_string
		self.random_seed_counter = random_seed_counter
		self.size_limit_of_new_cluster = size_limit_of_new_cluster
		self.condition_index = condition_index

		self.shuffle_option = self.value_combination[0]; self.info_type = self.value_combination[1]; self.clustering_param_0 = self.value_combination[2]	
		self.shuffle_option_string = self.shuffle_option[0] + str(self.shuffle_option[1])
		self.option_string = self.shuffle_option[0] + str(self.shuffle_option[1]) + "_Clustering" + str(self.clustering_param_0)
		if (original_clustering_param_0s_for_rotated_combination != -1.0):
			self.option_string = self.option_string + "_Original" + str(original_clustering_param_0s_for_rotated_combination)
		

		self.real_trial_indices_on_which_both_neurons_spiked = []
		if (self.shuffle_option[0] == "normal"):
			self.real_trial_indices_on_which_both_neurons_spiked = real_trial_indices_on_which_both_neurons_spiked

		self.do_use_clusters_in_analysis = False

		self.experiment_condition_units_string_no_dots = '%s' % self.experiment_condition_units_string
		self.experiment_condition_units_string_no_dots.replace('.', 'p')
		self.option_string_with_experiment_condition_units_string = self.option_string + "_" + self.experiment_condition_units_string


		self.number_of_trials_for_condition = 0
		self.stimulation_frequency_of_condition = 1
		self.neurophys_layer_strings_of_unit_pair = ['Temp', 'Temp']
		self.exc_inh_types_of_unit_pair = ['Temp', 'Temp']
		self.principle_condition_bool_for_each_unit = [True, True]

		self.non_singleton_single_cluster_holders = []
		self.non_intersecting_non_singleton_single_cluster_holders = []
		self.primary_single_cluster_holders = []
		self.secondary_single_cluster_holders = []

		if (scsuos != []):
			self.scsuo_0, self.scsuo_1 = scsuos

			self.number_of_trials_for_condition = self.scsuo_0.number_of_trials_for_condition
			self.stimulation_frequency_of_condition = self.scsuo_0.stimulation_frequency_of_condition
			self.neurophys_layer_strings_of_unit_pair = [self.scsuo_0.neurophys_corresponding_layer_string_of_unit, self.scsuo_1.neurophys_corresponding_layer_string_of_unit]
			self.exc_inh_types_of_unit_pair = [self.scsuo_0.exc_inh_type_of_unit, self.scsuo_1.exc_inh_type_of_unit]
			self.principle_condition_bool_for_each_unit = [self.scsuo_0.principle_condition_unit_combination, self.scsuo_1.principle_condition_unit_combination]
			self.barrels_of_unit_pair = [self.scsuo_0.barrel_of_unit, self.scsuo_1.barrel_of_unit]

		self.was_first_single_clustering_to_pass_for_pair = False


	def initial_statistics_and_clustering(self):

		self.calculate_initial_stats()
		self.specific_spike_shuffling_or_sampling_and_linear_regression()
		self.initial_clustering()



	def calculate_initial_stats(self):
		
		self.N0_not_only_conjunctive_original_FS_SD = self.scsuo_0.standard_deviation_of_first_spikes
		self.N1_not_only_conjunctive_original_FS_SD = self.scsuo_1.standard_deviation_of_first_spikes
		self.conjunctive_first_spike_0s_standard_deviation = np.std(self.first_spike_pairs_original[:,0])
		self.conjunctive_first_spike_1s_standard_deviation = np.std(self.first_spike_pairs_original[:,1])
		self.conjunctive_first_spike_differences_mean = np.mean(self.first_spike_pairs_original[:, 1] - self.first_spike_pairs_original[:,0])
		self.conjunctive_first_spike_differences_standard_deviation = np.std(self.first_spike_pairs_original[:, 1] - self.first_spike_pairs_original[:,0])
		self.conjunctive_first_spike_0s_mean = np.mean(self.first_spike_pairs_original[:,0])
		self.conjunctive_first_spike_1s_mean = np.mean(self.first_spike_pairs_original[:,1])
		self.conjunctive_range_0 = [np.min(self.first_spike_pairs_original[:,0]), np.max(self.first_spike_pairs_original[:,0])]
		self.conjunctive_range_1 = [np.min(self.first_spike_pairs_original[:,1]), np.max(self.first_spike_pairs_original[:,1])]



	def specific_spike_shuffling_or_sampling_and_linear_regression(self):

		first_spike_pairs_for_this_option_combination = np.copy(self.first_spike_pairs_original)
		np.random.seed(self.random_seed_counter)		
		number_of_pairs_to_generate = self.number_of_conjunctive_trials
		if (self.shuffle_option[1] != -1):
			number_of_pairs_to_generate = self.shuffle_option[1]

		if (self.shuffle_option[0] == 'shuffle'):
			random.Random(self.random_seed_counter).shuffle(first_spike_pairs_for_this_option_combination[:, 1])

		if (self.shuffle_option[0] == "sample"):
			a0 = np.random.normal(self.conjunctive_first_spike_0s_mean, self.conjunctive_first_spike_0s_standard_deviation, [number_of_pairs_to_generate])
			a1 = np.random.normal(self.conjunctive_first_spike_1s_mean, self.conjunctive_first_spike_1s_standard_deviation, [number_of_pairs_to_generate])
			first_spike_pairs_for_this_option_combination = np.stack((a0, a1), axis=1)


		self.list_of_outlier_first_spike_pairs = []
		self.spike_pair_distribution_spikes = first_spike_pairs_for_this_option_combination

		if (len(self.spike_pair_distribution_spikes) > 2):
			self.LR_all_conjunctive = linregress(self.spike_pair_distribution_spikes)

		else:
			self.LR_all_conjunctive = None


	

	def initial_clustering(self):

		did_cluster = False

		if (self.clustering_type == "Unclustered"):
			self.cluster_indices_for_non_outliers = np.ones((self.number_of_conjunctive_trials, 1), dtype=int)
			did_cluster = True

		if (self.clustering_type in ["GaussianMixtureInfCrit", "CustomGaussianMixtureInfCrit"]):

			if (self.number_of_conjunctive_trials > self.clustering_param_0):

				gmms = [GaussianMixture(n, covariance_type='full').fit(self.spike_pair_distribution_spikes) for n in range(1, self.clustering_param_0)]
				bics = [gmm.bic(self.spike_pair_distribution_spikes) for gmm in gmms]
				gmm = gmms[np.argmin(bics)]
				self.cluster_indices_for_non_outliers = gmm.predict(self.spike_pair_distribution_spikes)
				did_cluster = True

		if (self.clustering_type == "CustomDBSCAN"):
			db = DBSCAN(eps=self.clustering_param_0, min_samples=2).fit(self.spike_pair_distribution_spikes)
			num_clus = np.sum((np.unique(db.labels_) >= 0))

			self.cluster_indices_for_non_outliers = np.zeros((self.number_of_conjunctive_trials, 1), dtype=int)
			singleton_clus_index = num_clus
			for both_spiking_trial_index in range(self.number_of_conjunctive_trials):
				if (db.labels_[both_spiking_trial_index] == -1):
					self.cluster_indices_for_non_outliers[both_spiking_trial_index] = singleton_clus_index
					singleton_clus_index = singleton_clus_index + 1
				else:
					self.cluster_indices_for_non_outliers[both_spiking_trial_index] = db.labels_[both_spiking_trial_index]

			self.cluster_indices_for_non_outliers = self.cluster_indices_for_non_outliers + 1
			did_cluster = True

		if (did_cluster):
			self.unique_cluster_indices = np.unique(self.cluster_indices_for_non_outliers)
			self.number_of_original_clusters = len(self.unique_cluster_indices)
		else:
			self.number_of_original_clusters = 0



		
		# POST PROCESSING OF ORIGINAL CLUSTERS
		self.indices_for_each_original_cluster = []
		self.first_spike_pairs_for_each_original_cluster = []
		self.number_of_samples_for_each_original_cluster = []

		cmap = matplotlib.cm.get_cmap('magma')
		colors = ['g', 'b', 'r', 'c', 'y', 'm']
		self.rgbas_for_non_singleton_original_clusters  = []
		non_singleton_index = -1
		for cluster_index in self.unique_cluster_indices:

			original_cluster_indices = np.where(self.cluster_indices_for_non_outliers == cluster_index)[0]
			number_of_samples_in_original_cluster = len(original_cluster_indices)
			first_spike_pairs_for_original_cluster = self.spike_pair_distribution_spikes[original_cluster_indices]

			self.indices_for_each_original_cluster.append(original_cluster_indices)
			self.number_of_samples_for_each_original_cluster.append(number_of_samples_in_original_cluster)
			self.first_spike_pairs_for_each_original_cluster.append(first_spike_pairs_for_original_cluster)

			if (number_of_samples_in_original_cluster > 1):
				non_singleton_index = non_singleton_index + 1

			rgba = colors[non_singleton_index % 6]
			self.rgbas_for_non_singleton_original_clusters.append(matplotlib.colors.to_hex(rgba))

		self.number_of_non_singleton_original_clusters = non_singleton_index + 1



	def create_stage_1_clusters_and_optionally_pop_intersecting(self, pop_intersecting=False):

		non_singleton_index = -1
		
		for original_cluster_index in range(self.number_of_original_clusters):

			number_of_samples_in_original_cluster = self.number_of_samples_for_each_original_cluster[original_cluster_index]
			original_cluster_indices = self.indices_for_each_original_cluster[original_cluster_index]
			first_spike_pairs_for_original_cluster = self.first_spike_pairs_for_each_original_cluster[original_cluster_index]

			if (number_of_samples_in_original_cluster > 1):

				non_singleton_index = non_singleton_index + 1

				single_cluster_holder_options = sw.SingleClusterHolderOptions()				
				single_cluster_holder = sw.SingleClusterHolder(single_cluster_holder_options, 
																original_cluster_indices, 
																first_spike_pairs_for_original_cluster, 
																self.spike_pair_distribution_spikes, 
																self.number_of_conjunctive_trials, 
																self.number_of_trials_for_condition,
																non_singleton_index,
																self.scsuo_0.channel_unit_pair_index,
																self.scsuo_1.channel_unit_pair_index,
																self.condition_index,
																self.real_trial_indices_on_which_both_neurons_spiked,
																self.stimulation_frequency_of_condition,
																self.neurophys_layer_strings_of_unit_pair,
																self.exc_inh_types_of_unit_pair,
																self.barrels_of_unit_pair)


				
				single_cluster_holder.create_intersection_ellipse_and_stage_1_cluster(self.size_limit_of_new_cluster, with_custom_stage_1_cluster_expansion=self.custom_checks)

				self.non_singleton_single_cluster_holders.append(single_cluster_holder)
				self.non_intersecting_non_singleton_single_cluster_holders.append(single_cluster_holder)

		if (pop_intersecting):
			self.pop_single_cluster_holders0_which_intersect_with_non_singleton_single_cluster_holders(self.non_intersecting_non_singleton_single_cluster_holders, self.non_singleton_single_cluster_holders, option="flat_intersection_ellipse")



	def prepare_data_dictionaries_for_stage_1_and_2_and_calculate_stage_2_clusters_and_time_series_and_factor_analysis(self, total_angled_cluster_count, previous_single_clusterings, custom_stage_2_estimation=False, custom_stage_2_repeat_and_shared_pop=False):

		for non_intersecting_non_singleton_single_cluster_holder in self.non_intersecting_non_singleton_single_cluster_holders:
			for previous_single_clustering in previous_single_clusterings:
				for previous_single_clustering_primary_single_cluster_holder in previous_single_clustering.primary_single_cluster_holders:
					if ((non_intersecting_non_singleton_single_cluster_holder.stage_1_cluster != None) & (previous_single_clustering_primary_single_cluster_holder.stage_1_cluster != None)):
						if (np.array_equal(non_intersecting_non_singleton_single_cluster_holder.stage_1_cluster.cluster_indices, previous_single_clustering_primary_single_cluster_holder.stage_1_cluster.cluster_indices)):
							non_intersecting_non_singleton_single_cluster_holder.final_stage_1_cluster = False
							non_intersecting_non_singleton_single_cluster_holder.final_stage_2_cluster = False


		
		for non_intersecting_non_singleton_single_cluster_holder in self.non_intersecting_non_singleton_single_cluster_holders:
			if (non_intersecting_non_singleton_single_cluster_holder.final_stage_1_cluster):
				self.primary_single_cluster_holders.append(non_intersecting_non_singleton_single_cluster_holder)
			if (non_intersecting_non_singleton_single_cluster_holder.final_stage_2_cluster):
				self.secondary_single_cluster_holders.append(non_intersecting_non_singleton_single_cluster_holder)

		self.primary_data_dicts = {}
		self.secondary_data_dicts = {}

		for data_name in sw.list_of_first_stage_data_names:
			self.primary_data_dicts.update({data_name: []})
			self.secondary_data_dicts.update({data_name: []})
		
		for data_name in sw.list_of_second_stage_data_names:
			self.secondary_data_dicts.update({data_name: []})

		for single_cluster_holder in self.primary_single_cluster_holders:
			if (single_cluster_holder.final_stage_1_cluster):
				sw.add_single_cluster_holder_data_to_dict(self.primary_data_dicts, sw.list_of_first_stage_data_names, single_cluster_holder, self, "flat")


		for single_cluster_holder in self.secondary_single_cluster_holders:
			if (single_cluster_holder.final_stage_2_cluster):
				single_cluster_holder.calculate_stage_2_cluster_with_time_series_and_factor_analysis(custom_stage_2_estimation=custom_stage_2_estimation)
		
		# if (custom_stage_2):
		# 	self.pop_single_cluster_holders0_which_intersect_with_non_singleton_single_cluster_holders(self.secondary_single_cluster_holders, self.non_singleton_single_cluster_holders, option="PCA_BS_intersection_ellipse")

		number_of_secondary_single_cluster_holder_to_use = 0

		for single_cluster_holder in self.secondary_single_cluster_holders:

			for previous_single_clustering in previous_single_clusterings:

				if (custom_stage_2_repeat_and_shared_pop == True):

					if (single_cluster_holder.final_stage_2_cluster == False):
						break
					else:
						# Check if similar to any previous clusterings
						for previous_single_clustering_secondary_single_cluster_holder in previous_single_clustering.secondary_single_cluster_holders:

							if (previous_single_clustering_secondary_single_cluster_holder.final_stage_2_cluster):

								if ((single_cluster_holder.stage_2_cluster != None) & (previous_single_clustering_secondary_single_cluster_holder.stage_2_cluster != None)):

									max_number_of_samples_in_two_ellipses = max(single_cluster_holder.stage_2_cluster.number_of_samples_in_ellipse, previous_single_clustering_secondary_single_cluster_holder.stage_2_cluster.number_of_samples_in_ellipse)
									number_of_shared_samples = np.intersect1d(single_cluster_holder.stage_2_cluster.cluster_indices, previous_single_clustering_secondary_single_cluster_holder.stage_2_cluster.cluster_indices).shape[0]

									if (max_number_of_samples_in_two_ellipses == number_of_shared_samples):
										single_cluster_holder.final_stage_2_cluster = False

									if (number_of_shared_samples > 15):
										single_cluster_holder.final_stage_2_cluster = False


			if (single_cluster_holder.final_stage_2_cluster == True):
				number_of_secondary_single_cluster_holder_to_use = number_of_secondary_single_cluster_holder_to_use + 1
				if (total_angled_cluster_count == 0):
					self.was_first_single_clustering_to_pass_for_pair = True
					

			if (single_cluster_holder.final_stage_2_cluster):
				sw.add_single_cluster_holder_data_to_dict(self.secondary_data_dicts, sw.list_of_second_stage_data_names, single_cluster_holder, self, "angled")

		self.final_angled_cluster_count = number_of_secondary_single_cluster_holder_to_use


	def pop_single_cluster_holders0_which_intersect_with_non_singleton_single_cluster_holders(self, single_cluster_holders0, non_singleton_single_cluster_holders, option="PCA_BS_intersection_ellipse"):

		list_to_pop = []
		for single_cluster_holder_0_index, single_cluster_holder_0 in enumerate(single_cluster_holders0):
			not_yet_popped = True
			for non_singleton_index, non_singleton_single_cluster_holder in enumerate(non_singleton_single_cluster_holders):

				if (single_cluster_holder_0.non_singleton_index != non_singleton_index):

					if (option == "flat_intersection_ellipse"):
						ellipses_intersect = single_cluster_holder_0.flat_intersection_ellipse.intersects(non_singleton_single_cluster_holder.flat_intersection_ellipse)
						ellipses_within = single_cluster_holder_0.flat_intersection_ellipse.within(non_singleton_single_cluster_holder.flat_intersection_ellipse)


					if (option == "PCA_BS_intersection_ellipse"):
						
						ellipses_intersect = single_cluster_holder_0.stage_2_cluster.analysis_dicts['Original']['PCA_BS_intersection_ellipse'].intersects(non_singleton_single_cluster_holder.flat_intersection_ellipse)
						ellipses_within = single_cluster_holder_0.stage_2_cluster.analysis_dicts['Original']['PCA_BS_intersection_ellipse'].within(non_singleton_single_cluster_holder.flat_intersection_ellipse)
						

					if (ellipses_intersect):

						single_cluster_holder_0.final_stage_2_cluster = False

						if (not_yet_popped):

							list_to_pop.append(single_cluster_holder_0_index)
							not_yet_popped = False


		list_to_pop.reverse()
		for index in list_to_pop:
			single_cluster_holders0.pop(index)

