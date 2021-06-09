import numpy as np
import random



class SingleConditionSingleUnitObject(object):

	"""
	Class for storing and processing data of single condition single neuron combinations
	"""

	def __init__(self):

		self.spike_trains_for_condition_trials = None
		self.session_wise_lfp_principal_components_for_condition_trials = None
		self.condition_wise_lfp_principal_components_for_condition_trials = None
		self.condition_index = -1
		self.channel_unit_pair = None
		self.mean_of_first_spikes = -1.0
		self.standard_deviation_of_first_spikes = -1.0
		self.channel_unit_pair_index = -1
		self.true_if_spiked_atleast_once_on_trial = None
		self.stimulation_frequency_of_condition = 1
		self.lr_neighbouring_indices_vs_fs = None
		self.lr_actual_trial_indices_vs_fs = None
		self.lr_autocorrelation_lag_1 = None
		self.adfuller_dataframe = None

		self.neurophys_corresponding_layer_string_of_unit = "Temp"
		self.exc_inh_type_of_unit = "Temp"
		self.principle_condition_unit_combination = True



	def finalise(self, number_of_trials_for_condition=0):

		self.number_of_trials_for_condition = number_of_trials_for_condition
		

		if (self.session_wise_lfp_principal_components_for_condition_trials != None):
			if (self.session_wise_lfp_principal_components_for_condition_trials.any()):
				self.number_of_trials_for_condition = self.session_wise_lfp_principal_components_for_condition_trials.shape[0]
		if (self.condition_wise_lfp_principal_components_for_condition_trials != None):
			if (self.condition_wise_lfp_principal_components_for_condition_trials.any()):
				self.number_of_trials_for_condition = self.condition_wise_lfp_principal_components_for_condition_trials.shape[0]



		self.true_if_spiked_atleast_once_on_trial = np.zeros(self.number_of_trials_for_condition, dtype=bool)
		if (self.number_of_trials_for_condition > 0):

			self.were_spikes_on_all_trials = True
			self.first_spikes_for_condition_trials = np.zeros((self.number_of_trials_for_condition, 1))

			for sub_list_index in range(0, self.number_of_trials_for_condition):

				sub_list = self.spike_trains_for_condition_trials[sub_list_index]
				if sub_list == []:
					self.were_spikes_on_all_trials = False
				else:
					self.first_spikes_for_condition_trials[sub_list_index] = sub_list[0]



			self.true_if_spiked_atleast_once_on_trial = self.first_spikes_for_condition_trials > 0
			self.number_of_trials_with_atleast_one_spike = np.sum(self.true_if_spiked_atleast_once_on_trial)
			self.atleast_one_spike_reliability = float(self.number_of_trials_with_atleast_one_spike) / float(self.number_of_trials_for_condition)
			self.indices_of_trials_with_atleast_one_spike = np.argwhere(self.true_if_spiked_atleast_once_on_trial)
			self.indices_of_trials_with_atleast_one_spike = self.indices_of_trials_with_atleast_one_spike[:, 0]
			self.first_spikes_for_trials_which_spiked_atleast_once_on_trial = self.first_spikes_for_condition_trials[self.indices_of_trials_with_atleast_one_spike].flatten()
			if (self.number_of_trials_with_atleast_one_spike > 1):
				self.standard_deviation_of_first_spikes = np.std(self.first_spikes_for_trials_which_spiked_atleast_once_on_trial)
				self.mean_of_first_spikes = np.mean(self.first_spikes_for_trials_which_spiked_atleast_once_on_trial)
				if (self.number_of_trials_with_atleast_one_spike > 2):
					self.first_spikes_for_condition_trials_normalised = (self.first_spikes_for_condition_trials - self.mean_of_first_spikes) / self.standard_deviation_of_first_spikes
			if (self.condition_wise_lfp_principal_components_for_condition_trials != None):
				temp_bool_array = np.repeat(self.true_if_spiked_atleast_once_on_trial, self.condition_wise_lfp_principal_components_for_condition_trials.shape[1], axis=1)
				self.condition_wise_lfp_principal_components_for_trials_which_spiked_atleast_once_on_trial = self.condition_wise_lfp_principal_components_for_condition_trials[self.indices_of_trials_with_atleast_one_spike, :]


			# Spike counts 
			self.spike_counts_per_trial = np.zeros((self.number_of_trials_for_condition, 1))
			for trial_index in range(0, self.number_of_trials_for_condition):
				number_of_spikes = len(self.spike_trains_for_condition_trials[trial_index])
				self.spike_counts_per_trial[trial_index] = number_of_spikes


			self.fano_factor = -1
			if (self.number_of_trials_for_condition > 0):
				self.mean_trial_spike_count_for_condition = np.mean(self.spike_counts_per_trial)
				self.variance_trial_spike_count_for_condition = np.var(self.spike_counts_per_trial)
				if (self.mean_trial_spike_count_for_condition > 0):
					self.fano_factor = self.variance_trial_spike_count_for_condition / self.mean_trial_spike_count_for_condition

					if (self.number_of_trials_with_atleast_one_spike > 2):
						self.spike_counts_per_trial_normalised = (self.spike_counts_per_trial - self.mean_trial_spike_count_for_condition) / np.std(self.spike_counts_per_trial)

				spike_counts_per_trial_shuffled = np.copy(self.spike_counts_per_trial)
				random.Random().shuffle(spike_counts_per_trial_shuffled)

				added_spike_counts_per_trial = self.spike_counts_per_trial + spike_counts_per_trial_shuffled

				self.added_mean_trial_spike_count_for_condition = np.mean(added_spike_counts_per_trial)
				self.added_variance_trial_spike_count_for_condition = np.var(added_spike_counts_per_trial)
				if (self.added_mean_trial_spike_count_for_condition > 0):
					self.added_fano_factor = self.added_variance_trial_spike_count_for_condition / self.added_mean_trial_spike_count_for_condition

