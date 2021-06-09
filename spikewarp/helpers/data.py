import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from statsmodels.tsa.stattools import adfuller
import pandas as pd
from statsmodels.tsa.stattools import kpss
import warnings

import spikewarp as sw

"""
Helper class for loading spiking activity, condition info and neuron info from text files for single experiment
"""

class DataLoader(object):

	def __init__(self, original_data_root):
		
		self.original_data_root = original_data_root
		self.array_of_final_lfp_mean_values_for_EACH_condition = None


	def LOAD_SPIKES_FOR_ALL_UNITS_AND_CONDITIONS(self, spikes_directory_name='Spikes/', channel_unit_pairs_for_analysis_or_unit_list=None, channel_unit_pair_format=False):

		"""
		Method for loading spikes
		"""
		
		self.channel_unit_pairs_for_analysis_or_unit_list = channel_unit_pairs_for_analysis_or_unit_list
		self.channel_unit_pair_format = channel_unit_pair_format
		
		self.number_of_channel_unit_pairs = len(self.channel_unit_pairs_for_analysis_or_unit_list)
		
		self.all_channel_unit_spike_trains_for_all_trials = []
		self.all_cut_out_spike_trains = []
		self.all_spikes = []

		for channel_unit_or_unit_pair in self.channel_unit_pairs_for_analysis_or_unit_list:

			channel_unit_spike_trains_for_all_trials = []

			if (channel_unit_pair_format):
				spike_train_file = self.original_data_root + spikes_directory_name + "Channel" + str(channel_unit_or_unit_pair[0]) + "_ChannelUnit" + str(channel_unit_or_unit_pair[1]) + ".txt"
			else:
				spike_train_file = self.original_data_root + spikes_directory_name + "Unit" + str(channel_unit_or_unit_pair) + ".txt"

			with open(spike_train_file, 'r') as filehandle:  
				for trial_index, line in enumerate(filehandle):

					line = line.split() # to deal with blank 
					if line:
						line = [float(i) for i in line]
						channel_unit_spike_trains_for_all_trials.append(line)
						self.all_cut_out_spike_trains.append(line)
						self.all_spikes.extend(line)
					else:
						channel_unit_spike_trains_for_all_trials.append([])

			self.all_channel_unit_spike_trains_for_all_trials.append(channel_unit_spike_trains_for_all_trials)


	def LOAD_CONDITION_INFO(self, load_frequency_and_whisker_information=False, max_stimulation_frequency=1000.0):

		"""
		Method for loading condition info
		"""
		
		self.start_indices_for_each_condition_in_sorted_trial_arrays_ONE_INDEXED = np.loadtxt(self.original_data_root + "ConditionInformation/start_indices_for_each_condition_in_sorted_trial_arrays.txt", dtype='i', ndmin=1)
		self.number_of_trials_containing_mandatory_codes_for_each_sorted_condition = np.loadtxt(self.original_data_root + "ConditionInformation/number_of_trials_containing_mandatory_codes_for_each_sorted_condition.txt", dtype='i', ndmin=1)
		self.number_of_conditions = self.number_of_trials_containing_mandatory_codes_for_each_sorted_condition.shape[0]

		if (load_frequency_and_whisker_information):
			self.stimulated_whisker_for_each_condition = np.genfromtxt(self.original_data_root + "ConditionInformation/stimulated_whisker_for_each_condition.txt", dtype='str')
			self.stimulation_frequency_for_each_condition = np.loadtxt(self.original_data_root + "ConditionInformation/stimulation_frequency_for_each_condition.txt", dtype='f', ndmin=1)

		# self.condition_indices_to_process = condition_indices_to_process
		# if (self.condition_indices_to_process == None):
			# self.condition_indices_to_process = list(range(self.number_of_conditions))
		self.condition_indices_to_process = []
		for condition_index in range(self.number_of_conditions):
			if (self.stimulation_frequency_for_each_condition[condition_index] <= max_stimulation_frequency):
				self.condition_indices_to_process.append(condition_index)

	def LOAD_UNIT_INFO(self):

		"""
		Method for loading neuron info
		"""

		self.neurophys_layer_index_of_each_unit = np.loadtxt(self.original_data_root + "UNIT_LAYERS.txt", dtype='i', ndmin=1)
		self.barrel_of_each_unit = np.genfromtxt(self.original_data_root + "UNIT_BARRELS.txt", dtype='str')
		self.exc_inh_type_of_each_unit = np.genfromtxt(self.original_data_root + "UNIT_EXC_INH_CHARACTERISATIONS.txt", dtype='str')
		self.unique_barrels_some_removed = list(np.unique(self.barrel_of_each_unit))
		if ('S' in self.unique_barrels_some_removed):
			self.unique_barrels_some_removed.remove('S')
		if ('OUT' in self.unique_barrels_some_removed):
			self.unique_barrels_some_removed.remove('OUT')
		self.number_of_unique_barrels_some_removed = len(self.unique_barrels_some_removed)

		self.neurophys_layer_string_of_each_unit = []
		for neurophys_layer_index in self.neurophys_layer_index_of_each_unit:
			self.neurophys_layer_string_of_each_unit.append(sw.return_neurophys_layer_string_given_neurophys_layer_index(neurophys_layer_index))



	def create_scsuos(self, spike_train_start_time, spike_train_end_time, load_condition_wise_lfp=False, load_session_wise_lfp=False, condition_wise_lfp_principal_components_for_each_trial_directory="InputLFPs_TRIAL_MEANS_GroupingBy_CONDITION/ComponentsOfIndividualLFPsInGroupwisePCASpace/", custom_number_of_trials_for_condition=-1, do_single_unit_stationarity_tests=False):

		"""
		Method which creates SingleConditionSingleUnitObjects from data
		"""

		if (load_session_wise_lfp == True):
			session_wise_lfp_principal_components_for_each_trial = np.loadtxt(session_wise_lfp_principal_components_for_each_trial_file, dtype='f')

		self.list_of_all_scsuos = []
		for condition_index in self.condition_indices_to_process:

			start_trial_index_for_condition = self.start_indices_for_each_condition_in_sorted_trial_arrays_ONE_INDEXED[condition_index] - 1
			number_of_trials_for_condition = self.number_of_trials_containing_mandatory_codes_for_each_sorted_condition[condition_index]

			if (custom_number_of_trials_for_condition != -1):
				number_of_trials_for_condition = min(custom_number_of_trials_for_condition, number_of_trials_for_condition)

			if (number_of_trials_for_condition > 0):
				if (load_session_wise_lfp == True):
					session_wise_lfp_principal_components_for_condition_trials = session_wise_lfp_principal_components_for_each_trial[start_trial_index_for_condition:start_trial_index_for_condition + number_of_trials_for_condition, :]

				if (load_condition_wise_lfp == True):
					condition_wise_lfp_principal_components_for_each_trial_of_condition_file = condition_wise_lfp_principal_components_for_each_trial_directory + "Condition" + str(condition_index) + "_components_of_individual_lfps_in_groupwise_pca_space"
					condition_wise_lfp_principal_components_for_condition_trials = np.loadtxt(condition_wise_lfp_principal_components_for_each_trial_of_condition_file, dtype='f')


				for channel_unit_pair_index in range(self.number_of_channel_unit_pairs):

					scsuo = sw.SingleConditionSingleUnitObject()
					if (load_session_wise_lfp):
						scsuo.session_wise_lfp_principal_components_for_condition_trials = session_wise_lfp_principal_components_for_condition_trials
					if (load_condition_wise_lfp):
						scsuo.condition_wise_lfp_principal_components_for_condition_trials = condition_wise_lfp_principal_components_for_condition_trials
					if (self.array_of_final_lfp_mean_values_for_EACH_condition != None):
						scsuo.array_of_final_lfp_mean_values_for_each_trial_of_condition = self.array_of_final_lfp_mean_values_for_EACH_condition[condition_index]

					scsuo.condition_index = condition_index
					scsuo.channel_unit_pair = self.channel_unit_pairs_for_analysis_or_unit_list[channel_unit_pair_index]
					scsuo.channel_unit_pair_index = channel_unit_pair_index
					scsuo.stimulated_whisker_of_condition = self.stimulated_whisker_for_each_condition[condition_index]
					scsuo.stimulation_frequency_of_condition = self.stimulation_frequency_for_each_condition[condition_index]
					scsuo.cortical_layer_of_unit = self.neurophys_layer_index_of_each_unit[channel_unit_pair_index] # Should deprecate across code, as ambiguous
					scsuo.neurophys_layer_index_of_unit = self.neurophys_layer_index_of_each_unit[channel_unit_pair_index]
					scsuo.neurophys_corresponding_layer_string_of_unit = self.neurophys_layer_string_of_each_unit[channel_unit_pair_index]
					scsuo.barrel_of_unit = self.barrel_of_each_unit[channel_unit_pair_index]
					scsuo.exc_inh_type_of_unit = self.exc_inh_type_of_each_unit[channel_unit_pair_index]
					scsuo.principle_condition_unit_combination = False

					if (scsuo.barrel_of_unit == scsuo.stimulated_whisker_of_condition):
						scsuo.principle_condition_unit_combination = True

					scsuo.barrel_index_for_plot = -1
					if (scsuo.barrel_of_unit in self.unique_barrels_some_removed):
						scsuo.barrel_index_for_plot = self.unique_barrels_some_removed.index(scsuo.barrel_of_unit)

					if (scsuo.stimulated_whisker_of_condition in self.unique_barrels_some_removed):
						scsuo.stimulated_whisker_index_for_plot = self.unique_barrels_some_removed.index(scsuo.stimulated_whisker_of_condition)
					
					trimmed_trial_spike_trains_for_condition = []

					miss_first_spikes_buffer = 0

					for trial_index in range(start_trial_index_for_condition + miss_first_spikes_buffer, start_trial_index_for_condition + number_of_trials_for_condition):

						trial_spike_train = self.all_channel_unit_spike_trains_for_all_trials[channel_unit_pair_index][trial_index]
						trimmed_trial_spike_train = [spike for spike in trial_spike_train if (spike >= spike_train_start_time) and (spike < spike_train_end_time)]

						trimmed_trial_spike_trains_for_condition.append(trimmed_trial_spike_train)

					
					scsuo.spike_trains_for_condition_trials = trimmed_trial_spike_trains_for_condition
					scsuo.finalise(number_of_trials_for_condition - miss_first_spikes_buffer)

					if ((do_single_unit_stationarity_tests) & (scsuo.number_of_trials_with_atleast_one_spike > 20)):

						scsuo.lr_neighbouring_indices_vs_fs = linregress(list(range(scsuo.number_of_trials_with_atleast_one_spike)), scsuo.first_spikes_for_trials_which_spiked_atleast_once_on_trial)
						scsuo.lr_actual_trial_indices_vs_fs = linregress(scsuo.indices_of_trials_with_atleast_one_spike, scsuo.first_spikes_for_trials_which_spiked_atleast_once_on_trial)
						
						scsuo.lr_autocorrelation_lag_1 = linregress(scsuo.first_spikes_for_trials_which_spiked_atleast_once_on_trial[1:], scsuo.first_spikes_for_trials_which_spiked_atleast_once_on_trial[:-1])

						# See https://www.analyticsvidhya.com/blog/2018/09/non-stationary-time-series-python/
						dftest = adfuller(scsuo.first_spikes_for_trials_which_spiked_atleast_once_on_trial, autolag='AIC')
						scsuo.adfuller_dataframe = pd.Series(dftest[0:4], index=['Test Statistic','p-value','#Lags Used','Number of Observations Used'])

						with warnings.catch_warnings():
							warnings.filterwarnings("ignore")
							scsuo.kpss_pvalue = kpss(scsuo.first_spikes_for_trials_which_spiked_atleast_once_on_trial, nlags='auto')[1]


					self.list_of_all_scsuos.append(scsuo)


			