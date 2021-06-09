from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
from os.path import isfile, join
import seaborn as sns
sns.set()
from matplotlib.ticker import MultipleLocator
import math
from scipy.stats import linregress

import spikewarp as sw


def draw_spike_count_plots(single_neuron_analysis_figs_dir, data_loaders, investigation):
	"""
	"""
	dictionary_of_info_by_neuron_group = {'L23EXC':{},
										'L23INH':{},

										'L4EXC':{},
										'L4INH':{},

										'L5AEXC':{},
										'L5AINH':{},

										'L5BEXC':{},
										'L5BINH':{},
										}

	
	for neuron_group_dictionary in dictionary_of_info_by_neuron_group.values():
		neuron_group_dictionary['non_zero_spike_counts'] = []
	
	spike_count_means = []
	spike_count_variances = []
	atleast_one_spike_reliabilities = []
	all_non_zero_spike_counts = []

	for experiment_index in investigation["experiments_to_process"]:
		data_loader = data_loaders[investigation["experiments_to_process"].index(experiment_index)]
		for scsuo in data_loader.list_of_all_scsuos:

			neuron_group_key = scsuo.neurophys_corresponding_layer_string_of_unit + scsuo.exc_inh_type_of_unit

			dictionary_of_info_by_neuron_group[neuron_group_key]['non_zero_spike_counts'].extend(scsuo.spike_counts_per_trial[np.argwhere(scsuo.spike_counts_per_trial > 0)])
			spike_count_means.append(scsuo.mean_trial_spike_count_for_condition)
			spike_count_variances.append(scsuo.variance_trial_spike_count_for_condition)
			atleast_one_spike_reliabilities.append(scsuo.atleast_one_spike_reliability)
			all_non_zero_spike_counts.extend(scsuo.spike_counts_per_trial[np.argwhere(scsuo.spike_counts_per_trial > 0)])

	sw.draw_spike_count_histograms(dictionary_of_info_by_neuron_group, single_neuron_analysis_figs_dir + "NonZeroSpikeCountHistoByNeuronGroup.pdf")

	sw.normal_histo_plot([atleast_one_spike_reliabilities], single_neuron_analysis_figs_dir + "atleast_one_spike_reliabilities",	bins=20, density=False, x_axis_label='Did spike once reliability', y_axis_label='Frequency')

	n,_ = np.histogram(all_non_zero_spike_counts, bins=np.arange(1, 8) - 0.5, density=True)
	sw.draw_neighbouring_bar_chart([n], ('1', '2', '3', '4', '5', '6', '7', '8'), single_neuron_analysis_figs_dir + "NonZeroSpikeCountHisto.pdf", '', (''), 'Spike count', custom_y_tick_locators=[1, .2], optional_y_axis_label='Normalised frequency')


	sw.basic_x_y_plot([spike_count_means], [spike_count_variances], single_neuron_analysis_figs_dir + "SpikeCountMeansVsVariances.pdf", draw_y_equals_x=True, y_equals_x_max=3.0)


def single_unit_stationarity_tests(single_neuron_analysis_figs_dir, draw_individual_unit_plots, data_loaders, investigation):

	ALL_FIRST_SPIKES_FOR_HISTOGRAM = []
	ALL_lr_actual_trial_indices_vs_fs = []
	ALL_lr_autocorrelation_lag_1s = []
	ALL_adfuller_dataframes = []
	ALL_KPSS_pvalues = []
	acf_pvalues = []
	pacf_pvalues = []
	latencies_vs_trial_plots_directory = single_neuron_analysis_figs_dir + "LatenciesVsTrialPlots/"
	sw.mkdir(latencies_vs_trial_plots_directory)

	number_of_lags = 5


	dictionary_of_non_stationarity_dictionaries_by_neuron_group = {'L23EXC':{},
													'L23INH':{},

													'L4EXC':{},
													'L4INH':{},

													'L5AEXC':{},
													'L5AINH':{},

													'L5BEXC':{},
													'L5BINH':{},
					}


	freq_group_keys = ['0-0.2', '1.0', '2.0', '3.0', '4.0', '5.0', '6.0', '7.0', '8.0', '9.0', '10.0']
	number_of_blocks = 20

	for neuron_group_dictionary in dictionary_of_non_stationarity_dictionaries_by_neuron_group.values():
		# print(dictionary)
		neuron_group_dictionary['KPSS_pvalues'] = []
		neuron_group_dictionary['TI_LR_pvalues'] = []

		neuron_group_dictionary['TI_LR_pvalues_by_freq_group_and_block'] = {}
		neuron_group_dictionary['TI_LR_pvalues_by_freq_group'] = {}
		neuron_group_dictionary['spike_times_by_trial_index_by_freq_group'] = {}
		neuron_group_dictionary['spike_counts_by_trial_index_by_freq_group'] = {}

		for freq_group_key in freq_group_keys:
			neuron_group_dictionary['TI_LR_pvalues_by_freq_group_and_block'][freq_group_key] = {}
			neuron_group_dictionary['TI_LR_pvalues_by_freq_group'][freq_group_key] = []
			neuron_group_dictionary['spike_times_by_trial_index_by_freq_group'][freq_group_key] = {}
			neuron_group_dictionary['spike_counts_by_trial_index_by_freq_group'][freq_group_key] = {}

			for trial_index in range(100):
				neuron_group_dictionary['spike_times_by_trial_index_by_freq_group'][freq_group_key][str(trial_index)] = []
				neuron_group_dictionary['spike_counts_by_trial_index_by_freq_group'][freq_group_key][str(trial_index)] = []

			for block_index in range(number_of_blocks):
				neuron_group_dictionary['TI_LR_pvalues_by_freq_group_and_block'][freq_group_key][str(block_index)] = []


	TI_LR_pvalues_by_freq_group_NOT_neuron_group = {}
	for freq_group_key in freq_group_keys:
		TI_LR_pvalues_by_freq_group_NOT_neuron_group[freq_group_key] = []


	for experiment_index in investigation["experiments_to_process"]:

		experiment_latencies_vs_trial_plots_directory = latencies_vs_trial_plots_directory + str(experiment_index) + "/"
		sw.mkdir(experiment_latencies_vs_trial_plots_directory)

		data_loader = data_loaders[investigation["experiments_to_process"].index(experiment_index)]
		for scsuo in data_loader.list_of_all_scsuos:
			first_spikes_for_trials_which_spiked_atleast_once_on_trial_flattened = scsuo.first_spikes_for_trials_which_spiked_atleast_once_on_trial.flatten()

			if (scsuo.lr_actual_trial_indices_vs_fs != None):

				freq_group_key = str(scsuo.stimulation_frequency_of_condition)
				if (scsuo.stimulation_frequency_of_condition < 0.3):
					freq_group_key = '0-0.2'

				if (freq_group_key in freq_group_keys):

					ALL_lr_actual_trial_indices_vs_fs.append(scsuo.lr_actual_trial_indices_vs_fs)
					ALL_lr_autocorrelation_lag_1s.append(scsuo.lr_autocorrelation_lag_1)
					ALL_adfuller_dataframes.append(scsuo.adfuller_dataframe)
					ALL_KPSS_pvalues.append(scsuo.kpss_pvalue)
					
					acf_rvalues_0, acf_pvals_0, acf_positive_pvals_0, acf_negative_pvals_0 = sw.calculate_acf_pvalues_for_spikes(first_spikes_for_trials_which_spiked_atleast_once_on_trial_flattened, first_spikes_for_trials_which_spiked_atleast_once_on_trial_flattened, number_of_lags)
					pacf_rvalues, pacf_pvals, pacf_positive_pvals, pacf_negative_pvals = sw.calculate_pacf_pvalues_for_spikes(first_spikes_for_trials_which_spiked_atleast_once_on_trial_flattened, first_spikes_for_trials_which_spiked_atleast_once_on_trial_flattened, number_of_lags)
					acf_pvalues.append(acf_pvals_0)
					pacf_pvalues.append(pacf_pvals)

					neuron_group_key = scsuo.neurophys_corresponding_layer_string_of_unit + scsuo.exc_inh_type_of_unit

					for trial_index in range(100):
						if (trial_index < scsuo.number_of_trials_for_condition):
							if (scsuo.true_if_spiked_atleast_once_on_trial[trial_index]):
								dictionary_of_non_stationarity_dictionaries_by_neuron_group[neuron_group_key]['spike_times_by_trial_index_by_freq_group'][freq_group_key][str(trial_index)].append(scsuo.first_spikes_for_condition_trials_normalised[trial_index])
							dictionary_of_non_stationarity_dictionaries_by_neuron_group[neuron_group_key]['spike_counts_by_trial_index_by_freq_group'][freq_group_key][str(trial_index)].append(scsuo.spike_counts_per_trial_normalised[trial_index])

					if (scsuo.number_of_trials_with_atleast_one_spike > 20):


						
						dictionary_of_non_stationarity_dictionaries_by_neuron_group[neuron_group_key]['TI_LR_pvalues_by_freq_group'][freq_group_key].append(scsuo.lr_actual_trial_indices_vs_fs.pvalue)
						TI_LR_pvalues_by_freq_group_NOT_neuron_group[freq_group_key].append(scsuo.lr_actual_trial_indices_vs_fs.pvalue)

						dictionary_of_non_stationarity_dictionaries_by_neuron_group[neuron_group_key]['KPSS_pvalues'].append(scsuo.kpss_pvalue)
						dictionary_of_non_stationarity_dictionaries_by_neuron_group[neuron_group_key]['TI_LR_pvalues'].append(scsuo.lr_actual_trial_indices_vs_fs.pvalue)

						number_of_20_trial_spiking_blocks = math.floor(float(len(first_spikes_for_trials_which_spiked_atleast_once_on_trial_flattened)) / 20.0)

						for i in range(number_of_20_trial_spiking_blocks):
							if (i < 20):
							# print(i, number_of_20_trial_spiking_blocks)

								lr_pvalue_for_block = linregress(first_spikes_for_trials_which_spiked_atleast_once_on_trial_flattened[i * 20: i * 20 + 20], range(20)).pvalue

								if (freq_group_key in freq_group_keys):
									dictionary_of_non_stationarity_dictionaries_by_neuron_group[neuron_group_key]['TI_LR_pvalues_by_freq_group_and_block'][freq_group_key][str(i)].append(lr_pvalue_for_block)




						# print(scsuo.lr_actual_trial_indices_vs_fs.rvalue)
						# autocorrelation_shift_1 = np.corrcoef(scsuo.first_spikes_for_trials_which_spiked_atleast_once_on_trial[1:], scsuo.first_spikes_for_trials_which_spiked_atleast_once_on_trial[:-1])[0, 1]
						# ALL_autocorrelation_shift_1s.append(autocorrelation_shift_1)

						# if (scsuo.lr_actual_trial_indices_vs_fs.pvalue < 0.1):

						
						if (draw_individual_unit_plots):

							if (experiment_index == 1):

								sns.set()
								sns.set_style("ticks")
								plt.figure(figsize=(2, 4))
								plt.scatter(scsuo.first_spikes_for_trials_which_spiked_atleast_once_on_trial, scsuo.indices_of_trials_with_atleast_one_spike, s=2)
								plt.title(scsuo.neurophys_corresponding_layer_string_of_unit + " " + scsuo.exc_inh_type_of_unit + ", Stim Freq: " + str(scsuo.stimulation_frequency_of_condition) + '\nPrinciple Column: ' + str(scsuo.principle_condition_unit_combination) + '\nAdfuller p: ' + str(scsuo.adfuller_dataframe[1]) + '\nKPSS p-value: ' + str(scsuo.kpss_pvalue) + "\nLag-1 autocorrelation p-value: " + str(scsuo.lr_autocorrelation_lag_1.pvalue) + "\nTrial index vs spike time p-value: " + str(scsuo.lr_actual_trial_indices_vs_fs.pvalue))
								plt.gca().set_xlabel('ms', fontsize=24)
								plt.gca().set_ylabel('Trial', fontsize=24)
								plt.gca().set_ylim([0, plt.gca().get_ylim()[1]])
								plt.gca().set_xlim([5.0, 55.8])
								
								# plt.gca().xaxis.set_major_locator(FixedLocator([5.0, 30.0]))
								# plt.gca().yaxis.set_major_locator(FixedLocator([0, scsuo.number_of_trials_for_condition]))

								plt.gca().xaxis.set_minor_locator(MultipleLocator(5.0))
								# plt.gca().yaxis.set_minor_locator(MultipleLocator(20))

								for tick in plt.gca().xaxis.get_major_ticks():
									tick.label.set_fontsize(24) 
								for tick in plt.gca().yaxis.get_major_ticks():
									tick.label.set_fontsize(24) 

								if (scsuo.number_of_trials_for_condition == 100):
									plt.savefig(experiment_latencies_vs_trial_plots_directory + str(experiment_index) + "_" +  str(scsuo.channel_unit_pair_index) + "_" + str(scsuo.condition_index) + ".pdf", bbox_inches='tight')
								plt.close()


	
	count_or_time_keys = ['spike_times_by_trial_index_by_freq_group', 'spike_counts_by_trial_index_by_freq_group']
	count_or_time_y_axis_label = ['Normalised\nspike time', 'Normalised\nspike count']


	for count_or_time_index, count_or_time_key in enumerate(count_or_time_keys):

		block_dict = {}

		for neuron_group_key, neuron_group_dictionary in dictionary_of_non_stationarity_dictionaries_by_neuron_group.items():
			block_dict[neuron_group_key] = {}
			block_dict[neuron_group_key]['block_means'] = []
			block_dict[neuron_group_key]['block_pos_95_cis'] = []

			for block_index in range(10):

				block_spike_times = []
			
				for freq_group_index, freq_group_key in enumerate(freq_group_keys):

					for trial_index in range(10):
						trial_spike_times = neuron_group_dictionary[count_or_time_key][freq_group_key][str(block_index*10 + trial_index)]
						block_spike_times.extend(trial_spike_times)

				se = 0.0
				mean = 0.0
				if (len(block_spike_times) > 0):
					se = np.std(block_spike_times) / np.sqrt(len(block_spike_times))
					mean = np.mean(block_spike_times)
				block_dict[neuron_group_key]['block_means'].append(mean)
				block_dict[neuron_group_key]['block_pos_95_cis'].append(1.96 * se)


		plt.figure(figsize=(13, 14))
		fig, axes = plt.subplots(4, 1)

		for ei_type_index, (ei_type, ei_type_color) in enumerate(zip(['EXC', 'INH'], ['Reds', 'Blues'])):
			
			colour_map_string = ei_type_color
			my_cmap = plt.cm.get_cmap(colour_map_string, 4)
			xvals = np.arange(0, 10)
			for neuron_group_index, neuron_group_key in enumerate(dictionary_of_non_stationarity_dictionaries_by_neuron_group.keys()):

				if (ei_type in neuron_group_key):
					# print(6 - math.floor(neuron_group_index / 2) )
					# c = my_cmap(6 - math.floor(neuron_group_index / 2))
					c = my_cmap(2)
					ax_index = math.floor(neuron_group_index / 2)
					ax = axes[math.floor(neuron_group_index / 2)]
					ax.errorbar(xvals + ei_type_index * 0.1 - 0.05, block_dict[neuron_group_key]['block_means'], yerr=block_dict[neuron_group_key]['block_pos_95_cis'], c=c, label=neuron_group_key, ls='-', marker='o', markersize=2, lw=.3) # '.'
					ax.set_ylabel(count_or_time_y_axis_label[count_or_time_index], fontsize=10)

					# ax.set_ylim([-0.25, 0.6])
					ax.set_ylim([-0.4, 0.4])
					ax.yaxis.set_major_locator(MultipleLocator(0.4))	
					ax.yaxis.set_minor_locator(MultipleLocator(0.2))

					ax.grid(which='major', axis='y', linestyle='-', c='lightgrey', linewidth=0.5)
					ax.grid(which='minor', axis='y', linestyle='-', c='lightgrey', linewidth=0.5)

					ax.set_xticks(xvals)
					if (ax_index == 3):
						ax.set_xticklabels([str(block_index*10) + '-' + str(block_index*10 + 9) for block_index in range(10)])
						ax.set_xlabel('Trial group', fontsize=10)
					else:
						ax.set_xticklabels(['' for i in range(len(freq_group_keys))])

					for tick in ax.xaxis.get_major_ticks():
						tick.label.set_fontsize(10) 
					for tick in ax.yaxis.get_major_ticks():
						tick.label.set_fontsize(10) 

		plt.savefig(single_neuron_analysis_figs_dir + 'AllNeuronGroup' + count_or_time_key + '.pdf')
		plt.close()



		# meta_freq_groups = [['1.0', '3.0', '6.0', '10.0'], ['0-0.2', '2.0', '4.0'], ['5.0', '7.0', '8.0', '9.0']] # Mainz Paper
		meta_freq_groups = [['0-0.2', '1.0', '2.0', '3.0'], ['4.0', '5.0', '6.0', '7.0'], ['8.0', '9.0', '10.0']] # Thesis
		for meta_freq_group_index, meta_freq_group in enumerate(meta_freq_groups):
			for freq_group_index, freq_group_key in enumerate(freq_group_keys):

				block_dict[freq_group_key] = {}
				block_dict[freq_group_key]['block_means'] = []
				block_dict[freq_group_key]['block_pos_95_cis'] = []

				for block_index in range(10):

					block_spike_times = []
				
					for neuron_group_key, neuron_group_dictionary in dictionary_of_non_stationarity_dictionaries_by_neuron_group.items():

						for trial_index in range(10):
							trial_spike_times = neuron_group_dictionary[count_or_time_key][freq_group_key][str(block_index*10 + trial_index)]
							block_spike_times.extend(trial_spike_times)

					se = np.std(block_spike_times) / np.sqrt(len(block_spike_times))
					block_dict[freq_group_key]['block_means'].append(np.mean(block_spike_times))
					block_dict[freq_group_key]['block_pos_95_cis'].append(1.96 * se)




			plt.figure(figsize=(8, 3))
			colour_map_string = "Greens"
			my_cmap = plt.cm.get_cmap(colour_map_string, len(meta_freq_group) * 2)
			xvals = np.arange(0, 10)
			for freq_group_index, freq_group_key in enumerate(freq_group_keys):
				if (freq_group_key in meta_freq_group):

					index_in_meta_group = meta_freq_group.index(freq_group_key)

					c = my_cmap(index_in_meta_group + len(meta_freq_group))
					plt.errorbar(xvals + index_in_meta_group * 0.08 - .12, block_dict[freq_group_key]['block_means'], yerr=block_dict[freq_group_key]['block_pos_95_cis'], c=c, fmt='.', label=freq_group_key, ls='-', marker='o', markersize=2, lw=.3)


			ax = plt.gca()
			ax.set_ylabel(count_or_time_y_axis_label[count_or_time_index], fontsize=10)

			# plt.gca().set_ylim([-0.25, 0.65])
			ax.yaxis.set_major_locator(MultipleLocator(0.2))	

			# ax.set_ylim([-0.5, 0.4])
			# ax.yaxis.set_major_locator(MultipleLocator(0.4))
			# ax.yaxis.set_minor_locator(MultipleLocator(0.2))

			ax.grid(which='major', axis='y', linestyle='-', c='lightgrey', linewidth=0.5)
			ax.grid(which='minor', axis='y', linestyle='-', c='lightgrey', linewidth=0.5)

			ax.set_xticks(xvals)
			if (ax_index == 3):
				ax.set_xticklabels([str(block_index*10) + '-' + str(block_index*10 + 9) for block_index in range(10)])
				ax.set_xlabel('Trial group', fontsize=10)
			else:
				ax.set_xticklabels(['' for i in range(len(freq_group_keys))])

			for tick in ax.xaxis.get_major_ticks():
				tick.label.set_fontsize(10) 
			for tick in ax.yaxis.get_major_ticks():
				tick.label.set_fontsize(10) 

			plt.legend(edgecolor='k')
			plt.savefig(single_neuron_analysis_figs_dir + str(meta_freq_group_index) + 'AllByFreqGroup' + count_or_time_key + '.pdf')
			plt.close()

	







	for neuron_group_key, neuron_group_dictionary in dictionary_of_non_stationarity_dictionaries_by_neuron_group.items():
		neuron_group_KPSS_pvalues = np.asarray(neuron_group_dictionary['KPSS_pvalues'])
		neuron_group_LR_pvalues = np.asarray(neuron_group_dictionary['TI_LR_pvalues'])

		number_of_KPSS_samples_for_neuron_group = len(neuron_group_KPSS_pvalues)
		number_of_LR_samples_for_neuron_group = len(neuron_group_LR_pvalues)

		neuron_group_dictionary['Proportion_KPSS_pvalues_less_than_threshold'] = 0.0
		neuron_group_dictionary['Proportion_LR_pvalues_less_than_threshold'] = 0.0
		neuron_group_dictionary['Proportion_KPSS_pvalues_less_than_threshold_positive_confidence_interval'] = 0.0
		neuron_group_dictionary['Proportion_LR_pvalues_less_than_threshold_positive_confidence_interval'] = 0.0

		if (number_of_KPSS_samples_for_neuron_group > 0):
			neuron_group_dictionary['Proportion_KPSS_pvalues_less_than_threshold'] = np.sum(neuron_group_KPSS_pvalues < 0.05) / number_of_KPSS_samples_for_neuron_group
			neuron_group_dictionary['Proportion_KPSS_pvalues_less_than_threshold_positive_confidence_interval'] = 1.96 * math.sqrt((neuron_group_dictionary['Proportion_KPSS_pvalues_less_than_threshold'] * (1-neuron_group_dictionary['Proportion_KPSS_pvalues_less_than_threshold']) / number_of_KPSS_samples_for_neuron_group))
		if (number_of_LR_samples_for_neuron_group > 0):
			neuron_group_dictionary['Proportion_LR_pvalues_less_than_threshold'] = np.sum(neuron_group_LR_pvalues < 0.05) / number_of_LR_samples_for_neuron_group
			neuron_group_dictionary['Proportion_LR_pvalues_less_than_threshold_positive_confidence_interval'] = 1.96 * math.sqrt((neuron_group_dictionary['Proportion_LR_pvalues_less_than_threshold'] * (1-neuron_group_dictionary['Proportion_LR_pvalues_less_than_threshold']) / number_of_LR_samples_for_neuron_group))
		# Used https://www.dummies.com/education/math/statistics/how-to-determine-the-confidence-interval-for-a-population-proportion/
		
	

	# sw.draw_neighbouring_bar_chart(stimulation_frequency_of_condition_tests_passed_and_correlated_for_each_key, ('0-0.2', '1', '2', '3', '4', '5', '6', '7', '8'), dh.unit_type_analysis_directory + "stimulation_frequency_grouped.pdf", '', ('Stationary', 'ARIMA'), 'Hz', custom_y_tick_locators=[20, 5])
	sw.draw_neighbouring_bar_chart([[neuron_group_dictionary['Proportion_KPSS_pvalues_less_than_threshold'] for neuron_group_dictionary in dictionary_of_non_stationarity_dictionaries_by_neuron_group.values()], 
														[neuron_group_dictionary['Proportion_LR_pvalues_less_than_threshold'] for neuron_group_dictionary in dictionary_of_non_stationarity_dictionaries_by_neuron_group.values()]], 
														dictionary_of_non_stationarity_dictionaries_by_neuron_group.keys(), 
														single_neuron_analysis_figs_dir + "NonStationarityByNeuronGroup.pdf", 
														'', 
														('KPSS', 'Trial Index Linear Regression'), 
														"", 
														custom_y_tick_locators=[.3, .05], 
														rotate_x_labels=True, 
														y_tick_right=True, 
														optional_y_axis_label="Proportion",
														positive_confidence_intervals=[[neuron_group_dictionary['Proportion_KPSS_pvalues_less_than_threshold_positive_confidence_interval'] for neuron_group_dictionary in dictionary_of_non_stationarity_dictionaries_by_neuron_group.values()], 
														[neuron_group_dictionary['Proportion_LR_pvalues_less_than_threshold_positive_confidence_interval'] for neuron_group_dictionary in dictionary_of_non_stationarity_dictionaries_by_neuron_group.values()]],
														threshold_lin_value=0.05)


	proportion_TI_LR_pvalues_by_freq_group_NOT_neuron_group_less_than_threshold = {}
	proportion_TI_LR_pvalues_by_freq_group_NOT_neuron_group_less_than_threshold_positive_confidence_interval = {}
	for freq_group_key in freq_group_keys:

		TI_LR_pvalues_for_freq_group = np.asarray(TI_LR_pvalues_by_freq_group_NOT_neuron_group[freq_group_key])

		number_TI_LR_pvalues_by_freq_group_NOT_neuron_group = len(TI_LR_pvalues_for_freq_group)

		proportion_TI_LR_pvalues_by_freq_group_NOT_neuron_group_less_than_threshold[freq_group_key] = np.sum(TI_LR_pvalues_for_freq_group < 0.05) / number_TI_LR_pvalues_by_freq_group_NOT_neuron_group
		proportion_TI_LR_pvalues_by_freq_group_NOT_neuron_group_less_than_threshold_positive_confidence_interval[freq_group_key] = 1.96 * math.sqrt((proportion_TI_LR_pvalues_by_freq_group_NOT_neuron_group_less_than_threshold[freq_group_key] * (1-proportion_TI_LR_pvalues_by_freq_group_NOT_neuron_group_less_than_threshold[freq_group_key]) / number_TI_LR_pvalues_by_freq_group_NOT_neuron_group))


	sw.draw_neighbouring_bar_chart([[proportion_TI_LR_pvalues_by_freq_group_NOT_neuron_group_less_than_threshold[freq_group_key] for freq_group_key in freq_group_keys]], 
														freq_group_keys, 
														single_neuron_analysis_figs_dir + "NonStationarityByFrequency.pdf", 
														'', 
														('Trial Index Linear Regression'), 
														"", 
														custom_y_tick_locators=[.3, .05], 
														rotate_x_labels=True, 
														y_tick_right=True, 
														optional_y_axis_label="Proportion",
														positive_confidence_intervals=[[proportion_TI_LR_pvalues_by_freq_group_NOT_neuron_group_less_than_threshold[freq_group_key] for freq_group_key in freq_group_keys]],
														threshold_lin_value=0.05)


	acf_pvalues = np.asarray(acf_pvalues)
	pacf_pvalues = np.asarray(pacf_pvalues)

	sw.cumulative_histo_plot([acf_pvalues[:, lag_index_zeroed] for lag_index_zeroed in range(number_of_lags)], single_neuron_analysis_figs_dir + "_ACF_PVal_CumHist", bins=200, x_axis_label="p-value", y_axis_label="Normalised\ncumulative sum", custom_x_tick_locators=[1.0, 0.2], add_chi_squared_text=True)
	sw.cumulative_histo_plot([pacf_pvalues[:, lag_index_zeroed] for lag_index_zeroed in range(number_of_lags)], single_neuron_analysis_figs_dir + "_PACF_PVal_CumHist", bins=200, x_axis_label="p-value", y_axis_label="Normalised\ncumulative sum", custom_x_tick_locators=[1.0, 0.2], add_chi_squared_text=True)
	for lag_index_zeroed in range(number_of_lags):
		sw.normal_histo_plot([acf_pvalues[:, lag_index_zeroed]], single_neuron_analysis_figs_dir + "_ACFLag" + str(lag_index_zeroed + 1), bins=20, histo_range=[0.0, 1.0], x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)
		sw.normal_histo_plot([pacf_pvalues[:, lag_index_zeroed]], single_neuron_analysis_figs_dir + "_PACFLag" + str(lag_index_zeroed + 1), bins=20, histo_range=[0.0, 1.0], x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)

	ALL_pvalues_actual_trial_indices_vs_fs = [lr.pvalue for lr in ALL_lr_actual_trial_indices_vs_fs if lr != None]
	ALL_rvalues_actual_trial_indices_vs_fs = [lr.rvalue for lr in ALL_lr_actual_trial_indices_vs_fs if lr != None]
	ALL_rsquared_values_actual_trial_indices_vs_fs = [lr.rvalue**2 for lr in ALL_lr_actual_trial_indices_vs_fs if lr != None]
	PVALTHRESH_rvalues_actual_trial_indices_vs_fs = [lr.rvalue for lr in ALL_lr_actual_trial_indices_vs_fs if (lr != None) & (lr.pvalue < 0.005)]
	PVALTHRESH_rsquared_values_actual_trial_indices_vs_fs = [lr.rvalue**2 for lr in ALL_lr_actual_trial_indices_vs_fs if (lr != None) & (lr.pvalue < 0.005)]
	ALL_adfuller_pvalues = [adfuller_dataframe['p-value'] for adfuller_dataframe in ALL_adfuller_dataframes]
	PVALTHRESH_rvalues_actual_trial_indices_vs_fs = np.asarray(PVALTHRESH_rvalues_actual_trial_indices_vs_fs)
	ALL_KPSS_pvalues = np.asarray(ALL_KPSS_pvalues)
	ALL_LR_pvalues = np.asarray(ALL_pvalues_actual_trial_indices_vs_fs)

	print('Proportion of thresholded r-values positive: ' + str(np.sum(PVALTHRESH_rvalues_actual_trial_indices_vs_fs > 0)) + " of " + str(PVALTHRESH_rvalues_actual_trial_indices_vs_fs.shape[0]))
	print('Proportion of KPSS less than threshold: ' + str(np.sum(ALL_KPSS_pvalues < 0.05) / len(ALL_KPSS_pvalues)))
	print('Proportion of LRs p-values less than threshold: ' + str(np.sum(ALL_LR_pvalues < 0.05) / len(ALL_LR_pvalues)))

	sw.normal_histo_plot([ALL_pvalues_actual_trial_indices_vs_fs], single_neuron_analysis_figs_dir + "PVALUE_HISTO_trial_vs_fs", bins=20, x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[175, 175], alpha=0.78, add_chi_squared_text=True)
	sw.normal_histo_plot([ALL_rvalues_actual_trial_indices_vs_fs], single_neuron_analysis_figs_dir + "RVALUE_HISTO_trial_vs_fs", colors=['g'], bins=40, histo_range=[-1.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$r$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[100, 100], alpha=0.78, add_chi_squared_text=True)
	sw.normal_histo_plot([ALL_rsquared_values_actual_trial_indices_vs_fs], single_neuron_analysis_figs_dir + "RSQUARED_VALUE_HISTO _trial_vs_fs", colors=['b'], bins=20, histo_range=[0.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$r^2$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[650, 650], alpha=0.78)
	sw.normal_histo_plot([PVALTHRESH_rvalues_actual_trial_indices_vs_fs], single_neuron_analysis_figs_dir + "PVALTHRESH_RVALUE_HISTO_trial_vs_fs", colors=['g'], bins=40, histo_range=[-1.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$r$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[25, 5], alpha=0.78)
	sw.normal_histo_plot([PVALTHRESH_rsquared_values_actual_trial_indices_vs_fs], single_neuron_analysis_figs_dir + "PVALTHRESH_RSQUARED_VALUE_HISTO_trial_vs_fs", colors=['b'], bins=20, histo_range=[0.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$r^2$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[50, 50], alpha=0.78)
	sw.normal_histo_plot([PVALTHRESH_rsquared_values_actual_trial_indices_vs_fs, PVALTHRESH_rvalues_actual_trial_indices_vs_fs], single_neuron_analysis_figs_dir + "PVALTHRESH_R_AND_RSQUARED_VALUE_HISTO_trial_vs_fs", colors=['b', 'g'], bins=40, histo_range=[-1.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$r$, $r^2$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[16, 4], alpha=0.7)
	sw.normal_histo_plot([ALL_adfuller_pvalues], single_neuron_analysis_figs_dir + "PVALUE_HISTO_adfuller", bins=20, x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[120, 120], alpha=0.78, add_chi_squared_text=True)
	sw.normal_histo_plot([ALL_KPSS_pvalues], single_neuron_analysis_figs_dir + "PVALUE_HISTO_kpss", bins=6, histo_range=[0.0, 0.12], x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[0.12, 0.02], custom_y_tick_locators=[600, 50], alpha=0.78, add_chi_squared_text=True)





def load_spikes_calc_cort_onset_and_draw_PSTHs(investigation):

	random_seed_counter = 0

	# LOADING SPIKES AND CREATE PSTH
	data_loaders = []
	ALL_SPIKES_FOR_HISTOGRAM = []
	RELIABILITY_THRESHOLDED_SPIKES_FOR_HISTOGRAM = []
	for experiment_index in investigation["experiments_to_process"]:
		
		experiment_original_data_root = investigation['original_data_root'] + investigation["experiment_handles"][experiment_index] + '/'
		spikes_directory_name = 'CutOutSpikes/Spikes_' + investigation["original_spike_cut_string"] + '/'

		onlyfiles = [f for f in listdir(experiment_original_data_root + spikes_directory_name) if isfile(join(experiment_original_data_root + spikes_directory_name, f))]
		number_of_units = len(onlyfiles)
		channel_unit_pairs_or_units_for_analysis = [x for x in range(1, number_of_units + 1)]
		
		
		data_loader = sw.DataLoader(experiment_original_data_root)
		data_loader.LOAD_SPIKES_FOR_ALL_UNITS_AND_CONDITIONS(spikes_directory_name, channel_unit_pairs_or_units_for_analysis, investigation['channel_unit_pair_format'])
		data_loader.LOAD_CONDITION_INFO(load_frequency_and_whisker_information=True, max_stimulation_frequency=investigation['max_stimulation_frequency'])
		data_loader.LOAD_UNIT_INFO()
		data_loader.create_scsuos(spike_train_start_time=investigation["original_spike_train_start_time"], 
													spike_train_end_time=investigation["original_spike_train_end_time"], 
													custom_number_of_trials_for_condition=-1)	

		data_loaders.append(data_loader)

		for scsuo in data_loader.list_of_all_scsuos:
			for trial_index in range(scsuo.number_of_trials_for_condition):
				ALL_SPIKES_FOR_HISTOGRAM.extend([value for value in scsuo.spike_trains_for_condition_trials[trial_index]])
			# ALL_FIRST_SPIKES_FOR_HISTOGRAM.extend(scsuo.first_spikes_for_trials_which_spiked_atleast_once_on_trial.flatten())

			if (scsuo.number_of_trials_with_atleast_one_spike > 15):
				RELIABILITY_THRESHOLDED_SPIKES_FOR_HISTOGRAM.extend(scsuo.first_spikes_for_trials_which_spiked_atleast_once_on_trial.flatten())
				# for trial_index in range(scsuo.number_of_trials_for_condition):
					# RELIABILITY_THRESHOLDED_SPIKES_FOR_HISTOGRAM.extend([value for value in scsuo.spike_trains_for_condition_trials[trial_index]])


	sns.set()
	sns.set_style("ticks")
	plt.figure()
	n, _, _ = plt.hist(ALL_SPIKES_FOR_HISTOGRAM, bins=investigation["original_spike_cut_length"] * 5)
	
	sum_of_hists_over_all_trials_SMOOTHED = gaussian_filter(n, sigma=1*5)
	prestimulus_analysis_psth_start = -50
	prestimulus_analysis_psth_end = -1
	cortical_activity_search_start = 0
	zeroed_prestimulus_analysis_psth_end = 5 * (prestimulus_analysis_psth_end - investigation["original_spike_train_start_time"])
	zeroed_cortical_activity_search_start = 5 * (cortical_activity_search_start - investigation["original_spike_train_start_time"])
	cortical_activity_search_end = investigation["original_spike_train_end_time"]

	mean_to_use = np.mean(sum_of_hists_over_all_trials_SMOOTHED[:zeroed_prestimulus_analysis_psth_end])
	std_to_use = np.std(sum_of_hists_over_all_trials_SMOOTHED[:zeroed_prestimulus_analysis_psth_end])
	index_above_std = np.where(sum_of_hists_over_all_trials_SMOOTHED[zeroed_cortical_activity_search_start:] > (mean_to_use + 3*std_to_use))[0][0]
	x_vals = np.arange(float(prestimulus_analysis_psth_start), float(cortical_activity_search_end + 1.0), 1.0/5.0)
	cortical_onset_timepoint_above_threshold = (index_above_std / 5.0 + cortical_activity_search_start)
	determined_cortical_onset_timepoint = cortical_onset_timepoint_above_threshold

	plt.scatter([cortical_onset_timepoint_above_threshold], [mean_to_use + 3*std_to_use])
	plt.scatter([determined_cortical_onset_timepoint], [mean_to_use + 3*std_to_use], c='r')
	plt.savefig(investigation['full_spikes_experiment_figures_directory'] + "PSTH_AllStimuli")
	plt.close()


	cortical_activity_search_end = 100.0
	for extra_string, x_axis_lims, xlocator in zip(['Smaller', 'Larger'], [[cortical_activity_search_start, 25.0], [prestimulus_analysis_psth_start, cortical_activity_search_end]], [5, 50]):

		plt.figure()
		plt.plot(x_vals, sum_of_hists_over_all_trials_SMOOTHED)
		plt.scatter([cortical_onset_timepoint_above_threshold], [mean_to_use + 3*std_to_use])
		plt.scatter([determined_cortical_onset_timepoint], [mean_to_use + 3*std_to_use], c='r', s=100)
		plt.gca().set_xlim(x_axis_lims)
		plt.gca().set_xlabel('ms', fontsize=24)
		plt.gca().set_ylabel('Count', fontsize=24)

		plt.gca().xaxis.set_major_locator(MultipleLocator(xlocator))

		for tick in plt.gca().xaxis.get_major_ticks():
			tick.label.set_fontsize(10) 
		for tick in plt.gca().yaxis.get_major_ticks():
			tick.label.set_fontsize(10) 

		plt.savefig(investigation['full_spikes_experiment_figures_directory'] + "StimulusSpikePSTHOverAllTrials" + extra_string + '.pdf', bbox_inches='tight')
		plt.close()



	
	n, _, _ = plt.hist(RELIABILITY_THRESHOLDED_SPIKES_FOR_HISTOGRAM, bins=investigation["original_spike_cut_length"] * 5)
	plt.close()
	plt.figure()
	RELIABILITY_THRESHOLDED_sum_of_hists_over_all_trials_SMOOTHED = gaussian_filter(n, sigma=3)
	plt.plot(x_vals, RELIABILITY_THRESHOLDED_sum_of_hists_over_all_trials_SMOOTHED)
	plt.gca().set_xlim(x_axis_lims)
	plt.savefig(investigation['full_spikes_experiment_figures_directory'] + "FSReliabilityThresholdedPSTH" + extra_string)
	plt.close()


	plt.figure()
	plt.hist(ALL_SPIKES_FOR_HISTOGRAM, bins=investigation["original_spike_cut_length"] * 5, density=True)
	plt.savefig(investigation['full_spikes_experiment_figures_directory'] + "Normalised_PSTH_AllStimuli")
	plt.close()

	plt.figure()
	plt.hist(ALL_SPIKES_FOR_HISTOGRAM, bins=investigation["original_spike_cut_length"] * 5, density=True, cumulative=True)
	plt.savefig(investigation['full_spikes_experiment_figures_directory'] + "CumulativeNormalised_PSTH_AllStimuli")
	plt.close()



	return data_loaders, determined_cortical_onset_timepoint
