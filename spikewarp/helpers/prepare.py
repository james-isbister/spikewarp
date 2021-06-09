import pickle
import numpy as np
import itertools

import spikewarp as sw

"""
Helpers for preparing analysis and investigation dictionaries for main analyses.
"""


class OptionHelper(object):

	"""
	Class for managing shuffle options
	"""

	def __init__(self, shuffle_options, info_types, clustering_param_0s):

		self.shuffle_options = shuffle_options
		self.shuffle_option_strings = [shuffle_option[0] + str(shuffle_option[1]) for shuffle_option in shuffle_options]
		self.info_types = info_types
		self.clustering_param_0s = clustering_param_0s
		
		self.number_of_shuffle_options = len(self.shuffle_options)
		self.number_of_info_types_options = len(self.info_types)
		self.number_of_clustering_param_0_options = len(self.clustering_param_0s)

		self.all_value_combinations = []
		self.all_value_combinations_as_shuffle_type_organised_list = []

		for shuffle_option in self.shuffle_options:
			
			array_of_variable_value_arrays_to_combine = [[shuffle_option], self.info_types, self.clustering_param_0s]
			all_value_combinations_for_shuffle_option = list(itertools.product(*array_of_variable_value_arrays_to_combine))
			self.all_value_combinations.extend(all_value_combinations_for_shuffle_option)
			self.all_value_combinations_as_shuffle_type_organised_list.append(all_value_combinations_for_shuffle_option)

		self.total_number_of_combinations = len(self.all_value_combinations)



def prepare_investigation(investigation):

	"""
	Function for preparing investigatioin dictionary
	"""

	investigation["original_spike_cut_string"] = str(investigation["original_spike_cut"][0]) + "-" + str(investigation["original_spike_cut"][1]);
	investigation["original_spike_cut_string_first_minus_replaced"] = investigation["original_spike_cut_string"].replace('-', 'minus', 1)
	investigation["original_spike_cut_length"] = investigation["original_spike_cut"][1] - investigation["original_spike_cut"][0] + 1; 
	investigation["original_spike_train_start_time"] = investigation["original_spike_cut"][0]
	investigation["original_spike_train_end_time"] = investigation["original_spike_cut"][1]
	investigation["channel_unit_pair_format"] = False
	investigation['full_spikes_experiment_figures_directory'] = investigation['figures_root_directory'] + investigation["original_spike_cut_string_first_minus_replaced"] + "/"

	sw.makedirs(investigation['figures_root_directory'])
	sw.makedirs(investigation['full_spikes_experiment_figures_directory'])


def create_time_window_dependent_dirs(analysis, investigation):

	"""
	Function for creating directory and main subdirectories for analysis of particular time window
	"""

	analysis['time_span_figures_directory'] = investigation['figures_root_directory'] + str(analysis['main_analysis_spike_train_start_time']) + '-' + str(analysis['main_analysis_spike_train_end_time']) + "/"
	analysis['single_neuron_analysis_figs_dir'] = analysis['time_span_figures_directory'] + "SingleNeuronAnalysis/"
	analysis['neuron_pair_figs_dir'] = analysis['time_span_figures_directory'] + "PairwiseAnalysis/"
	sw.makedirs(analysis['time_span_figures_directory'])
	sw.makedirs(analysis['single_neuron_analysis_figs_dir'])
	sw.makedirs(analysis['neuron_pair_figs_dir'])


def prepare_investigation_analysis_recipes(investigation):

	"""
	Function for preparing each analysis dictionary
	"""

	for analysis_name in investigation['analyses_to_process']:

		analysis = investigation['analyses'][analysis_name]

		analysis['number_of_processors'] = investigation['number_of_processors']

		analysis['calulate_time_span'] = False
		analysis['draw_spike_count_plots'] = False
		analysis['do_single_unit_stationarity_tests'] = False
		analysis['run_pairwise_analysis'] = False
		analysis['run_triplet_analysis'] = False
		analysis['do_plot_super_time_warp_cluster_busters'] = False
		analysis['do_plot_super_time_warp_cluster_buster_videos_ALL_PAIRS'] = False
		analysis['do_plot_super_time_warp_cluster_buster_videos_SINGLE_PAIR'] = False
		analysis['do_plot_basic_spike_pair_rasters'] = False
		analysis['do_long_meta_analyses'] = False
		analysis['custom_number_of_trials_for_condition'] = -1
		analysis['main_analysis_spike_train_start_time'] = 0.0
		analysis['main_analysis_spike_train_end_time'] = 0.0


		if (analysis['analysis_type'] == "Mainz_CalculatingTimeSpan"):

			analysis['calulate_time_span'] = True
			analysis['run_pairwise_analysis'] = True
			analysis['main_analysis_spike_train_start_time'] = investigation['determined_cortical_onset_timepoint']
			analysis['main_analysis_spike_train_end_time'] = analysis['main_analysis_spike_train_start_time'] + analysis['window_length']
			create_time_window_dependent_dirs(analysis, investigation)

			shuffle_options = [['normal', -1]]

		if (analysis['analysis_type'] == "SpikeCountPlots"):

			analysis['draw_spike_count_plots'] = True
			analysis['main_analysis_spike_train_start_time'] = investigation['determined_cortical_onset_timepoint']
			# analysis['main_analysis_spike_train_start_time'] = 0.0
			analysis['main_analysis_spike_train_end_time'] = analysis['main_analysis_spike_train_start_time'] + analysis['window_length']
			analysis['main_analysis_spike_train_length'] = analysis['main_analysis_spike_train_end_time'] - analysis['main_analysis_spike_train_start_time']
			create_time_window_dependent_dirs(analysis, investigation)

			shuffle_options = [['normal', -1]]

		if (analysis['analysis_type'] == "SingleUnitStationarity"):

			analysis['do_single_unit_stationarity_tests'] = True
			analysis['main_analysis_spike_train_start_time'] = investigation['determined_cortical_onset_timepoint']
			analysis['main_analysis_spike_train_end_time'] = analysis['main_analysis_spike_train_start_time'] + analysis['window_length']
			analysis['main_analysis_spike_train_length'] = analysis['main_analysis_spike_train_end_time'] - analysis['main_analysis_spike_train_start_time']
			create_time_window_dependent_dirs(analysis, investigation)

			shuffle_options = [['normal', -1]]


		if (analysis['analysis_type'] == "PairwiseAnalysis"):

			analysis['run_pairwise_analysis'] = True
			analysis['do_long_meta_analyses'] = True
			analysis['calulate_time_span'] = True
			analysis['do_plot_super_time_warp_cluster_busters'] = True
			analysis['do_plot_super_time_warp_cluster_buster_videos_ALL_PAIRS'] = True
			# do_plot_super_time_warp_cluster_busters = False
			# do_plot_super_time_warp_cluster_buster_videos_ALL_PAIRS = False

			analysis['main_analysis_spike_train_start_time'] = investigation['determined_cortical_onset_timepoint']
			analysis['main_analysis_spike_train_end_time'] = analysis['main_analysis_spike_train_start_time'] + analysis['window_length']
			analysis['main_analysis_spike_train_length'] = analysis['main_analysis_spike_train_end_time'] - analysis['main_analysis_spike_train_start_time']
			create_time_window_dependent_dirs(analysis, investigation)

			if (analysis['clustering_type'] == "Unclustered"):
				analysis['clustering_param_0s'] = np.arange(analysis['n_components_start'], analysis['n_components_end'], analysis['n_components_step'])
			if (analysis['clustering_type'] == "CustomDBSCAN"):
				analysis['clustering_param_0s'] = np.arange(analysis['epsilon_start'], analysis['epsilon_end'], analysis['epsilon_step'])
			if (analysis['clustering_type'] == "GaussianMixtureInfCrit"):
				analysis['clustering_param_0s'] = np.arange(analysis['n_components_start'], analysis['n_components_end'], analysis['n_components_step'])
			if (analysis['clustering_type'] == "CustomGaussianMixtureInfCrit"):
				analysis['clustering_param_0s'] = np.arange(analysis['n_components_start'], analysis['n_components_end'], analysis['n_components_step'])

			if (analysis['clustering_type'] == "DBSCANDensityPeaks"):
				analysis['clustering_param_0s'] = np.arange(analysis['epsilon_start'], analysis['epsilon_end'], analysis['epsilon_step'])
			if (analysis['clustering_type'] == "GaussianMixture"):
				analysis['clustering_param_0s'] = np.arange(analysis['n_components_start'], analysis['n_components_end'], analysis['n_components_step'])

			analysis['option_helper'] = sw.OptionHelper(analysis['shuffle_options'], [''], analysis['clustering_param_0s'])


			analysis['super_time_warp_cluster_buster_option_objects'] = [
				{
					'type_string': 'type_0',
					'do_plot_3sd_conjunctive_ring': True,
					'do_plot_4sd_conjunctive_ring': True,
					'do_plot_3sd_flat_bounding_ellipses': True,
					'do_plot_4sd_flat_bounding_ellipses': True,
					'do_plot_3sd_flat_bounding_ellipses_dashed': True,
					'do_plot_3sd_angled_bounding_ellipses': True,
					'do_plot_4sd_angled_bounding_ellipses': True,
					'do_plot_3sd_angled_bounding_ellipses_dashed': True,
					'do_plot_fa_line_and_3sd_ellipse': True,
					'do_plot_fa_line_and_4sd_ellipse': True,
					'do_plot_extra': False,
					'is_mainz': True,
					'create_pdf': True,
					'create_video': False,
					'optional_fixed_lims': [5.0, 5.0 + analysis['window_length']]

				},
				{
					'type_string': 'type_1',
					'do_plot_3sd_conjunctive_ring': False,
					'do_plot_4sd_conjunctive_ring': True,
					'do_plot_3sd_flat_bounding_ellipses': False,
					'do_plot_4sd_flat_bounding_ellipses': True,
					'do_plot_3sd_flat_bounding_ellipses_dashed': False,
					'do_plot_3sd_angled_bounding_ellipses': False,
					'do_plot_4sd_angled_bounding_ellipses': True,
					'do_plot_3sd_angled_bounding_ellipses_dashed': False,
					'do_plot_fa_line_and_3sd_ellipse': False,
					'do_plot_fa_line_and_4sd_ellipse': True,
					'do_plot_extra': False,
					'is_mainz': True,
					'create_pdf': False,
					'create_video': True,
					'optional_fixed_lims': [5.0, 5.0 + analysis['window_length']]

				},
				{
					'type_string': 'type_2',
					'do_plot_3sd_conjunctive_ring': False,
					'do_plot_4sd_conjunctive_ring': True,
					'do_plot_3sd_flat_bounding_ellipses': False,
					'do_plot_4sd_flat_bounding_ellipses': True,
					'do_plot_3sd_flat_bounding_ellipses_dashed': False,
					'do_plot_3sd_angled_bounding_ellipses': False,
					'do_plot_4sd_angled_bounding_ellipses': True,
					'do_plot_3sd_angled_bounding_ellipses_dashed': False,
					'do_plot_fa_line_and_3sd_ellipse': False,
					'do_plot_fa_line_and_4sd_ellipse': True,
					'do_plot_extra': True,
					'is_mainz': True,
					'create_pdf': True,
					'create_video': True,
					'optional_fixed_lims': [5.0, 5.0 + analysis['window_length']]

				},
			]


			analysis['all_tw_files_by_shuffle_type_and_stwcboo'] = {}
			for shuffle_option in analysis['shuffle_options']:
				shuffle_option_string = shuffle_option[0] + str(shuffle_option[1])
				analysis['all_tw_files_by_shuffle_type_and_stwcboo'][shuffle_option_string] = [[] for stwcboo in analysis['super_time_warp_cluster_buster_option_objects']]

			analysis['directory_holder'] = sw.DirectoriesHolder()
			analysis['directory_holder'].collated_root_output_directory = analysis['neuron_pair_figs_dir'] + analysis_name + "/"
			analysis['directory_holder'].default_initialisation()
			analysis['directory_holder'].create_directories()

			with open(analysis['directory_holder'].collated_root_output_directory + 'directory_holder.pickle', 'wb') as handle:
				pickle.dump(analysis['directory_holder'], handle, protocol=pickle.HIGHEST_PROTOCOL)


class DirectoriesHolder(object):

	"""
	Class for managing and creating main Pairwise analysis subdirectories
	"""

	def __init__(self):

		self.collated_root_output_directory = ''
		self.basic_scatter_plot_directory = 'BasicScatterPlots/'
		self.pair_cluster_plots_directory = 'PairClusterPlots/'
		self.pair_cluster_stats_comparison_directory = 'PairClusterStatsComparison/'
		self.standard_error_comparison_directory = 'StandardErrorComparisons/'
		# self.fano_factor_directory = 'FanoFactor/'
		self.pairwise_reliabilities_directory = 'PairwiseReliabilities/'
		self.unit_type_analysis_directory = 'UnitTypeAnalysis/'
		self.angle_analysis_directory = 'AngleAnalysis/'
		self.clus_non_stationarity_dir = 'ClusterNonStationarity/'
		self.clus_time_spans_dir = 'ClusterTimeSpans/'
		self.prim_clus_corr_dir = 'PrimaryClusterCorrelations/'
		self.sec_clus_corr_dir = 'SecondaryClusterCorrelations/'
		self.clus_pair_differences_dir = 'ClusterPairDifferencesAnalysis/'
		self.fa_varaince_comparisons_directory = "FAVarianceComparisons/"
		self.cluster_reliabilities_dir = "ClusterReliabilities/"


		self.current_shuffle_option = ''


	def default_initialisation(self):

		self.basic_scatter_plot_directory = self.collated_root_output_directory + self.basic_scatter_plot_directory
		self.pair_cluster_plots_directory = self.collated_root_output_directory + self.pair_cluster_plots_directory
		self.pair_cluster_stats_comparison_directory = self.collated_root_output_directory + self.pair_cluster_stats_comparison_directory
		self.pairwise_reliabilities_directory = self.collated_root_output_directory + self.pairwise_reliabilities_directory
		self.unit_type_analysis_directory = self.collated_root_output_directory + self.unit_type_analysis_directory
		self.angle_analysis_directory = self.collated_root_output_directory + self.angle_analysis_directory
		self.clus_non_stationarity_dir = self.collated_root_output_directory + self.clus_non_stationarity_dir
		self.clus_time_spans_dir = self.collated_root_output_directory + self.clus_time_spans_dir
		self.prim_clus_corr_dir = self.collated_root_output_directory + self.prim_clus_corr_dir
		self.sec_clus_corr_dir = self.collated_root_output_directory + self.sec_clus_corr_dir
		self.clus_pair_differences_dir = self.collated_root_output_directory + self.clus_pair_differences_dir
		self.fa_varaince_comparisons_directory = self.collated_root_output_directory + self.fa_varaince_comparisons_directory
		self.cluster_reliabilities_dir = self.collated_root_output_directory + self.cluster_reliabilities_dir



	def create_directories(self):

		sw.makedirs(self.collated_root_output_directory)
		sw.makedirs(self.pair_cluster_plots_directory)
		sw.makedirs(self.basic_scatter_plot_directory)
		sw.makedirs(self.pair_cluster_stats_comparison_directory)
		sw.makedirs(self.unit_type_analysis_directory)
		sw.makedirs(self.angle_analysis_directory)
		sw.makedirs(self.clus_non_stationarity_dir)
		sw.makedirs(self.clus_time_spans_dir)
		sw.makedirs(self.prim_clus_corr_dir)
		sw.makedirs(self.sec_clus_corr_dir)
		sw.makedirs(self.clus_pair_differences_dir)
		sw.makedirs(self.fa_varaince_comparisons_directory)
		sw.makedirs(self.cluster_reliabilities_dir)


		if (self.pairwise_reliabilities_directory != ''):
			self.EXPERIMENT_pairwise_reliabilities_directory = self.pairwise_reliabilities_directory + "ByExperiment/"
			sw.makedirs(self.pairwise_reliabilities_directory)
			sw.makedirs(self.EXPERIMENT_pairwise_reliabilities_directory)



