'''
1_single_and_pairwise_analysis.py

Script for producing all figures up to Fig. 7a (exlcuding Figs. 3d & 4c) and all supplementary figures.



### Configuration

1. Update the variable "original_data_root" on line 2 of "recipes/investigation_recipe.json" to the path at which the in vivo data is stored. For example:
> "original_data_root": "/Users/james/IsbisterEtAl-SciRep2021-Data/",

2. Update the variable "figures_root_directory" on line 3 of "recipes/investigation_recipe.json" to the directory you would like the output of the analyses to be stored. For example:
> "figures_root_directory": "/Users/james/IsbisterEtAl-SciRep2021-Figures/",

3. Update the variable "number_of_processors" on line 4 of "recipes/investigation_recipe.json" to the number of threads you would like to use (initially set to 8).



### Options

1. The two main analyses (CustomDBSCANExtra and Unclustered), which are neccessary for the 2nd script to be run, are initially selected on line 142 of "recipes/investigation_recipe.json". 

2. For additional speed, the clustering controls can be removed by changing the value of the "shuffle_options" parameters on lines 74 and 132 of "recipes/investigation_recipe.json" to [["normal", -1]]. 

3. All additional analyses (including the alternative clustering algorithms) can be run by changing the value of the "analyses_to_process" parameter on line 142 to include keys from the line above (141).

4. To run the analyses for a subset of the data (i.e. to test that the end-to-end analysis pipeline works), a smaller subset of experimental sessions can be selected on line 8 of "recipes/investigation_recipe.json". Experiment 0 produces a significant number of clusters for example.




### Run

The script can be run with the following command:
python 1_single_and_pairwise_analysis.py

'''


import pickle
import sys
import multiprocessing

import spikewarp as sw


print(multiprocessing.cpu_count())

if (len(sys.argv) > 1):
	investigation = sw.load_json_recipe(sys.argv[1])
else:
	investigation = sw.load_json_recipe("recipes/investigation_recipe.json")

sw.prepare_investigation(investigation)
data_loaders, investigation['determined_cortical_onset_timepoint'] = sw.load_spikes_calc_cort_onset_and_draw_PSTHs(investigation)
sw.prepare_investigation_analysis_recipes(investigation)


for analysis_name in investigation['analyses_to_process']:

	print(analysis_name)

	analysis = investigation['analyses'][analysis_name]

	for exp_ind in investigation["experiments_to_process"]:
		data_loaders[investigation["experiments_to_process"].index(exp_ind)].create_scsuos(spike_train_start_time=analysis['main_analysis_spike_train_start_time'],
																							spike_train_end_time=analysis['main_analysis_spike_train_end_time'], 
																							custom_number_of_trials_for_condition=analysis['custom_number_of_trials_for_condition'],
																							do_single_unit_stationarity_tests=analysis['do_single_unit_stationarity_tests'])



	if (analysis['draw_spike_count_plots']):
		sw.draw_spike_count_plots(analysis['single_neuron_analysis_figs_dir'], data_loaders, investigation)


	if (analysis['do_single_unit_stationarity_tests']):
		sw.single_unit_stationarity_tests(analysis['single_neuron_analysis_figs_dir'], analysis['draw_individual_unit_plots'], data_loaders, investigation)


	if (analysis['run_pairwise_analysis']):

		mcah_for_each_shuffle_type = {shuffle_option_string:sw.MetaClusterAnalysisHolder(shuffle_option_string) for shuffle_option_string in analysis['option_helper'].shuffle_option_strings}
		random_seed_counter = 0
		single_cluster_holders_by_exp_cond_and_shuffle_option = {}

		for exp_ind in investigation["experiments_to_process"]:

			# Analyse pairwise distributions
			mcah_for_each_shuffle_type_single_experiment, \
			random_seed_counter, \
			tw_files_by_stwcboo_and_shuffle_option, \
			single_cluster_holders_by_condition_and_shuffle_option = sw.run_pairwise_analysis(random_seed_counter,
																						str(exp_ind) + "_" + investigation["experiment_codes"][exp_ind],
																						data_loaders[investigation["experiments_to_process"].index(exp_ind)], 
																						analysis)
																						# [6])
		
			# Pickle for later analysis
			single_cluster_holders_by_exp_cond_and_shuffle_option[str(exp_ind)] = single_cluster_holders_by_condition_and_shuffle_option
			sw.add_appropriate_pairwise_response_dist_files(analysis, tw_files_by_stwcboo_and_shuffle_option)
			
							
			for shuffle_option_string in analysis['option_helper'].shuffle_option_strings:

				mcah_for_shuffle_type = mcah_for_each_shuffle_type[shuffle_option_string]
				mcah_for_shuffle_type.extend_standard_cluster_arrays_using_another_mcah(mcah_for_each_shuffle_type_single_experiment[shuffle_option_string])
				
				if ((shuffle_option_string == 'normal-1') & (analysis['do_long_meta_analyses'])):
					mcah_for_shuffle_type.extend_pairwise_reliabilities(mcah_for_each_shuffle_type_single_experiment[shuffle_option_string].all_both_spiking_reliabilities, 
																		mcah_for_each_shuffle_type_single_experiment[shuffle_option_string].all_number_of_conjunctive_trials)



		sw.pickle_clusters(analysis, single_cluster_holders_by_exp_cond_and_shuffle_option); sw.pickle_directory_holders(analysis)
		sw.create_pairwise_response_dist_gifs_when_appropriate(analysis)
		
		for shuffle_option_string in analysis['option_helper'].shuffle_option_strings:

			mcah_for_shuffle_type = mcah_for_each_shuffle_type[shuffle_option_string]

			if (shuffle_option_string == 'normal-1'):

				sw.pickle_mcah(analysis, mcah_for_shuffle_type, shuffle_option_string)
				
				# Run meta analyses
				if analysis['calulate_time_span']:
					mcah_for_shuffle_type.calculate_time_span_info_and_plots(analysis['directory_holder'], investigation['determined_cortical_onset_timepoint'], 50.0, analysis['main_analysis_spike_train_end_time'])
				if analysis['do_long_meta_analyses']:
					mcah_for_shuffle_type.plot_p_value_histos(analysis['directory_holder'], do_extra_plots=False)
					mcah_for_shuffle_type.create_comparison_to_standard_error_plots(analysis['directory_holder'])
					mcah_for_shuffle_type.plot_pairwise_reliability_plots(analysis['directory_holder'].pairwise_reliabilities_directory, 'ALL')

			if (shuffle_option_string in ['rotate_correlatd_clusters_to_45_degrees-1', 'sample_correlated_cluster_pca_ellipse_rotated_to_45-1', '45cluster_with_second_cluster-1']):
				mcah_for_shuffle_type.plot_p_value_histos(analysis['directory_holder'], do_extra_plots=False)

			print(shuffle_option_string + " analysis complete.")

		# Meta meta analyses
		sw.plot_multi_shuffle_type_plots(mcah_for_each_shuffle_type, analysis)
