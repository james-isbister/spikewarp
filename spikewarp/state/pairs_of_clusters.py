import itertools
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from sklearn.decomposition import PCA
from factor_analyzer import FactorAnalyzer
from sklearn.decomposition import FactorAnalysis
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from scipy.stats import t
import math

import spikewarp as sw


def process_sch_combinations(schs_for_cond, exp_cond_units_holder, exp_cond_str, meta_dict, dh, plot_opt):

	"""
	Function for running analyses of pairs of clusters for a single experimental condition
	"""

	combinations_of_condition_test_schs = list(itertools.combinations(schs_for_cond, r=2))
	
	for sch_index, sch in enumerate(schs_for_cond):
		sch.sch_index = sch_index

	placeholder_1 = None
	if hasattr(sw, 'Placeholder1'):
		placeholder_1 = sw.Placeholder1(exp_cond_units_holder.num_units, schs_for_cond)

	for [sch_0, sch_1] in combinations_of_condition_test_schs:

		if (sch_0.unit_indices != sch_1.unit_indices):

			process(sch_0, sch_1, exp_cond_str, meta_dict, dh, placeholder_1, make_single_quad_vid=plot_opt['make_single_quad_vid'], plot_cluster_pair_multi_stats=plot_opt['plot_cluster_pair_multi_stats'], plot_pairwise_dewarp=['plot_pairwise_dewarp'])
			
			if (placeholder_1 != None):
				placeholder_1.update(sch_0, sch_1)

	if (placeholder_1 != None):
		placeholder_2(exp_cond_units_holder, combinations_of_condition_test_schs, schs_for_cond, meta_dict, plot_opt,  dh, exp_cond_str, state_diff_mat, state_angle_mat, placeholder_1)


def process(sch_0, sch_1, exp_cond_str, meta_dict, dh, placeholder_1, iteration="Original", states_to_use_string="normalised_predicted_states", make_single_quad_vid=False, plot_cluster_pair_multi_stats=False, plot_pairwise_dewarp=False):

	"""
	Main function for analysis of a pair of clusters
	"""

	exp_cond_and_sch_pair_units_str = exp_cond_str + "_" + sch_0.unit_pair_str + "_" + sch_1.unit_pair_str

	conj_indices_0 = np.nonzero(np.in1d(sch_0.trials, sch_1.trials))[0]
	conj_indices_1 = np.nonzero(np.in1d(sch_1.trials, sch_0.trials))[0]

	num_conj_trials = len(conj_indices_0)

	quad_draw_dict = {}

	if (iteration=="Original"):
		meta_dict["all_clusters"]['proportion_of_pairwise_conj_trials_by_single_unit_ALL'].extend([sch_0.proportion_of_conjunctive_trials_in_cluster, sch_1.proportion_of_conjunctive_trials_in_cluster])

	if (num_conj_trials > 15):

		conj_state_0 = getattr(sch_0, states_to_use_string)[conj_indices_0]
		conj_state_1 = getattr(sch_1, states_to_use_string)[conj_indices_1]
		spikes_to_use_0_chosen = sch_0.spikes[conj_indices_0, :]
		spikes_to_use_1_chosen = sch_1.spikes[conj_indices_1, :]
		differences_0 = spikes_to_use_0_chosen[:, 0] - spikes_to_use_0_chosen[:, 1]
		differences_1 = spikes_to_use_1_chosen[:, 0] - spikes_to_use_1_chosen[:, 1]
		
		cluster_pair_type = "quadruplet"
		if ((sch_0.unit_indices[0] in sch_1.unit_indices) | (sch_0.unit_indices[1] in sch_1.unit_indices)):
			cluster_pair_type = "triplet"

		state_lr = linregress(conj_state_0, conj_state_1)
		
		ld = {}
		if (placeholder_1 != None):
			ld = placeholder_1.update_1(spikes_to_use_0_chosen, spikes_to_use_1_chosen, conj_state_0, conj_state_1, cluster_pair_type, meta_dict, sch_0, sch_1)

		if (plot_cluster_pair_multi_stats):
			cluster_pair_multi_stat_plot(sch_0, sch_1, conj_state_0, conj_state_1, state_lr, num_conj_trials, ld, spikes_to_use_0_chosen, spikes_to_use_1_chosen, dh, cluster_pair_type, differences_0, differences_1, iteration, exp_cond_and_sch_pair_units_str)

			

		if (iteration=="Original"):
			meta_dict[cluster_pair_type]['cluster_pairs_state_corr_pvals'].append(state_lr.pvalue)
			meta_dict[cluster_pair_type]['cluster_pairs_state_corr_rvals'].append(state_lr.rvalue)
			meta_dict[cluster_pair_type]['cluster_pairs_state_corr_r2vals'].append(state_lr.rvalue**2)
			meta_dict[cluster_pair_type]['proportion_of_pairwise_conj_trials_by_single_unit_ENOUGH_QUAD_CONJ'].extend([sch_0.proportion_of_conjunctive_trials_in_cluster, sch_1.proportion_of_conjunctive_trials_in_cluster])
			meta_dict[cluster_pair_type]['quadruplet_keys_ENOUGH_QUAD_CONJ'].append(str(min([sch_0.unique_pair_str, sch_1.unique_pair_str])) + "_" + str(max([sch_0.unique_pair_str, sch_1.unique_pair_str])))

			if (state_lr.pvalue < 0.05):
				meta_dict[cluster_pair_type]['proportion_of_pairwise_conj_trials_by_single_unit_ENOUGH_QUAD_CONJ_CORR'].extend([sch_0.proportion_of_conjunctive_trials_in_cluster, sch_1.proportion_of_conjunctive_trials_in_cluster])
				meta_dict[cluster_pair_type]['quadruplet_keys_ENOUGH_QUAD_CONJ_CORR'].append(str(min([sch_0.unique_pair_str, sch_1.unique_pair_str])) + "_" + str(max([sch_0.unique_pair_str, sch_1.unique_pair_str])))

			if ((cluster_pair_type == "quadruplet") & (hasattr(sw, 'placeholder_3'))):
				sw.placeholder_3(conj_state_0, conj_state_1, spikes_to_use_0_chosen, spikes_to_use_1_chosen, dh, exp_cond_and_sch_pair_units_str, sch_0, sch_1, num_conj_trials, make_single_quad_vid)


		if (iteration=="Original"):						
			if (cluster_pair_type == 'triplet'):

				sch_0_unit_ind_to_use = 0
				if (sch_0.unit_indices[0] in sch_1.unit_indices):
					sch_0_unit_ind_to_use = 1
				all_3_spikes = np.hstack((np.asarray([spikes_to_use_0_chosen[:,sch_0_unit_ind_to_use].tolist()]).T, spikes_to_use_1_chosen))
				pca = PCA(n_components=3); pca.fit_transform(all_3_spikes)

			if (cluster_pair_type == 'quadruplet'):
			
				all_4_spikes = np.hstack((spikes_to_use_0_chosen, spikes_to_use_1_chosen))
				pca = PCA(n_components=4); pca.fit_transform(all_4_spikes)				

			meta_dict[cluster_pair_type]['explained_variance_ratio_by_PCA_component'].append(pca.explained_variance_ratio_) 



		if ((iteration=="Original") & (cluster_pair_type == 'quadruplet')):
			# Differences
			ds = []
			if (sch_0.is_PCA_BS_empirical_pvalue_different_from_45):
				ds.append(sw.times_state_prediction(conj_state_1, differences_0, FA_predicted_sd=sch_0.predicted_state_dependent_differences_sd, conj_orig_sd=sch_0.Conjunctive_original_FS_diff_SD, plot_dewarp_and_state_dependence=plot_pairwise_dewarp, plot_path=dh.single_pairwise_diff_dewarp_fig_dir + exp_cond_and_sch_pair_units_str + "_sch0"))
			if (sch_1.is_PCA_BS_empirical_pvalue_different_from_45):
				ds.append(sw.times_state_prediction(conj_state_0, differences_1, FA_predicted_sd=sch_1.predicted_state_dependent_differences_sd, conj_orig_sd=sch_1.Conjunctive_original_FS_diff_SD, plot_dewarp_and_state_dependence=plot_pairwise_dewarp, plot_path=dh.single_pairwise_diff_dewarp_fig_dir + exp_cond_and_sch_pair_units_str+ "_sch1"))
			for d in ds:

				meta_dict[cluster_pair_type]['difference_prediction_from_state_pvalues'].append(d['p'])
				meta_dict[cluster_pair_type]['difference_prediction_from_state_rvalues'].append(d['r'])
				meta_dict[cluster_pair_type]['difference_prediction_from_state_rsquared_values'].append(d['rsquared'])
				if (d['p'] < 0.05):
					meta_dict[cluster_pair_type]['difference_prediction_from_state_FA_predicted_sds'].append(d['FA_predicted_sd'])
					meta_dict[cluster_pair_type]['difference_prediction_from_state_conj_orig_sds'].append(d['conj_orig_sd'])
					meta_dict[cluster_pair_type]['difference_prediction_from_state_cluster_orig_sds'].append(d['orig_sd'])
					meta_dict[cluster_pair_type]['difference_prediction_from_state_state_predic_error_sds'].append(d['state_predic_error_sd'])


			# Spike times
			ds = []
			for unit_index in range(2):
				ds.append(sw.times_state_prediction(conj_state_1, spikes_to_use_0_chosen[:, unit_index], sch=sch_0, FA_predicted_sd=sch_0.FA_N_sds[unit_index], conj_orig_sd=sch_0.Conj_original_N_sds[unit_index], plot_dewarp_and_state_dependence=plot_pairwise_dewarp, plot_path=dh.single_pairwise_st_dewarp_fig_dir + exp_cond_and_sch_pair_units_str + "_sch0_u" + str(unit_index)))
				ds.append(sw.times_state_prediction(conj_state_0, spikes_to_use_1_chosen[:, unit_index], sch=sch_1, FA_predicted_sd=sch_1.FA_N_sds[unit_index], conj_orig_sd=sch_1.Conj_original_N_sds[unit_index], plot_dewarp_and_state_dependence=plot_pairwise_dewarp, plot_path=dh.single_pairwise_st_dewarp_fig_dir + exp_cond_and_sch_pair_units_str + "_sch1_u" + str(unit_index)))

				for d in ds:
					meta_dict[cluster_pair_type]['sch_single_unit_prediction_from_state_pvalues'].append(d['p'])
					
					if (d['p'] < 0.05):

						meta_dict[cluster_pair_type]['sch_single_unit_prediction_from_state_rvalues'].append(d['r'])
						meta_dict[cluster_pair_type]['sch_single_unit_prediction_from_state_rsquared_values'].append(d['rsquared'])
						meta_dict[cluster_pair_type]['sch_single_unit_prediction_from_state_FA_predicted_sds'].append(d['FA_predicted_sd'])
						meta_dict[cluster_pair_type]['sch_single_unit_prediction_from_state_conj_orig_sds'].append(d['conj_orig_sd'])
						meta_dict[cluster_pair_type]['sch_single_unit_prediction_from_state_cluster_orig_sds'].append(d['orig_sd'])
						meta_dict[cluster_pair_type]['sch_single_unit_prediction_from_state_state_predic_error_sds'].append(d['state_predic_error_sd'])



def cluster_pair_multi_stat_plot(sch_0, sch_1, conj_state_0, conj_state_1, state_lr, num_conj_trials, ld, spikes_to_use_0_chosen, spikes_to_use_1_chosen, dh, cluster_pair_type, differences_0, differences_1, iteration, exp_cond_and_sch_pair_units_str):

	"""
	Function for plotting analyses of a pair of clusters
	"""

	fig, axes = plt.subplots(1, 8, figsize=(20, 3)) 
	sns.set()
	sns.set_style("ticks")

	axes[0].scatter(sch_0.spikes[:, 0], sch_0.spikes[:, 1], c=sch_0.predicted_states)
	# lims = [np.min(sch_0.spikes) - 2.0, np.max(sch_0.spikes) + 2.0]
	lims = [5.0, 50.0]
	axes[0].set_aspect('equal', 'box')
	ticks = np.arange(5.0, 55.0, 5.0)
	axes[0].set_yticks(ticks)
	axes[0].set_xticks(ticks)
	axes[0].set_xlim(lims); axes[0].set_ylim(lims)
	axes[0].set_xlabel('Time (ms)'); axes[0].set_ylabel('Time (ms)')
	
	axes[1].scatter(sch_1.spikes[:, 0], sch_1.spikes[:, 1], c=sch_1.predicted_states)
	# lims = [np.min(sch_1.spikes) - 2.0, np.max(sch_1.spikes) + 2.0]
	lims = [5.0, 50.0]
	axes[1].set_aspect('equal', 'box')
	ticks = np.arange(5.0, 55.0, 5.0)
	axes[1].set_yticks(ticks)
	axes[1].set_xticks(ticks)
	axes[1].set_xlim(lims); axes[1].set_ylim(lims)
	axes[1].set_xlabel('Time (ms)'); axes[1].set_ylabel('Time (ms)')
	
	axes[2].scatter(conj_state_0, conj_state_1)
	lims = [np.min(np.hstack((conj_state_0, conj_state_1))) - 2.0, np.max(np.hstack((conj_state_0, conj_state_1))) + 2.0]
	# axes[2].set_xlim(lims); axes[2].set_ylim(lims)
	# print(math.isnan(state_lr.slope))
	subtitle = 'State corr - p: ' + str(sw.round_to_2(state_lr.pvalue)) + ', r: ' + str(sw.round_to_2(state_lr.rvalue)) + ', num: ' + str(num_conj_trials)
	if (ld != {}):
		subtitle += f", PCA Angle: {ld['BS_PCA_mean_angle']:.6f}"

	plt.suptitle(subtitle)
	axes[2].set_aspect('equal', 'box')
	ticks = np.arange(-3.0, 3.5, 1.0)
	axes[2].set_yticks(ticks)
	axes[2].set_xticks(ticks)
	axes[2].set_xlim([-3.5, 3.5]); axes[2].set_ylim([-3.5, 3.5])
	x = np.linspace(-3,3,10)
	axes[2].plot(x, x, '--k', linewidth=2)
	axes[2].set_xlabel('Predicted state\ncluster 1'); axes[2].set_ylabel('Predicted state\ncluster 2')

	axes[3].scatter(differences_0, differences_1)

	axes[4].scatter(spikes_to_use_0_chosen[:, 0], spikes_to_use_1_chosen[:, 0])
	axes[4].set_title(str(sch_0.unit_indices[0]) + ', ' + str(sch_1.unit_indices[0]))
	axes[5].scatter(spikes_to_use_0_chosen[:, 0], spikes_to_use_1_chosen[:, 1])
	axes[5].set_title(str(sch_0.unit_indices[0]) + ', ' + str(sch_1.unit_indices[1]))
	axes[6].scatter(spikes_to_use_0_chosen[:, 1], spikes_to_use_1_chosen[:, 0])
	axes[6].set_title(str(sch_0.unit_indices[1]) + ', ' + str(sch_1.unit_indices[0]))
	axes[7].scatter(spikes_to_use_0_chosen[:, 1], spikes_to_use_1_chosen[:, 1])
	axes[7].set_title(str(sch_0.unit_indices[1]) + ', ' + str(sch_1.unit_indices[1]))

	# fig.tight_layout()
	ftypes = ['pdf', 'png']
	for ftype in ftypes:
		plt.savefig(dh.pairs_of_clusters_state_comp_fig_dir + cluster_pair_type + "TimesAndState" + iteration + exp_cond_and_sch_pair_units_str + '.' + ftype, bbox_inches='tight') #
	plt.close()


