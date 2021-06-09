import matplotlib.pyplot as plt
import numpy as np

import spikewarp as sw

def cluster_pair_plots(meta_dict, dh):

	"""
	Function for plotting meta analyses of cluster pair analyses
	"""


	data_to_plot = [meta_dict['quadruplet']['proportion_of_pairwise_conj_trials_by_single_unit_ENOUGH_QUAD_CONJ'], meta_dict['all_clusters']['proportion_of_pairwise_conj_trials_by_single_unit_ALL']]
	sw.normal_histo_plot(data_to_plot, dh.state_corr_meta_stats + 'proportion_of_pairwise_conj_trials_by_single_unit_ENOUGH_QUAD_CONJ_FOR_TEST_AND_ALL', bins=20, histo_range=[0.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="Proportion", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[5, 5], alpha=0.78, labels=['Quadruplets', "All"])

	data_to_plot = [meta_dict['quadruplet']['proportion_of_pairwise_conj_trials_by_single_unit_ENOUGH_QUAD_CONJ'], meta_dict['quadruplet']['proportion_of_pairwise_conj_trials_by_single_unit_ENOUGH_QUAD_CONJ_CORR']]
	sw.normal_histo_plot(data_to_plot, dh.state_corr_meta_stats + 'proportion_of_pairwise_conj_trials_by_single_unit_ENOUGH_QUAD_CONJ_ENOUGH_QUAD_CONJ_FOR_TEST_AND_CORR', bins=20, histo_range=[0.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="Proportion", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[5, 5], alpha=0.78, labels=['Quadruplets enough conj for test', 'Quadruplets enough conj for test corr'])
	


	############## EXPLAINED VARIANCE ##############
	# TRIPLETS
	exp_var_ratios = np.asarray(meta_dict['triplet']['explained_variance_ratio_by_PCA_component'])
	if (exp_var_ratios != []):
		plt.figure()
		plt.scatter(exp_var_ratios[:, 0], exp_var_ratios[:, 1])
		plt.gca().set_xlim([0.0, 1.0]); plt.gca().set_ylim([0.0, 1.0])
		plt.savefig(dh.pca_explained_var_meta_stats + 'exp_var_ratios_triplets_pca_comp_0_and_1.pdf')
		plt.close()

	plt.figure()
	plt.bar(['1st', '2nd', '3rd'], np.mean(exp_var_ratios, axis=0))
	plt.savefig(dh.pca_explained_var_meta_stats + 'mean_exp_var_ratios_triplets.pdf')
	plt.close()

	# QUADRUPLETS
	exp_var_ratios = np.asarray(meta_dict['quadruplet']['explained_variance_ratio_by_PCA_component'])
	if (exp_var_ratios != []):
		plt.figure()
		plt.scatter(exp_var_ratios[:, 0], exp_var_ratios[:, 1])
		plt.gca().set_xlim([0.0, 1.0]); plt.gca().set_ylim([0.0, 1.0])
		plt.savefig(dh.pca_explained_var_meta_stats + 'exp_var_ratios_quadruplets_pca_comp_0_and_1.pdf')
		plt.close()

		plt.figure()
		plt.scatter(exp_var_ratios[:, 2], exp_var_ratios[:, 3])
		plt.gca().set_xlim([0.0, 1.0]); plt.gca().set_ylim([0.0, 1.0])
		plt.savefig(dh.pca_explained_var_meta_stats + 'exp_var_ratios_quadruplets_pca_comp_2_and_3.pdf')
		plt.close()

	plt.figure()
	plt.bar(['1st', '2nd', '3rd', '4th'], np.mean(exp_var_ratios, axis=0))
	plt.gca().set_xlabel('PCA component'); plt.gca().set_ylabel('Proportion of\nexplained variance');
	plt.savefig(dh.pca_explained_var_meta_stats + 'mean_exp_var_ratios_quadruplets.pdf')
	plt.close()


	############## STATE CORR ##############
	sw.normal_histo_plot([meta_dict['quadruplet']['cluster_pairs_state_corr_pvals']], dh.state_corr_meta_stats + 'cluster_pairs_state_corr_pvals', bins=20, histo_range=[0.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$p value$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[5, 5], alpha=0.78)
	sw.normal_histo_plot([np.asarray(meta_dict['quadruplet']['cluster_pairs_state_corr_rvals'])[np.argwhere(np.asarray(meta_dict['quadruplet']['cluster_pairs_state_corr_pvals']) < 0.05).flatten()]], dh.state_corr_meta_stats + 'cluster_pairs_state_corr_rvals_quadruplets', bins=40, histo_range=[-1.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$r$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[5, 5], alpha=0.78)
	sw.normal_histo_plot([np.asarray(meta_dict['quadruplet']['cluster_pairs_state_corr_r2vals'])[np.argwhere(np.asarray(meta_dict['quadruplet']['cluster_pairs_state_corr_pvals']) < 0.05).flatten()]], dh.state_corr_meta_stats + 'cluster_pairs_state_corr_r2vals_quadruplets', bins=20, histo_range=[0.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$r^2$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[5, 5], alpha=0.78)

	############## DIFFERENCES PREDICTION AND CORR ##############
	sw.normal_histo_plot([meta_dict['quadruplet']['difference_prediction_from_state_pvalues']], dh.diff_corr_meta_stats + 'difference_prediction_from_state_pvalues', bins=20, histo_range=[0.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$p value$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[5, 5], alpha=0.78)
	sw.normal_histo_plot([np.asarray(meta_dict['quadruplet']['difference_prediction_from_state_rvalues'])[np.argwhere(np.asarray(meta_dict['quadruplet']['difference_prediction_from_state_pvalues']) < 0.05).flatten()]], dh.diff_corr_meta_stats + 'difference_prediction_from_state_rvalues', bins=40, histo_range=[-1.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$r$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[5, 5], alpha=0.78)
	sw.normal_histo_plot([np.asarray(meta_dict['quadruplet']['difference_prediction_from_state_rsquared_values'])[np.argwhere(np.asarray(meta_dict['quadruplet']['difference_prediction_from_state_pvalues']) < 0.05).flatten()]], dh.diff_corr_meta_stats + 'difference_prediction_from_state_rsquared_values', bins=20, histo_range=[0.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$r^2$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[5, 5], alpha=0.78)

	############## SINGLE UNIT PREDICTION AND CORR ##############
	# Already p-value thresholded
	sw.normal_histo_plot([meta_dict['quadruplet']['sch_single_unit_prediction_from_state_pvalues']], dh.state_corr_meta_stats + 'sch_single_unit_prediction_from_state_pvalues', bins=20, histo_range=[0.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$p value$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[100, 100], alpha=0.78)
	sw.normal_histo_plot([np.asarray(meta_dict['quadruplet']['sch_single_unit_prediction_from_state_rvalues'])], dh.state_corr_meta_stats + 'sch_single_unit_prediction_from_state_rvalues', bins=40, histo_range=[-1.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$r$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[5, 5], alpha=0.78)
	sw.normal_histo_plot([np.asarray(meta_dict['quadruplet']['sch_single_unit_prediction_from_state_rsquared_values'])], dh.state_corr_meta_stats + 'sch_single_unit_prediction_from_state_rsquared_values', bins=20, histo_range=[0.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$r^2$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[5, 5], alpha=0.78)



	scatter_point_color_groups = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

	sw.plot_fa_comparison([meta_dict['quadruplet']['difference_prediction_from_state_cluster_orig_sds']], 
						[meta_dict['quadruplet']['difference_prediction_from_state_state_predic_error_sds']], 
						dh.diff_corr_meta_stats + "_ClusterSD_vs_StatePredicErrorSD",
						y_equals_x_max=30, 
						optional_y_max=30,
						scatter_point_color_groups=scatter_point_color_groups)

	sw.plot_fa_comparison([meta_dict['quadruplet']['difference_prediction_from_state_FA_predicted_sds']], 
						[meta_dict['quadruplet']['difference_prediction_from_state_state_predic_error_sds']], 
						dh.diff_corr_meta_stats + "_FAPredicSD_vs_StatePredicErrorSD",
						y_equals_x_max=14, 
						optional_y_max=14,
						scatter_point_color_groups=scatter_point_color_groups)

	sw.plot_fa_comparison([meta_dict['quadruplet']['difference_prediction_from_state_conj_orig_sds']], 
						[meta_dict['quadruplet']['difference_prediction_from_state_state_predic_error_sds']], 
						dh.diff_corr_meta_stats + "_ConjSD_vs_StatePredicErrorSD",
						y_equals_x_max=25, 
						optional_y_max=25,
						scatter_point_color_groups=scatter_point_color_groups)

	sw.plot_fa_comparison([meta_dict['quadruplet']['sch_single_unit_prediction_from_state_cluster_orig_sds']], 
						[meta_dict['quadruplet']['sch_single_unit_prediction_from_state_state_predic_error_sds']], 
						dh.state_corr_meta_stats + "_ClusterSD_vs_StatePredicErrorSD",
						y_equals_x_max=30, 
						optional_y_max=30,
						scatter_point_color_groups=scatter_point_color_groups)

	sw.plot_fa_comparison([meta_dict['quadruplet']['sch_single_unit_prediction_from_state_FA_predicted_sds']], 
						[meta_dict['quadruplet']['sch_single_unit_prediction_from_state_state_predic_error_sds']], 
						dh.state_corr_meta_stats + "_FAPredicSD_vs_StatePredicErrorSD",
						y_equals_x_max=14, 
						optional_y_max=14,
						scatter_point_color_groups=scatter_point_color_groups)

	sw.plot_fa_comparison([meta_dict['quadruplet']['sch_single_unit_prediction_from_state_conj_orig_sds']], 
						[meta_dict['quadruplet']['sch_single_unit_prediction_from_state_state_predic_error_sds']], 
						dh.state_corr_meta_stats + "_ConjSD_vs_StatePredicErrorSD",
						y_equals_x_max=25, 
						optional_y_max=25,
						scatter_point_color_groups=scatter_point_color_groups)


	
