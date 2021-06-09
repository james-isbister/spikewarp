import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import spikewarp as sw

"""
Class and helpers for main clustering meta analyses
"""

class MetaClusterAnalysisHolder(object):

	def __init__(self, shuffle_option_string, is_mainz=True):

		self.shuffle_option_string = shuffle_option_string
		self.suf = "_" + shuffle_option_string
		self.is_mainz = is_mainz

		self.pdds = {}
		self.sdds = {}
		for data_name in sw.list_of_first_stage_data_names:
			self.pdds.update({data_name: []})
		for data_name in sw.list_of_second_stage_data_names:
			self.sdds.update({data_name: []})

		self.final_angled_cluster_count = 0
		self.did_contribute_atleast_one_final_angled_cluster_count = 0

		self.all_both_spiking_reliabilities = []; self.all_both_spiking_reliabilities_0s_removed = []
		self.all_number_of_conjunctive_trials = []; self.all_number_of_conjunctive_trials_0s_removed = []



	def extend_standard_cluster_arrays(self, single_clustering):
		
		if (single_clustering.do_use_clusters_in_analysis):

			self.final_angled_cluster_count += single_clustering.final_angled_cluster_count
			self.did_contribute_atleast_one_final_angled_cluster_count += single_clustering.was_first_single_clustering_to_pass_for_pair

			for key in single_clustering.primary_data_dicts.keys():
				if (key not in self.pdds.keys()):
					self.pdds[key] = []
				self.pdds[key].extend(single_clustering.primary_data_dicts[key])

			for key in single_clustering.secondary_data_dicts.keys():
				if (key not in self.sdds.keys()):
					self.sdds[key] = []
				self.sdds[key].extend(single_clustering.secondary_data_dicts[key])
				


	def extend_standard_cluster_arrays_using_another_mcah(self, mcah):

		self.final_angled_cluster_count += mcah.final_angled_cluster_count
		self.did_contribute_atleast_one_final_angled_cluster_count += mcah.did_contribute_atleast_one_final_angled_cluster_count

		for key in mcah.pdds.keys():
			if (key not in self.pdds.keys()):
				self.pdds[key] = []
			self.pdds[key].extend(mcah.pdds[key])

		for key in mcah.sdds.keys():
			if (key not in self.sdds.keys()):
				self.sdds[key] = []
			self.sdds[key].extend(mcah.sdds[key])


	def calculate_time_span_info_and_plots(self, directory_holder, cortical_onset, time_window_following_cortical_onset, end_of_spiking_activity):

		sdds = self.sdds
		pdds = self.pdds
		dh = directory_holder
		suf = self.suf

		tex_tag_file_name = dh.collated_root_output_directory + "AnalysisOutputLatexTimeSpan.tex"
		with open(tex_tag_file_name, "w") as tex_file:
			print(f"", file=tex_file)


		# Cluster Time Spans
		sw.basic_x_y_plot([pdds['FlatClusterStats_FlatCluster_FS_Mean0']], [pdds['FlatClusterStats_FlatCluster_FS_Mean1']], dh.clus_time_spans_dir + "PrimaryClusterMeans" + suf, s=4, draw_y_equals_x=True, y_equals_x_max=100, x_axis_label='ms', y_axis_label='ms', scatter_point_color_groups=['g'], custom_x_tick_locators=[50, 10])
		sw.basic_x_y_plot([sdds['FlatClusterStats_FlatCluster_FS_Mean0']], [sdds['FlatClusterStats_FlatCluster_FS_Mean1']], dh.clus_time_spans_dir + "SecondaryClusterMeans" + suf, s=4, draw_y_equals_x=True, y_equals_x_max=100, x_axis_label='ms', y_axis_label='ms', scatter_point_color_groups=['g'], custom_x_tick_locators=[50, 10])
		sw.basic_x_y_plot([2.0*np.hstack((pdds['FlatClusterStats_FlatCluster_N0_FS_SD'], pdds['FlatClusterStats_FlatCluster_N1_FS_SD']))], [np.hstack((pdds['FlatClusterStats_FlatCluster_FS_Mean0'], pdds['FlatClusterStats_FlatCluster_FS_Mean1']))], dh.clus_time_spans_dir + "PrimaryClusterMeans_VS_2sds" + suf, s=4, x_axis_label='ms', y_axis_label='ms', scatter_point_color_groups=['g'], custom_x_tick_locators=[50, 10], opt_x_and_y_max=[40.0, 100.0], y_axis_on_right=False)
		sw.basic_x_y_plot([2.0*np.hstack((sdds['FlatClusterStats_FlatCluster_N0_FS_SD'], sdds['FlatClusterStats_FlatCluster_N1_FS_SD']))], [np.hstack((sdds['FlatClusterStats_FlatCluster_FS_Mean0'], sdds['FlatClusterStats_FlatCluster_FS_Mean1']))], dh.clus_time_spans_dir + "SecondaryClusterMeans_VS_2sds" + suf, s=4, x_axis_label='ms', y_axis_label='ms', scatter_point_color_groups=['g'], custom_x_tick_locators=[50, 10], opt_x_and_y_max=[40.0, 100.0], y_axis_on_right=False)
		secondary_flat_cluster_means = np.hstack((sdds['FlatClusterStats_FlatCluster_FS_Mean0'], sdds['FlatClusterStats_FlatCluster_FS_Mean1']))
		secondary_flat_cluster_pre_limits = secondary_flat_cluster_means - 4.0 * np.hstack((sdds['FlatClusterStats_FlatCluster_N0_FS_SD'], sdds['FlatClusterStats_FlatCluster_N1_FS_SD']))
		secondary_flat_cluster_post_limits = secondary_flat_cluster_means + 4.0 * np.hstack((sdds['FlatClusterStats_FlatCluster_N0_FS_SD'], sdds['FlatClusterStats_FlatCluster_N1_FS_SD']))
		sw.normal_histo_plot([secondary_flat_cluster_post_limits], dh.clus_time_spans_dir + "LimitsOfFlatClustersForAngledClustersOnly" + suf, bins=20, histo_range=[0.0, 100.0], x_axis_label="ms", y_axis_label="Frequency", custom_x_tick_locators=[100.0, 10.0], custom_y_tick_locators=[10.0, 10.0], alpha=0.78, add_chi_squared_text=True)
		
		time_threshold = cortical_onset + time_window_following_cortical_onset
		num_before = np.sum(secondary_flat_cluster_post_limits < time_threshold)
		num_after = np.sum(secondary_flat_cluster_post_limits > time_threshold)
		percent_before = 100.0 * float(num_before) / float(num_after + num_before)
		percent_before_string = "{:.{}f}".format(percent_before, 1)
		

		data_part = percent_before_string + "\\%"
		cluster_time_span_string = "As " + data_part + " of Stage 2 clusters extracted over 90ms following cortical activation onset lied within " + str(int(time_window_following_cortical_onset)) + "ms following onset (Supplementary Fig. 12), analysis was constrained to spikes in the first " + str(int(time_window_following_cortical_onset)) + "ms following activation onset. "
		sw.append_new_tag(data_part, "ClusterTimeSpanSummaryNum", tex_tag_file_name)
		sw.append_new_tag(cluster_time_span_string, "ClusterTimeSpanSummary", tex_tag_file_name)



	def plot_p_value_histos(self, directory_holder, do_extra_plots=False):

		sdds = self.sdds
		pdds = self.pdds
		dh = directory_holder
		suf = self.suf

		plot_all_lag_histograms = False
		if (do_extra_plots):
			plot_all_lag_histograms = True

		tex_tag_file_name = dh.collated_root_output_directory + suf + "AnalysisOutputLatex.tex"
		with open(tex_tag_file_name, "w") as tex_file:
			print(f"", file=tex_file)


		specific_prim_clus_corr_dir = dh.prim_clus_corr_dir + suf + "/"; sw.makedirs(specific_prim_clus_corr_dir)
		specific_sec_clus_corr_dir = dh.sec_clus_corr_dir + suf + "/"; sw.makedirs(specific_sec_clus_corr_dir)
		

		# Cluster Correlations Primary
		sw.normal_histo_plot([pdds['FlatClusterStats_FlatCluster_LR_pvalue']], specific_prim_clus_corr_dir + "PVal_ZoomHist", bins=20, histo_range=[0.0, 0.1], x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[0.1, 0.01], custom_y_tick_locators=[30, 30], alpha=0.78, add_chi_squared_text=True)
		flat_cluster_correlations_chi_squared_table_strings_array = sw.cumulative_histo_plot([pdds['FlatClusterStats_FlatCluster_LR_pvalue']], specific_prim_clus_corr_dir + "PVal_CumHist", bins=200, x_axis_label="p-value", y_axis_label="Normalised\ncumulative sum", custom_x_tick_locators=[1.0, 0.2], add_chi_squared_text=True)
		sw.normal_histo_plot([pdds['FlatClusterStats_FlatCluster_LR_pvalue']], specific_prim_clus_corr_dir + "PVal_LowResHist", bins=40, x_axis_label="p-value", y_axis_label="Frequency", custom_y_tick_locators=[100, 100], alpha=0.78, add_chi_squared_text=True)
		sw.cumulative_histo_plot([pdds['FlatClusterStats_FlatCluster_LR_pvalue']], specific_prim_clus_corr_dir + "LowRes_LowResCumHist", bins=20, x_axis_label="p-value", y_axis_label="Normalised\ncumulative sum", add_chi_squared_text=True)

		if ('FlatClusterStats_FlatCluster_LR_rsquared' in sdds.keys()):

			# Cluster Correlations Secondary
			sw.normal_histo_plot([sdds['FlatClusterStats_FlatCluster_LR_rsquared'], sdds['FlatClusterStats_FlatCluster_LR_rvalue']], specific_sec_clus_corr_dir + "RVal_Hist", bins=40, histo_range=[-1.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$r$, $r^2$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[50, 10], alpha=0.78)
			sw.normal_histo_plot([sdds['FlatClusterStats_FlatCluster_LR_rsquared']], specific_sec_clus_corr_dir + "R^2_Hist", colors=['g'], bins=20, x_axis_left_buffer=0.01, x_axis_label="r^2-value", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[20, 20])
			
			cluster_p_minus_unclustered_conj_p = np.asarray(sdds['FlatClusterStats_FlatCluster_LR_pvalue']) - np.asarray(sdds['Unclustered_Conj_LR_pvalue'])
			num_improved_by_clustering = np.sum(cluster_p_minus_unclustered_conj_p < 0.0)
			num_not_improved_by_clustering = np.sum(cluster_p_minus_unclustered_conj_p >= 0.0)
			percent_improved_by_clustering = 100.0 * float(num_improved_by_clustering) / float(num_improved_by_clustering + num_not_improved_by_clustering)
			percent_improved_by_clustering_string = "{:.{}f}".format(percent_improved_by_clustering, 1)

			num_non_significant_before_clustering = np.sum(np.asarray(sdds['Unclustered_Conj_LR_pvalue']) > 0.05)
			num_sdd_clusters = len(sdds['Unclustered_Conj_LR_pvalue'])
			percent_non_significant_before_clustering = 100.0*(num_non_significant_before_clustering/num_sdd_clusters)
			percent_non_significant_before_clustering_string = "{:.{}f}".format(percent_non_significant_before_clustering, 1)



			sw.basic_x_y_plot([sdds['Unclustered_Conj_LR_pvalue']], [sdds['FlatClusterStats_FlatCluster_LR_pvalue']], specific_sec_clus_corr_dir + "NonConjPVal_Vs_ClusPVal", draw_y_equals_x=True, y_equals_x_max=1.0, x_axis_label='p-value', y_axis_label='p-value', scatter_point_color_groups=['b'], custom_x_tick_locators=[1.0, 0.2], dashes=(8, 2))
			sw.normal_histo_plot([sdds['Unclustered_Conj_LR_pvalue']], specific_sec_clus_corr_dir + "ConjPVal_Vs_ClusPVal", bins=20, histo_range=[0.0, 1.0], x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[10, 10], alpha=0.78)
			sw.normal_histo_plot([np.asarray(sdds['FlatClusterStats_FlatCluster_LR_pvalue']) - np.asarray(sdds['Unclustered_Conj_LR_pvalue'])], specific_sec_clus_corr_dir + "ClusPVal_Minus_ConjPVal_Hist", bins=21, histo_range=[-1.0, 0.05], x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[10, 10], alpha=0.78)

			# Cluster Differences Correlations
			sw.normal_histo_plot([sdds['FlatClusterStats_FlatCluster_Diff_LR_pvalue']], dh.clus_pair_differences_dir + "FS0_Vs_Diff_LR_PVal_ZoomHist" + suf, bins=20, histo_range=[0.0, 0.1], x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[0.1, 0.01], custom_y_tick_locators=[200, 200], alpha=0.78, add_chi_squared_text=True)
			differences_chi_squared_table_strings_array = sw.cumulative_histo_plot([sdds['FlatClusterStats_FlatCluster_Diff_LR_pvalue']], dh.clus_pair_differences_dir + "FS0_Vs_Diff_LR_PVal_CumHist" + suf, bins=200, x_axis_label="p-value", y_axis_label="Normalised\ncumulative sum", custom_x_tick_locators=[1.0, 0.2], add_chi_squared_text=True)
			sw.normal_histo_plot([sdds['FlatClusterStats_FlatCluster_Diff_LR_pvalue']], dh.clus_pair_differences_dir + "FS0_Vs_Diff_LR_PVal_LowResHist" + suf, bins=20, x_axis_label="p-value", y_axis_label="Frequency", custom_y_tick_locators=[100, 20], alpha=0.78, add_chi_squared_text=True)

			
			# Cluster Correlation Summary Latex
			sw.append_new_tag(str(len(pdds['FlatClusterStats_FlatCluster_LR_pvalue'])) + " Stage 1 clusters were extracted", "NumStage1ClustersFullString", tex_tag_file_name)
			sw.append_new_tag(str(len(pdds['FlatClusterStats_FlatCluster_LR_pvalue'])), "NumStage1ClustersData", tex_tag_file_name)

			
			cluster_correlation_string0 = "Spike pairs within Stage 1 cluster ellipses were linearly correlated above chance levels (Fisher's method: " + flat_cluster_correlations_chi_squared_table_strings_array[0] + ")"
			sw.append_new_tag(cluster_correlation_string0, "Stage1ClusterFisherFullString", tex_tag_file_name)
			sw.append_new_tag(flat_cluster_correlations_chi_squared_table_strings_array[0], "Stage1ClusterFisherData", tex_tag_file_name)

			cluster_correlation_string0p1 = "spike pair differences were correlated with the spike time of the first neuron in the pair for Stage 2 clusters (Fisher's method: " + differences_chi_squared_table_strings_array[0] + "; Fig. 3g), shows that correlations are not explained by a model of the form $s_1 = s_0 + d + independent\\_noise$ where $d$ is a fixed difference."
			sw.append_new_tag(cluster_correlation_string0p1, "ClusterCorrelationSummary0p1", tex_tag_file_name)

			num_greaterthan = np.sum(np.asarray(sdds['FlatClusterStats_FlatCluster_LR_rvalue']) > 0.0)
			data_part = sw.percent_and_frac_string(num_greaterthan, self.final_angled_cluster_count)
			cluster_correlation_string1 = data_part + " of Stage 2 clusters were positively correlated "
			sw.append_new_tag(cluster_correlation_string1, "Stage2PositivelyCorrelatedFullString", tex_tag_file_name)
			sw.append_new_tag(data_part, "Stage2PositivelyCorrelatedNum", tex_tag_file_name)
			
			
			cluster_correlation_string2 = percent_improved_by_clustering_string + "\\% (" + str(num_improved_by_clustering) + "/" + str(num_improved_by_clustering + num_not_improved_by_clustering) + ") of the Stage 2 clusters had correlations of higher significance than correlations calculated for all unclustered first spike pairs in the originating response distribution (Fig. 3h). Moreover, " + percent_non_significant_before_clustering_string + "\\% (" + str(num_non_significant_before_clustering) + '/' + str(num_sdd_clusters) + ") of the original response distributions from which Stage 2 clusters were extracted were not correlated significantly (p>0.05) (Fig. 3h). "
			sw.append_new_tag(cluster_correlation_string2, "ClusterCorrelationSummary2", tex_tag_file_name)

			angled_clusters_unique_pairs_summary_string = "A total of " + str(self.final_angled_cluster_count) + " unique Stage 2 clusters were extracted from " + str(self.did_contribute_atleast_one_final_angled_cluster_count) + " unique response distributions." #, confirming that there were no repeated or similar clusters."
			sw.append_new_tag(angled_clusters_unique_pairs_summary_string, "AngledClustersUniquePairsSummary", tex_tag_file_name)

			
			# Angle Comparisons
			sw.basic_x_y_plot([sdds["Original" + '_BS_PCA_mean_angle']], [sdds["SelectivelyDifferencedBoxJenkins" + '_FA_angle_BS_mean']], dh.angle_analysis_directory + "BS_PCA_VS_SelectivelyDifferencedBoxJenkins_FA_Angles" + suf, draw_y_equals_x=True, y_equals_x_max=90, x_axis_label='Degrees', y_axis_label='Degrees', s=4, scatter_point_color_groups=['g'], custom_x_tick_locators=[90, 10])
			
			# Cluster Reliabilities
			sw.plot_cluster_reliability_plots(sdds['PCA_ellipse_overall_reliability'], sdds['PCA_ellipse_conj_reliability'], dh.cluster_reliabilities_dir, suf)


			analysis_dict_keys= ['Original', 'OriginalTestsPassed', "SelectivelyDifferenced", "SelectivelyDifferencedTestsPassedActuallyDifferenced", "SelectivelyDifferencedBoxJenkins", "SelectivelyDifferencedBoxJenkinsTestsPassed"]
			if ('analysis_dict_member_keys' in sdds.keys()):

				analysis_dict_member_keys = sdds['analysis_dict_member_keys']
				
				for analysis_dict_key in analysis_dict_keys:

					# Directories
					specific_angle_analysis_dir = dh.angle_analysis_directory + analysis_dict_key + "/" + suf + "/"; sw.makedirs(specific_angle_analysis_dir)
					specific_nonstationarity_dir = dh.clus_non_stationarity_dir + analysis_dict_key + "/" + suf + "/"; sw.makedirs(specific_nonstationarity_dir)
					sharipo_normality_specific_nonstationarity_dir = specific_nonstationarity_dir + "SharipoNormality/"; sw.makedirs(sharipo_normality_specific_nonstationarity_dir)
					KPSS_stationarity_specific_nonstationarity_dir = specific_nonstationarity_dir + "KPSSStationarity/"; sw.makedirs(KPSS_stationarity_specific_nonstationarity_dir)
					ADF_stationarity_specific_nonstationarity_dir = specific_nonstationarity_dir + "ADFStationarity/"; sw.makedirs(ADF_stationarity_specific_nonstationarity_dir)
					LR_specific_nonstationarity_dir = specific_nonstationarity_dir + "LRStationarity/"; sw.makedirs(LR_specific_nonstationarity_dir)
					HZ_specific_nonstationarity_dir = specific_nonstationarity_dir + "HZStationarity/"; sw.makedirs(HZ_specific_nonstationarity_dir)
					bartlett_specific_nonstationarity_dir = specific_nonstationarity_dir + "BartlettSphericity/"; sw.makedirs(bartlett_specific_nonstationarity_dir)
					specific_lag_pvals_nonstationary_dir = specific_nonstationarity_dir + "LagPVals/"; sw.makedirs(specific_lag_pvals_nonstationary_dir)
					LR_correlation_specific_nonstationarity_dir = specific_nonstationarity_dir + "LRCorrelation/"; sw.makedirs(LR_correlation_specific_nonstationarity_dir)

					true_where_tests_not_passed_ORIGINAL = np.asarray(sdds['Original_tests_passed'])
					num_tests_not_passed_ORIGINAL = np.sum(true_where_tests_not_passed_ORIGINAL == False)


					if (analysis_dict_key in ["Original", "SelectivelyDifferencedBoxJenkins", "SelectivelyDifferenced"]):

						num_for_type = np.sum(np.bitwise_not(np.asarray(sdds[analysis_dict_key + '_is_empty'])))

						true_where_normal = np.asarray(sdds[analysis_dict_key + '_normal'])
						num_normal = np.sum(true_where_normal)
						where_normal = np.where(true_where_normal)

						true_where_tests_passed = np.asarray(sdds[analysis_dict_key + '_tests_passed'])
						num_tests_passed = np.sum(true_where_tests_passed)
						where_tests_passed = np.where(true_where_tests_passed)

						true_where_tests_not_passed = np.asarray(sdds[analysis_dict_key + '_tests_passed'])
						num_tests_not_passed = np.sum(true_where_tests_not_passed == False)

						true_where_tests_passed_and_normal = np.asarray(sdds[analysis_dict_key + '_tests_passed_and_normal'])
						num_tests_passed_and_normal = np.sum(true_where_tests_passed_and_normal)
						where_tests_passed_and_normal = np.where(true_where_tests_passed_and_normal)

						true_where_correlated = np.asarray(sdds[analysis_dict_key + '_is_still_correlated'])
						number_correlated = np.sum(true_where_correlated)
						where_correlated = np.where(true_where_correlated)

						true_where_tests_passed_and_correlated = np.logical_and(true_where_correlated, true_where_tests_passed)
						num_tests_passed_and_correlated = np.sum(true_where_tests_passed_and_correlated)
						where_tests_passed_and_correlated = np.where(true_where_tests_passed_and_correlated)

						where_different_from_45 = np.logical_and(np.asarray(sdds[analysis_dict_key + '_is_PCA_BS_empirical_pvalue_different_from_45']), np.asarray(sdds[analysis_dict_key + '_is_PCA_BS_empirical_pvalue_different_from_0']))
						num_different_from_45 = np.sum(where_different_from_45)

						true_where_correlated_and_different_from_45 = np.logical_and(true_where_correlated, np.asarray(sdds[analysis_dict_key + '_is_PCA_BS_empirical_pvalue_different_from_45']))
						num_correlated_and_different_from_45 = np.sum(true_where_correlated_and_different_from_45)
						where_correlated_and_different_from_45 = np.where(true_where_correlated_and_different_from_45)

						true_where_correlated_and_different_from_45_tests_passed = np.logical_and(true_where_correlated_and_different_from_45, true_where_tests_passed)
						num_correlated_and_different_from_45_tests_passed = np.sum(true_where_correlated_and_different_from_45_tests_passed)
						where_correlated_and_different_from_45_tests_passed = np.where(true_where_correlated_and_different_from_45_tests_passed)

						true_where_correlated_and_different_from_45_tests_passed_and_normal = np.logical_and(true_where_correlated_and_different_from_45, true_where_tests_passed_and_normal)
						num_correlated_and_different_from_45_tests_passed_and_normal = np.sum(true_where_correlated_and_different_from_45_tests_passed_and_normal)
						where_correlated_and_different_from_45_tests_passed_and_normal = np.where(true_where_correlated_and_different_from_45_tests_passed_and_normal)

						true_where_correlated_and_different_from_45_and_different_from_0 = np.logical_and(true_where_correlated_and_different_from_45, np.asarray(sdds[analysis_dict_key + '_is_PCA_BS_empirical_pvalue_different_from_0']))
						num_correlated_and_different_from_45_and_different_from_0 = np.sum(true_where_correlated_and_different_from_45_and_different_from_0)
						where_correlated_and_different_from_45_and_different_from_0 = np.where(true_where_correlated_and_different_from_45_and_different_from_0)

						true_where_correlated_and_different_from_45_and_different_from_0_tests_passed = np.logical_and(true_where_correlated_and_different_from_45_and_different_from_0, true_where_tests_passed)
						num_correlated_and_different_from_45_and_different_from_0_tests_passed = np.sum(true_where_correlated_and_different_from_45_and_different_from_0_tests_passed)
						where_correlated_and_different_from_45_and_different_from_0_tests_passed = np.where(true_where_correlated_and_different_from_45_and_different_from_0_tests_passed)

						true_where_correlated_and_different_from_45_and_different_from_0_tests_passed_and_normal = np.logical_and(true_where_correlated_and_different_from_45_and_different_from_0, true_where_tests_passed_and_normal)
						num_correlated_and_different_from_45_and_different_from_0_tests_passed_and_normal = np.sum(true_where_correlated_and_different_from_45_and_different_from_0_tests_passed_and_normal)
						where_correlated_and_different_from_45_and_different_from_0_tests_passed_and_normal = np.where(true_where_correlated_and_different_from_45_and_different_from_0_tests_passed_and_normal)

						num_correlated_diferent_from_45_but_not_different_from_0 = num_correlated_and_different_from_45 - num_correlated_and_different_from_45_and_different_from_0
						num_correlated_diferent_from_45_but_not_different_from_0_tests_passed = num_correlated_and_different_from_45_tests_passed - num_correlated_and_different_from_45_and_different_from_0_tests_passed
						num_correlated_diferent_from_45_but_not_different_from_0_tests_passed_and_normal = num_correlated_and_different_from_45_tests_passed_and_normal - num_correlated_and_different_from_45_and_different_from_0_tests_passed_and_normal

						data_part = sw.percent_and_frac_string(num_different_from_45, self.final_angled_cluster_count)
						ps_0 = data_part + " of $\\theta_{45}$ angles were between and significantly different from $0^{\circ}$ and $45^{\circ}$ for " + analysis_dict_key
						sw.append_new_tag(ps_0, analysis_dict_key + "AnglesDifferentFrom45FullText", tex_tag_file_name)
						sw.append_new_tag(data_part, analysis_dict_key + "AnglesDifferentFrom45Num", tex_tag_file_name)

						data_part = sw.percent_and_frac_string(num_correlated_and_different_from_45_and_different_from_0_tests_passed, num_tests_passed)
						ps_1_1 = data_part + analysis_dict_key + " had $\\theta_{45}$ angles between and significantly different from $45^{\circ}$ (p<0.025) and $0^{\circ}$ (p<0.025). "
						sw.append_new_tag(ps_1_1, analysis_dict_key + "TestsPassedAngleSummaryFullString", tex_tag_file_name)
						sw.append_new_tag(data_part, analysis_dict_key + "TestsPassedAngleSummaryNum", tex_tag_file_name)

						data_part = sw.percent_and_frac_string(num_for_type - num_normal, num_for_type)
						ps_new_normality_str = "It is important to note that the Henze-Zirkler null hypothesis of normality was rejected (p < 0.05) for " + data_part + " " + analysis_dict_key
						sw.append_new_tag(ps_new_normality_str, analysis_dict_key + "NormalityFullString", tex_tag_file_name)
						sw.append_new_tag(data_part, analysis_dict_key + "NormalityNum", tex_tag_file_name)

						data_part = sw.percent_and_frac_string(num_tests_passed - num_tests_passed_and_normal, num_tests_passed)
						ps_new_tests_passed_normality_str = "It is important to note that the Henze-Zirkler null hypothesis (p < 0.05) of normality was rejected for " + data_part + " " + analysis_dict_key + "TestsPassed"
						sw.append_new_tag(ps_new_tests_passed_normality_str, analysis_dict_key + "TestsPassedNormalityFullString", tex_tag_file_name)
						sw.append_new_tag(data_part, analysis_dict_key + "TestsPassedNormalityNum", tex_tag_file_name)

						if (analysis_dict_key == "Original"):

							num_stage_2_clusters_string = str(self.final_angled_cluster_count) + " unique Stage 2 clusters were extracted" 
							sw.append_new_tag(num_stage_2_clusters_string, "NumStage2ClustersFullString", tex_tag_file_name)
							sw.append_new_tag(str(self.final_angled_cluster_count), "NumStage2ClustersNum", tex_tag_file_name)


							data_part = sw.percent_and_frac_string(number_correlated, self.final_angled_cluster_count)
							ps_p1 = data_part + " of which remained correlated above chance (p < 0.005)."
							sw.append_new_tag(ps_p1, "NumStage2ClustersCorrelatedFullString", tex_tag_file_name)
							sw.append_new_tag(data_part, "NumStage2ClustersCorrelatedNum", tex_tag_file_name)

							ps_p2 = str(self.final_angled_cluster_count - number_correlated) + " Stage 2 clusters were not significantly correlated (p >= 0.005) and were removed from further analysis, leaving " + str(number_correlated) + " Stage 2 clusters."
							sw.append_new_tag(ps_p2, "OriginalStage2Uncorrelated", tex_tag_file_name)

							data_part = sw.percent_and_frac_string(num_tests_passed, self.final_angled_cluster_count)
							ps_1_0 = data_part + " of Stage 2 clusters were determined as stationary "
							sw.append_new_tag(ps_1_0, "OriginalTestsPassedSummaryFullString", tex_tag_file_name)
							sw.append_new_tag(data_part, "OriginalTestsPassedSummaryNum", tex_tag_file_name)

							ps_1_1 = str(num_tests_passed_and_normal) + " stationary Stage 2 clusters were determined as normal (Henze-Zirkler p > 0.05), of which " + str(num_correlated_and_different_from_45_and_different_from_0_tests_passed_and_normal) + " had $\\theta_{45}$ angles between and significantly different from $45^{\circ}$ (p-value < 0.025) and $0^{\circ}$ (p-value < 0.025)."
							sw.append_new_tag(ps_1_1, "OriginalTestsPassedAndNormalAngleSummary", tex_tag_file_name)

							ps_1_2 = "To test whether trial-to-trial excitability fluctuations also modulate the remaining " + str(num_tests_not_passed) + " non-stationary Stage 2 clusters, "
							sw.append_new_tag(ps_1_2, "NumNonStationary", tex_tag_file_name)
							

						if (analysis_dict_key in ["SelectivelyDifferenced", "SelectivelyDifferencedBoxJenkins"]):

							true_where_undifferenced = np.asarray(sdds[analysis_dict_key + '_selective_differences_undifferenced'])
							where_undifferenced = np.where(true_where_undifferenced)
							num_undifferenced = np.sum(true_where_undifferenced)

							true_where_differenced = np.asarray(sdds[analysis_dict_key + '_selective_differences_differenced'])
							where_differenced = np.where(true_where_differenced)
							num_differenced = np.sum(true_where_differenced)

							sw.draw_neighbouring_bar_chart([[num_undifferenced], [num_differenced]], ('Clusters'), specific_nonstationarity_dir + "ClustersDifferenced.pdf", '', ('Undifferenced', 'Differenced'), "")


							true_where_differenced_and_undifferenced = np.logical_or(true_where_undifferenced, true_where_differenced)

							number_of_differences = np.hstack((np.asarray(sdds[analysis_dict_key + '_index_to_use_0'])[np.where(true_where_differenced_and_undifferenced)], np.asarray(sdds[analysis_dict_key + '_index_to_use_1'])[np.where(true_where_differenced_and_undifferenced)]))
							hist_number_of_selective_differences, _ = np.histogram(number_of_differences, range=[0,2], bins=3)


						if (analysis_dict_key in ["SelectivelyDifferenced"]):

							data_part = sw.percent_and_frac_string(num_differenced, num_tests_not_passed_ORIGINAL)
							ps_2 = data_part + " of non-stationary clusters had differencing applied to at least one neuron. "
							sw.append_new_tag(ps_2, "SelectivelyDifferencedSummaryFullString", tex_tag_file_name)
							sw.append_new_tag(data_part, "SelectivelyDifferencedSummaryNum", tex_tag_file_name)

							ps_3 = "In total, " + str(number_correlated) + " of the " + str(num_for_type) + " criteria fulfilling clusters were linearly correlated following differencing (p-value < 0.05) (Fig. 6h). "
							ps_3 += str(num_correlated_and_different_from_45_and_different_from_0_tests_passed) + " of the correlated criteria fulfilling clusters had $\\theta_{45}$ angles between and significantly different from $0^{\circ}$ (p-value < 0.025) and less than $45^{\circ}$ (p-value < 0.025). "
							ps_3 += "Of these, normality was not rejected for " + str(num_correlated_and_different_from_45_tests_passed_and_normal) + " clusters (Henze-Zirkler p > 0.05). "
							ps_3 += "Fig. 6j shows that the distribution of estimated angles for the successfully differenced clusters. Fig. 8? shows the difference between original and differenced estimated clusters?"
							sw.append_new_tag(ps_3, "s", tex_tag_file_name)

							
							data_part = sw.percent_and_frac_string(num_tests_passed_and_correlated, num_differenced)
							ps_5 = data_part + " of these clusters fulfilled the stationarity and correlation criteria (see Methods), of which "
							sw.append_new_tag(ps_5, "SelectivelyDifferencedTestsPassedNewSummaryFullString1", tex_tag_file_name)
							sw.append_new_tag(data_part, "SelectivelyDifferencedTestsPassedNewSummaryNum1", tex_tag_file_name)

						if (analysis_dict_key in ["SelectivelyDifferencedBoxJenkins"]):

		
							true_where_was_box_jenkinsed = np.asarray(sdds[analysis_dict_key + '_was_box_jeckinsed'])
							true_where_differenced_and_boxed = np.logical_and(true_where_differenced, true_where_was_box_jenkinsed)
							num_true_where_differenced_and_boxed = np.sum(true_where_differenced_and_boxed)
							true_where_tests_passed_differenced_and_boxed = np.logical_and(true_where_tests_passed, true_where_differenced_and_boxed)
							true_where_tests_passed_and_normal_differenced_and_boxed = np.logical_and(true_where_tests_passed_and_normal, true_where_differenced_and_boxed)
							num_true_where_tests_passed_differenced_and_boxed = np.sum(true_where_tests_passed_differenced_and_boxed)
							num_true_where_tests_passed_and_normal_differenced_and_boxed = np.sum(true_where_tests_passed_and_normal_differenced_and_boxed)

							data_part = sw.percent_and_frac_string(num_true_where_differenced_and_boxed, num_differenced)
							ps_8 = data_part + " of these differenced clusters had an AR and/or MA model applied to at least one neuron. "
							sw.append_new_tag(ps_8, "SelectivelyDifferencedBoxJenkinsDifferencedSummary", tex_tag_file_name)
							sw.append_new_tag(data_part, "SelectivelyDifferencedBoxJenkinsDifferencedSummaryNum", tex_tag_file_name)

							data_part_1 = sw.percent_and_frac_string(num_tests_passed, num_true_where_differenced_and_boxed)
							ps_9 = data_part_1 + " of ARIMA modelled clusters fulfilled the criteria, of which " 
							sw.append_new_tag(ps_9, "SelectivelyDifferencedBoxJenkinsCorrelatedSummary1", tex_tag_file_name)
							sw.append_new_tag(data_part_1, "SelectivelyDifferencedBoxJenkinsCorrelatedSummaryNum1", tex_tag_file_name)

							differenced_boxed_ARs = np.hstack((np.asarray(sdds[analysis_dict_key + '_AR_p_0'])[np.where(true_where_differenced_and_boxed)], np.asarray(sdds[analysis_dict_key + '_AR_p_1'])[np.where(true_where_differenced_and_boxed)]))
							differenced_boxed_MAs = np.hstack((np.asarray(sdds[analysis_dict_key + '_MA_q_0'])[np.where(true_where_differenced_and_boxed)], np.asarray(sdds[analysis_dict_key + '_MA_q_1'])[np.where(true_where_differenced_and_boxed)]))
							differenced_boxed_ARs_histo_counts = np.histogram(differenced_boxed_ARs, range=[0, 5], bins=5)[0].tolist()
							differenced_boxed_MAs_histo_counts = np.histogram(differenced_boxed_MAs, range=[0, 5], bins=5)[0].tolist()
							sw.draw_neighbouring_bar_chart([differenced_boxed_ARs_histo_counts, differenced_boxed_MAs_histo_counts], ('0', '1', '2', '3', '4'), specific_nonstationarity_dir + "Differenced_ARIMA.pdf", '', ('AR', 'MA'), 'Order', custom_y_tick_locators=[220, 20])


					# Angles
					sw.plot_angle_vs_reliability_plots(sdds[analysis_dict_key + '_BS_PCA_mean_angle_up_to_45'], sdds['PCA_ellipse_overall_reliability'], sdds['PCA_ellipse_conj_reliability'], specific_angle_analysis_dir, "Angle_Vs_Reliability_BS_PCA")
					sw.normal_histo_plot([sdds[analysis_dict_key + '_BS_PCA_different_from_45_sd_method']], specific_angle_analysis_dir + "Diff45_BS_Sd_PVal_Hist", bins=20, x_axis_label="p-value", y_axis_label="Frequency", alpha=0.78, add_chi_squared_text=True)
					sw.normal_histo_plot([sdds[analysis_dict_key + '_PCA_BS_empirical_pvalue_different_from_45']], specific_angle_analysis_dir + "Diff45_BS_PCA_empirical_pvalue_hist", bins=40, x_axis_label="p-value", y_axis_label="Frequency", alpha=0.78, add_chi_squared_text=True)
					sw.normal_histo_plot([sdds[analysis_dict_key + '_PCA_BS_empirical_pvalue_different_from_0']], specific_angle_analysis_dir + "Diff0_BS_PCA_empirical_pvalue_hist", bins=40, x_axis_label="p-value", y_axis_label="Frequency", alpha=0.78, add_chi_squared_text=True)				
					sw.plot_angle_confidence_bound_plots(sdds[analysis_dict_key + '_BS_PCA_mean_angle_up_to_45'], sdds[analysis_dict_key + '_PCA_BS_empirical_CI_lower'], sdds[analysis_dict_key + '_PCA_BS_empirical_CI_upper'], sdds[analysis_dict_key + '_is_still_correlated'], sdds[analysis_dict_key + '_is_PCA_BS_empirical_pvalue_different_from_45'], specific_angle_analysis_dir, "AngleCI_BS_PCA")
					sw.basic_x_y_plot([sdds["Original" + '_BS_PCA_mean_angle']], [sdds[analysis_dict_key + '_BS_PCA_mean_angle']], specific_angle_analysis_dir + "OriginalAngle_Vs_Angle_BS_PCA", draw_y_equals_x=True, y_equals_x_max=90, x_axis_label='Degrees', y_axis_label='Degrees', s=4, scatter_point_color_groups=['g'], custom_x_tick_locators=[90, 10])
					sw.basic_x_y_plot([sdds[analysis_dict_key + '_BS_PCA_mean_angle']], [sdds[analysis_dict_key + '_FA_angle_BS_mean']], specific_angle_analysis_dir + "Angle_BS_PCA_Vs_FA", draw_y_equals_x=True, y_equals_x_max=90, x_axis_label='Degrees', y_axis_label='Degrees', s=4, scatter_point_color_groups=['g'], custom_x_tick_locators=[90, 10])

					# Non Stationarity
					acf_rvalues = np.vstack((np.asarray(sdds[analysis_dict_key + '_STs_acf_rvalues_0']), np.asarray(sdds[analysis_dict_key + '_STs_acf_rvalues_1'])))
					acf_pvalues = np.vstack((np.asarray(sdds[analysis_dict_key + '_STs_acf_pvalues_0']), np.asarray(sdds[analysis_dict_key + '_STs_acf_pvalues_1'])))
					acf_positive_pvalues = np.vstack((np.asarray(sdds[analysis_dict_key + '_STs_acf_positive_pvalues_0']), np.asarray(sdds[analysis_dict_key + '_STs_acf_positive_pvalues_1'])))
					acf_negative_pvalues = np.vstack((np.asarray(sdds[analysis_dict_key + '_STs_acf_negative_pvalues_0']), np.asarray(sdds[analysis_dict_key + '_STs_acf_negative_pvalues_1'])))

					ccf_pvalues = np.vstack((np.asarray(sdds[analysis_dict_key + '_STs_ccf_pvalues_0']), np.asarray(sdds[analysis_dict_key + '_STs_ccf_pvalues_1'])))
					ccf_positive_pvalues = np.vstack((np.asarray(sdds[analysis_dict_key + '_STs_ccf_positive_pvalues_0']), np.asarray(sdds[analysis_dict_key + '_STs_ccf_positive_pvalues_1']))) 
					ccf_negative_pvalues = np.vstack((np.asarray(sdds[analysis_dict_key + '_STs_ccf_negative_pvalues_0']), np.asarray(sdds[analysis_dict_key + '_STs_ccf_negative_pvalues_1']))) 

					pccf_pvalues = np.vstack((np.asarray(sdds[analysis_dict_key + '_STs_pccf_pvalues_0']), np.asarray(sdds[analysis_dict_key + '_STs_pccf_pvalues_1'])))
					pccf_positive_pvalues = np.vstack((np.asarray(sdds[analysis_dict_key + '_STs_pccf_positive_pvalues_0']), np.asarray(sdds[analysis_dict_key + '_STs_pccf_positive_pvalues_1']))) 
					pccf_negative_pvalues = np.vstack((np.asarray(sdds[analysis_dict_key + '_STs_pccf_negative_pvalues_0']), np.asarray(sdds[analysis_dict_key + '_STs_pccf_negative_pvalues_1']))) 
					
					pacf_pvalues = np.vstack((np.asarray(sdds[analysis_dict_key + '_STs_pacf_pvalues_0']), np.asarray(sdds[analysis_dict_key + '_STs_pacf_pvalues_1'])))
					pacf_positive_pvalues = np.vstack((np.asarray(sdds[analysis_dict_key + '_STs_pacf_positive_pvalues_0']), np.asarray(sdds[analysis_dict_key + '_STs_pacf_positive_pvalues_1'])))
					pacf_negative_pvalues = np.vstack((np.asarray(sdds[analysis_dict_key + '_STs_pacf_negative_pvalues_0']), np.asarray(sdds[analysis_dict_key + '_STs_pacf_negative_pvalues_1'])))

					sharipo_normality_pvalues = np.hstack((np.asarray(sdds[analysis_dict_key + '_sharipo_normality_p_0']), np.asarray(sdds[analysis_dict_key + '_sharipo_normality_p_1'])))
					kpss_stationarity_pvalues = np.hstack((np.asarray(sdds[analysis_dict_key + '_KPSS_STs_0_pvalue']), np.asarray(sdds[analysis_dict_key + '_KPSS_STs_1_pvalue'])))
					ADFuller_STs_0_and_1_pvalue = np.hstack((sdds[analysis_dict_key + '_ADFuller_STs_0_pvalue'], sdds[analysis_dict_key + '_ADFuller_STs_1_pvalue']))
					TI_Vs_STs_LR_0_and_1_pvalue = np.hstack((sdds[analysis_dict_key + '_TI_Vs_STs_LR_0_pvalue'], sdds[analysis_dict_key + '_TI_Vs_STs_LR_1_pvalue']))

					(number_of_examples, number_of_lags) = acf_pvalues.shape

					for number_of_lags_for_cumulative_plot in [5, number_of_lags]:

						particular_number_of_lags_specific_lag_pvals_nonstationary_dir = specific_lag_pvals_nonstationary_dir + str(number_of_lags_for_cumulative_plot) + "Lags/"
						sw.mkdir(particular_number_of_lags_specific_lag_pvals_nonstationary_dir)

						sw.cumulative_histo_plot([acf_pvalues[:, lag_index_zeroed] for lag_index_zeroed in range(number_of_lags_for_cumulative_plot)], particular_number_of_lags_specific_lag_pvals_nonstationary_dir + analysis_dict_key + "_ACF_PVal_CumHist" + suf, bins=200, x_axis_label="p-value", y_axis_label="Normalised\ncumulative sum", custom_x_tick_locators=[1.0, 0.2], add_chi_squared_text=True)
						sw.cumulative_histo_plot([acf_positive_pvalues[:, lag_index_zeroed] for lag_index_zeroed in range(number_of_lags_for_cumulative_plot)], particular_number_of_lags_specific_lag_pvals_nonstationary_dir + analysis_dict_key + "_ACF_PVal_positive_CumHist" + suf, bins=200, x_axis_label="p-value", y_axis_label="Normalised\ncumulative sum", custom_x_tick_locators=[1.0, 0.2], add_chi_squared_text=True)
						sw.cumulative_histo_plot([acf_negative_pvalues[:, lag_index_zeroed] for lag_index_zeroed in range(number_of_lags_for_cumulative_plot)], particular_number_of_lags_specific_lag_pvals_nonstationary_dir + analysis_dict_key + "_ACF_PVal_negative_CumHist" + suf, bins=200, x_axis_label="p-value", y_axis_label="Normalised\ncumulative sum", custom_x_tick_locators=[1.0, 0.2], add_chi_squared_text=True)

						sw.cumulative_histo_plot([ccf_pvalues[:, lag_index_zeroed] for lag_index_zeroed in range(number_of_lags_for_cumulative_plot)], particular_number_of_lags_specific_lag_pvals_nonstationary_dir + analysis_dict_key + "_CCF_PVal_CumHist" + suf, bins=200, x_axis_label="p-value", y_axis_label="Normalised\ncumulative sum", custom_x_tick_locators=[1.0, 0.2], add_chi_squared_text=True)
						sw.cumulative_histo_plot([ccf_positive_pvalues[:, lag_index_zeroed] for lag_index_zeroed in range(number_of_lags_for_cumulative_plot)], particular_number_of_lags_specific_lag_pvals_nonstationary_dir + analysis_dict_key + "_CCF_PVal_positive_CumHist" + suf, bins=200, x_axis_label="p-value", y_axis_label="Normalised\ncumulative sum", custom_x_tick_locators=[1.0, 0.2], add_chi_squared_text=True)
						sw.cumulative_histo_plot([ccf_negative_pvalues[:, lag_index_zeroed] for lag_index_zeroed in range(number_of_lags_for_cumulative_plot)], particular_number_of_lags_specific_lag_pvals_nonstationary_dir + analysis_dict_key + "_CCF_PVal_negative_CumHist" + suf, bins=200, x_axis_label="p-value", y_axis_label="Normalised\ncumulative sum", custom_x_tick_locators=[1.0, 0.2], add_chi_squared_text=True)

						sw.cumulative_histo_plot([pacf_pvalues[:, lag_index_zeroed] for lag_index_zeroed in range(number_of_lags_for_cumulative_plot)], particular_number_of_lags_specific_lag_pvals_nonstationary_dir + analysis_dict_key + "_PACF_PVal_CumHist" + suf, bins=200, x_axis_label="p-value", y_axis_label="Normalised\ncumulative sum", custom_x_tick_locators=[1.0, 0.2], add_chi_squared_text=True)
						sw.cumulative_histo_plot([pacf_positive_pvalues[:, lag_index_zeroed] for lag_index_zeroed in range(number_of_lags_for_cumulative_plot)], particular_number_of_lags_specific_lag_pvals_nonstationary_dir + analysis_dict_key + "_PACF_PVal_positive_CumHist" + suf, bins=200, x_axis_label="p-value", y_axis_label="Normalised\ncumulative sum", custom_x_tick_locators=[1.0, 0.2], add_chi_squared_text=True)
						sw.cumulative_histo_plot([pacf_negative_pvalues[:, lag_index_zeroed] for lag_index_zeroed in range(number_of_lags_for_cumulative_plot)], particular_number_of_lags_specific_lag_pvals_nonstationary_dir + analysis_dict_key + "_PACF_PVal_negative_CumHist" + suf, bins=200, x_axis_label="p-value", y_axis_label="Normalised\ncumulative sum", custom_x_tick_locators=[1.0, 0.2], add_chi_squared_text=True)

						sw.cumulative_histo_plot([pccf_pvalues[:, lag_index_zeroed] for lag_index_zeroed in range(number_of_lags_for_cumulative_plot)], particular_number_of_lags_specific_lag_pvals_nonstationary_dir + analysis_dict_key + "_PCCF_PVal_CumHist" + suf, bins=200, x_axis_label="p-value", y_axis_label="Normalised\ncumulative sum", custom_x_tick_locators=[1.0, 0.2], add_chi_squared_text=True)
						sw.cumulative_histo_plot([pccf_positive_pvalues[:, lag_index_zeroed] for lag_index_zeroed in range(number_of_lags_for_cumulative_plot)], particular_number_of_lags_specific_lag_pvals_nonstationary_dir + analysis_dict_key + "_PCCF_PVal_positive_CumHist" + suf, bins=200, x_axis_label="p-value", y_axis_label="Normalised\ncumulative sum", custom_x_tick_locators=[1.0, 0.2], add_chi_squared_text=True)
						sw.cumulative_histo_plot([pccf_negative_pvalues[:, lag_index_zeroed] for lag_index_zeroed in range(number_of_lags_for_cumulative_plot)], particular_number_of_lags_specific_lag_pvals_nonstationary_dir + analysis_dict_key + "_PCCF_PVal_negative_CumHist" + suf, bins=200, x_axis_label="p-value", y_axis_label="Normalised\ncumulative sum", custom_x_tick_locators=[1.0, 0.2], add_chi_squared_text=True)


					if (plot_all_lag_histograms):

						acf_dir = specific_lag_pvals_nonstationary_dir + "ACFs/"; sw.mkdir(acf_dir)
						acf_positive_dir = specific_lag_pvals_nonstationary_dir + "ACFs_Positive/"; sw.mkdir(acf_positive_dir)
						acf_negative_dir = specific_lag_pvals_nonstationary_dir + "ACFs_Negative/"; sw.mkdir(acf_negative_dir)
						ccf_dir = specific_lag_pvals_nonstationary_dir + "CCFs/"; sw.mkdir(ccf_dir)
						ccf_positive_dir = specific_lag_pvals_nonstationary_dir + "CCFs_Positive/"; sw.mkdir(ccf_positive_dir)
						ccf_negative_dir = specific_lag_pvals_nonstationary_dir + "CCFs_Negative/"; sw.mkdir(ccf_negative_dir)
						pccf_dir = specific_lag_pvals_nonstationary_dir + "PCCFs/"; sw.mkdir(pccf_dir)
						pccf_positive_dir = specific_lag_pvals_nonstationary_dir + "PCCFs_Positive/"; sw.mkdir(pccf_positive_dir)
						pccf_negative_dir = specific_lag_pvals_nonstationary_dir + "PCCFs_Negative/"; sw.mkdir(pccf_negative_dir)
						pacf_dir = specific_lag_pvals_nonstationary_dir + "PACFs/"; sw.mkdir(pacf_dir)
						pacf_positive_dir = specific_lag_pvals_nonstationary_dir + "PACFs_Positive/"; sw.mkdir(pacf_positive_dir)
						pacf_negative_dir = specific_lag_pvals_nonstationary_dir + "PACFs_Negative/"; sw.mkdir(pacf_negative_dir)
						
						
						for lag_index_zeroed in range(number_of_lags):
							sw.normal_histo_plot([acf_rvalues[:, lag_index_zeroed]], acf_dir + "RValues_ACFLag" + str(lag_index_zeroed + 1) + suf, bins=40, histo_range=[-1.0, 1.0], x_axis_label="r-value", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[20, 20], alpha=0.78)

						for number_of_bins, histo_range in zip([20, 20], [[0.0, 1.0], [0.0, 0.1]]):

							for lag_index_zeroed in range(number_of_lags):

								sw.normal_histo_plot([acf_pvalues[:, lag_index_zeroed]], acf_dir + str(histo_range[1]) + "_ACFLag" + str(lag_index_zeroed + 1) + suf, bins=number_of_bins, histo_range=histo_range, x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[histo_range[1], 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)
								sw.normal_histo_plot([acf_positive_pvalues[:, lag_index_zeroed]], acf_positive_dir + str(histo_range[1]) + "_ACFPositiveLag" + str(lag_index_zeroed + 1) + suf, bins=number_of_bins, histo_range=histo_range, x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[histo_range[1], 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)
								sw.normal_histo_plot([acf_negative_pvalues[:, lag_index_zeroed]], acf_negative_dir + str(histo_range[1]) + "_ACFNegativeLag" + str(lag_index_zeroed + 1) + suf, bins=number_of_bins, histo_range=histo_range, x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[histo_range[1], 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)

								sw.normal_histo_plot([ccf_pvalues[:, lag_index_zeroed]], ccf_dir + str(histo_range[1]) + "_CCFLag" + str(lag_index_zeroed + 1) + suf, bins=number_of_bins, histo_range=histo_range, x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[histo_range[1], 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)
								sw.normal_histo_plot([ccf_positive_pvalues[:, lag_index_zeroed]], ccf_positive_dir + str(histo_range[1]) + "_CCFPositiveLag" + str(lag_index_zeroed + 1) + suf, bins=number_of_bins, histo_range=histo_range, x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[histo_range[1], 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)
								sw.normal_histo_plot([ccf_negative_pvalues[:, lag_index_zeroed]], ccf_negative_dir + str(histo_range[1]) + "_CCFNegativeLag" + str(lag_index_zeroed + 1) + suf, bins=number_of_bins, histo_range=histo_range, x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[histo_range[1], 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)

								sw.normal_histo_plot([pacf_pvalues[:, lag_index_zeroed]], pacf_dir + str(histo_range[1]) + "_PACFLag" + str(lag_index_zeroed + 1) + suf, bins=number_of_bins, histo_range=histo_range, x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[histo_range[1], 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)
								sw.normal_histo_plot([pacf_positive_pvalues[:, lag_index_zeroed]], pacf_positive_dir + str(histo_range[1]) + "_PACFPositiveLag" + str(lag_index_zeroed + 1) + suf, bins=number_of_bins, histo_range=histo_range, x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[histo_range[1], 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)
								sw.normal_histo_plot([pacf_negative_pvalues[:, lag_index_zeroed]], pacf_negative_dir + str(histo_range[1]) + "_PACFNegativeLag" + str(lag_index_zeroed + 1) + suf, bins=number_of_bins, histo_range=histo_range, x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[histo_range[1], 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)

								sw.normal_histo_plot([pccf_pvalues[:, lag_index_zeroed]], pccf_dir + str(histo_range[1]) + "_PCCFLag" + str(lag_index_zeroed + 1) + suf, bins=number_of_bins, histo_range=histo_range, x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[histo_range[1], 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)
								sw.normal_histo_plot([pccf_positive_pvalues[:, lag_index_zeroed]], pccf_positive_dir + str(histo_range[1]) + "_PCCFPositiveLag" + str(lag_index_zeroed + 1) + suf, bins=number_of_bins, histo_range=histo_range, x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[histo_range[1], 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)
								sw.normal_histo_plot([pccf_negative_pvalues[:, lag_index_zeroed]], pccf_negative_dir + str(histo_range[1]) + "_PCCFNegativeLag" + str(lag_index_zeroed + 1) + suf, bins=number_of_bins, histo_range=histo_range, x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[histo_range[1], 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)


					# Bartlett's Sphericity
					sw.normal_histo_plot([sdds[analysis_dict_key + '_bartlett_spherecity_p_value']], bartlett_specific_nonstationarity_dir + "BartlettsSphericity_pvalues" + suf, bins=20, histo_range=[0.0, 1.0], x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)

					# Henze-Zirkler
					sw.normal_histo_plot([sdds[analysis_dict_key + '_henze-zirkler_multivariate_normality_p']], HZ_specific_nonstationarity_dir + "Henze-Zirkler_pvalues" + suf, bins=20, histo_range=[0.0, 1.0], x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)
							
					# Sharipo
					sw.normal_histo_plot([sharipo_normality_pvalues], sharipo_normality_specific_nonstationarity_dir + "Sharipo_pvalues" + suf, bins=20, histo_range=[0.0, 1.0], x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)

					# KPSS
					sw.normal_histo_plot([kpss_stationarity_pvalues], KPSS_stationarity_specific_nonstationarity_dir + "KPSS_pvalues" + suf, bins=6, histo_range=[0.0, 0.12], x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[0.12, 0.02], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)
					if (analysis_dict_key == "Original"):
						kpss_num_less_than_pvalue = np.sum(kpss_stationarity_pvalues < 0.05)
						kpss_string = "For the " + str(num_for_type) + " angled clusters, " + str(kpss_num_less_than_pvalue) + " of the corresponding " + str(2 * num_for_type) + " ($=2*" + str(num_for_type) + "$) single neuron cluster first spike sequences were determined as KPSS non-stationary (p < 0.05; Fig. 4d)."
						sw.append_new_tag(kpss_string, "ClusterSingleUnitKPSSStationaritySummary", tex_tag_file_name)


					# LR
					cluster_single_unit_trial_index_LR_chi_squared_table_strings_array = sw.normal_histo_plot([TI_Vs_STs_LR_0_and_1_pvalue], LR_specific_nonstationarity_dir + "LR_PVal_Hist_TIVsSpikeTimes" + suf, bins=20, histo_range=[0.0, 1.0], x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)
					sw.basic_x_y_plot([sdds[analysis_dict_key + '_TI_Vs_STs_LR_0_pvalue']], [sdds[analysis_dict_key + '_TI_Vs_STs_LR_1_pvalue']], LR_specific_nonstationarity_dir + "LR_PVal_U0_Vs_U1_TIVsSpikeTimes" + suf, draw_y_equals_x=True, y_equals_x_max=1.0, x_axis_label='p-value', y_axis_label='p-value', scatter_point_color_groups=['b'], custom_x_tick_locators=[1.0, 0.2], dashes=(8, 2), opt_min_lim_buffer=0.01)
					sw.normal_histo_plot([sdds[analysis_dict_key + '_TI_Vs_STs_LR_multiple_output_rsquared']], LR_specific_nonstationarity_dir + "MultipleOutputLR_R^2_Hist_TIVsSpikePair" + suf, bins=20, histo_range=[0.0, 1.0], x_axis_label="$r^2$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[5, 5], alpha=0.78)


					if (analysis_dict_key == "Original"):
						cluster_single_unit_trial_index_LR_string = "The cluster spike times of either 0, 1 or 2 of the neurons in a cluster pair were correlated with trial index (Fisher's method: " + cluster_single_unit_trial_index_LR_chi_squared_table_strings_array[0] + "; Fig. 4b,c)."
						sw.append_new_tag(cluster_single_unit_trial_index_LR_string, "ClusterSingleUnitTrialIndexLRSummary", tex_tag_file_name)

					# LR Correlation
					sw.normal_histo_plot([sdds[analysis_dict_key + '_LR_pvalue']], LR_correlation_specific_nonstationarity_dir + "PVal_LowResHist" + suf, bins=40, x_axis_label="p-value", y_axis_label="Frequency", custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)
					sw.normal_histo_plot([sdds[analysis_dict_key + '_LR_rsquared']], LR_correlation_specific_nonstationarity_dir + "RSquared_Hist" + suf, bins=20, histo_range=[0.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$r^2$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[5, 5], alpha=0.78)
					
					#JIRev
					if (analysis_dict_key in ["OriginalTestsPassed", "SelectivelyDifferencedBoxJenkinsTestsPassed", "SelectivelyDifferencedTestsPassedActuallyDifferenced"]):
						sw.normal_histo_plot([np.asarray(sdds['Unclustered_Conj_LR_rsquaredvalue'])[where_tests_passed_and_correlated], sdds[analysis_dict_key + '_LR_rsquared']], LR_correlation_specific_nonstationarity_dir + "RSquared_Hist_with_unclustered" + suf, bins=20, histo_range=[0.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$r^2$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[10, 10], alpha=0.78)
					if (analysis_dict_key in ["Original"]):
						sw.normal_histo_plot([sdds['Unclustered_Conj_LR_rsquaredvalue'], sdds[analysis_dict_key + '_LR_rsquared']], LR_correlation_specific_nonstationarity_dir + "RSquared_Hist_with_unclustered" + suf, bins=20, histo_range=[0.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$r^2$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[10, 10], alpha=0.78)

					
					sw.normal_histo_plot([sdds[analysis_dict_key + '_LR_rsquared'], sdds[analysis_dict_key + '_TI_Vs_STs_LR_multiple_output_rsquared']], LR_correlation_specific_nonstationarity_dir + "RSquared_Hist_AND_MultipleOutputLR_R^2_Hist_TIVsSpikePair" + suf, bins=20, histo_range=[0.0, 1.0], x_axis_left_buffer=0.01, x_axis_label="$r^2$", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[5, 5], alpha=0.78)


					# ADFuller
					sw.basic_x_y_plot([ADFuller_STs_0_and_1_pvalue], [TI_Vs_STs_LR_0_and_1_pvalue], ADF_stationarity_specific_nonstationarity_dir + "ADFPVals_Vs_TISpikeTimeLRPVals_Original" + suf, draw_y_equals_x=True, y_equals_x_max=1, x_axis_label='p-value', y_axis_label='p-value', s=4, scatter_point_color_groups=['b'], custom_x_tick_locators=[1.0, 0.2], opt_min_lim_buffer=0.01)
					if (analysis_dict_key != "Original"):
						Original_ADFuller_STs_0_and_1_pvalue = np.hstack((sdds["Original" + '_ADFuller_STs_0_pvalue'], sdds["Original" + '_ADFuller_STs_1_pvalue']))
						sw.normal_histo_plot([Original_ADFuller_STs_0_and_1_pvalue, ADFuller_STs_0_and_1_pvalue], ADF_stationarity_specific_nonstationarity_dir + "ADFullerPVals_OriginalAnd" + analysis_dict_key + suf, bins=20, histo_range=[0.0, 1.0], x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=True)
						sw.basic_x_y_plot([Original_ADFuller_STs_0_and_1_pvalue], [ADFuller_STs_0_and_1_pvalue], ADF_stationarity_specific_nonstationarity_dir + "ADFullerPVals_OriginalVs" + analysis_dict_key + suf, draw_y_equals_x=True, y_equals_x_max=1, x_axis_label='p-value', y_axis_label='p-value', s=4, scatter_point_color_groups=['g'], custom_x_tick_locators=[1.0, 0.2], opt_min_lim_buffer=0.01)


		
		# Unit Type Analysis
		if (self.is_mainz):

			self.count_for_each_principle_column_bool_combination = np.zeros((4), dtype=int)
			self.count_for_each_neurophys_layer_ei_type_combination = np.zeros((sw.number_of_combinations_of_all_pairs_of_neurophys_corresponding_layer_name_and_ei_types), dtype=int)

			for secondary_cluster_index in range(len(sdds['principle_condition_bool_for_each_unit'])):

				principle_condition_bools = sdds['principle_condition_bool_for_each_unit'][secondary_cluster_index]
				principle_column_bool_index = sw.index_of_principle_column_options_pair(principle_condition_bools)
				self.count_for_each_principle_column_bool_combination[principle_column_bool_index] = self.count_for_each_principle_column_bool_combination[principle_column_bool_index] + 1				

				layer_strings_for_unit_pair = sdds["neurophys_layer_strings_of_unit_pair"][secondary_cluster_index]
				exc_inh_types = sdds["exc_inh_types_of_unit_pair"][secondary_cluster_index]
				combination_index = sw.index_of_neurophys_layer_string_pair_ei_type_in_combination_list(layer_strings_for_unit_pair, exc_inh_types)
				self.count_for_each_neurophys_layer_ei_type_combination[combination_index] = self.count_for_each_neurophys_layer_ei_type_combination[combination_index] + 1

			plt.figure(figsize=(10, 5))
			plt.barh(sw.strings_for_combinations_of_all_pairs_of_neurophys_corresponding_layer_name_and_ei_types(), self.count_for_each_neurophys_layer_ei_type_combination)
			plt.subplots_adjust(left=0.2)
			plt.savefig(dh.unit_type_analysis_directory + "NeuronGroupCombinationsHist_CorrelatedClusters" + suf)
			plt.close()


			plt.figure()
			plt.bar(list(range(4)), self.count_for_each_principle_column_bool_combination)
			plt.savefig(dh.unit_type_analysis_directory + "ColumnCombinationsHist_CorrelatedClusters" + suf)
			plt.close()


			new_keys = ["OriginalTestsPassed", "SelectivelyDifferencedBoxJenkinsTestsPassed"]
			stimulation_frequency_of_condition_tests_passed_and_correlated_for_each_key = []
			dictionaries_to_merge = []
			cum_atleast_one_in_septum_count = 0
			for key in new_keys:

				key_dir = dh.unit_type_analysis_directory + key + "/"; sw.mkdir(key_dir)
				
				true_where_tests_passed = np.asarray(sdds[key + '_tests_passed'])
				true_where_correlated = np.asarray(sdds[key + '_is_still_correlated'])
				true_where_tests_passed_and_correlated = np.logical_and(true_where_correlated, true_where_tests_passed)
				where_tests_passed_and_correlated = np.where(true_where_tests_passed_and_correlated)

				stimulation_frequency_of_condition_tests_passed_and_correlated = np.asarray(sdds['stimulation_frequency_of_condition'])[where_tests_passed_and_correlated]
				sw.draw_stimulus_frequency_histo(stimulation_frequency_of_condition_tests_passed_and_correlated, key_dir + key + "_StimulationFrequencyHist_CorrelatedClusters" + suf)


				dictionary_of_counts, atleast_one_in_septum_count = sw.layer_column_neuron_type_preprocessing(dh.unit_type_analysis_directory, ['', ''], key, self.sdds, self.pdds, self.suf)

				cum_atleast_one_in_septum_count += atleast_one_in_septum_count
				
				sw.draw_cortex_spatial_seperation_plot(dictionary_of_counts, key_dir + key + "ExcCortexPlot_CorrelatedClusters" + suf + '.pdf', atleast_one_in_septum_count)

				if (key in ["OriginalTestsPassed", "SelectivelyDifferencedBoxJenkinsTestsPassed"]):
					dictionaries_to_merge.append(dictionary_of_counts)
					stimulation_frequency_of_condition_tests_passed_and_correlated_for_each_key.append(np.histogram(stimulation_frequency_of_condition_tests_passed_and_correlated, range=[0, 10], bins=11)[0].tolist())


			merged_dictionary_of_counts = {k: dictionaries_to_merge[0].get(k, 0) + dictionaries_to_merge[1].get(k, 0) for k in set(dictionaries_to_merge[0]) | set(dictionaries_to_merge[1])}
			sw.draw_cortex_spatial_seperation_plot(merged_dictionary_of_counts,  dh.unit_type_analysis_directory + "MergedOrigAndBoxJenkCortexPlot_CorrelatedClusters" + suf + '.pdf', cum_atleast_one_in_septum_count)	

			sw.draw_neighbouring_bar_chart(stimulation_frequency_of_condition_tests_passed_and_correlated_for_each_key, ('0-0.2', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10'), dh.unit_type_analysis_directory + "stimulation_frequency_grouped" + suf + ".pdf", '', ('Stationary', 'ARIMA'), 'Hz', custom_y_tick_locators=[45, 5])

	

	def create_comparison_to_standard_error_plots(self, directory_holder):

		sdds = self.sdds
		pdds = self.pdds
		dh = directory_holder
		suf = self.suf

		y_equals_x_max = 12

		# analysis_dict_keys = ['SelectivelyDifferencedBoxJenkinsTestsPassedAndNormal', 'SelectivelyDifferencedBoxJenkinsTestsPassed', 'SelectivelyDifferencedTestsPassedAndNormal', 'SelectivelyDifferencedTestsPassed', 'OriginalTestsPassedAndNormal', 'OriginalTestsPassed']
		analysis_dict_keys = ['SelectivelyDifferencedBoxJenkinsTestsPassed', 'SelectivelyDifferencedTestsPassed', 'SelectivelyDifferencedTestsPassedActuallyDifferenced', 'OriginalTestsPassed']

		list_fa_1sds = []; list_flat_1sds = []; list_angled_1sds = []; list_conj_1sds = []; list_not_only_conj_1sds = []; list_fa_diffs = []; list_flat_diffs = []; list_angled_diffs = []; list_conj_diffs = []; legend_labels = []

		if ('analysis_dict_member_keys' in sdds.keys()):

			analysis_dict_member_keys = sdds['analysis_dict_member_keys']
			
			for key in analysis_dict_keys:

				if (key + '_tests_passed' in sdds.keys()):

					true_where_tests_passed = np.asarray(sdds[key + '_tests_passed']); num_tests_passed = np.sum(true_where_tests_passed); where_tests_passed = np.where(true_where_tests_passed)
					true_where_tests_passed_but_not_normal = np.asarray(sdds[key + '_tests_passed_but_not_normal']); num_tests_passed_but_not_normal = np.sum(true_where_tests_passed_but_not_normal); where_tests_passed_but_not_normal = np.where(true_where_tests_passed_but_not_normal)
					true_where_tests_passed_and_normal = np.asarray(sdds[key + '_tests_passed_and_normal']); num_tests_passed_and_normal = np.sum(true_where_tests_passed_and_normal); where_tests_passed_and_normal = np.where(true_where_tests_passed_and_normal)

					specif_fa_dir = dh.fa_varaince_comparisons_directory + key + "/"; sw.mkdir(specif_fa_dir)
					one_d_specif_fa_dir = specif_fa_dir + "1D/"; sw.mkdir(one_d_specif_fa_dir)
					two_d_specif_fa_dir = specif_fa_dir + "2D/"; sw.mkdir(two_d_specif_fa_dir)
					diff_specif_fa_dir = specif_fa_dir + "Diff/"; sw.mkdir(diff_specif_fa_dir)
					
					BS_PCA_mean_component_sds_tests_passed_but_not_normal = np.hstack((np.asarray(sdds[key + '_BS_PCA_mean_component_0_sd'])[where_tests_passed_but_not_normal], np.asarray(sdds[key + '_BS_PCA_mean_component_1_sd'])[where_tests_passed_but_not_normal]))
					BS_PCA_mean_component_sds_tests_passed_and_normal = np.hstack((np.asarray(sdds[key + '_BS_PCA_mean_component_0_sd'])[where_tests_passed_and_normal], np.asarray(sdds[key + '_BS_PCA_mean_component_1_sd'])[where_tests_passed_and_normal]))
					FlatCluster_FS_SDs_tests_passed_but_not_normal = np.hstack((np.asarray(sdds['FlatClusterStats_FlatCluster_N0_FS_SD'])[where_tests_passed_but_not_normal], np.asarray(sdds['FlatClusterStats_FlatCluster_N1_FS_SD'])[where_tests_passed_but_not_normal]))
					FlatCluster_FS_SDs_tests_passed_and_normal = np.hstack((np.asarray(sdds['FlatClusterStats_FlatCluster_N0_FS_SD'])[where_tests_passed_and_normal], np.asarray(sdds['FlatClusterStats_FlatCluster_N1_FS_SD'])[where_tests_passed_and_normal]))
					Conjunctive_FS_SDs_tests_passed_but_not_normal = np.hstack((np.asarray(sdds['Conjunctive_FS_SDs_N0'])[where_tests_passed_but_not_normal], np.asarray(sdds['Conjunctive_FS_SDs_N1'])[where_tests_passed_but_not_normal]))
					Conjunctive_FS_SDs_tests_passed_and_normal = np.hstack((np.asarray(sdds['Conjunctive_FS_SDs_N0'])[where_tests_passed_and_normal], np.asarray(sdds['Conjunctive_FS_SDs_N1'])[where_tests_passed_and_normal]))
					NotOnlyConjunctive_FS_SDs_tests_passed_but_not_normal = np.hstack((np.asarray(sdds['NotOnlyConjunctive_FS_SDs_N0'])[where_tests_passed_but_not_normal], np.asarray(sdds['NotOnlyConjunctive_FS_SDs_N1'])[where_tests_passed_but_not_normal]))
					NotOnlyConjunctive_FS_SDs_tests_passed_and_normal = np.hstack((np.asarray(sdds['NotOnlyConjunctive_FS_SDs_N0'])[where_tests_passed_and_normal], np.asarray(sdds['NotOnlyConjunctive_FS_SDs_N1'])[where_tests_passed_and_normal]))
					FA_sds_tests_passed_but_not_normal = np.hstack((np.asarray(sdds[key + '_FA_N0_sd_BS_mean'])[where_tests_passed_but_not_normal], np.asarray(sdds[key + '_FA_N1_sd_BS_mean'])[where_tests_passed_but_not_normal]))
					FA_sds_tests_passed_and_normal = np.hstack((np.asarray(sdds[key + '_FA_N0_sd_BS_mean'])[where_tests_passed_and_normal], np.asarray(sdds[key + '_FA_N1_sd_BS_mean'])[where_tests_passed_and_normal]))



					# 1SDs
					sw.plot_fa_comparison([BS_PCA_mean_component_sds_tests_passed_but_not_normal, BS_PCA_mean_component_sds_tests_passed_and_normal], 
											[FA_sds_tests_passed_but_not_normal, FA_sds_tests_passed_and_normal], one_d_specif_fa_dir + key + "_AngledCluster_V_ClusterFA__1SD_" + suf, 
											y_equals_x_max=12, optional_y_max=6)
					
					sw.plot_fa_comparison([FlatCluster_FS_SDs_tests_passed_but_not_normal, FlatCluster_FS_SDs_tests_passed_and_normal], 
											[FA_sds_tests_passed_but_not_normal, FA_sds_tests_passed_and_normal], one_d_specif_fa_dir + key + "_FlatCluster_V_ClusterFA__1SD_" + suf, 
											y_equals_x_max=12, optional_y_max=6)
					
					sw.plot_fa_comparison([Conjunctive_FS_SDs_tests_passed_but_not_normal, Conjunctive_FS_SDs_tests_passed_and_normal], 
											[FA_sds_tests_passed_but_not_normal, FA_sds_tests_passed_and_normal], one_d_specif_fa_dir + key + "_Conjunctive_V_ClusterFA__1SD_" + suf, 
											y_equals_x_max=18, optional_y_max=6)

					sw.plot_fa_comparison([NotOnlyConjunctive_FS_SDs_tests_passed_but_not_normal, NotOnlyConjunctive_FS_SDs_tests_passed_and_normal], 	
											[FA_sds_tests_passed_but_not_normal, FA_sds_tests_passed_and_normal], one_d_specif_fa_dir + key + "_NotOnlyConjunctive_V_ClusterFA__1SD_" + suf, 
											y_equals_x_max=18, optional_y_max=6)

					# 1SD Area
					sw.plot_fa_comparison([np.asarray(sdds['Original_BS_PCA_mean_sd_area'])[where_tests_passed], np.asarray(sdds['Original_BS_PCA_mean_sd_area'])[where_tests_passed_and_normal]], 
											[np.asarray(sdds[key + '_FA_1sd_area_BS_mean'])[where_tests_passed], np.asarray(sdds[key + '_FA_1sd_area_BS_mean'])[where_tests_passed_and_normal]], 
											two_d_specif_fa_dir + key + "_AngledCluster_V_ClusterFA__1SDArea_" + suf, y_equals_x_max=30, optional_y_max=6)

					sw.plot_fa_comparison([np.asarray(sdds['FlatClusterStats_FlatCluster_1sd_area'])[where_tests_passed], np.asarray(sdds['FlatClusterStats_FlatCluster_1sd_area'])[where_tests_passed_and_normal]], 
											[np.asarray(sdds[key + '_FA_1sd_area_BS_mean'])[where_tests_passed], np.asarray(sdds[key + '_FA_1sd_area_BS_mean'])[where_tests_passed_and_normal]], 
											two_d_specif_fa_dir + key + "_FlatCluster_V_ClusterFA__1SDArea_" + suf, y_equals_x_max=30, optional_y_max=6)				
					
					sw.plot_fa_comparison([np.asarray(sdds['Conjunctive_original_1sd_area'])[where_tests_passed], np.asarray(sdds['Conjunctive_original_1sd_area'])[where_tests_passed_and_normal]], 
											[np.asarray(sdds[key + '_FA_1sd_area_BS_mean'])[where_tests_passed], np.asarray(sdds[key + '_FA_1sd_area_BS_mean'])[where_tests_passed_and_normal]], 
											two_d_specif_fa_dir + key + "_Conjunctive_V_ClusterFA__1SDArea_" + suf, y_equals_x_max=60, optional_y_max=6)

					sw.plot_fa_comparison([np.asarray(sdds['NotOnlyConjunctive_original_1sd_area'])[where_tests_passed], np.asarray(sdds['NotOnlyConjunctive_original_1sd_area'])[where_tests_passed_and_normal]], 
											[np.asarray(sdds[key + '_FA_1sd_area_BS_mean'])[where_tests_passed], np.asarray(sdds[key + '_FA_1sd_area_BS_mean'])[where_tests_passed_and_normal]], 
											two_d_specif_fa_dir + key + "_NotOnlyConjunctive_V_ClusterFA__1SDArea_" + suf, y_equals_x_max=60, optional_y_max=6)


					# 1SD Differences
					sw.plot_fa_comparison([np.asarray(sdds['angled_cluster_spike_pairs_differences_sd'])[where_tests_passed], np.asarray(sdds['angled_cluster_spike_pairs_differences_sd'])[where_tests_passed_and_normal]], 
											[np.asarray(sdds[key + '_FA_1sd_estimated_difference_area_BS_mean'])[where_tests_passed], np.asarray(sdds[key + '_FA_1sd_estimated_difference_area_BS_mean'])[where_tests_passed_and_normal]], 
											diff_specif_fa_dir + key + "_AngledClusterDiff_V_ClusterFADiff__1SD_" + suf, y_equals_x_max=12, optional_y_max=6)

					sw.plot_fa_comparison([np.asarray(sdds['FlatClusterStats_FlatCluster_FS_diff_SD'])[where_tests_passed], np.asarray(sdds['FlatClusterStats_FlatCluster_FS_diff_SD'])[where_tests_passed_and_normal]], 
											[np.asarray(sdds[key + '_FA_1sd_estimated_difference_area_BS_mean'])[where_tests_passed], np.asarray(sdds[key + '_FA_1sd_estimated_difference_area_BS_mean'])[where_tests_passed_and_normal]], 
											diff_specif_fa_dir + key + "_FlatClusterDiff_V_ClusterFADiff__1SD_" + suf, y_equals_x_max=12, optional_y_max=6)

					sw.plot_fa_comparison([np.asarray(sdds['Conjunctive_original_FS_diff_SD'])[where_tests_passed], np.asarray(sdds['Conjunctive_original_FS_diff_SD'])[where_tests_passed_and_normal]], 
											[np.asarray(sdds[key + '_FA_1sd_estimated_difference_area_BS_mean'])[where_tests_passed], np.asarray(sdds[key + '_FA_1sd_estimated_difference_area_BS_mean'])[where_tests_passed_and_normal]], 
											diff_specif_fa_dir + key + "_ConjunctiveDiff_V_ClusterFADiff__1SD_" + suf, y_equals_x_max=18, optional_y_max=6)

					if (key in ["SelectivelyDifferencedBoxJenkinsTestsPassed", "SelectivelyDifferencedTestsPassedActuallyDifferenced", "OriginalTestsPassed"]):

						legend_labels.append(key)
						list_fa_1sds.append(np.hstack((np.asarray(sdds[key + '_FA_N0_sd_BS_mean'])[where_tests_passed], np.asarray(sdds[key + '_FA_N1_sd_BS_mean'])[where_tests_passed])))
						list_flat_1sds.append(np.hstack((np.asarray(sdds['FlatClusterStats_FlatCluster_N0_FS_SD'])[where_tests_passed], np.asarray(sdds['FlatClusterStats_FlatCluster_N1_FS_SD'])[where_tests_passed])))
						list_angled_1sds.append(np.hstack((np.asarray(sdds[key + '_BS_PCA_mean_component_0_sd'])[where_tests_passed], np.asarray(sdds[key + '_BS_PCA_mean_component_1_sd'])[where_tests_passed])))
						list_conj_1sds.append(np.hstack((np.asarray(sdds['Conjunctive_FS_SDs_N0'])[where_tests_passed], np.asarray(sdds['Conjunctive_FS_SDs_N1'])[where_tests_passed])))
						list_not_only_conj_1sds.append(np.hstack((np.asarray(sdds['NotOnlyConjunctive_FS_SDs_N0'])[where_tests_passed], np.asarray(sdds['NotOnlyConjunctive_FS_SDs_N1'])[where_tests_passed])))
						list_fa_diffs.append(np.asarray(sdds[key + '_FA_1sd_estimated_difference_area_BS_mean'])[where_tests_passed])
						list_flat_diffs.append(np.asarray(sdds['FlatClusterStats_FlatCluster_FS_diff_SD'])[where_tests_passed])
						list_angled_diffs.append(np.asarray(sdds['angled_cluster_spike_pairs_differences_sd'])[where_tests_passed])
						list_conj_diffs.append(np.asarray(sdds['Conjunctive_original_FS_diff_SD'])[where_tests_passed])
						


					# Factor Correlation Matrix Eigen Values
					factor_correlation_matrix_eigen_value_0s = np.asarray(sdds[key + '_factor_correlation_matrix_eigen_value_0'])
					factor_correlation_matrix_eigen_value_1s = np.asarray(sdds[key + '_factor_correlation_matrix_eigen_value_1'])

					factor_correlation_matrix_ev0s_greater_equal_1 = np.sum((factor_correlation_matrix_eigen_value_0s >=0.0) & (factor_correlation_matrix_eigen_value_0s >= 1.0))
					factor_correlation_matrix_ev0s_less_than_1 = np.sum((factor_correlation_matrix_eigen_value_0s >=0.0) & (factor_correlation_matrix_eigen_value_0s < 1.0))
					factor_correlation_matrix_ev1s_greater_equal_1 = np.sum((factor_correlation_matrix_eigen_value_1s >= 0.0) & (factor_correlation_matrix_eigen_value_1s >= 1.0))
					factor_correlation_matrix_ev1s_less_than_1 = np.sum((factor_correlation_matrix_eigen_value_1s >= 0.0) & (factor_correlation_matrix_eigen_value_1s < 1.0))

					with open(specif_fa_dir + key + "_factor_correlation_matrix_evs.txt", "w") as text_file:
						print(f"factor_correlation_matrix_ev0s_greater_equal_1: {factor_correlation_matrix_ev0s_greater_equal_1}\nfactor_correlation_matrix_ev0s_less_than_1: {factor_correlation_matrix_ev0s_less_than_1}\nfactor_correlation_matrix_ev1s_greater_equal_1: {factor_correlation_matrix_ev1s_greater_equal_1}\nfactor_correlation_matrix_ev1s_less_than_1: {factor_correlation_matrix_ev1s_less_than_1}\n", file=text_file)


			# 1SD
			sw.plot_fa_comparison(list_angled_1sds, list_fa_1sds, dh.fa_varaince_comparisons_directory + "Multi__AngledCluster_V_ClusterFA__1SD_" + suf, 
								y_equals_x_max=12, optional_y_max=6)

			sw.plot_fa_comparison(list_flat_1sds, list_fa_1sds, dh.fa_varaince_comparisons_directory + "Multi__FlatCluster_V_ClusterFA__1SD_" + suf, 
								y_equals_x_max=12, optional_y_max=6)

			sw.plot_fa_comparison(list_conj_1sds, list_fa_1sds, dh.fa_varaince_comparisons_directory + "Multi__Conjunctive_V_ClusterFA__1SD_" + suf, 
								y_equals_x_max=18, optional_y_max=6)

			sw.plot_fa_comparison(list_not_only_conj_1sds, list_fa_1sds, dh.fa_varaince_comparisons_directory + "Multi__NotOnlyConjunctive_V_ClusterFA__1SD_" + suf, 
								y_equals_x_max=18, optional_y_max=6)


			# 1SD Differences
			sw.plot_fa_comparison(list_angled_diffs, list_fa_diffs, dh.fa_varaince_comparisons_directory + "Multi_AngledClusterDiff_V_ClusterFADiff__1SD_" + suf, 
								y_equals_x_max=12, optional_y_max=6, legend_labels=legend_labels)

			sw.plot_fa_comparison(list_flat_diffs, list_fa_diffs, dh.fa_varaince_comparisons_directory + "Multi_FlatClusterDiff_V_ClusterFADiff__1SD_" + suf, 
								y_equals_x_max=12, optional_y_max=6, legend_labels=legend_labels)
			
			sw.plot_fa_comparison(list_conj_diffs, list_fa_diffs, dh.fa_varaince_comparisons_directory + "Multi_ConjunctiveDiff_V_ClusterFADiff__1SD_" + suf, 
								y_equals_x_max=18, optional_y_max=6, legend_labels=legend_labels)
			
			sw.plot_fa_comparison(list_conj_diffs, list_angled_diffs, dh.fa_varaince_comparisons_directory + "Multi_ConjunctiveDiff_V_AngledClusterDiff__1SD_" + suf, 
								y_equals_x_max=18, optional_y_max=12, legend_labels=legend_labels)



	def extend_pairwise_reliabilities_and_plot_pairwise_reliabilities(self, both_spiking_reliabilities, number_of_conjunctive_trials, pairwise_reliabilities_directory, experiment_code):

		self.extend_pairwise_reliabilities(both_spiking_reliabilities, number_of_conjunctive_trials)
		self.plot_pairwise_reliability_plots(pairwise_reliabilities_directory, experiment_code)


	def extend_pairwise_reliabilities(self, both_spiking_reliabilities, number_of_conjunctive_trials):

		self.all_both_spiking_reliabilities.extend(both_spiking_reliabilities)
		self.all_number_of_conjunctive_trials.extend(number_of_conjunctive_trials)


	def plot_pairwise_reliability_plots(self, pairwise_reliabilities_directory, experiment_code):

		sw.normal_histo_plot([self.all_both_spiking_reliabilities], pairwise_reliabilities_directory + 'all_both_spiking_reliabilities_HIST_' + experiment_code, bins=20, histo_range=[0.0, 1.0], x_axis_label="Pairwise Reliability", y_axis_label="Count")
		sw.normal_histo_plot([self.all_number_of_conjunctive_trials], pairwise_reliabilities_directory + 'all_number_of_conjunctive_trials_HIST_' + experiment_code, bins=20, histo_range=[0.0, 1.0], x_axis_label="Number of trials", y_axis_label="Count")

		self.all_both_spiking_reliabilities_0s_removed = [value for value in self.all_both_spiking_reliabilities if value != 0.0]
		self.all_number_of_conjunctive_trials_0s_removed =  [value for value in self.all_number_of_conjunctive_trials if value != 0]

		sw.normal_histo_plot([self.all_both_spiking_reliabilities_0s_removed], pairwise_reliabilities_directory + 'all_both_spiking_reliabilities_0s_removed_HIST_' + experiment_code, bins=20, histo_range=[0.0, 1.0])
		sw.normal_histo_plot([self.all_number_of_conjunctive_trials_0s_removed], pairwise_reliabilities_directory + 'all_number_of_conjunctive_trials_0s_removed_HIST_' + experiment_code, bins=20, histo_range=[0.0, 1.0])

		sw.basic_x_y_plot([self.all_both_spiking_reliabilities], [self.all_number_of_conjunctive_trials], pairwise_reliabilities_directory + 'Reliability_Vs_NumberSpikingTrials_SCATTER_' + experiment_code, s=2)



def plot_multi_shuffle_type_plots(mcah_for_each_shuffle_type, analysis):

	shuffle_options_for_plot = ['normal-1', 'shuffle-1', 'sample-1']		
	p_values_for_each_shuffle_type = [results_holder.pdds['FlatClusterStats_FlatCluster_LR_pvalue'] for results_holder in mcah_for_each_shuffle_type.values() if results_holder.shuffle_option_string in shuffle_options_for_plot]
	sw.cumulative_histo_plot(p_values_for_each_shuffle_type, analysis['directory_holder'].prim_clus_corr_dir + "COLLATED_CUM_HISTO_", bins=200, x_axis_label="p-value", y_axis_label="Normalised\ncumulative sum", custom_x_tick_locators=[1.0, 0.2], add_chi_squared_text=True)

	shuffle_options_for_plot = ['normal-1', 'sample_correlated_cluster_pca_ellipse_rotated_to_45-1']		
	p_values_for_each_shuffle_type = [results_holder.sdds['Original_PCA_BS_empirical_pvalue_different_from_45'] for results_holder in mcah_for_each_shuffle_type.values() if results_holder.shuffle_option_string in shuffle_options_for_plot]
	sw.cumulative_histo_plot(p_values_for_each_shuffle_type, analysis['directory_holder'].angle_analysis_directory + "PCA_not_45_bootstrap_pvalue_from_sd_COLLATED_CUM_HISTO_", bins=20, x_axis_label="p-value", y_axis_label="Normalised\ncumulative sum", custom_x_tick_locators=[1.0, 0.2])
	sw.normal_histo_plot(p_values_for_each_shuffle_type, analysis['directory_holder'].angle_analysis_directory + "PCA_not_45_bootstrap_pvalue_from_sd_COLLATED_HISTO", bins=20, histo_range=[0.0, 1.0], x_axis_label="p-value", y_axis_label="Frequency", custom_x_tick_locators=[1.0, 0.2], custom_y_tick_locators=[20, 20], alpha=0.78, add_chi_squared_text=False)


	