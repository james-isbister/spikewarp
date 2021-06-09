'''
3_order_figures.py

As the analysis scripts produce a large number of figures, this script is provided to copy the figures that are included 
in the paper and supplementary figures into an easy to navigate directory. Moreover, the figure names of these copied figures 
include their figure number (i.e. Fig. 3a). The script can be run at any time although warnings will be given for figures 
which are yet to be generated. The script provides a table (below) that displays the correspondence between figures produced 
by the main analyses (left) and the figure numbers of the paper (right). This provides a useful way of understanding how the 
initial analysis scripts relate to the figures of the paper.



### Configuration
1. Update the variable "figs_root" on line 27 of this script to the path originally specified for the variable "figures_root_directory" in "recipes/investigation_recipe.json".



### Run
The script can be run with the following command:
python 3_order_figures.py 

'''

######################################
# CONFIGURATION
######################################
figs_root = [PATH_TO_BE_SET_BY_USER]


import os
import glob

import spikewarp as sw

def test_exists_and_move(original_fig_path, new_fig_path):
	if (os.path.exists(original_fig_path)):
		os.system("cp " + original_fig_path + " " + new_fig_path)
	else:
		print(original_fig_path + " doesn't exist")




######################################
# SETUP
######################################
pairwise_anal_figs_root = figs_root + "5.4-55.4/PairwiseAnalysis/"
timespan_pairwise_anal_figs_root = figs_root + "5.4-95.4/PairwiseAnalysis/"
cellwise_anal_figs_root = figs_root + "5.4-55.4/SingleNeuronAnalysis/"
orig_state_anal_figs_root = figs_root + "5.4-55.4/StateAnalysis/"

new_fig_root = figs_root + "OrderedFigs/"
fig_2_dir = new_fig_root + "1-Figure2/"
fig_3_dir = new_fig_root + "2-Figure3/"
fig_4_dir = new_fig_root + "3-Figure4/"
fig_5_dir = new_fig_root + "4-Figure5/"
fig_6_dir = new_fig_root + "5-Figure6/"
fig_7_dir = new_fig_root + "6-Figure7/"
fig_8_dir = new_fig_root + "7-Figure8/"
subfigs_fig_root = new_fig_root + "8-SuppFigs/"
normality_analysis_fig_dir = new_fig_root + "11-NormalityAnalysis/"


sw.makedirs(new_fig_root)
sw.makedirs(fig_2_dir)
sw.makedirs(fig_3_dir)
sw.makedirs(fig_4_dir)
sw.makedirs(fig_5_dir)
sw.makedirs(fig_6_dir)
sw.makedirs(fig_7_dir)
sw.makedirs(fig_8_dir)
sw.makedirs(subfigs_fig_root)
sw.makedirs(normality_analysis_fig_dir)

if hasattr(sw, "placeholder_9"):
	sw.placeholder_9(new_fig_root)



src_dest_pairs = [

###### Latex statistics
[pairwise_anal_figs_root + "CustomDBSCANExtra/_normal-1AnalysisOutputLatex.tex", 																												new_fig_root + "CustomDBSCANExtraAnalysisOutputLatex.tex"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/_sample_correlated_cluster_pca_ellipse_rotated_to_45-1AnalysisOutputLatex.tex", 																	new_fig_root + "CustomDBSCANExtra_sample_correlated_cluster_pca_ellipse_rotated_to_45-1AnalysisOutputLatex.tex"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/_45cluster_with_second_cluster-1AnalysisOutputLatex.tex", 																						new_fig_root + "CustomDBSCANExtra_45cluster_with_second_cluster-1AnalysisOutputLatex.tex"],
[pairwise_anal_figs_root + "Unclustered/_normal-1AnalysisOutputLatex.tex", 																														new_fig_root + "Unclustered_normal-1AnalysisOutputLatex.tex"],
[pairwise_anal_figs_root + "CustomGaussianMixtureInfCrit/_normal-1AnalysisOutputLatex.tex", 																									new_fig_root + "CustomGaussianMixtureInfCrit_normal-1AnalysisOutputLatex.tex"],
[pairwise_anal_figs_root + "CustomGaussianMixtureInfCrit/_sample_correlated_cluster_pca_ellipse_rotated_to_45-1AnalysisOutputLatex.tex", 														new_fig_root + "CustomGaussianMixtureInfCrit_sample_correlated_cluster_pca_ellipse_rotated_to_45-1AnalysisOutputLatex.tex"],
[pairwise_anal_figs_root + "CustomGaussianMixtureInfCrit/_45cluster_with_second_cluster-1AnalysisOutputLatex.tex", 																				new_fig_root + "CustomGaussianMixtureInfCrit_45cluster_with_second_cluster-1AnalysisOutputLatex.tex"],
[timespan_pairwise_anal_figs_root + "CustomDBSCANExtraTimeSpanCheck/AnalysisOutputLatexTimeSpan.tex",																							new_fig_root + "ClusterTimeSpan.tex"],


###### Supplementary Videos
[pairwise_anal_figs_root + "CustomDBSCANExtra/PairClusterPlots/normal-1/type_1_all_super_time_warp_cluster_buster_files.mp4", 																	subfigs_fig_root + "SupplementaryVideo1.mp4"],
[pairwise_anal_figs_root + "Unclustered/PairClusterPlots/normal-1/type_1_all_super_time_warp_cluster_buster_files.mp4", 																		subfigs_fig_root + "SupplementaryVideo2.mp4"],


###### FIGURE 2
[cellwise_anal_figs_root + "NonZeroSpikeCountHisto.pdf",																																		fig_2_dir + "Fig2B_NonZeroSpikeCountHisto_50ms_post_cort_onset.pdf"],
[cellwise_anal_figs_root + "atleast_one_spike_reliabilities.pdf",																																fig_2_dir + "Fig2C_atleast_one_spike_reliabilities_50ms_post_cort_onset.pdf"],


###### FIGURE 3
[pairwise_anal_figs_root + "CustomDBSCANExtra/PairClusterPlots/normal-1/normal-1_Clustering1.3499999999999996/type_2_STWCB_0_S001E026G001_6_U13_U17_normal-1_Clustering1.3499999999999996.pdf",	fig_3_dir + "Fig3A_ClusteringExample.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/PrimaryClusterCorrelations/_normal-1/PVal_LowResHist.pdf",																						fig_3_dir + "Fig3C1_CusDB_P_Stage1.pdf"],
[pairwise_anal_figs_root + "Unclustered/PrimaryClusterCorrelations/_normal-1/PVal_LowResHist.pdf", 																								fig_3_dir + "Fig3C2_Unclus_P_Stage1.pdf"],
[orig_state_anal_figs_root + "including_non_stationary_ClustVenn.pdf",																															fig_3_dir + "Fig3D1_including_non_stationary_ClustVenn.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/SecondaryClusterCorrelations/_normal-1/RVal_Hist.pdf",																							fig_3_dir + "Fig3D2_CusDB_R_Stage2.pdf"],
[pairwise_anal_figs_root + "Unclustered/SecondaryClusterCorrelations/_normal-1/RVal_Hist.pdf",																									fig_3_dir + "Fig3D3_Unclus_R_Stage2.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/Original/_normal-1/LRCorrelation/RSquared_Hist_with_unclustered_normal-1.pdf", 											fig_3_dir + "Fig3E1_R2_original_vs_same_unclustered.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/SecondaryClusterCorrelations/_normal-1/ClusPVal_Minus_ConjPVal_Hist.pdf", 																		fig_3_dir + "Fig3E2_R2_original_pval_decrease_by_clustering.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterReliabilities/_normal-1OverallAndConjHist.pdf", 																							fig_3_dir + "Fig3F_CLUSTER_RELIABILITY_original.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/Original/_normal-1/LRCorrelation/RSquared_Hist_AND_MultipleOutputLR_R^2_Hist_TIVsSpikePair_normal-1.pdf",					fig_3_dir + "Fig3G1_CusDB_R2_original_and_explained_by_trial.pdf"],
[pairwise_anal_figs_root + "Unclustered/ClusterNonStationarity/Original/_normal-1/LRCorrelation/RSquared_Hist_AND_MultipleOutputLR_R^2_Hist_TIVsSpikePair_normal-1.pdf",						fig_3_dir + "Fig3G2_Unclus_R2_original_and_explained_by_trial.pdf"],


###### FIGURE 4
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterPairDifferencesAnalysis/FS0_Vs_Diff_LR_PVal_LowResHist_normal-1.pdf",																		fig_4_dir + "Fig4A1_CusDB_P_Stage2_Difference_Vs_FirstSpike.pdf"],
[pairwise_anal_figs_root + "Unclustered/ClusterPairDifferencesAnalysis/FS0_Vs_Diff_LR_PVal_LowResHist_normal-1.pdf",																			fig_4_dir + "Fig4A2_Unclus_P_Stage2_Difference_Vs_FirstSpike.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/AngleAnalysis/Original/_normal-1/AngleCI_BS_PCA.pdf", 																							fig_4_dir + "Fig4B1_CusDB_ANGLE_original.pdf"],
[pairwise_anal_figs_root + "Unclustered/AngleAnalysis/Original/_normal-1/AngleCI_BS_PCA.pdf", 																									fig_4_dir + "Fig4B2_Unclus_ANGLE_original.pdf"],
[orig_state_anal_figs_root + "Confidence_interval_widths.pdf",																																	fig_4_dir + "Fig4C_confidence_interval_widths.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/OriginalTestsPassed/_normal-1/LRCorrelation/RSquared_Hist_with_unclustered_normal-1.pdf", 									fig_4_dir + "Fig4E1_R2_original_tests_passed_vs_same_unclustered.pdf"],
[pairwise_anal_figs_root + "Unclustered/ClusterNonStationarity/OriginalTestsPassed/_normal-1/LRCorrelation/RSquared_Hist_normal-1.pdf", 														fig_4_dir + "Fig4E2_unclus_R2_original_tests_passed.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/AngleAnalysis/OriginalTestsPassed/_normal-1/AngleCI_BS_PCA.pdf", 																					fig_4_dir + "Fig4F1_CusDB_ANGLE_original_tests_passed.pdf"],
[pairwise_anal_figs_root + "Unclustered/AngleAnalysis/OriginalTestsPassed/_normal-1/AngleCI_BS_PCA.pdf", 																						fig_4_dir + "Fig4F2_Unclus_ANGLE_original_tests_passed.pdf"],


###### FIGURE 5
[pairwise_anal_figs_root + "CustomDBSCANExtra/PairClusterPlots/normal-1/normal-1_Clustering3.499999999999999/type_2_STWCB_7_S001E045G001_3_U8_U14_normal-1_Clustering3.499999999999999.pdf", 	fig_5_dir + "Fig5B_EXAMPLE.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/FAVarianceComparisons/Multi_ConjunctiveDiff_V_AngledClusterDiff__1SD__normal-1.pdf", 																fig_5_dir + "Fig5C1_CusDB_ConjunctiveDiff_V_AngledClusterDiff__1SD.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/FAVarianceComparisons/Multi_AngledClusterDiff_V_ClusterFADiff__1SD__normal-1.pdf", 																fig_5_dir + "Fig5C2_CusDB_AngledClusterDiff_V_ClusterFADiff__1SD.pdf"],
[pairwise_anal_figs_root + "Unclustered/FAVarianceComparisons/Multi_AngledClusterDiff_V_ClusterFADiff__1SD__normal-1.pdf", 																		fig_5_dir + "Fig5D_Unclus_AngledClusterDiff_V_ClusterFADiff__1SD.pdf"],


###### FIGURE 6
[pairwise_anal_figs_root + "CustomDBSCANExtra/AngleAnalysis/SelectivelyDifferencedTestsPassedActuallyDifferenced/_normal-1/AngleCI_BS_PCA.pdf", 												fig_6_dir + "Fig6B1_CusDB_ANGLE_seldif_tests_passed.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/AngleAnalysis/SelectivelyDifferencedBoxJenkinsTestsPassed/_normal-1/AngleCI_BS_PCA.pdf", 															fig_6_dir + "Fig6B2_CusDB_ANGLE_boxjenk_tests_passed.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/SelectivelyDifferencedBoxJenkins/_normal-1/Differenced_ARIMA.pdf",															fig_6_dir + "Fig6C_CusDB_DifferencedAndARIMA_TestsPassedCounts.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/SelectivelyDifferencedTestsPassedActuallyDifferenced/_normal-1/LRCorrelation/RSquared_Hist_normal-1.pdf",					fig_6_dir + "Fig6D1_CusDB_Differenced_tests_passed_explained_var.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/SelectivelyDifferencedBoxJenkinsTestsPassed/_normal-1/LRCorrelation/RSquared_Hist_normal-1.pdf",							fig_6_dir + "Fig6D2_CusDB_ARIMA_tests_passed_explained_var.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/UnitTypeAnalysis/stimulation_frequency_grouped_normal-1.pdf", 																					fig_6_dir + "Fig6E_CusDB_ClusterStimFreq.pdf"],


###### FIGURE 7
[pairwise_anal_figs_root + "CustomDBSCANExtra/UnitTypeAnalysis/MergedOrigAndBoxJenkCortexPlot_CorrelatedClusters_normal-1.pdf",																	fig_7_dir + "Fig7A_CusDB_ClusterPairLocations.pdf"],
[orig_state_anal_figs_root + "Conj/PairsOfClusters/StateComparison/quadrupletTimesAndStateOriginalE2_C3_Conj_137_1819.pdf",																		fig_7_dir + "Fig7B-pair_of_clusters_example_E2_C3_Conj_137_1819.pdf"],
[orig_state_anal_figs_root + "conj_dbscan_and_unclustered_cluster_pairs_state_corr_pvals.pdf",																									fig_7_dir + "Fig7C1-dbscan_and_unclustered_cluster_pairs_state_corr_pvals.pdf"],
[orig_state_anal_figs_root + "Conj/MetaStats/StateCorr/cluster_pairs_state_corr_rvals_quadruplets.pdf",																							fig_7_dir + "Fig7C2-conj_cluster_pairs_state_corr_rvals.pdf"],
[orig_state_anal_figs_root + "Conj/MetaStats/PCA/mean_exp_var_ratios_quadruplets.pdf", 																											fig_7_dir + "Fig7D-PCA_ExpVar_Pool.pdf"],


###### FIGURE 8
[orig_state_anal_figs_root + "Conj/MetaStats/StateCorr/sch_single_unit_prediction_from_state_pvalues.pdf",																						fig_8_dir + "8A1-sch_single_unit_prediction_from_state_pvalues.pdf"],
[orig_state_anal_figs_root + "Conj/MetaStats/StateCorr/sch_single_unit_prediction_from_state_rsquared_values.pdf",																				fig_8_dir + "8A2-sch_single_unit_prediction_from_state_rsquared_values.pdf"],
[orig_state_anal_figs_root + "Conj/MetaStats/DiffCorr/difference_prediction_from_state_pvalues.pdf",																							fig_8_dir + "8B1-NON45_ConjOnly_difference_prediction_from_state_pvalues.pdf"],
[orig_state_anal_figs_root + "Conj/MetaStats/DiffCorr/difference_prediction_from_state_rsquared_values.pdf",																					fig_8_dir + "8B2-NON45_ConjOnly_difference_prediction_from_state_rsquared_values.pdf"],
[orig_state_anal_figs_root + "Conj/MetaStats/StateCorr/_ClusterSD_vs_StatePredicErrorSD.pdf",																									fig_8_dir + "8C1-SCHUnit_ClusterSD_vs_StatePredicErrorSD.pdf"],
[orig_state_anal_figs_root + "Conj/MetaStats/StateCorr/_FAPredicSD_vs_StatePredicErrorSD.pdf",																									fig_8_dir + "8C2-SCHUnit_FAPredicSD_vs_StatePredicErrorSD.pdf"],
[orig_state_anal_figs_root + "Conj/MetaStats/DiffCorr/_ClusterSD_vs_StatePredicErrorSD.pdf",																									fig_8_dir + "8D1-Diff_ClusterSD_vs_StatePredicErrorSD.pdf"],
[orig_state_anal_figs_root + "Conj/MetaStats/DiffCorr/_FAPredicSD_vs_StatePredicErrorSD.pdf",																									fig_8_dir + "8D2-Diff_FAPredicSD_vs_StatePredicErrorSD.pdf"],


###### SUP FIG 1
[cellwise_anal_figs_root + "0AllByFreqGroupspike_counts_by_trial_index_by_freq_group.pdf",																										subfigs_fig_root + "S1_0AllByFreqGroupspike_counts_by_trial_index_by_freq_group.pdf"],
[cellwise_anal_figs_root + "1AllByFreqGroupspike_counts_by_trial_index_by_freq_group.pdf",																										subfigs_fig_root + "S1_1AllByFreqGroupspike_counts_by_trial_index_by_freq_group.pdf"],
[cellwise_anal_figs_root + "2AllByFreqGroupspike_counts_by_trial_index_by_freq_group.pdf",																										subfigs_fig_root + "S1_2AllByFreqGroupspike_counts_by_trial_index_by_freq_group.pdf"],
[cellwise_anal_figs_root + "0AllByFreqGroupspike_times_by_trial_index_by_freq_group.pdf",																										subfigs_fig_root + "S1_0AllByFreqGroupspike_spikes_by_trial_index_by_freq_group.pdf"],
[cellwise_anal_figs_root + "1AllByFreqGroupspike_times_by_trial_index_by_freq_group.pdf",																										subfigs_fig_root + "S1_1AllByFreqGroupspike_spikes_by_trial_index_by_freq_group.pdf"],
[cellwise_anal_figs_root + "2AllByFreqGroupspike_times_by_trial_index_by_freq_group.pdf",																										subfigs_fig_root + "S1_2AllByFreqGroupspike_spikes_by_trial_index_by_freq_group.pdf"],


###### SUP FIG 2
[cellwise_anal_figs_root + "AllNeuronGroupspike_times_by_trial_index_by_freq_group.pdf",																										subfigs_fig_root + "S2A_All_NeuronGroupspike_times_by_trial_index_by_freq_group.pdf"],
[cellwise_anal_figs_root + "AllNeuronGroupspike_counts_by_trial_index_by_freq_group.pdf",																										subfigs_fig_root + "S2B_All_NeuronGroupspike_counts_by_trial_index_by_freq_group.pdf"],


###### SUP FIG 3
[cellwise_anal_figs_root + "LatenciesVsTrialPlots/1/1_22_1.pdf",																																subfigs_fig_root + "S3A_SingleUnitSingleConditionExample1.pdf"],
[cellwise_anal_figs_root + "LatenciesVsTrialPlots/1/1_22_5.pdf",																																subfigs_fig_root + "S3B_SingleUnitSingleConditionExample2.pdf"],
[cellwise_anal_figs_root + "LatenciesVsTrialPlots/1/1_48_5.pdf",																																subfigs_fig_root + "S3C_SingleUnitSingleConditionExample3.pdf"],
[cellwise_anal_figs_root + "LatenciesVsTrialPlots/1/1_35_6.pdf",																																subfigs_fig_root + "S3D_SingleUnitSingleConditionExample4.pdf"],
[cellwise_anal_figs_root + "LatenciesVsTrialPlots/1/1_4_13.pdf",																																subfigs_fig_root + "S3E_SingleUnitSingleConditionExample5.pdf"],
[cellwise_anal_figs_root + "LatenciesVsTrialPlots/1/1_29_13.pdf",																																subfigs_fig_root + "S3F_SingleUnitSingleConditionExample6.pdf"],
[cellwise_anal_figs_root + "LatenciesVsTrialPlots/1/1_20_17.pdf",																																subfigs_fig_root + "S3G_SingleUnitSingleConditionExample7.pdf"],
[cellwise_anal_figs_root + "LatenciesVsTrialPlots/1/1_20_18.pdf",																																subfigs_fig_root + "S3H_SingleUnitSingleConditionExample8.pdf"],


###### SUP FIG 4
[cellwise_anal_figs_root + "_ACF_PVal_CumHist.pdf", 																																			subfigs_fig_root + "S4A_SingleUnitACF.pdf"],
[cellwise_anal_figs_root + "_PACF_PVal_CumHist.pdf", 																																			subfigs_fig_root + "S4B_SingleUnitPACF.pdf"],


###### SUP FIG 5
[pairwise_anal_figs_root + "CustomDBSCANExtra/PrimaryClusterCorrelations/COLLATED_CUM_HISTO_.pdf", 																								subfigs_fig_root + "S5-1-A4_CUM_P_PrimaryCorrelations_CustomDBSCANExtra.pdf"],
[pairwise_anal_figs_root + "Unclustered/PrimaryClusterCorrelations/COLLATED_CUM_HISTO_.pdf", 																									subfigs_fig_root + "S5_2-A_CUM_P_PrimaryCorrelations_Unclustered.pdf"],
[pairwise_anal_figs_root + "SciRepGMMInfCrit/PrimaryClusterCorrelations/COLLATED_CUM_HISTO_.pdf", 																								subfigs_fig_root + "S5_2-B_CUM_P_PrimaryCorrelations_SciRepGMMInfCrit.pdf"],
[pairwise_anal_figs_root + "SciRepGMMInfCrit/PairClusterPlots/shuffle-1/shuffle-1_Clustering10/type_2_STWCB_4_S001E050G001_14_U5_U15_shuffle-1_Clustering10.pdf",								subfigs_fig_root + "S5_2-C_GMMExampleDetectionOfCorrelatedClusterThroughFragmentation.pdf"],
[pairwise_anal_figs_root + "CustomGaussianMixtureInfCrit/PrimaryClusterCorrelations/COLLATED_CUM_HISTO_.pdf", 																					subfigs_fig_root + "S5_3-A_CUM_P_PrimaryCorrelations_CustomGaussianMixtureInfCrit.pdf"],
[pairwise_anal_figs_root + "CustomGaussianMixtureInfCrit/ClusterReliabilities/_normal-1OverallAndConjHist.pdf",																					subfigs_fig_root + "S5_3-B_CustomGMM_CLUSTER_RELIABILITY_original.pdf"],
[pairwise_anal_figs_root + "CustomGaussianMixtureInfCrit/AngleAnalysis/Original/_normal-1/AngleCI_BS_PCA.pdf", 																					subfigs_fig_root + "S5_3-C1_CustomGMM_ANGLE_original.pdf"],
[pairwise_anal_figs_root + "CustomGaussianMixtureInfCrit/AngleAnalysis/OriginalTestsPassed/_normal-1/AngleCI_BS_PCA.pdf",																		subfigs_fig_root + "S5_3-C2_CustomGMM_original_tests_passed.pdf"],


###### SUP FIG 6
[pairwise_anal_figs_root + "CustomDBSCANExtra/PairClusterPlots/normal-1/normal-1_Clustering1.1999999999999997/type_0_STWCB_2_S001E075F009_0_U43_U48_normal-1_Clustering1.1999999999999997.pdf", subfigs_fig_root + "S6_1-A1_45_pca_clus_example_BEFORE.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/PairClusterPlots/sample_correlated_cluster_pca_ellipse_rotated_to_45-1/sample_correlated_cluster_pca_ellipse_rotated_to_45-1_Clustering1.1999999999999997_Original1.1999999999999997/type_0_STWCB_2_S001E075F009_0_U43_U48_sample_correlated_cluster_pca_ellipse_rotated_to_45-1_Clustering1.1999999999999997_Original1.1999999999999997.pdf", subfigs_fig_root + "S6_1-A2_45_pca_clus_example_AFTER.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/AngleAnalysis/Original/_sample_correlated_cluster_pca_ellipse_rotated_to_45-1/AngleCI_BS_PCA.pdf",												subfigs_fig_root + "S6_1-B_ANGLE_sample_correlated_cluster_pca_ellipse_rotated_to_45.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/AngleAnalysis/Original/_sample_correlated_cluster_pca_ellipse_rotated_to_45-1/Diff45_BS_Sd_PVal_Hist.pdf",										subfigs_fig_root + "S6_1-C_ANGLE_DIFF_FROM_45_PVAL_sample_correlated_cluster_pca_ellipse_rotated_to_45-1.pdf"],

[pairwise_anal_figs_root + "CustomDBSCANExtra/PairClusterPlots/45cluster_with_second_cluster-1/45cluster_with_second_cluster-1_Clustering0.75_Original1.15/type_0_STWCB_4_S001E050G001_17_U15_U37_45cluster_with_second_cluster-1_Clustering0.75_Original1.15.pdf", subfigs_fig_root + "S6_2-A_double_clus_example.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/AngleAnalysis/Original/_45cluster_with_second_cluster-1/AngleCI_BS_PCA.pdf",																		subfigs_fig_root + "S6_2-B_ANGLE_45cluster_with_second_cluster.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/AngleAnalysis/Original/_45cluster_with_second_cluster-1/Diff45_BS_Sd_PVal_Hist.pdf",																subfigs_fig_root + "S6_2-C_ANGLE_DIFF_FROM_45_PVAL_45cluster_with_second_cluster.pdf"],
[pairwise_anal_figs_root + "Unclustered/AngleAnalysis/Original/_45cluster_with_second_cluster-1/AngleCI_BS_PCA.pdf",																			subfigs_fig_root + "S6_2-D_Unclus_ANGLE_45cluster_with_second_cluster.pdf"],

[pairwise_anal_figs_root + "CustomGaussianMixtureInfCrit/AngleAnalysis/Original/_sample_correlated_cluster_pca_ellipse_rotated_to_45-1/AngleCI_BS_PCA.pdf",										subfigs_fig_root + "S6-3-A_CustomGMM_ANGLE_sample_correlated_cluster_pca_ellipse_rotated_to_45.pdf"],
[pairwise_anal_figs_root + "CustomGaussianMixtureInfCrit/AngleAnalysis/Original/_45cluster_with_second_cluster-1/AngleCI_BS_PCA.pdf",															subfigs_fig_root + "S6-3-B_CustomGMM_ANGLE_45cluster_with_second_cluster.pdf"],


###### SUP FIG 7
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/Original/_normal-1/LRStationarity/LR_PVal_Hist_TIVsSpikeTimes_normal-1.pdf",												subfigs_fig_root + "S7Atop_clus_P_LR_ClusterSingleUnitSpikeTimesWithTrialIndex.pdf"],
[pairwise_anal_figs_root + "Unclustered/ClusterNonStationarity/Original/_normal-1/LRStationarity/LR_PVal_Hist_TIVsSpikeTimes_normal-1.pdf",														subfigs_fig_root + "S7Abot_unclus_P_LR_ClusterSingleUnitSpikeTimesWithTrialIndex.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/Original/_normal-1/LRStationarity/LR_PVal_U0_Vs_U1_TIVsSpikeTimes_normal-1.pdf",											subfigs_fig_root + "S7Btop_clus_P_LR_ClusterSingleUnitVsSingleUnitSpikeTimes.pdf"],
[pairwise_anal_figs_root + "Unclustered/ClusterNonStationarity/Original/_normal-1/LRStationarity/LR_PVal_U0_Vs_U1_TIVsSpikeTimes_normal-1.pdf",													subfigs_fig_root + "S7Bbot_unclus_P_LR_ClusterSingleUnitVsSingleUnitSpikeTimes.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/Original/_normal-1/KPSSStationarity/KPSS_pvalues_normal-1.pdf",															subfigs_fig_root + "S7Ctop_P_clus_KPSS_ClusterSingleUnitSpikeTimes.pdf"],
[pairwise_anal_figs_root + "Unclustered/ClusterNonStationarity/Original/_normal-1/KPSSStationarity/KPSS_pvalues_normal-1.pdf",																	subfigs_fig_root + "S7Cbot_P_unclus_KPSS_ClusterSingleUnitSpikeTimes.pdf"],


###### SUP FIG 8
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/OriginalTestsPassed/_normal-1/LagPVals/5Lags/OriginalTestsPassed_ACF_PVal_CumHist_normal-1.pdf",							subfigs_fig_root + "S8-1-A_ACF_P_original_tests_passed_cluster_units.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/OriginalTestsPassed/_normal-1/LagPVals/5Lags/OriginalTestsPassed_PACF_PVal_CumHist_normal-1.pdf",							subfigs_fig_root + "S8-1-B_PACF_P_original_tests_passed_cluster_units.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/OriginalTestsPassed/_normal-1/LagPVals/5Lags/OriginalTestsPassed_CCF_PVal_CumHist_normal-1.pdf",							subfigs_fig_root + "S8-1-C_CACF_P_original_tests_passed_cluster_units.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/OriginalTestsPassed/_normal-1/LagPVals/5Lags/OriginalTestsPassed_PCCF_PVal_CumHist_normal-1.pdf",							subfigs_fig_root + "S8-1-D_PCCF_P_original_tests_passed_cluster_units.pdf"],

[pairwise_anal_figs_root + "Unclustered/ClusterNonStationarity/OriginalTestsPassed/_normal-1/LagPVals/5Lags/OriginalTestsPassed_ACF_PVal_CumHist_normal-1.pdf",									subfigs_fig_root + "S8-2-A_Unclustered_ACF_P_original_tests_passed_cluster_units_.pdf"],
[pairwise_anal_figs_root + "Unclustered/ClusterNonStationarity/OriginalTestsPassed/_normal-1/LagPVals/5Lags/OriginalTestsPassed_PACF_PVal_CumHist_normal-1.pdf",								subfigs_fig_root + "S8-2-B_Unclustered_PACF_P_original_tests_passed_cluster_units.pdf"],
[pairwise_anal_figs_root + "Unclustered/ClusterNonStationarity/OriginalTestsPassed/_normal-1/LagPVals/5Lags/OriginalTestsPassed_CCF_PVal_CumHist_normal-1.pdf",									subfigs_fig_root + "S8-2-C_Unclustered_CACF_P_original_tests_passed_cluster_units.pdf"],
[pairwise_anal_figs_root + "Unclustered/ClusterNonStationarity/OriginalTestsPassed/_normal-1/LagPVals/5Lags/OriginalTestsPassed_PCCF_PVal_CumHist_normal-1.pdf",								subfigs_fig_root + "S8-2-D_Unclustered_PCCF_P_original_tests_passed_cluster_units.pdf"],

[pairwise_anal_figs_root + "CustomDBSCANExtra/PairClusterPlots/normal-1/normal-1_Clustering0.4/type_2_STWCB_6_S001E067G001_12_U14_U20_normal-1_Clustering0.4.pdf",								subfigs_fig_root + "S8-3_DifferencingNeccessaryExample.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/PairClusterPlots/normal-1/normal-1_Clustering1.0499999999999998/type_2_STWCB_6_S001E067G001_12_U19_U20_normal-1_Clustering1.0499999999999998.pdf",subfigs_fig_root + "S8-3_DifferencingNeccessaryExample2.pdf"],


###### SUP FIG 9
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/SelectivelyDifferencedTestsPassedActuallyDifferenced/_normal-1/LagPVals/5Lags/SelectivelyDifferencedTestsPassedActuallyDifferenced_ACF_PVal_CumHist_normal-1.pdf",  		subfigs_fig_root + "S9-1-A_ACF_P_differenced_tests_passed.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/SelectivelyDifferencedTestsPassedActuallyDifferenced/_normal-1/LagPVals/5Lags/SelectivelyDifferencedTestsPassedActuallyDifferenced_PACF_PVal_CumHist_normal-1.pdf", 		subfigs_fig_root + "S9-1-B_PACF_P_differenced_tests_passed.pdf"],

[pairwise_anal_figs_root + "Unclustered/ClusterNonStationarity/SelectivelyDifferencedTestsPassedActuallyDifferenced/_normal-1/LagPVals/5Lags/SelectivelyDifferencedTestsPassedActuallyDifferenced_ACF_PVal_CumHist_normal-1.pdf",  		subfigs_fig_root + "S9-2-A_Unclustered_ACF_P_differenced_tests_passed.pdf"],
[pairwise_anal_figs_root + "Unclustered/ClusterNonStationarity/SelectivelyDifferencedTestsPassedActuallyDifferenced/_normal-1/LagPVals/5Lags/SelectivelyDifferencedTestsPassedActuallyDifferenced_PACF_PVal_CumHist_normal-1.pdf", 		subfigs_fig_root + "S9-2-B_Unclustered_PACF_P_differenced_tests_passed.pdf"],

[pairwise_anal_figs_root + "Unclustered/AngleAnalysis/SelectivelyDifferencedTestsPassedActuallyDifferenced/_normal-1/AngleCI_BS_PCA.pdf", 														subfigs_fig_root + "S9-2-A1_equiv6B1_Unclus_ANGLE_seldif_tests_passed.pdf"],
[pairwise_anal_figs_root + "Unclustered/AngleAnalysis/SelectivelyDifferencedBoxJenkinsTestsPassed/_normal-1/AngleCI_BS_PCA.pdf", 																subfigs_fig_root + "S9-2-A2_equiv6B2_Unclus_ANGLE_boxjenk_tests_passed.pdf"],
[pairwise_anal_figs_root + "Unclustered/ClusterNonStationarity/SelectivelyDifferencedBoxJenkins/_normal-1/Differenced_ARIMA.pdf",																subfigs_fig_root + "S9-2-B_equiv6C_Unclus_DifferencedAndARIMA_TestsPassedCounts.pdf"],
[pairwise_anal_figs_root + "Unclustered/ClusterNonStationarity/SelectivelyDifferencedTestsPassedActuallyDifferenced/_normal-1/LRCorrelation/RSquared_Hist_normal-1.pdf",						subfigs_fig_root + "S9-2-C1_equiv6D1_Unclus_Differenced_tests_passed_explained_var.pdf"],
[pairwise_anal_figs_root + "Unclustered/ClusterNonStationarity/SelectivelyDifferencedBoxJenkinsTestsPassed/_normal-1/LRCorrelation/RSquared_Hist_normal-1.pdf",									subfigs_fig_root + "S9-2-C2_equiv6D2_Unclus_ARIMA_tests_passed_explained_var.pdf"],
[pairwise_anal_figs_root + "Unclustered/UnitTypeAnalysis/stimulation_frequency_grouped_normal-1.pdf", 																							subfigs_fig_root + "S9-2-D_equiv6E_Unclus_ClusterStimFreq.pdf"],



###### SUP FIG 10
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/SelectivelyDifferencedBoxJenkinsTestsPassed/_normal-1/LagPVals/5Lags/SelectivelyDifferencedBoxJenkinsTestsPassed_ACF_PVal_CumHist_normal-1.pdf",  subfigs_fig_root + "S10-1-A_ACF_P_arima_tests_passed.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/SelectivelyDifferencedBoxJenkinsTestsPassed/_normal-1/LagPVals/5Lags/SelectivelyDifferencedBoxJenkinsTestsPassed_CCF_PVal_CumHist_normal-1.pdf", subfigs_fig_root + "S10-1-B_CCF_P_arima_tests_passed.pdf"],

[pairwise_anal_figs_root + "Unclustered/ClusterNonStationarity/SelectivelyDifferencedBoxJenkinsTestsPassed/_normal-1/LagPVals/5Lags/SelectivelyDifferencedBoxJenkinsTestsPassed_ACF_PVal_CumHist_normal-1.pdf",  subfigs_fig_root + "S10-2-A_Unclustered_ACF_P_arima_tests_passed.pdf"],
[pairwise_anal_figs_root + "Unclustered/ClusterNonStationarity/SelectivelyDifferencedBoxJenkinsTestsPassed/_normal-1/LagPVals/5Lags/SelectivelyDifferencedBoxJenkinsTestsPassed_CCF_PVal_CumHist_normal-1.pdf", subfigs_fig_root + "S10-2-B_Unclustered_CCF_P_arima_tests_passed.pdf"],


###### SUP FIG 11
[pairwise_anal_figs_root + "CustomDBSCANExtra/AngleAnalysis/OriginalTestsPassed/_normal-1/Angle_BS_PCA_Vs_FA.pdf",																				subfigs_fig_root + "S11-1-A_ANGLE_PCA_vs_FA_original_tests_passed.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/AngleAnalysis/SelectivelyDifferencedTestsPassedActuallyDifferenced/_normal-1/OriginalAngle_Vs_Angle_BS_PCA.pdf", 									subfigs_fig_root + "S11-1-B_ANGLE_PCA_before_and_after_differencing_for_actually_differenced_tests_passed.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/AngleAnalysis/SelectivelyDifferencedTestsPassedActuallyDifferenced/_normal-1/Angle_BS_PCA_Vs_FA.pdf",												subfigs_fig_root + "S11-1-C_ANGLE_PCA_vs_FA_actually_differenced_tests_passed.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/AngleAnalysis/SelectivelyDifferencedBoxJenkinsTestsPassed/_normal-1/OriginalAngle_Vs_Angle_BS_PCA.pdf", 											subfigs_fig_root + "S11-1-D_ANGLE_PCA_before_and_after_differencing_for_arima_tests_passed.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/AngleAnalysis/SelectivelyDifferencedBoxJenkinsTestsPassed/_normal-1/Angle_BS_PCA_Vs_FA.pdf",														subfigs_fig_root + "S11-1-E_ANGLE_PCA_vs_FA_arima_tests_passed.pdf"],


###### SUP FIG 11-2
[pairwise_anal_figs_root + "Unclustered/UnitTypeAnalysis/MergedOrigAndBoxJenkCortexPlot_CorrelatedClusters_normal-1.pdf",																		subfigs_fig_root + "S11-2_equiv7A_Unclus_ClusterPairLocations.pdf"],


###### SUP FIG 12
[figs_root + "minus50-150/StimulusSpikePSTHOverAllTrialsLarger.pdf",																															subfigs_fig_root + "S12-1_StimulusSpikePSTHOverAllTrialsLarger.pdf"],
[figs_root + "minus50-150/StimulusSpikePSTHOverAllTrialsSmaller.pdf",																															subfigs_fig_root + "S12-2_StimulusSpikePSTHOverAllTrialsLarger.pdf"],


###### SUP FIG 13
[pairwise_anal_figs_root + "CustomDBSCANExtra/PairClusterPlots/normal-1/normal-1_Clustering1.3499999999999996/type_0_STWCB_0_S001E026G001_6_U13_U17_normal-1_Clustering1.3499999999999996.pdf",	subfigs_fig_root + "S13_ClusteringExample.pdf"],


###### SUP FIG 14
[timespan_pairwise_anal_figs_root + "CustomDBSCANExtraTimeSpanCheck/ClusterTimeSpans/LimitsOfFlatClustersForAngledClustersOnly_normal-1.pdf",													subfigs_fig_root + "S14_ClusterTimeSpans.pdf"],


###### Normality analysis
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/Original/_normal-1/HZStationarity/Henze-Zirkler_pvalues_normal-1.pdf", 													normality_analysis_fig_dir + "FigXa_CusDB_Henze-Zirkler_Original.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/OriginalTestsPassed/_normal-1/HZStationarity/Henze-Zirkler_pvalues_normal-1.pdf", 											normality_analysis_fig_dir + "FigXb_CusDB_Henze-Zirkler_OriginalTestsPassed.pdf"],
[pairwise_anal_figs_root + "CustomDBSCANExtra/ClusterNonStationarity/SelectivelyDifferencedBoxJenkinsTestsPassed/_normal-1/HZStationarity/Henze-Zirkler_pvalues_normal-1.pdf", 					normality_analysis_fig_dir + "FigXc_CusDB_Henze-Zirkler_SelectivelyDifferencedBoxJenkinsTestsPassed.pdf"],

[pairwise_anal_figs_root + "Unclustered/ClusterNonStationarity/Original/_normal-1/HZStationarity/Henze-Zirkler_pvalues_normal-1.pdf", 															normality_analysis_fig_dir + "FigXd_Unclus_Henze-Zirkler_Original.pdf"],
[pairwise_anal_figs_root + "Unclustered/ClusterNonStationarity/OriginalTestsPassed/_normal-1/HZStationarity/Henze-Zirkler_pvalues_normal-1.pdf", 												normality_analysis_fig_dir + "FigXe_Unclus_Henze-Zirkler_OriginalTestsPassed.pdf"],
[pairwise_anal_figs_root + "Unclustered/ClusterNonStationarity/SelectivelyDifferencedBoxJenkinsTestsPassed/_normal-1/HZStationarity/Henze-Zirkler_pvalues_normal-1.pdf", 						normality_analysis_fig_dir + "FigXf_Unclus_Henze-Zirkler_SelectivelyDifferencedBoxJenkinsTestsPassed.pdf"],




]


if hasattr(sw, "placeholder_list"):
	src_dest_pairs + sw.placeholder_list


for [src, dest] in src_dest_pairs:
	test_exists_and_move(src, dest)
