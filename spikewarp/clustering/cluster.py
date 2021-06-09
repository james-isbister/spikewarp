import sys
import math
import warnings

from sklearn.decomposition import PCA
from sklearn.decomposition import FactorAnalysis
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score # Temp used to avoid deprecation warning

from scipy.stats import multivariate_normal
from scipy import stats
from scipy.stats import linregress
from scipy.stats import shapiro
from scipy.spatial.distance import cdist

from statsmodels.tsa.stattools import adfuller
from statsmodels.tsa.stattools import kpss
from statsmodels.tsa.arima_model import ARIMA

import numpy as np

import pingouin as pg

from factor_analyzer import FactorAnalyzer
from factor_analyzer.factor_analyzer import calculate_bartlett_sphericity
from factor_analyzer.factor_analyzer import calculate_kmo


import spikewarp as sw




"""
Class for defining cluster and calculating cluster stats
"""

class Cluster(object):

	def __init__(self):
		
		self.analysis_dicts = {}


	def create_from_pairs(self, single_cluster_holder_options, original_spike_pairs, spike_pair_distribution_spikes, number_trials_both_spiking, number_of_trials_for_condition, real_trial_indices_on_which_both_neurons_spiked):

		mean = [np.mean(original_spike_pairs[:, 0]), np.mean(original_spike_pairs[:, 1])]

		flat_gauss_cov_mat = np.zeros((2, 2), dtype=float)
		flat_gauss_cov_mat[0, 0] = np.var(original_spike_pairs[:, 0])
		flat_gauss_cov_mat[1, 1] = np.var(original_spike_pairs[:, 1])

		self.create_from_covariance_matrix(single_cluster_holder_options, mean, flat_gauss_cov_mat, spike_pair_distribution_spikes, number_trials_both_spiking, number_of_trials_for_condition, real_trial_indices_on_which_both_neurons_spiked)


	def create_nonflat_cluster_from_cluster_spikes(self, single_cluster_holder_options, original_spike_pairs, original_cluster_indices, number_trials_both_spiking, number_of_trials_for_condition, real_trial_indices_on_which_both_neurons_spiked, make_differencing_calculations=False):

		pca = PCA(n_components=2)
		pca.fit_transform(original_spike_pairs)
		covariance_matrix = pca.get_covariance()
		mean = [np.mean(original_spike_pairs[:, 0]), np.mean(original_spike_pairs[:, 1])]

		self.create_from_final_spikes_and_covariance(single_cluster_holder_options, original_spike_pairs, original_cluster_indices, mean, covariance_matrix, number_trials_both_spiking, number_of_trials_for_condition, real_trial_indices_on_which_both_neurons_spiked, make_differencing_calculations=make_differencing_calculations)


	def create_from_covariance_matrix(self, single_cluster_holder_options, mean, covariance_matrix, spike_pair_distribution_spikes, number_trials_both_spiking, number_of_trials_for_condition, real_trial_indices_on_which_both_neurons_spiked, make_differencing_calculations=False):

		if ((covariance_matrix[0][0] == 0.0) | (covariance_matrix[1][1] == 0.0)):
			return False
		
		mahalanobis_distances = cdist(spike_pair_distribution_spikes, [mean], metric='mahalanobis', VI=np.linalg.inv(covariance_matrix))	
		indices_within_outer_bound = np.argwhere(mahalanobis_distances <= single_cluster_holder_options.outlier_bound).T[0]
		cluster_indices = indices_within_outer_bound
		cluster_spike_pairs = spike_pair_distribution_spikes[cluster_indices, :]

		self.create_from_final_spikes_and_covariance(single_cluster_holder_options, cluster_spike_pairs, cluster_indices, mean, covariance_matrix, number_trials_both_spiking, number_of_trials_for_condition, real_trial_indices_on_which_both_neurons_spiked, make_differencing_calculations=make_differencing_calculations)


	def create_from_final_spikes_and_covariance(self, single_cluster_holder_options, cluster_spike_pairs, cluster_indices, mean, covariance_matrix, number_trials_both_spiking, number_of_trials_for_condition, real_trial_indices_on_which_both_neurons_spiked, make_differencing_calculations=False):

		self.cluster_spike_pairs = cluster_spike_pairs
		self.cluster_indices = cluster_indices
		self.real_trial_indices_of_cluster = []
		if (real_trial_indices_on_which_both_neurons_spiked != []):
			self.real_trial_indices_of_cluster = real_trial_indices_on_which_both_neurons_spiked[self.cluster_indices]
		self.number_of_samples_in_ellipse = len(cluster_indices)
		self.mean = np.asarray(mean)
		self.covariance_matrix = covariance_matrix

		self.cluster_spike_pairs_differences = self.cluster_spike_pairs[:, 1] - self.cluster_spike_pairs[:, 0]
		self.cluster_spike_pairs_differences_sd = np.std(self.cluster_spike_pairs_differences)

		self.ellipse_conj_reliability = float(self.number_of_samples_in_ellipse) / float(number_trials_both_spiking)
		self.ellipse_overall_reliability = float(self.number_of_samples_in_ellipse) / float(number_of_trials_for_condition)


		self.unangled_ellipses = []
		for ring in single_cluster_holder_options.rings:
			ellr = sw.create_ellipse(self.mean, math.sqrt(self.covariance_matrix[0][0])*ring, math.sqrt(self.covariance_matrix[1][1])*ring, 0.0)
			self.unangled_ellipses.append(ellr)


		self.cluster_spike_pairs_N0_sd = np.std(self.cluster_spike_pairs[:, 0])
		self.cluster_spike_pairs_N1_sd = np.std(self.cluster_spike_pairs[:, 1])


		if (make_differencing_calculations):
			# Useful links:
				# https://www.machinelearningplus.com/time-series/arima-model-time-series-forecasting-python/
				# http://people.duke.edu/~rnau/Notes_on_nonseasonal_ARIMA_models--Robert_Nau.pdf
				# https://www.machinelearningplus.com/time-series/arima-model-time-series-forecasting-python/

			self.number_of_lags = 12
			self.use_selective_differences = True; self.was_box_jeckinsed = False; self.was_auto_arimad = False
			self.AR_p_0 = -1; self.AR_p_1 = -1
			self.MA_q_0 = -1; self.MA_q_1 = -1
			self.index_to_use_0 = -1; self.index_to_use_1 = -1

			self.analysis_dicts['Original'] = self.create_analysis_dictionary_for_spike_pairs(self.cluster_spike_pairs, self.cluster_indices, 'Original', create_additional_dicts=True)
			self.pairs_differenced = (self.cluster_spike_pairs[1:, :] - self.cluster_spike_pairs[:-1, :])
			self.analysis_dicts['Differenced'] = self.create_analysis_dictionary_for_spike_pairs(self.pairs_differenced, self.cluster_indices[1:], 'Differenced')
		
			self.index_to_use_0 = sw.decisions_for_selective_differencing(self.analysis_dicts['Original']['TI_Vs_STs_LR_0_pvalue'], self.analysis_dicts['Original']['KPSS_STs_0_pvalue'], self.analysis_dicts['Original']['STs_acf_pvalues_0'][0], self.analysis_dicts['Differenced']['TI_Vs_STs_LR_0_pvalue'], self.analysis_dicts['Differenced']['KPSS_STs_0_pvalue'])
			self.index_to_use_1 = sw.decisions_for_selective_differencing(self.analysis_dicts['Original']['TI_Vs_STs_LR_1_pvalue'], self.analysis_dicts['Original']['KPSS_STs_1_pvalue'], self.analysis_dicts['Original']['STs_acf_pvalues_1'][0], self.analysis_dicts['Differenced']['TI_Vs_STs_LR_1_pvalue'], self.analysis_dicts['Differenced']['KPSS_STs_1_pvalue'])
			selective_difference_indices = [self.index_to_use_0, self.index_to_use_1]
			
			if ((-1 in selective_difference_indices) | (1 not in selective_difference_indices)):
				self.use_selective_differences = False

				keys_for_empty_analysis_dict = ['SelectivelyDifferenced', 'SelectivelyDifferencedTestsPassedAndNormal', 'SelectivelyDifferencedTestsPassed', 'SelectivelyDifferencedTestsPassedActuallyDifferenced', 'SelectivelyDifferencedBoxJenkins', 'SelectivelyDifferencedBoxJenkinsTestsPassedAndNormal', 'SelectivelyDifferencedBoxJenkinsTestsPassed']
				
				for key in keys_for_empty_analysis_dict:
					self.analysis_dicts[key] = self.create_empty_analysis_dictionary_for_spike_pairs(self.number_of_lags)

				
			if (self.use_selective_differences):
				self.set_spike_pairs_for_selective_differences()
				self.analysis_dicts['SelectivelyDifferenced'] = self.create_analysis_dictionary_for_spike_pairs(self.selective_differences_scaled_down, self.cluster_indices[self.starting_index_for_selective_difference_indices:], 'SelectivelyDifferenced', create_additional_dicts=True, selective_difference_indices=selective_difference_indices)
				
				self.box_jenkins_selective_differences()
				self.analysis_dicts['SelectivelyDifferencedBoxJenkins'] = self.create_analysis_dictionary_for_spike_pairs(self.selective_differences_arima_residuals, self.cluster_indices[self.starting_index_for_selective_difference_indices:], 'SelectivelyDifferencedBoxJenkins', create_additional_dicts=True, selective_difference_indices=selective_difference_indices)

		return True





	def box_jenkins_selective_differences(self):

		# BOX JENKINS SELCTIVE DIFFERENCES
		self.AR_p_0, self.MA_q_0, self.use_arima_model_0 = sw.calculate_arima_variables_for_acf_and_pacf_postive_arrays(self.analysis_dicts['SelectivelyDifferenced']['STs_acf_rvalues_0'], self.analysis_dicts['SelectivelyDifferenced']['STs_acf_pvalues_0'], self.analysis_dicts['SelectivelyDifferenced']['STs_pacf_pvalues_0'], self.analysis_dicts['SelectivelyDifferenced']['STs_acf_positive_pvalues_0'], self.analysis_dicts['SelectivelyDifferenced']['STs_pacf_positive_pvalues_0'])
		self.AR_p_1, self.MA_q_1, self.use_arima_model_1 = sw.calculate_arima_variables_for_acf_and_pacf_postive_arrays(self.analysis_dicts['SelectivelyDifferenced']['STs_acf_rvalues_1'], self.analysis_dicts['SelectivelyDifferenced']['STs_acf_pvalues_1'], self.analysis_dicts['SelectivelyDifferenced']['STs_pacf_pvalues_1'], self.analysis_dicts['SelectivelyDifferenced']['STs_acf_positive_pvalues_1'], self.analysis_dicts['SelectivelyDifferenced']['STs_pacf_positive_pvalues_1'])
		self.selective_differences_arima_residuals = np.copy(self.selective_differences)


		if (self.use_arima_model_0):
			self.was_box_jeckinsed = True
			model_0 = ARIMA(self.selective_differences[:, 0], order=(self.AR_p_0, 0, self.MA_q_0))
			with warnings.catch_warnings():
				warnings.filterwarnings("ignore")
				try:
					model_fit_0 = model_0.fit(disp=0, transparams=False)
					self.selective_differences_arima_residuals[:, 0] = model_fit_0.resid
				except:
					print((self.AR_p_0, 0, self.MA_q_0))
			

		if (self.use_arima_model_1):
			self.was_box_jeckinsed = True
			with warnings.catch_warnings():
				warnings.filterwarnings("ignore")
				model_1 = ARIMA(self.selective_differences[:, 1], order=(self.AR_p_1, 0, self.MA_q_1))
				try:
					model_fit_1 = model_1.fit(disp=0, transparams=False)
					self.selective_differences_arima_residuals[:, 1] = model_fit_1.resid
				except:
					print((self.AR_p_1, 0, self.MA_q_1))
				



	def set_spike_pairs_for_selective_differences(self):

		if ((self.index_to_use_0 == 0) & (self.index_to_use_1 == 0)):
			self.selective_differences = self.cluster_spike_pairs
			self.selective_differences_scaled_down = self.cluster_spike_pairs
			self.number_of_samples_for_selective_differences = self.number_of_samples_in_ellipse
		

		elif ((self.index_to_use_0 <= 1) & (self.index_to_use_1 <= 1)):
			self.selective_differences = np.copy(self.cluster_spike_pairs[1:, :])
			self.selective_differences_scaled_down = np.copy(self.cluster_spike_pairs[1:, :])
			self.number_of_samples_for_selective_differences = self.number_of_samples_in_ellipse - 1

			if (self.index_to_use_0 == 1):
				self.selective_differences[:, 0] = self.pairs_differenced[:, 0]
				# self.selective_differences_scaled_down[:, 0] = self.pairs_differenced[:, 0] / 2.0
				self.selective_differences_scaled_down[:, 0] = self.pairs_differenced[:, 0] / math.sqrt(2.0)

			if (self.index_to_use_1 == 1):
				self.selective_differences[:, 1] = self.pairs_differenced[:, 1]
				# self.selective_differences_scaled_down[:, 1] = self.pairs_differenced[:, 1] / 2.0
				self.selective_differences_scaled_down[:, 1] = self.pairs_differenced[:, 1] / math.sqrt(2.0)


		self.starting_index_for_selective_difference_indices = self.number_of_samples_in_ellipse - self.number_of_samples_for_selective_differences


	def create_stage_1_cluster_statistics_dictionary(self):

		FlatCluster_LR_object = linregress(self.cluster_spike_pairs[:, 0], self.cluster_spike_pairs[:, 1])

		FlatCluster_N0_FS_SD = np.std(self.cluster_spike_pairs[:, 0])
		FlatCluster_N1_FS_SD = np.std(self.cluster_spike_pairs[:, 1])

		Cluster_Diff_LR_object = linregress(self.cluster_spike_pairs[:, 0], self.cluster_spike_pairs[:, 1] - self.cluster_spike_pairs[:, 0])		

		self.analysis_dicts['FlatClusterStats'] = {}
		self.analysis_dicts['FlatClusterStats']['FlatCluster_LR_object'] = FlatCluster_LR_object
		self.analysis_dicts['FlatClusterStats']['FlatCluster_LR_pvalue'] = FlatCluster_LR_object.pvalue
		self.analysis_dicts['FlatClusterStats']['FlatCluster_LR_rsquared'] = FlatCluster_LR_object.rvalue**2
		self.analysis_dicts['FlatClusterStats']['FlatCluster_LR_rvalue'] = FlatCluster_LR_object.rvalue
		self.analysis_dicts['FlatClusterStats']['FlatCluster_N0_FS_SD'] = FlatCluster_N0_FS_SD
		self.analysis_dicts['FlatClusterStats']['FlatCluster_N1_FS_SD'] = FlatCluster_N1_FS_SD
		self.analysis_dicts['FlatClusterStats']['FlatCluster_1sd_area'] = FlatCluster_N0_FS_SD * FlatCluster_N1_FS_SD
		self.analysis_dicts['FlatClusterStats']['FlatCluster_FS_diff_SD'] = np.std(self.cluster_spike_pairs[:, 1] - self.cluster_spike_pairs[:, 0])
		self.analysis_dicts['FlatClusterStats']['FlatCluster_FS_Mean0'] = np.mean(self.cluster_spike_pairs[:, 0])
		self.analysis_dicts['FlatClusterStats']['FlatCluster_FS_Mean1'] = np.mean(self.cluster_spike_pairs[:, 1])
		self.analysis_dicts['FlatClusterStats']['FlatCluster_Diff_LR_pvalue'] = Cluster_Diff_LR_object.pvalue


	def create_empty_analysis_dictionary_for_spike_pairs(self, number_of_lags):

		analysis_dict = {}

		analysis_dict['is_empty'] = True

		analysis_dict['TI_Vs_STs_LR_0'] = None
		analysis_dict['TI_Vs_STs_LR_1'] = None
		analysis_dict['TI_Vs_STs_LR_0_pvalue'] = -1.0
		analysis_dict['TI_Vs_STs_LR_1_pvalue'] = -1.0

		analysis_dict['TI_Vs_STs_LR_multiple_output'] = None
		analysis_dict['TI_Vs_STs_LR_multiple_output_rsquared'] = -1.0

		analysis_dict['ADFuller_STs_0_pvalue'] = -1.0
		analysis_dict['ADFuller_STs_1_pvalue'] = -1.0

		analysis_dict['STs_acf_rvalues_0'], analysis_dict['STs_acf_pvalues_0'], analysis_dict['STs_acf_positive_pvalues_0'], analysis_dict['STs_acf_negative_pvalues_0'] = [-2.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)]
		analysis_dict['STs_acf_rvalues_1'], analysis_dict['STs_acf_pvalues_1'], analysis_dict['STs_acf_positive_pvalues_1'], analysis_dict['STs_acf_negative_pvalues_1'] = [-2.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)]
		
		analysis_dict['STs_ccf_rvalues_0'], analysis_dict['STs_ccf_pvalues_0'], analysis_dict['STs_ccf_positive_pvalues_0'], analysis_dict['STs_ccf_negative_pvalues_0'] = [-2.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)]
		analysis_dict['STs_ccf_rvalues_1'], analysis_dict['STs_ccf_pvalues_1'], analysis_dict['STs_ccf_positive_pvalues_1'], analysis_dict['STs_ccf_negative_pvalues_1'] = [-2.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)]

		analysis_dict['STs_pacf_rvalues_0'], analysis_dict['STs_pacf_pvalues_0'], analysis_dict['STs_pacf_positive_pvalues_0'], analysis_dict['STs_pacf_negative_pvalues_0'] = [-2.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)]
		analysis_dict['STs_pacf_rvalues_1'], analysis_dict['STs_pacf_pvalues_1'], analysis_dict['STs_pacf_positive_pvalues_1'], analysis_dict['STs_pacf_negative_pvalues_1'] = [-2.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)]

		analysis_dict['STs_pccf_rvalues_0'], analysis_dict['STs_pccf_pvalues_0'], analysis_dict['STs_pccf_positive_pvalues_0'], analysis_dict['STs_pccf_negative_pvalues_0'] = [-2.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)]
		analysis_dict['STs_pccf_rvalues_1'], analysis_dict['STs_pccf_pvalues_1'], analysis_dict['STs_pccf_positive_pvalues_1'], analysis_dict['STs_pccf_negative_pvalues_1'] = [-2.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)], [-1.0 for i in range(number_of_lags)]

		analysis_dict['distances_from_factor_line'] = []
		analysis_dict['PCA_predicited_state_for_cluster_trials'] = []
		analysis_dict['BS_PCA_mean_angle_up_to_45'] = -1.0
		analysis_dict['BS_PCA_mean_angle'] = -1.0
		analysis_dict['PCA_BS_empirical_CI_lower'] = -1.0
		analysis_dict['PCA_BS_empirical_CI_upper'] = -1.0
		analysis_dict['PCA_BS_empirical_pvalue_different_from_45'] = -1.0
		analysis_dict['is_PCA_BS_empirical_pvalue_different_from_45'] = False
		analysis_dict['PCA_BS_empirical_pvalue_different_from_0'] = -1.0
		analysis_dict['is_PCA_BS_empirical_pvalue_different_from_0'] = False
		analysis_dict['BS_PCA_different_from_45_sd_method'] = -1.0
		analysis_dict['is_BS_PCA_different_from_45_sd_method'] = False
		analysis_dict['BS_PCA_mean_component_0_sd'] = -1.0
		analysis_dict['BS_PCA_mean_component_1_sd'] = -1.0
		analysis_dict['BS_PCA_mean_of_mean0'] = -1.0
		analysis_dict['BS_PCA_mean_of_mean1'] = -1.0

		analysis_dict['FA_predicited_state_for_cluster_trials'] = []
		analysis_dict['FA_angle_BS_mean'] = -1.0
		analysis_dict['FA_N0_sd_BS_mean'] = -1.0
		analysis_dict['FA_N1_sd_BS_mean'] = -1.0
		analysis_dict['FA_1sd_area_BS_mean'] = -1.0
		analysis_dict['FA_1sd_estimated_difference_area_BS_mean'] = -1.0
		analysis_dict['FA_mean_0_mean'] = -1.0
		analysis_dict['FA_mean_1_mean'] = -1.0
		analysis_dict['FA_component_0_mean'] = -1.0
		analysis_dict['FA_component_1_mean'] = -1.0
		analysis_dict['FA_noise_variance_0_mean'] = -1.0
		analysis_dict['FA_noise_variance_1_mean'] = -1.0
		analysis_dict['FA_angle_BS_mean_45d'] = -1.0
		analysis_dict['FA_angle_BS_empirical_CI'] = -1.0

		analysis_dict['sharipo_normality_p_0'] = -1.0
		analysis_dict['sharipo_normality_p_1'] = -1.0

		analysis_dict['KPSS_STs_0_pvalue'] = -1.0
		analysis_dict['KPSS_STs_1_pvalue'] = -1.0

		analysis_dict['henze-zirkler_multivariate_normality_p'] = -1.0

		analysis_dict['bartlett_spherecity_p_value'] = -1.0

		analysis_dict['LR'] = None
		analysis_dict['LR_pvalue'] = -1.0
		analysis_dict['LR_rvalue'] = -1.1
		analysis_dict['LR_rsquared'] = -1.1
		analysis_dict['is_still_correlated'] = False

		analysis_dict['factor_correlation_matrix_eigen_value_0'] = -1.0
		analysis_dict['factor_correlation_matrix_eigen_value_1'] = -1.0

		analysis_dict['tests_passed'] = False
		analysis_dict['tests_passed_and_normal'] = False
		analysis_dict['tests_passed_but_not_normal'] = False
		analysis_dict['selective_differences_undifferenced'] = False
		analysis_dict['selective_differences_differenced'] = False

		analysis_dict['index_to_use_0'] = -1
		analysis_dict['index_to_use_1'] = -1

		analysis_dict['was_box_jeckinsed'] = False
		analysis_dict['was_not_box_jeckinsed'] = False
		analysis_dict['AR_p_0'] = -1.0
		analysis_dict['MA_q_0'] = -1.0
		analysis_dict['AR_p_1'] = -1.0
		analysis_dict['MA_q_1'] = -1.0


		return analysis_dict

	def create_analysis_dictionary_for_spike_pairs(self, spike_pairs_to_use, indices_to_use, analysis_key, create_additional_dicts=False, selective_difference_indices=[]):

		analysis_dict = {}

		analysis_dict['is_empty'] = False

		analysis_dict['Mean_0'] = np.mean(spike_pairs_to_use[:, 0])
		analysis_dict['Mean_1'] = np.mean(spike_pairs_to_use[:, 1])


		analysis_dict['TI_Vs_STs_LR_0'] = linregress(indices_to_use, spike_pairs_to_use[: , 0])
		analysis_dict['TI_Vs_STs_LR_1'] = linregress(indices_to_use, spike_pairs_to_use[: , 1])


		
		
		analysis_dict['TI_Vs_STs_LR_0_pvalue'] = analysis_dict['TI_Vs_STs_LR_0'].pvalue
		analysis_dict['TI_Vs_STs_LR_1_pvalue'] = analysis_dict['TI_Vs_STs_LR_1'].pvalue

		analysis_dict['TI_Vs_STs_LR_multiple_output'] = LinearRegression().fit(indices_to_use.reshape(-1, 1), spike_pairs_to_use)
		analysis_dict['TI_Vs_STs_LR_multiple_output_rsquared'] = r2_score(spike_pairs_to_use, analysis_dict['TI_Vs_STs_LR_multiple_output'].predict(indices_to_use.reshape(-1, 1)), multioutput='variance_weighted')

		adfuller_0 = adfuller(spike_pairs_to_use[:, 0])
		adfuller_1 = adfuller(spike_pairs_to_use[:, 1])
		analysis_dict['ADFuller_STs_0_pvalue'] = adfuller_0[1]
		analysis_dict['ADFuller_STs_1_pvalue'] = adfuller_1[1]

		with warnings.catch_warnings():
			warnings.filterwarnings("ignore")

			kpss_0 = kpss(spike_pairs_to_use[:, 0], nlags='auto')
			kpss_1 = kpss(spike_pairs_to_use[:, 1], nlags='auto')
			analysis_dict['KPSS_STs_0_pvalue'] = kpss_0[1]
			analysis_dict['KPSS_STs_1_pvalue'] = kpss_1[1]



		analysis_dict['STs_acf_rvalues_0'], analysis_dict['STs_acf_pvalues_0'], analysis_dict['STs_acf_positive_pvalues_0'], analysis_dict['STs_acf_negative_pvalues_0'] = sw.calculate_acf_pvalues_for_spikes(spike_pairs_to_use[:, 0], spike_pairs_to_use[:, 0], self.number_of_lags)
		analysis_dict['STs_acf_rvalues_1'], analysis_dict['STs_acf_pvalues_1'], analysis_dict['STs_acf_positive_pvalues_1'], analysis_dict['STs_acf_negative_pvalues_1'] = sw.calculate_acf_pvalues_for_spikes(spike_pairs_to_use[:, 1], spike_pairs_to_use[:, 1], self.number_of_lags)

		analysis_dict['STs_ccf_rvalues_0'], analysis_dict['STs_ccf_pvalues_0'], analysis_dict['STs_ccf_positive_pvalues_0'], analysis_dict['STs_ccf_negative_pvalues_0'] = sw.calculate_acf_pvalues_for_spikes(spike_pairs_to_use[:, 0], spike_pairs_to_use[:, 1], self.number_of_lags)
		analysis_dict['STs_ccf_rvalues_1'], analysis_dict['STs_ccf_pvalues_1'], analysis_dict['STs_ccf_positive_pvalues_1'], analysis_dict['STs_ccf_negative_pvalues_1'] = sw.calculate_acf_pvalues_for_spikes(spike_pairs_to_use[:, 1], spike_pairs_to_use[:, 0], self.number_of_lags)

		analysis_dict['STs_pccf_rvalues_0'], analysis_dict['STs_pccf_pvalues_0'], analysis_dict['STs_pccf_positive_pvalues_0'], analysis_dict['STs_pccf_negative_pvalues_0'] = sw.calculate_pacf_pvalues_for_spikes(spike_pairs_to_use[:, 0], spike_pairs_to_use[:, 1], self.number_of_lags)
		analysis_dict['STs_pccf_rvalues_1'], analysis_dict['STs_pccf_pvalues_1'], analysis_dict['STs_pccf_positive_pvalues_1'], analysis_dict['STs_pccf_negative_pvalues_1'] = sw.calculate_pacf_pvalues_for_spikes(spike_pairs_to_use[:, 1], spike_pairs_to_use[:, 0], self.number_of_lags)

		analysis_dict['STs_pacf_rvalues_0'], analysis_dict['STs_pacf_pvalues_0'], analysis_dict['STs_pacf_positive_pvalues_0'], analysis_dict['STs_pacf_negative_pvalues_0'] = sw.calculate_pacf_pvalues_for_spikes(spike_pairs_to_use[:, 0], spike_pairs_to_use[:, 0], self.number_of_lags)
		analysis_dict['STs_pacf_rvalues_1'], analysis_dict['STs_pacf_pvalues_1'], analysis_dict['STs_pacf_positive_pvalues_1'], analysis_dict['STs_pacf_negative_pvalues_1'] = sw.calculate_pacf_pvalues_for_spikes(spike_pairs_to_use[:, 1], spike_pairs_to_use[:, 1], self.number_of_lags)

		analysis_dict['BS_PCA_mean_angle_up_to_45'] = -1.0
		analysis_dict['PCA_BS_empirical_CI_lower'] = -1.0
		analysis_dict['PCA_BS_empirical_CI_upper'] = -1.0
		analysis_dict['PCA_BS_empirical_pvalue_different_from_45'] = -1.0
		analysis_dict['is_PCA_BS_empirical_pvalue_different_from_45'] = False
		analysis_dict['PCA_BS_empirical_pvalue_different_from_0'] = -1.0
		analysis_dict['is_PCA_BS_empirical_pvalue_different_from_0'] = False
		analysis_dict['BS_PCA_different_from_45_sd_method'] = -1.0
		analysis_dict['is_BS_PCA_different_from_45_sd_method'] = False
		analysis_dict['BS_PCA_mean_of_mean0'] = -1.0
		analysis_dict['BS_PCA_mean_of_mean1'] = -1.0
		

		stat, analysis_dict['sharipo_normality_p_0'] = shapiro(spike_pairs_to_use[: , 0])
		stat, analysis_dict['sharipo_normality_p_1'] = shapiro(spike_pairs_to_use[: , 1])

		k2, analysis_dict['henze-zirkler_multivariate_normality_p'], _ = pg.multivariate_normality(spike_pairs_to_use, alpha=.05)

		analysis_dict['normal'] = False
		if (analysis_dict['henze-zirkler_multivariate_normality_p'] >= 0.05):
			analysis_dict['normal'] = True
			


		analysis_dict['LR'] = linregress(spike_pairs_to_use[:, 0], spike_pairs_to_use[:, 1])
		
		
		analysis_dict['LR_pvalue'] = analysis_dict['LR'].pvalue
		analysis_dict['LR_rvalue'] = analysis_dict['LR'].rvalue
		analysis_dict['LR_rsquared'] = analysis_dict['LR'].rvalue**2
		analysis_dict['is_still_correlated'] = analysis_dict['LR_pvalue'] < 0.005

		chi_square_value, analysis_dict['bartlett_spherecity_p_value'] = calculate_bartlett_sphericity(spike_pairs_to_use)

		fa = FactorAnalyzer()
		fa.set_params(n_factors=2, rotation=None)
		fa.fit(self.cluster_spike_pairs)
		factor_correlation_matrix_eigen_values, v = fa.get_eigenvalues()
		analysis_dict['factor_correlation_matrix_eigen_value_0'] = factor_correlation_matrix_eigen_values[0]
		analysis_dict['factor_correlation_matrix_eigen_value_1'] = factor_correlation_matrix_eigen_values[1]

		analysis_dict['index_to_use_0'] = self.index_to_use_0
		analysis_dict['index_to_use_1'] = self.index_to_use_1

		analysis_dict['was_box_jeckinsed'] = self.was_box_jeckinsed
		analysis_dict['was_not_box_jeckinsed'] = not self.was_box_jeckinsed

		analysis_dict['AR_p_0'] = self.AR_p_0
		analysis_dict['MA_q_0'] = self.MA_q_0
		analysis_dict['AR_p_1'] = self.AR_p_1
		analysis_dict['MA_q_1'] = self.MA_q_1

		analysis_dict['tests_passed'] = False
		analysis_dict['tests_passed_and_normal'] = False
		analysis_dict['tests_passed_but_not_normal'] = False

		autocorrelation_criteria_pvalue_threshold = 0.05

		if (analysis_key == "SelectivelyDifferenced"):
			autocorrelation_criteria_pvalue_threshold = 0.0
		if ((analysis_dict['LR_pvalue'] < 0.005) & (analysis_dict['LR'].rvalue > 0.3) & (analysis_dict['bartlett_spherecity_p_value'] < 0.05) & (analysis_dict['TI_Vs_STs_LR_multiple_output_rsquared'] < 0.05) & (analysis_dict['TI_Vs_STs_LR_0_pvalue'] > 0.05) & (analysis_dict['TI_Vs_STs_LR_1_pvalue'] > 0.05) & (analysis_dict['KPSS_STs_0_pvalue'] > 0.05) & (analysis_dict['KPSS_STs_1_pvalue'] > 0.05) & (analysis_dict['STs_acf_pvalues_0'][0] > autocorrelation_criteria_pvalue_threshold) & (analysis_dict['STs_acf_pvalues_1'][0] > autocorrelation_criteria_pvalue_threshold)):
		
			analysis_dict['tests_passed'] = True
		
			if (analysis_dict['henze-zirkler_multivariate_normality_p'] >= 0.05):
				analysis_dict['tests_passed_and_normal'] = True
			else:
				analysis_dict['tests_passed_but_not_normal'] = True



		analysis_dict['selective_differences_undifferenced'] = False
		analysis_dict['selective_differences_differenced'] = False		
		if (len(selective_difference_indices) > 0):
			if ((selective_difference_indices[0] == 0) & (selective_difference_indices[1] == 0)):
				analysis_dict['selective_differences_undifferenced'] = True
			else:
				analysis_dict['selective_differences_differenced'] = True


		if (create_additional_dicts):
			
			if (analysis_dict['tests_passed']):
				self.analysis_dicts[analysis_key + 'TestsPassed'] = analysis_dict.copy()
				if (analysis_key == "SelectivelyDifferenced"):
					if (analysis_dict['selective_differences_differenced']):
						self.analysis_dicts[analysis_key + 'TestsPassedActuallyDifferenced'] = analysis_dict.copy()
					else:
						self.analysis_dicts[analysis_key + 'TestsPassedActuallyDifferenced'] = self.create_empty_analysis_dictionary_for_spike_pairs(self.number_of_lags)

					
			else: 
				self.analysis_dicts[analysis_key + 'TestsPassed'] = self.create_empty_analysis_dictionary_for_spike_pairs(self.number_of_lags)
				if (analysis_key == "SelectivelyDifferenced"):
					self.analysis_dicts[analysis_key + 'TestsPassedActuallyDifferenced'] = self.create_empty_analysis_dictionary_for_spike_pairs(self.number_of_lags)

			if (analysis_dict['tests_passed_and_normal']):
				self.analysis_dicts[analysis_key + 'TestsPassedAndNormal'] = analysis_dict.copy()
			else: 
				self.analysis_dicts[analysis_key + 'TestsPassedAndNormal'] = self.create_empty_analysis_dictionary_for_spike_pairs(self.number_of_lags)
				



		return analysis_dict

				