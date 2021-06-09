import warnings

import spikewarp as sw


"""
ANGLE HELPERS
Helper functions for calculating bootstrap correlation angles and confidence intervals
"""

import math
import numpy as np
from sklearn.decomposition import PCA

def run_pca_and_get_positive_first_component_angle(spikes):

	pca = PCA(n_components=2)
	pca.fit_transform(spikes)
	pca_angle_unadjusted = np.arctan2([pca.components_[0][1]], [pca.components_[0][0]]) * 180 / np.pi
	pca_angle = pca_angle_unadjusted
	if (pca_angle < 0.0):
		pca_angle  = pca_angle + 180.0


	pca_angle_half = pca_angle
	pca_angle_half_expanded_to_360 = pca_angle_half * 2.0
	pca_angle_half_expanded_to_360_radians = math.radians(pca_angle_half_expanded_to_360)
	pca_angle_half_expanded_to_360_x_y = [math.cos(pca_angle_half_expanded_to_360_radians), math.sin(pca_angle_half_expanded_to_360_radians)]

	return pca, pca_angle[0], pca_angle_unadjusted[0], pca_angle_half_expanded_to_360_x_y

def empirical_confidence_bound_new(all_angles):

    alpha = 0.95

    p_lower = ((1.0-alpha)/2.0) * 100
    lower = np.percentile(all_angles, p_lower)
    
    p_upper = (alpha+((1.0-alpha)/2.0)) * 100
    upper = np.percentile(all_angles, p_upper)

    return lower, upper


def empirical_confidence_bound(all_angles):

    alpha = 0.95

    p_lower = ((1.0-alpha)/2.0) * 100
    lower = np.percentile(all_angles, p_lower)
    
    p_upper = (alpha+((1.0-alpha)/2.0)) * 100
    upper = np.percentile(all_angles, p_upper)

    return lower - np.mean(all_angles), upper - np.mean(all_angles)


def convert_expanded_360_xy_to_180_angle_degrees(expanded_points_x, expanded_points_y):

	angle = convert_expanded_360_xy_to_360_angle_degrees(expanded_points_x, expanded_points_y)
	angle = angle / 2.0

	return angle

def convert_expanded_360_xy_to_360_angle_degrees(expanded_points_x, expanded_points_y):

	angle = convert_xy_to_degrees_basic(expanded_points_x, expanded_points_y)
	where_angle_less_than_0 = np.where(angle < 0.0)
	angle[where_angle_less_than_0] = angle[where_angle_less_than_0] + 360.0

	# if (angle < 0.0):
	# 	angle = angle + 360.0
	return angle

def convert_xy_to_degrees_basic(expanded_points_x, expanded_points_y):
	angle = np.arctan2(expanded_points_y, expanded_points_x) * 180 / np.pi
	return angle


def dotproduct(v1, v2):
	return sum((a*b) for a, b in zip(v1, v2))

def length(v):
	return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
	return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

def angle_degrees(v1, v2):
	return math.degrees(angle(v1, v2))

def signed_angle_degrees(v1, v2):
	return math.degrees(math.atan2(v2[1], v2[0]) - math.atan2(v1[1], v1[0]))



"""
BOOTSTRAP PCA 
Helper function for bootstrap PCA
"""

import numpy as np
from scipy import stats
import math
from shapely.geometry.point import Point
from shapely import affinity

def create_ellipse(center, width, height, angle):
   	
	circ = Point(center).buffer(1)
	ell = affinity.scale(circ, width, height)
	ellr = affinity.rotate(ell, angle)
	return ellr

def bootstrap_PCA_body(spikes_to_use, number_of_bootstrap_iterations, single_cluster_holder_options=None):

	number_of_samples_per_bootstrap = spikes_to_use.shape[0]

	BS_PCA_angle_half_expanded_to_360_x_ys = []
	BS_PCA_component_0s = []
	BS_PCA_component_1s = []
	BS_PCA_component_0_sds = []
	BS_PCA_component_1_sds = []
	BS_PCA_mean_0s = []
	BS_PCA_mean_1s = []
	BS_PCA_covariance_matrices = []

	for bootstrap_iteration in range(number_of_bootstrap_iterations):
		bootstrap_sample_indices = np.random.choice(number_of_samples_per_bootstrap, number_of_samples_per_bootstrap)
		bootstrap_sample = spikes_to_use[bootstrap_sample_indices]

		pca, _, _, pca_angle_half_expanded_to_360_x_y = sw.run_pca_and_get_positive_first_component_angle(bootstrap_sample)

		BS_PCA_angle_half_expanded_to_360_x_ys.append(pca_angle_half_expanded_to_360_x_y)
		BS_PCA_component_0_sds.append(np.sqrt(pca.explained_variance_[0]))
		BS_PCA_component_1_sds.append(np.sqrt(pca.explained_variance_[1]))
		BS_PCA_mean_0s.append(pca.mean_[0])
		BS_PCA_mean_1s.append(pca.mean_[1])
		BS_PCA_covariance_matrices.append(pca.get_covariance())



	BS_PCA_angle_half_expanded_to_360_x_y_mean = np.mean(BS_PCA_angle_half_expanded_to_360_x_ys, axis=0)
	BS_PCA_angle_back_to_180_mean = sw.convert_expanded_360_xy_to_180_angle_degrees([BS_PCA_angle_half_expanded_to_360_x_y_mean[0]], [BS_PCA_angle_half_expanded_to_360_x_y_mean[1]])[0]
	BS_PCA_angle_half_expanded_to_360_x_ys = np.asarray(BS_PCA_angle_half_expanded_to_360_x_ys)

	angle_differences_to_45_which_is_90_expanded = []
	for i in range(number_of_bootstrap_iterations):
		angle_differences_to_45_which_is_90_expanded.append(sw.signed_angle_degrees(BS_PCA_angle_half_expanded_to_360_x_ys[i, :], [0.0, 1.0]))

	angle_differences_to_0_which_is_0_expanded = []
	for i in range(number_of_bootstrap_iterations):
		angle_differences_to_0_which_is_0_expanded.append(sw.signed_angle_degrees(BS_PCA_angle_half_expanded_to_360_x_ys[i, :], [0.0, 0.0]))


	angles_from_mean_half_expanded_to_360 = []
	for i in range(number_of_bootstrap_iterations):
		angle_from_mean_half_expanded_to_360 = sw.signed_angle_degrees(BS_PCA_angle_half_expanded_to_360_x_y_mean, BS_PCA_angle_half_expanded_to_360_x_ys[i,:])
		angles_from_mean_half_expanded_to_360.append(angle_from_mean_half_expanded_to_360)
	std_of_angles_from_mean_half_expanded_to_360 = np.std(angles_from_mean_half_expanded_to_360)
	std_of_angles_from_mean_back_to_half = std_of_angles_from_mean_half_expanded_to_360 / 2.0

	######################################################################################################

	BS_PCA_mean_component_0_sd = np.mean(BS_PCA_component_0_sds)
	BS_PCA_mean_component_1_sd = np.mean(BS_PCA_component_1_sds)
	BS_PCA_mean_of_mean0 = np.mean(BS_PCA_mean_0s)
	BS_PCA_mean_of_mean1 = np.mean(BS_PCA_mean_1s)
	BS_PCA_mean_sd_area = BS_PCA_mean_component_0_sd * BS_PCA_mean_component_1_sd
	BS_PCA_mean_covariance_matrix = np.mean(BS_PCA_covariance_matrices, axis=0)


	PCA_BS_empirical_CI_lower_expanded, PCA_BS_empirical_CI_upper_expanded  = sw.empirical_confidence_bound_new(angles_from_mean_half_expanded_to_360)
	PCA_BS_empirical_CI_lower = abs(PCA_BS_empirical_CI_lower_expanded / 2.0)
	PCA_BS_empirical_CI_upper = abs(PCA_BS_empirical_CI_upper_expanded / 2.0)

	PCA_BS_empirical_CI_lower_original = PCA_BS_empirical_CI_lower
	PCA_BS_empirical_CI_upper_original = PCA_BS_empirical_CI_upper

	# ANGLE CALCULATION
	BS_PCA_angle_sd = std_of_angles_from_mean_back_to_half
	BS_PCA_mean_angle = BS_PCA_angle_back_to_180_mean
	BS_PCA_mean_angle_up_to_45 = BS_PCA_mean_angle

	if (BS_PCA_mean_angle_up_to_45 > 45.0):
		BS_PCA_mean_angle_up_to_45 = 90.0 - BS_PCA_mean_angle_up_to_45

		temp = PCA_BS_empirical_CI_lower
		PCA_BS_empirical_CI_lower = PCA_BS_empirical_CI_upper
		PCA_BS_empirical_CI_upper = temp


	BS_PCA_mean_angle_difference_from_45 = 45.0 - BS_PCA_mean_angle
	
	
	# Empirical p-value calculation
	number_less_than_45 = np.sum(np.asarray(angle_differences_to_45_which_is_90_expanded) < 0.0)
	pvalue_1 = float(number_less_than_45) / float(number_of_bootstrap_iterations)
	number_greater_than_45 = np.sum(np.asarray(angle_differences_to_45_which_is_90_expanded) > 0.0)
	pvalue_2 = float(number_greater_than_45) / float(number_of_bootstrap_iterations)
	empirical_pvalue_different_from_45 = 2.0*np.min([pvalue_1, pvalue_2])


	# Empirical p-value calculation
	number_less_than_0 = np.sum(np.asarray(angle_differences_to_0_which_is_0_expanded) < 0.0)
	pvalue_1 = float(number_less_than_0) / float(number_of_bootstrap_iterations)
	number_greater_than_0 = np.sum(np.asarray(angle_differences_to_0_which_is_0_expanded) > 0.0)
	pvalue_2 = float(number_greater_than_0) / float(number_of_bootstrap_iterations)
	empirical_pvalue_different_from_0 = 2.0*np.min([pvalue_1, pvalue_2])


	# P-VALUES DIFFERENT FROM 45
	new_temp = stats.norm(np.mean(angle_differences_to_45_which_is_90_expanded), np.std(angle_differences_to_45_which_is_90_expanded)).sf(0.0)
	if (new_temp > 0.5):
		new_temp = 1.0 - new_temp
	BS_PCA_different_from_45_sd_method = 2.0*new_temp


	local_data_dict = {}
	local_data_dict['BS_PCA_mean_angle_up_to_45'] = BS_PCA_mean_angle_up_to_45
	local_data_dict['PCA_BS_empirical_CI_lower'] = PCA_BS_empirical_CI_lower
	local_data_dict['PCA_BS_empirical_CI_upper'] = PCA_BS_empirical_CI_upper
	local_data_dict['PCA_BS_empirical_CI_lower_original'] = PCA_BS_empirical_CI_lower_original
	local_data_dict['PCA_BS_empirical_CI_upper_original'] = PCA_BS_empirical_CI_upper_original
	local_data_dict['PCA_BS_empirical_pvalue_different_from_45'] = empirical_pvalue_different_from_45
	local_data_dict['is_PCA_BS_empirical_pvalue_different_from_45'] = empirical_pvalue_different_from_45 <= 0.025
	local_data_dict['PCA_BS_empirical_pvalue_different_from_0'] = empirical_pvalue_different_from_0
	local_data_dict['is_PCA_BS_empirical_pvalue_different_from_0'] = empirical_pvalue_different_from_0 <= 0.025
	local_data_dict['BS_PCA_different_from_45_sd_method'] = BS_PCA_different_from_45_sd_method
	local_data_dict['is_BS_PCA_different_from_45_sd_method'] = BS_PCA_different_from_45_sd_method < 0.05
	local_data_dict['BS_PCA_mean_component_0_sd'] = BS_PCA_mean_component_0_sd
	local_data_dict['BS_PCA_mean_component_1_sd'] = BS_PCA_mean_component_1_sd
	local_data_dict['BS_PCA_mean_of_mean0'] = BS_PCA_mean_of_mean0
	local_data_dict['BS_PCA_mean_of_mean1'] = BS_PCA_mean_of_mean1
	local_data_dict['BS_PCA_mean_sd_area'] = BS_PCA_mean_sd_area
	local_data_dict['BS_PCA_mean_covariance_matrix'] = BS_PCA_mean_covariance_matrix
	local_data_dict['BS_PCA_angle_sd'] = BS_PCA_angle_sd
	local_data_dict['BS_PCA_mean_angle'] = BS_PCA_mean_angle
	# local_data_dict['BS_PCA_mean_angle_right_half'] = BS_PCA_mean_angle_right_half
	local_data_dict['BS_PCA_mean_angle_difference_from_45'] = BS_PCA_mean_angle_difference_from_45
	local_data_dict['BS_PCA_different_from_45_sd_method'] = BS_PCA_different_from_45_sd_method
	local_data_dict['BS_PCA_different_from_45_generalised_method'] = -1.0
	local_data_dict['BS_PCA_ellipses'] = []
	local_data_dict['BS_PCA_ellipses_mean_zeroed'] = []
	local_data_dict['BS_PCA_ellipses_ORIGINAL_ELLIPSE_AT_NEW_MEAN'] = []


	if (single_cluster_holder_options != None):
		for ring in single_cluster_holder_options.rings:
			pca_ellr = create_ellipse([BS_PCA_mean_of_mean0, BS_PCA_mean_of_mean1], BS_PCA_mean_component_0_sd*ring, BS_PCA_mean_component_1_sd*ring, BS_PCA_mean_angle)
			local_data_dict['BS_PCA_ellipses'].append(pca_ellr)

			pca_ellr_mean_zeroed = create_ellipse([0, 0], BS_PCA_mean_component_0_sd*ring, BS_PCA_mean_component_1_sd*ring, BS_PCA_mean_angle)
			local_data_dict['BS_PCA_ellipses_mean_zeroed'].append(pca_ellr_mean_zeroed)

		local_data_dict['PCA_BS_intersection_ellipse'] = create_ellipse([BS_PCA_mean_of_mean0, BS_PCA_mean_of_mean1], BS_PCA_mean_component_0_sd*single_cluster_holder_options.intersection_dist, BS_PCA_mean_component_1_sd*single_cluster_holder_options.intersection_dist, BS_PCA_mean_angle)


	components_ = np.asarray([1.0, math.tan(math.radians(BS_PCA_mean_angle))])

	# print(BS_PCA_mean_angle, components_)

	means_ = np.asarray([local_data_dict['BS_PCA_mean_of_mean0'], local_data_dict['BS_PCA_mean_of_mean1']])
	local_data_dict['PCA_predicited_state_for_cluster_trials'] = np.dot(spikes_to_use - means_, components_)

	factor_line_m = components_[1] / components_[0]
	factor_line_c = means_[1] - factor_line_m * means_[0]
	perpendicular_m = -1.0 / factor_line_m
	perpendicular_c = spikes_to_use[:,1] - perpendicular_m * spikes_to_use[:,0]
	intersect_x =  (perpendicular_c - factor_line_c) / (factor_line_m - perpendicular_m)
	intersect_y = factor_line_m * intersect_x + factor_line_c

	local_data_dict['distances_from_factor_line'] = np.sqrt((intersect_x - spikes_to_use[:,0])**2 + (intersect_y - spikes_to_use[:,1])**2) / BS_PCA_mean_component_1_sd


	return local_data_dict
	

"""
ELLIPSE HELPER
"""
from scipy.spatial.distance import cdist

def indices_within_x_mahalanobis_distances(mean, covariance_matrix, spike_pair_distribution_spikes, distance_criterion):
	mahalanobis_distances = cdist(spike_pair_distribution_spikes, [mean], metric='mahalanobis', VI=np.linalg.inv(covariance_matrix))
	indices_within_bound = np.argwhere(mahalanobis_distances <= distance_criterion).T[0]
	return indices_within_bound



"""
AUTOCORRELATION HELPERS
Helper functions for calculating autocorrelations and partial autocorrelation p values
"""

from scipy.stats import linregress
from scipy import stats

def calculate_acf_pvalues_for_spikes(spikes, spikes2, number_of_lags):

	acf_rvalues_0 = []
	acf_pvals_0 = []
	acf_positive_pvals_0 = []
	acf_negative_pvals_0 = []

	for shift in range(1,number_of_lags + 1):	

		lr_0 = linregress(spikes[:-shift], spikes2[shift:]);

		acf_rvalues_0.append(lr_0.rvalue)
		acf_pvals_0.append(lr_0.pvalue)

		df = len(spikes[:-shift]) - 1

		t_statistic = (lr_0.slope - 0.0) / lr_0.stderr
		p_greater_than_0 = 1 - stats.t.cdf(t_statistic,df=df)
		acf_positive_pvals_0.append(p_greater_than_0)

		t_statistic = (lr_0.slope) / lr_0.stderr
		p_less_than_0 = stats.t.cdf(t_statistic,df=df)
		acf_negative_pvals_0.append(p_less_than_0)

	return acf_rvalues_0, acf_pvals_0, acf_positive_pvals_0, acf_negative_pvals_0


def calculate_pacf_pvalues_for_spikes(spikes, residuals, number_of_lags):

	pacf_rvalues = []
	pacf_pvals = []
	pacf_positive_pvals = []
	pacf_negative_pvals = []

	for shift in range(1,number_of_lags + 1):			

		lr = linregress(spikes[:-shift], residuals[1:]); slope = lr.slope; intercept = lr.intercept
		
		pacf_rvalues.append(lr.rvalue)
		pacf_pvals.append(lr.pvalue)
		
		df = len(spikes[:-shift]) - 1

		t_statistic = (lr.slope - 0.0) / lr.stderr
		p_greater_than = 1 - stats.t.cdf(t_statistic,df=df)
		pacf_positive_pvals.append(p_greater_than)

		t_statistic = (lr.slope) / lr.stderr
		p_less_than = stats.t.cdf(t_statistic,df=df)
		pacf_negative_pvals.append(p_less_than)


		estimate = intercept + slope * spikes[:-shift]
		residuals = residuals[1:] - estimate

	return pacf_rvalues, pacf_pvals, pacf_positive_pvals, pacf_negative_pvals



"""
DIFFERENCING & ARIMA HELPERS
Helper functions for making differencing decision & calculating ARIMA orders
"""

def calculate_arima_variables_for_acf_and_pacf_postive_arrays(acf_rvalues, acf_pvalues, pacf_pvalues, acf_positive_pvalues, pacf_positive_pvalues):

	AR_p = 0
	MA_q = 0
	use_arima_model = False

	# From https://people.duke.edu/~rnau/411arim2.htm (Rest of the differencing rules below)
	# Rule 6: If the PACF of the differenced series displays a sharp cutoff and/or the lag-1 autocorrelation is positive--i.e., if the series appears slightly "underdifferenced"--then consider adding an AR term to the model. The lag at which the PACF cuts off is the indicated number of AR terms.
	# Rule 7: If the ACF of the differenced series displays a sharp cutoff and/or the lag-1 autocorrelation is negative--i.e., if the series appears slightly "overdifferenced"--then consider adding an MA term to the model. The lag at which the ACF cuts off is the indicated number of MA terms.


	if (acf_rvalues[0] < -0.1): # If negative first autocorrelation
		
		last_pval = 0.0
		for acf_pval_index, acf_pval in enumerate(acf_pvalues[:4]):
			if (acf_pval < 0.05):
				last_pval = acf_pval
			if (acf_pval > 0.05):
				difference_to_last = acf_pval - last_pval
				if (difference_to_last > 0.15):
					MA_q = acf_pval_index
					use_arima_model = True
					break

	if (acf_rvalues[0] > 0.1): # If positive first autocorrelation
		
		last_pval = 0.0
		for pacf_pval_index, pacf_pval in enumerate(pacf_pvalues[:4]):
			if (pacf_pval < 0.05):
				last_pval = pacf_pval
			if (pacf_pval > 0.05):
				difference_to_last = pacf_pval - last_pval
				if (difference_to_last > 0.15):
					AR_p = pacf_pval_index
					use_arima_model = True
					break

	return AR_p, MA_q, use_arima_model



def decisions_for_selective_differencing(TI_Vs_SpikeTimes_LR_pvalue, KPSS_SpikeTimes_pvalue, AutoCorr1_SpikeTimes_pvalue, TI_Vs_PairsDifferenced_LR_pvalue, KPSS_PairsDifferenced_pvalue):

	# https://people.duke.edu/~rnau/411arim2.htm
	# Rule 1: If the series has positive autocorrelations out to a high number of lags, then it probably needs a higher order of differencing.
	# Rule 2: If the lag-1 autocorrelation is zero or negative, or the autocorrelations are all small and patternless, then the series does not need a higher order of  differencing. If the lag-1 autocorrelation is -0.5 or more negative, the series may be overdifferenced.  BEWARE OF OVERDIFFERENCING!!
	# Rule 3: The optimal order of differencing is often the order of differencing at which the standard deviation is lowest.
	# Rule 4: A model with no orders of differencing assumes that the original series is stationary (mean-reverting). A model with one order of differencing assumes that the original series has a constant average trend (e.g. a random walk or SES-type model, with or without growth). A model with two orders of total differencing assumes that the original series has a time-varying trend (e.g. a random trend or LES-type model).
	# Rule 5: A model with no orders of differencing normally includes a constant term (which allows for a non-zero mean value). A model with two orders of total differencing normally does not include a constant term. In a model with one order of total differencing, a constant term should be included if the series has a non-zero average trend.
	# Rule 6: If the PACF of the differenced series displays a sharp cutoff and/or the lag-1 autocorrelation is positive--i.e., if the series appears slightly "underdifferenced"--then consider adding an AR term to the model. The lag at which the PACF cuts off is the indicated number of AR terms.
	# Rule 7: If the ACF of the differenced series displays a sharp cutoff and/or the lag-1 autocorrelation is negative--i.e., if the series appears slightly "overdifferenced"--then consider adding an MA term to the model. The lag at which the ACF cuts off is the indicated number of MA terms.

	UNCORR_ORIG = False
	UNCORR_DIFF = False
	STAT_ORIG = False
	STAT_DIFF = False
	AUTO_ORIG = False
	# AUTO_DIFF = False
	


	threshold_1 = 0.05
	if (TI_Vs_SpikeTimes_LR_pvalue >= threshold_1):
		UNCORR_ORIG = True
	if (TI_Vs_PairsDifferenced_LR_pvalue >= threshold_1):
		UNCORR_DIFF = True

	threshold_2 = 0.05
	if (KPSS_SpikeTimes_pvalue >= threshold_2):
		STAT_ORIG = True  
	if (KPSS_PairsDifferenced_pvalue >= threshold_2):
		STAT_DIFF = True

	threshold_3 = 0.05
	if (AutoCorr1_SpikeTimes_pvalue >= threshold_3):
		AUTO_ORIG = True  
	# if (AutoCorr1_PairsDifferenced_pvalue >= threshold_2):
	# 	AUTO_DIFF = True

	index_to_use = -1

	original_needs_differencing = False
	if ((not UNCORR_ORIG) | (not STAT_ORIG) | (not AUTO_ORIG)):
		original_needs_differencing = True
	else:
		index_to_use = 0

	
	if (original_needs_differencing):
		differencining_fixed_original = False
		if (UNCORR_DIFF & STAT_DIFF):
			differencining_fixed_original = True
			index_to_use = 1


	return index_to_use




"""
BOOTSTRAP FACTOR ANALYSIS
Helper functions for bootstrap factor analysis
"""

from sklearn.decomposition import FactorAnalysis
import numpy as np
import math
from scipy import linalg

def bootstrap_factor_analysis(cluster, number_of_bootstrap_iterations, analysis_dict_key):

	if (analysis_dict_key == "Original"):
		spikes_to_use = cluster.cluster_spike_pairs
		number_of_samples_per_bootstrap = cluster.number_of_samples_in_ellipse

	elif (analysis_dict_key == "SelectivelyDifferenced"):
		if (cluster.use_selective_differences):
			spikes_to_use = cluster.selective_differences_scaled_down
			number_of_samples_per_bootstrap = cluster.number_of_samples_for_selective_differences
		else:
			return

	elif (analysis_dict_key == "SelectivelyDifferencedBoxJenkins"):
		if (cluster.use_selective_differences):
			spikes_to_use = cluster.selective_differences_arima_residuals
			number_of_samples_per_bootstrap = cluster.number_of_samples_for_selective_differences
		else:
			return

	angle_half_expanded_to_360_x_ys = []
	bootstrapped_angles = []
	bootstrapped_FA_N0_sds = []
	bootstrapped_FA_N1_sds = []
	bootstrapped_FA_1sd_areas = []
	bootstrapped_FA_1sd_estimated_difference_area = []
	bootstrapped_FA_mean_0s = []
	bootstrapped_FA_mean_1s = []
	bootstrapped_FA_components_0s = []
	bootstrapped_FA_components_1s = []
	bootstrapped_FA_noise_variances_0 = []
	bootstrapped_FA_noise_variances_1 = []

	bootstrapped_FA_predicited_states_for_cluster_trials = []

	for bootstrap_iteration in range(number_of_bootstrap_iterations):
		bootstrap_sample_indices = np.random.choice(number_of_samples_per_bootstrap, number_of_samples_per_bootstrap)
		bootstrap_sample = spikes_to_use[bootstrap_sample_indices]

		factor_analysis = FactorAnalysis(n_components=1, random_state=bootstrap_iteration)
		factor_analysis.fit(bootstrap_sample)
		

		FA_N0_sd = np.sqrt(factor_analysis.noise_variance_[0])
		FA_N1_sd = np.sqrt(factor_analysis.noise_variance_[1])
		FA_1sd_area = FA_N0_sd * FA_N1_sd
		FA_1sd_estimated_difference_area = np.sqrt(factor_analysis.noise_variance_[0] + factor_analysis.noise_variance_[1])

		factor_analysis_angle = (np.arctan2([factor_analysis.components_[0][1]], [factor_analysis.components_[0][0]]) * 180 / np.pi)[0]

		if (factor_analysis_angle < 0.0):
			factor_analysis_angle  = factor_analysis_angle + 180.0
			bootstrapped_FA_components_0s.append(-factor_analysis.components_[0][0])
			bootstrapped_FA_components_1s.append(-factor_analysis.components_[0][1])
		else:
			bootstrapped_FA_components_0s.append(factor_analysis.components_[0][0])
			bootstrapped_FA_components_1s.append(factor_analysis.components_[0][1])


		angle_half_expanded_to_360 = factor_analysis_angle * 2.0
		angle_half_expanded_to_360_radians = math.radians(angle_half_expanded_to_360)
		angle_half_expanded_to_360_x_y = [math.cos(angle_half_expanded_to_360_radians), math.sin(angle_half_expanded_to_360_radians)]
		angle_half_expanded_to_360_x_ys.append(angle_half_expanded_to_360_x_y)


		bootstrapped_FA_predicited_states_for_cluster_trials.append(factor_analysis.transform(bootstrap_sample))
		bootstrapped_angles.append(factor_analysis_angle)
		bootstrapped_FA_N0_sds.append(FA_N0_sd)
		bootstrapped_FA_N1_sds.append(FA_N1_sd)
		bootstrapped_FA_1sd_areas.append(FA_1sd_area)
		bootstrapped_FA_1sd_estimated_difference_area.append(FA_1sd_estimated_difference_area)
		bootstrapped_FA_mean_0s.append(factor_analysis.mean_[0])
		bootstrapped_FA_mean_1s.append(factor_analysis.mean_[1])
		bootstrapped_FA_noise_variances_0.append(factor_analysis.noise_variance_[0])
		bootstrapped_FA_noise_variances_1.append(factor_analysis.noise_variance_[1])

	angle_half_expanded_to_360_x_y_mean = np.mean(angle_half_expanded_to_360_x_ys, axis=0)
	angle_back_to_180_mean = sw.convert_expanded_360_xy_to_180_angle_degrees([angle_half_expanded_to_360_x_y_mean[0]], [angle_half_expanded_to_360_x_y_mean[1]])[0]

	FA_angle_BS_mean_45d = angle_back_to_180_mean
	if (FA_angle_BS_mean_45d > 45.0):
		FA_angle_BS_mean_45d = 90.0 - FA_angle_BS_mean_45d
	

	array_of_extra_analysis_keys = [analysis_dict_key]
	if (analysis_dict_key in ["Original", "SelectivelyDifferenced", "SelectivelyDifferencedBoxJenkins"]):

		sublist_suffixes = ["TestsPassed", "TestsPassedAndNormal"]
		for sublist_suffix in sublist_suffixes:

			if (not cluster.analysis_dicts[analysis_dict_key + sublist_suffix]['is_empty']):
				array_of_extra_analysis_keys.append(analysis_dict_key + sublist_suffix)

		if (analysis_dict_key == "SelectivelyDifferenced"):
			if (not cluster.analysis_dicts[analysis_dict_key + "TestsPassedActuallyDifferenced"]['is_empty']):
				array_of_extra_analysis_keys.append(analysis_dict_key + "TestsPassedActuallyDifferenced")


	for extra_analysis_dict_key in array_of_extra_analysis_keys:

		bootstrapped_FA_predicited_states_for_cluster_trials = np.asarray(bootstrapped_FA_predicited_states_for_cluster_trials)

		cluster.analysis_dicts[extra_analysis_dict_key]['FA_predicited_state_for_cluster_trials'] = np.mean(bootstrapped_FA_predicited_states_for_cluster_trials, axis=0)

		cluster.analysis_dicts[extra_analysis_dict_key]['FA_angle_BS_mean'] = angle_back_to_180_mean
		cluster.analysis_dicts[extra_analysis_dict_key]['FA_N0_sd_BS_mean'] = np.mean(bootstrapped_FA_N0_sds)
		cluster.analysis_dicts[extra_analysis_dict_key]['FA_N1_sd_BS_mean'] = np.mean(bootstrapped_FA_N1_sds)
		cluster.analysis_dicts[extra_analysis_dict_key]['FA_1sd_area_BS_mean'] = np.mean(bootstrapped_FA_1sd_areas)
		cluster.analysis_dicts[extra_analysis_dict_key]['FA_1sd_estimated_difference_area_BS_mean'] = np.mean(bootstrapped_FA_1sd_estimated_difference_area)
		cluster.analysis_dicts[extra_analysis_dict_key]['FA_mean_0_mean'] = np.mean(bootstrapped_FA_mean_0s)
		cluster.analysis_dicts[extra_analysis_dict_key]['FA_mean_1_mean'] = np.mean(bootstrapped_FA_mean_1s)
		cluster.analysis_dicts[extra_analysis_dict_key]['FA_component_0_mean'] = np.mean(bootstrapped_FA_components_0s)
		cluster.analysis_dicts[extra_analysis_dict_key]['FA_component_1_mean'] = np.mean(bootstrapped_FA_components_1s)
		cluster.analysis_dicts[extra_analysis_dict_key]['FA_noise_variance_0_mean'] = np.mean(bootstrapped_FA_noise_variances_0)
		cluster.analysis_dicts[extra_analysis_dict_key]['FA_noise_variance_1_mean'] = np.mean(bootstrapped_FA_noise_variances_1)
		cluster.analysis_dicts[extra_analysis_dict_key]['FA_angle_BS_mean_45d'] = FA_angle_BS_mean_45d
		cluster.analysis_dicts[extra_analysis_dict_key]['FA_angle_BS_empirical_CI'] = sw.empirical_confidence_bound(bootstrapped_angles)


		mean_ = np.asarray([cluster.analysis_dicts[extra_analysis_dict_key]['FA_mean_0_mean'], cluster.analysis_dicts[extra_analysis_dict_key]['FA_mean_1_mean']])
		components_ = np.asarray([[cluster.analysis_dicts[extra_analysis_dict_key]['FA_component_0_mean'], cluster.analysis_dicts[extra_analysis_dict_key]['FA_component_1_mean']]])
		noise_variance_ = np.asarray([cluster.analysis_dicts[extra_analysis_dict_key]['FA_noise_variance_0_mean'], cluster.analysis_dicts[extra_analysis_dict_key]['FA_noise_variance_1_mean']])

		cluster.analysis_dicts[extra_analysis_dict_key]['FA_predicited_state_for_cluster_trials'] = mean_fa_transform(spikes_to_use, mean_, components_, noise_variance_)


def mean_fa_transform(X, mean_, components_, noise_variance_):
	# JI: Adapted from FactorAnalysis 'transform' function to take arguments
	"""Apply dimensionality reduction to X using the model.
	Compute the expected mean of the latent variables.
	See Barber, 21.2.33 (or Bishop, 12.66).
	Parameters
	"""

	Ih = np.eye(len(components_))

	X_transformed = X - mean_

	Wpsi = components_ / noise_variance_
	cov_z = linalg.inv(Ih + np.dot(Wpsi, components_.T))
	tmp = np.dot(X_transformed, Wpsi.T)
	X_transformed = np.dot(tmp, cov_z)

	return X_transformed



"""
CHI SQUARED (FISHER STATISTIC)
Helper function for chi squared (Fisher) statistic
"""

from scipy.stats import combine_pvalues

def calculate_chi_squared_string_for_values(values_for_histogram_in_range):

	number_of_values = len(values_for_histogram_in_range)
	with warnings.catch_warnings():
		warnings.filterwarnings("ignore")
		chi_sqaured_statistic, chi_sqaured_p_value = combine_pvalues(values_for_histogram_in_range)

	chi_squared_latex_string = ""
	part_of_table_string = ""
	if (not math.isnan(chi_sqaured_p_value)):
		chi_squared_latex_string = "$\\chi^{2}(" + str(2 * number_of_values) + ", N = " + str(number_of_values) + ") = " + str(round(chi_sqaured_statistic, 1)) + ", $ $" + str(sw.string_with_pval_lessthan_to_lower_limit(chi_sqaured_p_value)) + "$"
		part_of_table_string = "df=" + str(2 * number_of_values) + ", N=" + str(number_of_values) + ", " + str(round(chi_sqaured_statistic, 1)) + ", " + str(chi_sqaured_p_value)

	return number_of_values, chi_sqaured_statistic, chi_sqaured_p_value, chi_squared_latex_string, part_of_table_string


"""
Matrix Helpers
"""

def is_pos_def(x):
	return np.all(np.linalg.eigvals(x) > 0)


