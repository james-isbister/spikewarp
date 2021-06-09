from sklearn.linear_model import LinearRegression
from scipy.stats import linregress
import numpy as np

import spikewarp as sw

def times_state_prediction(states, times, num_states_averaged=-1, sch=None, FA_predicted_sd=-1.0, conj_orig_sd=-1.0, plot_dewarp_and_state_dependence=False, plot_path='', only_plot_if_corr=True):

	# Helper function for linear prediction of one quantity from another.

	LR = LinearRegression().fit(states.reshape(-1, 1), times.reshape(-1, 1))
	lr = linregress(states, times)
	linear_state_prediction_errors = times.reshape(-1, 1) - LR.predict(states.reshape(-1, 1))

	d = {'orig_sd': np.std(times),
		'state_predic_error_sd': np.std(linear_state_prediction_errors),
		'r': lr.rvalue,
		'rsquared': lr.rvalue**2,
		'p': lr.pvalue,
		'FA_predicted_sd': FA_predicted_sd,
		'conj_orig_sd': conj_orig_sd,
		'num_states_averaged':num_states_averaged}


	if hasattr(sw, "placeholder_6"):
		if (plot_dewarp_and_state_dependence):
			if ((not only_plot_if_corr) | (only_plot_if_corr & (d['p'] < 0.05))):
				sw.placeholder_6(d, times, num_states_averaged, FA_predicted_sd, sch, states, linear_state_prediction_errors)

	return d

