import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

import numpy as np
import cv2
import sys
import math
from scipy.ndimage import gaussian_filter
from pylab import *
from matplotlib.ticker import MultipleLocator
from matplotlib import rc
import pandas as pd

import spikewarp as sw


"""
VIDEO HELPERS
Helper functions for creating videos
"""

def create_gif_from_list_of_image_file_names(list_of_files, 
											output_file_name, 
											frames_per_second=1):

	if (list_of_files != []):

		img_array = []
		for filename in list_of_files:
			img = cv2.imread(filename)
			height, width, layers = img.shape
			size = (width,height)
			img_array.append(img)

		out = cv2.VideoWriter(output_file_name, cv2.VideoWriter_fourcc(*'mp4v'), frames_per_second, size)
		for i in range(len(img_array)):
			frame = cv2.resize(img_array[i], size)
			out.write(frame)
		out.release()



def add_appropriate_pairwise_response_dist_files(analysis, 
												tw_files_by_stwcboo_and_shuffle_option):

	"""
	Function adds the appropriate pairwise response distribution files to list for later gif creation
	"""

	if (analysis['do_plot_super_time_warp_cluster_busters'] & analysis['do_plot_super_time_warp_cluster_buster_videos_ALL_PAIRS']):
		for shuffle_option_string in analysis['option_helper'].shuffle_option_strings:
			for stwcboo_index, stwcboo in enumerate(analysis['super_time_warp_cluster_buster_option_objects']):
				if stwcboo['create_video']:
					analysis['all_tw_files_by_shuffle_type_and_stwcboo'][shuffle_option_string][stwcboo_index].extend(tw_files_by_stwcboo_and_shuffle_option[shuffle_option_string][stwcboo_index])



def create_pairwise_response_dist_gifs_when_appropriate(analysis):

	"""
	Function creates response distribution gifs when appropriate
	"""

	for shuffle_option_string in analysis['option_helper'].shuffle_option_strings:
		if (analysis['do_plot_super_time_warp_cluster_busters'] & analysis['do_plot_super_time_warp_cluster_buster_videos_ALL_PAIRS']):
			for stwcboo_index, stwcboo in enumerate(analysis['super_time_warp_cluster_buster_option_objects']):
				if (analysis['all_tw_files_by_shuffle_type_and_stwcboo'][shuffle_option_string][stwcboo_index] != []):
					if (stwcboo['create_video']):
						sw.create_gif_from_list_of_image_file_names(analysis['all_tw_files_by_shuffle_type_and_stwcboo'][shuffle_option_string][stwcboo_index], 
																	analysis['directory_holder'].pair_cluster_plots_directory + shuffle_option_string + "/" + stwcboo['type_string'] + '_all_super_time_warp_cluster_buster_files.mp4')



"""
GENERAL PLOTTING FUNCTIONS
Collection of multi purpose general plotting functions
"""

def draw_neighbouring_bar_chart(data, 
								x_labels, 
								dir_and_filename, 
								title, 
								bar_legend_labels, 
								xaxis_label, 
								custom_y_tick_locators=[-1, -1], 
								rotate_x_labels=False, 
								y_tick_right=False, 
								optional_y_axis_label='',
								positive_confidence_intervals=[],
								threshold_lin_value=0.0,
								width=0.35):

	"""
	Custom general function for plotting multiple neighbouring bar charts
	"""

	N_per_xpos = len(data)
	N_xpos = len(data[0])
	ind = np.arange(N_xpos)  # the x locations for the groups
	# width = 0.35       # the width of the bars

	fig, ax = plt.subplots()

	rects = []
	colours = ['b', 'g', 'r']

	rect_tupple_for_legend = []
	for single_data_index, single_data in enumerate(data):

		confidence_intervals=None
		if (positive_confidence_intervals != []):
			confidence_intervals = positive_confidence_intervals[single_data_index]
		
		rects1 = ax.bar(ind + single_data_index * width, single_data, width, color=colours[single_data_index], yerr=confidence_intervals)
		rects.append(rects1)
		rect_tupple_for_legend = rect_tupple_for_legend + [rects1[0]]


	if (y_tick_right):
		ax.yaxis.tick_right()

	# add some text for labels, title and axes ticks
	ax.set_xlabel(xaxis_label, fontsize=24)
	y_axis_label = 'Count'
	if (optional_y_axis_label != ""):
		y_axis_label = optional_y_axis_label
	ax.set_ylabel(y_axis_label, fontsize=24)
	ax.set_title(title)
	ax.set_xticks(ind + width / 2)
	ax.set_xticklabels(x_labels)
	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(24) 
	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(24) 

	if (custom_y_tick_locators[0] != -1):
		plt.gca().yaxis.set_major_locator(MultipleLocator(custom_y_tick_locators[0]))	
		plt.gca().yaxis.set_minor_locator(MultipleLocator(custom_y_tick_locators[1]))
		plt.gca().set_ylim([0, custom_y_tick_locators[0]])

	if (rotate_x_labels):
		plt.xticks(rotation=90)

	if (threshold_lin_value != 0.0):
		plt.axhline(y=threshold_lin_value, xmin=0.0, xmax=ind[-1], c='k', linestyle='--')


	ax.legend(rect_tupple_for_legend, bar_legend_labels, prop={'size': 20})
	plt.savefig(dir_and_filename, bbox_inches='tight')
	plt.close()


def normal_histo_plot(values_for_histogram, 
						file_name, 
						colors=['b', 'g', 'r'], 
						bins=200, 
						histo_range=[0.0, 1.0], 
						x_axis_left_buffer=0.0, 
						x_axis_label='', 
						y_axis_label='', 
						custom_x_tick_locators=[-1, -1], 
						custom_y_tick_locators=[-1, -1], 
						alpha=1.0, 
						add_chi_squared_text=False,
						density=False,
						labels=[]):

	"""
	Custom general function for plotting histogram
	"""

	sns.set()
	sns.set_style("ticks")
	sns.despine()

	num_hists = len(values_for_histogram)
	rwidth=0.99
	if (num_hists > 1):
		rwidth = 0.75

		
	plt.figure()
	plt.hist(values_for_histogram, bins=bins, range=histo_range, color=colors[:num_hists], rwidth=rwidth, alpha=alpha, density=density, label=labels)
	# for index, hist_values in enumerate(values_for_histogram):
	# 	if (labels != []):
			# plt.hist(hist_values, bins=bins, range=histo_range, color=colors[index], rwidth=0.99, alpha=alpha, density=density, label=labels[index])
		# else:
			# plt.hist(hist_values, bins=bins, range=histo_range, color=colors[index], rwidth=0.99, alpha=alpha, density=density)

	if (x_axis_label != ''):
		plt.gca().set_xlabel(x_axis_label, fontsize=24)
	if (y_axis_label != ''):
		plt.gca().set_ylabel(y_axis_label, fontsize=24)


	for tick in plt.gca().xaxis.get_major_ticks():
		tick.label.set_fontsize(24) 
	for tick in plt.gca().yaxis.get_major_ticks():
		tick.label.set_fontsize(24) 

	plt.gca().set_xlim([histo_range[0] - x_axis_left_buffer, histo_range[1]])

	if (custom_x_tick_locators[0] != -1):
		plt.gca().xaxis.set_major_locator(MultipleLocator(custom_x_tick_locators[0]))
		plt.gca().xaxis.set_minor_locator(MultipleLocator(custom_x_tick_locators[1]))

	if (not density):
		if (custom_y_tick_locators[0] != -1):
			plt.gca().yaxis.set_major_locator(MultipleLocator(custom_y_tick_locators[0]))	
			plt.gca().yaxis.set_minor_locator(MultipleLocator(custom_y_tick_locators[1]))

	
	if (add_chi_squared_text):
		chi_squared_table_strings_array = []
		chi_squared_latex_strings = ''
		chi_squared_table_strings = ''
		
		values_for_histogram = np.asarray(values_for_histogram)
		values_for_histogram_in_range = values_for_histogram[np.where(np.logical_and(values_for_histogram > 0.0, values_for_histogram < 1.0))]
		number_of_values, chi_sqaured_statistic, chi_sqaured_p_value, chi_squared_latex_string, part_of_table_string = sw.calculate_chi_squared_string_for_values(values_for_histogram_in_range)
		
		plt.text(0.5, 0.5, part_of_table_string, fontsize=9, horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)

		chi_squared_table_strings_array.append(chi_squared_latex_string)


	if (labels != []):
		plt.gca().legend()

	plt.savefig(file_name + '.pdf', bbox_inches='tight')
	plt.close()

	if (add_chi_squared_text):
		return chi_squared_table_strings_array
	else:
		return




def cumulative_histo_plot(list_of_values_for_histogram, 
						file_name, 
						bins=200, 
						histo_range=[0.0, 1.0], 
						x_axis_label='', 
						y_axis_label='', 
						custom_x_tick_locators=[-1, -1], 
						xaxis_lims=[-1.0, -1.0], 
						add_chi_squared_text=False):

	"""
	Custom general function for plotting multiple cumulative histograms
	"""

	sns.set()
	sns.set_style("ticks")
	plt.figure()

	chi_squared_table_strings_array = []
	if (len(list_of_values_for_histogram) > 0):

		x = np.linspace(0,1,5)
		y = x
		plt.plot(x, y, '--k', linewidth=10)

		colors = ['b', 'g', 'r', 'c', 'orange', 'k', 'm', 'y', 'slateblue', 'hotpink', 'b', 'g']

		chi_squared_latex_strings = ''
		chi_squared_table_strings = ''
		

		one_to_draw = False
		for index, values_for_histogram in enumerate(list_of_values_for_histogram):
			values_for_histogram = np.asarray(values_for_histogram)
			values_for_histogram = values_for_histogram[values_for_histogram>=0.0]
			
			if (add_chi_squared_text):
				number_of_values, chi_sqaured_statistic, chi_sqaured_p_value, chi_squared_latex_string, part_of_table_string = sw.calculate_chi_squared_string_for_values(values_for_histogram)
				chi_squared_table_strings += str(index + 1) + ", " + colors[index] + ", " + part_of_table_string + "\n"
				chi_squared_table_strings_array.append(chi_squared_latex_string)

			if (values_for_histogram.shape[0] > 0):

				hist, bin_edges = np.histogram(values_for_histogram, bins=bins, range=histo_range)
				if (np.sum(hist) > 0):
					one_to_draw = True
					plt.step(bin_edges[:-1], np.cumsum(hist).astype(float) / np.sum(hist).astype(float), linewidth=3, color=colors[index])
			
		if (add_chi_squared_text):
			plt.text(0.1, 0.1, chi_squared_table_strings, fontsize=9)
			with open(file_name + "_TABLE.txt", "w") as tex_file:
				# print(file_name)
				print(chi_squared_table_strings, file=tex_file)

		if (not one_to_draw):
			plt.close()
			return


		if (x_axis_label != ''):
			plt.gca().set_xlabel(x_axis_label, fontsize=24)
		if (y_axis_label != ''):
			plt.gca().set_ylabel(y_axis_label, fontsize=24)

		for tick in plt.gca().xaxis.get_major_ticks():
			tick.label.set_fontsize(24) 
		for tick in plt.gca().yaxis.get_major_ticks():
			tick.label.set_fontsize(24) 

		if (custom_x_tick_locators[0] != -1):
			plt.gca().xaxis.set_major_locator(MultipleLocator(custom_x_tick_locators[0]))
			plt.gca().yaxis.set_major_locator(MultipleLocator(custom_x_tick_locators[0]))

			plt.gca().xaxis.set_minor_locator(MultipleLocator(custom_x_tick_locators[1]))
			plt.gca().yaxis.set_minor_locator(MultipleLocator(custom_x_tick_locators[1]))


		if (xaxis_lims[0] == -1.0):
			plt.gca().set_xlim(histo_range)
		else:
			plt.gca().set_xlim(xaxis_lims)
		plt.gca().set_ylim([0.0, 1.0]) 

		# sns.set_style("ticks")
		sns.despine()
			
		plt.savefig(file_name + '.pdf', bbox_inches='tight')
	plt.close()

	return chi_squared_table_strings_array



def basic_x_y_plot(x_point_groups, 
					y_point_groups, 
					file_name, 
					draw_y_equals_x=False, 
					y_equals_x_max=10, 
					s=2, 
					x_axis_label='', 
					y_axis_label='', 
					title='', 
					scatter_point_color_groups=['b'], 
					custom_x_tick_locators=[-1, -1], 
					dashes=(8, 2), 
					x_error_groups=[], 
					y_error_groups=[], 
					opt_x_and_y_max=[], 
					y_axis_on_right=False, 
					opt_min_lim_buffer=0.0, 
					optional_y_max=-1.0, 
					axis_to_apply_to=None, 
					hide_right_and_top_splines=False, 
					legend_labels=[]):

	"""
	Custom general scatter plot function
	"""

	sns.set()
	sns.set_style("ticks")

	if (axis_to_apply_to):
		ax = axis_to_apply_to
	else:
		plt.figure()
		ax = plt.gca()
		
	if (x_axis_label != ''):
		ax.set_xlabel(x_axis_label, fontsize=24)
	if (y_axis_label != ''):
		ax.set_ylabel(y_axis_label, fontsize=24)
	if (title != ''):
		ax.set_title(title)

	sns.despine()

	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(24) 
	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(24) 
	
	if (x_error_groups != []):
		for x_points, y_points, scatter_point_color, x_errors, y_errors in zip(x_point_groups, y_point_groups, scatter_point_color_groups, x_error_groups, y_error_groups):
			ax.errorbar(x_points, y_points, xerr=x_errors, yerr=y_errors, c=scatter_point_color, fmt='o', elinewidth=0.1, ms=0.5, zorder=2)
	else:
		index=0
		for x_points, y_points, scatter_point_color in zip(x_point_groups, y_point_groups, scatter_point_color_groups):

			if (legend_labels != []):
				ax.scatter(x_points, y_points, c=scatter_point_color, s=s, zorder=2, label=legend_labels[index]) 
			else:
				ax.scatter(x_points, y_points, c=scatter_point_color, s=s, zorder=2) 
			index += 1

	if (custom_x_tick_locators[0] != -1):
		ax.xaxis.set_major_locator(MultipleLocator(custom_x_tick_locators[0]))
		ax.yaxis.set_major_locator(MultipleLocator(custom_x_tick_locators[0]))

		ax.xaxis.set_minor_locator(MultipleLocator(custom_x_tick_locators[1]))
		ax.yaxis.set_minor_locator(MultipleLocator(custom_x_tick_locators[1]))


	if (draw_y_equals_x):
		x = np.linspace(0,int(y_equals_x_max), int(y_equals_x_max+1))
		y = x
		if (dashes != None):
			ax.plot(x, y, '--k', linewidth=2, dashes=dashes, zorder=1)
		else:
			ax.plot(x, y, '--k', zorder=1)


	if (draw_y_equals_x):
		ax.set_xlim([0 - opt_min_lim_buffer, y_equals_x_max])
		ax.set_ylim([0 - opt_min_lim_buffer, y_equals_x_max])
		ax.set_aspect('equal', adjustable='box')
	elif (opt_x_and_y_max != []):
		ax.set_aspect('equal', adjustable='box')
		ax.set_xlim([0 - opt_min_lim_buffer, opt_x_and_y_max[0]])
		ax.set_ylim([0 - opt_min_lim_buffer, opt_x_and_y_max[1]])

	if (y_axis_on_right):
		ax.yaxis.tick_right()


	if (optional_y_max != -1.0):
		ax.set_ylim([0.0, optional_y_max])

	if (hide_right_and_top_splines):
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)

	if (legend_labels != []):
		ax.legend()
		

	if (not axis_to_apply_to):
		plt.savefig(file_name + '.pdf', bbox_inches='tight')
		plt.close()



"""
ANGLE HELPERS
Helper functions for creating plots related to correlation angles
"""

def plot_angle_vs_reliability_plots(means, 
									pca_overall_reliabilities, 
									pca_conj_reliabilities, 
									directory, 
									filename):

	"""
	Function for plotting correlation angles vs reliability
	"""

	pca_conj_reliabilities = np.asarray(pca_conj_reliabilities)
	pca_overall_reliabilities = np.asarray(pca_overall_reliabilities)
	
	plt.figure()
	plt.scatter(means, pca_conj_reliabilities)
	plt.gca().set_xlim([0.0, 45.0])
	plt.savefig(directory + filename + "_vs_conjreliability.pdf")
	plt.close()

	plt.figure()
	plt.scatter(means, pca_overall_reliabilities)
	plt.gca().set_xlim([0.0, 45.0])
	plt.savefig(directory + filename + "_vs_overallreliability.pdf")
	plt.close()



def plot_angle_confidence_bound_plots(means, 
									confidence_bounds_lower, 
									confidence_bounds_upper, 
									is_still_correlated, 
									is_PCA_BS_empirical_pvalue_different_from_45, 
									directory, 
									filename, 
									ylim_upper=46.0):

	"""
	Function for plotting correlation angles with confidence intervals
	"""

	means = np.asarray(means)
	confidence_bounds_lower = np.asarray(confidence_bounds_lower)
	confidence_bounds_upper = np.asarray(confidence_bounds_upper)
	confidence_bounds_lower_plus_means = means - confidence_bounds_lower
	confidence_bounds_upper_plus_means = means + confidence_bounds_upper

	is_still_correlated = np.asarray(is_still_correlated)
	is_PCA_BS_empirical_pvalue_different_from_45 = np.asarray(is_PCA_BS_empirical_pvalue_different_from_45)
	number_of_points = np.sum(is_still_correlated)

	indices_with_confidence_bounds_between_0_and_45 = np.argwhere(((confidence_bounds_lower_plus_means > 0.0) & (confidence_bounds_upper_plus_means < 45.0)) & (is_still_correlated)).flatten()
	indices_with_confidence_bounds_outside_0_and_45 = np.argwhere(((confidence_bounds_lower_plus_means < 0.0) | (confidence_bounds_upper_plus_means > 45.0)) & (is_still_correlated)).flatten()
	number_of_indices_with_confidence_bound_between_0_and_45 = len(indices_with_confidence_bounds_between_0_and_45)
	number_of_indices_with_confidence_bound_outside_0_and_45 = len(indices_with_confidence_bounds_outside_0_and_45)

	means_with_confidence_bounds_between_0_and_45 = means[indices_with_confidence_bounds_between_0_and_45]
	lower_confidence_bounds_with_confidence_bounds_between_0_and_45 = confidence_bounds_lower[indices_with_confidence_bounds_between_0_and_45]
	upper_confidence_bounds_with_confidence_bounds_between_0_and_45 = confidence_bounds_upper[indices_with_confidence_bounds_between_0_and_45]

	means_with_confidence_bounds_outside_0_and_45 = means[indices_with_confidence_bounds_outside_0_and_45]
	lower_confidence_bounds_with_confidence_bounds_outside_0_and_45 = confidence_bounds_lower[indices_with_confidence_bounds_outside_0_and_45]
	upper_confidence_bounds_with_confidence_bounds_outside_0_and_45 = confidence_bounds_upper[indices_with_confidence_bounds_outside_0_and_45]

	means_sorted_indices_between_45 = np.argsort(means_with_confidence_bounds_between_0_and_45)
	means_sorted_indices_outside_45 = np.argsort(means_with_confidence_bounds_outside_0_and_45)
	sorted_means_sorted_indices_between_45 = means_with_confidence_bounds_between_0_and_45[means_sorted_indices_between_45]
	sorted_means_sorted_indices_outside_45 = means_with_confidence_bounds_outside_0_and_45[means_sorted_indices_outside_45]

	sorted_lower_confidence_bounds_with_confidence_bounds_between_0_and_45 = lower_confidence_bounds_with_confidence_bounds_between_0_and_45[means_sorted_indices_between_45]
	sorted_upper_confidence_bounds_with_confidence_bounds_between_0_and_45 = upper_confidence_bounds_with_confidence_bounds_between_0_and_45[means_sorted_indices_between_45]

	sorted_lower_confidence_bounds_with_confidence_bounds_outside_0_and_45 = lower_confidence_bounds_with_confidence_bounds_outside_0_and_45[means_sorted_indices_outside_45]
	sorted_upper_confidence_bounds_with_confidence_bounds_outside_0_and_45 = upper_confidence_bounds_with_confidence_bounds_outside_0_and_45[means_sorted_indices_outside_45]


	plt.figure(figsize=(8, 3))
	sns.set()
	sns.set_style("ticks")

	if (number_of_points > 0):
		for tick in plt.gca().xaxis.get_major_ticks():
			tick.label.set_fontsize(14) 
		for tick in plt.gca().yaxis.get_major_ticks():
			tick.label.set_fontsize(14) 

		plt.gca().set_xlabel('Cluster index', fontsize=14)
		plt.gca().set_ylabel('Degrees', fontsize=14)

		x = np.linspace(0,number_of_points,number_of_points)
		y = np.linspace(45,45,number_of_points)
		plt.plot(x, y, '--k', linewidth=2, dashes=(8, 2))

		placement_counter_1 = 0
		for i in range(number_of_indices_with_confidence_bound_between_0_and_45):
			yerr = np.asarray([[sorted_lower_confidence_bounds_with_confidence_bounds_between_0_and_45[i], sorted_upper_confidence_bounds_with_confidence_bounds_between_0_and_45[i]]]).T
			plt.errorbar(i + 1, sorted_means_sorted_indices_between_45[i], yerr=yerr, color='b', fmt='o', elinewidth=0.4, ms=0.5)

		for i in range(number_of_indices_with_confidence_bound_outside_0_and_45):
			mean = sorted_means_sorted_indices_outside_45[i]

			lower = sorted_lower_confidence_bounds_with_confidence_bounds_outside_0_and_45[i]
			if (mean - sorted_lower_confidence_bounds_with_confidence_bounds_outside_0_and_45[i] < 0.0):
				lower = mean

			upper = sorted_upper_confidence_bounds_with_confidence_bounds_outside_0_and_45[i]
			if (mean + sorted_upper_confidence_bounds_with_confidence_bounds_outside_0_and_45[i] > 45.0):
				upper = 45.0 - mean

			yerr = np.asarray([[lower, upper]]).T
			plt.errorbar(number_of_indices_with_confidence_bound_between_0_and_45 + 1 + i, mean, yerr=yerr, color='r', fmt='o', elinewidth=0.4, ms=1.0)

		plt.gca().set_xlim([0, number_of_points])		
		plt.gca().set_ylim([0, ylim_upper])

		plt.gca().xaxis.set_major_locator(MultipleLocator(number_of_points))
		plt.gca().xaxis.set_minor_locator(MultipleLocator(number_of_points))
		plt.gca().yaxis.set_major_locator(MultipleLocator(45))
		plt.gca().yaxis.set_minor_locator(MultipleLocator(5))
		plt.gca().spines['right'].set_visible(False)
		plt.gca().spines['top'].set_visible(False)

	plt.savefig(directory + filename + '.pdf', bbox_inches='tight')
	plt.close()



"""
PAIR CORTEX SEPERATION
Plotting and preprocessing functions
"""

def draw_cortex_spatial_seperation_plot(dictionary_of_counts, 
										final_file_name_and_directory, 
										atleast_one_in_septum_counts):

	"""
	Function for generating correlated pairs cortex spatial seperation plot
	"""

	figure_background_colour = 'silver'
	linewidth_multiple = 1.5
	# colour_map_string = 'YlGnBu'
	colour_map_string = "Greens"
	layer_strings = ['L5B', 'L5A', 'L4', 'L23']

	plt.figure(figsize=(14, 12), facecolor=figure_background_colour)
	sns.set()
	sns.set_style("ticks")
	
	location_dictionary = {}


	for column_index in range(4):
		for layer_index, layer_string in enumerate(layer_strings):

			y_pos = float(layer_index)*3.0 + 1.5

			location_dictionary[str(column_index) + layer_string + 'EXC'] = {}
			location_dictionary[str(column_index) + layer_string + 'EXC']['column_index'] = column_index
			location_dictionary[str(column_index) + layer_string + 'EXC']['x'] = column_index * 60.0 + 10.0
			location_dictionary[str(column_index) + layer_string + 'EXC']['y'] = y_pos
			location_dictionary[str(column_index) + layer_string + 'EXC']['c'] = 'r'

			location_dictionary[str(column_index) + layer_string + 'INH'] = {}
			location_dictionary[str(column_index) + layer_string + 'INH']['column_index'] = column_index
			location_dictionary[str(column_index) + layer_string + 'INH']['x'] = column_index * 60.0 + 30.0
			location_dictionary[str(column_index) + layer_string + 'INH']['y'] = y_pos
			location_dictionary[str(column_index) + layer_string + 'INH']['c'] = 'b'

		

	list_of_counts = list(dictionary_of_counts.values())
	if (list_of_counts != []):

		max_count = np.max(list_of_counts)
		my_cmap = cm.get_cmap(colour_map_string, max_count + 1)

		for count_key, count in dictionary_of_counts.items():
			count_keys = count_key.split('_')

			count_key_0 = count_keys[0]
			count_key_1 = count_keys[1]

			if ((count_key_0[0] != 'S') & (count_key_1[0] != 'S')):

				if (count > 0):

					c = my_cmap(count)

					plt.plot([location_dictionary[count_key_0]['x'], location_dictionary[count_key_1]['x']], 
							[location_dictionary[count_key_0]['y'], location_dictionary[count_key_1]['y']], 
							linestyle='-', color=c, lw=linewidth_multiple*count, zorder=0, label=count)

					plt.scatter([location_dictionary[count_key_0]['x']], [location_dictionary[count_key_0]['y']], c=location_dictionary[count_key_0]['c'], s=1000, zorder=2)
					plt.scatter([location_dictionary[count_key_1]['x']], [location_dictionary[count_key_1]['y']], c=location_dictionary[count_key_1]['c'], s=1000, zorder=2)


				if (count_key_0 == count_key_1):

					plt.gca().annotate(str(count), (location_dictionary[count_key_0]['x'], location_dictionary[count_key_0]['y']), size=15, c='w')


		plt.gca().set_ylim([0.0, 12.0])
		plt.gca().set_xlim([0, 216])

		start_buff = 0.0185

		weird_width = 0.147
		weird_gap = 0.131
		for y_val in [0.0, 3.0, 6.0, 9.0, 12.0]:
			plt.axhline(y=y_val, xmin=start_buff, xmax=start_buff + weird_width, c='k')
			plt.axhline(y=y_val, xmin=start_buff + weird_width + weird_gap, xmax=start_buff + 2.0*weird_width + weird_gap, c='k')
			plt.axhline(y=y_val, xmin=start_buff + 2.0*weird_width + 2.0*weird_gap, xmax=start_buff + 3.0*weird_width + 2.0*weird_gap, c='k')
			plt.axhline(y=y_val, xmin=start_buff + 3.0*weird_width + 3.0*weird_gap, xmax=start_buff + 4.0*weird_width + 3.0*weird_gap, c='k')

		plt.axvline(x=4, ymin=0, ymax=1, c='k')
		plt.axvline(x=36, ymin=0, ymax=1, c='k')
		plt.axvline(x=64, ymin=0, ymax=1, c='k')
		plt.axvline(x=96, ymin=0, ymax=1, c='k')
		plt.axvline(x=124, ymin=0, ymax=1, c='k')
		plt.axvline(x=156, ymin=0, ymax=1, c='k')
		plt.axvline(x=184, ymin=0, ymax=1, c='k')
		plt.axvline(x=216, ymin=0, ymax=1, c='k')

		plt.gca().axis('off')

		hand, labels = plt.gca().get_legend_handles_labels()
		labels = [int(i) for i in labels]
		unique_labels = np.unique(labels)
		unique_labels = np.sort(unique_labels)
		list_of_unique_handles = []
		for unique_label in unique_labels:
			for label_index, label in enumerate(labels):
				if (label == unique_label):
					list_of_unique_handles.append(hand[label_index])
					break
		plt.legend(list_of_unique_handles, unique_labels, facecolor=figure_background_colour, edgecolor='k')


		plt.title(str(atleast_one_in_septum_counts))


	plt.savefig(final_file_name_and_directory, facecolor=figure_background_colour)
	plt.close()


def layer_column_neuron_type_preprocessing(file_dir_and_name, unit_pair_type, key, sdds, pdds, suf):

	"""
	Preprocessing for correlated pairs cortex spatial seperation plot
	"""

	dictionary_of_counts = {}
	atleast_one_in_septum_count = 0

	for secondary_cluster_index in range(len(sdds['principle_condition_bool_for_each_unit'])):

		if (sdds[key + '_tests_passed'][secondary_cluster_index]):

			means = [sdds['Original_BS_PCA_mean_of_mean0'][secondary_cluster_index], sdds['Original_BS_PCA_mean_of_mean1'][secondary_cluster_index]]
			principle_condition_bools = sdds['principle_condition_bool_for_each_unit'][secondary_cluster_index]
			layer_strings_for_unit_pair = sdds["neurophys_layer_strings_of_unit_pair"][secondary_cluster_index]
			exc_inh_types = sdds["exc_inh_types_of_unit_pair"][secondary_cluster_index]
			single_clustering = sdds['single_clustering'][secondary_cluster_index]
						
			layer_string_0 = layer_strings_for_unit_pair[0]
			layer_string_1 = layer_strings_for_unit_pair[1]
			layer_0_index = sw.neurophys_layer_group_strings.index(layer_string_0)
			layer_1_index = sw.neurophys_layer_group_strings.index(layer_string_1)

			exc_inh_type_0 = exc_inh_types[0]
			exc_inh_type_1 = exc_inh_types[1]
			exc_inh_type_index_0 = sw.exc_inh_types.index(exc_inh_type_0)
			exc_inh_type_index_1 = sw.exc_inh_types.index(exc_inh_type_1)

			num_units_in_principle_column = np.sum(principle_condition_bools)

			if ((single_clustering.scsuo_0.barrel_index_for_plot > -1) & (single_clustering.scsuo_1.barrel_index_for_plot > -1)):

				# If both in principle column, GOOD
				if (num_units_in_principle_column == 2):
					column_index_for_plot_0 = 1
					column_index_for_plot_1 = 1

				# If just one in principle column, GOOD
				if (num_units_in_principle_column == 1):
					if (single_clustering.scsuo_0.barrel_index_for_plot == single_clustering.scsuo_0.stimulated_whisker_index_for_plot):
						column_index_for_plot_0 = 1
						column_index_for_plot_1 = abs(single_clustering.scsuo_0.barrel_index_for_plot - single_clustering.scsuo_1.barrel_index_for_plot) + 1

					if (single_clustering.scsuo_1.barrel_index_for_plot == single_clustering.scsuo_0.stimulated_whisker_index_for_plot):
						column_index_for_plot_0 = abs(single_clustering.scsuo_0.barrel_index_for_plot - single_clustering.scsuo_1.barrel_index_for_plot) + 1
						column_index_for_plot_1 = 1
						

				# If neither in principle column
				if (num_units_in_principle_column == 0):

					# If in same barrel, GOOD
					if (single_clustering.scsuo_0.barrel_index_for_plot == single_clustering.scsuo_1.barrel_index_for_plot):

						if (abs(single_clustering.scsuo_0.barrel_index_for_plot - single_clustering.scsuo_0.stimulated_whisker_index_for_plot) == 1):
							column_index_for_plot_0 = 2
							column_index_for_plot_1 = 2

						elif (abs(single_clustering.scsuo_0.barrel_index_for_plot - single_clustering.scsuo_0.stimulated_whisker_index_for_plot) == 2):
							column_index_for_plot_0 = 3
							column_index_for_plot_1 = 3

					# If neighbouring each other, unit 0 neighbouring principal, GOOD
					elif (single_clustering.scsuo_1.barrel_index_for_plot - single_clustering.scsuo_0.barrel_index_for_plot == 1):
						column_index_for_plot_0 = 2
						column_index_for_plot_1 = 3

					# If neighbouring each other, unit 1 in neighbouring principal, GOOD
					elif (single_clustering.scsuo_0.barrel_index_for_plot - single_clustering.scsuo_1.barrel_index_for_plot == 1):
						column_index_for_plot_0 = 3
						column_index_for_plot_1 = 2

					# If seperated by two 
					elif (abs(single_clustering.scsuo_0.barrel_index_for_plot - single_clustering.scsuo_1.barrel_index_for_plot) == 2):
						# If either side of the principle column
						if ((single_clustering.scsuo_0.barrel_index_for_plot + 1 == single_clustering.scsuo_0.stimulated_whisker_index_for_plot) & (single_clustering.scsuo_1.barrel_index_for_plot - 1 == single_clustering.scsuo_0.stimulated_whisker_index_for_plot)):
							column_index_for_plot_0 = 0
							column_index_for_plot_1 = 2

						if ((single_clustering.scsuo_0.barrel_index_for_plot - 1 == single_clustering.scsuo_0.stimulated_whisker_index_for_plot) & (single_clustering.scsuo_1.barrel_index_for_plot + 1 == single_clustering.scsuo_0.stimulated_whisker_index_for_plot)):
							column_index_for_plot_0 = 2
							column_index_for_plot_1 = 0

					# If seperated by 3
					elif (abs(single_clustering.scsuo_0.barrel_index_for_plot - single_clustering.scsuo_1.barrel_index_for_plot) == 3):
						
						if (abs(single_clustering.scsuo_0.barrel_index_for_plot - single_clustering.scsuo_0.stimulated_whisker_index_for_plot) == 1):
							column_index_for_plot_0 = 0
							column_index_for_plot_1 = 3

						elif (abs(single_clustering.scsuo_1.barrel_index_for_plot - single_clustering.scsuo_0.stimulated_whisker_index_for_plot) == 1):
							column_index_for_plot_0 = 3
							column_index_for_plot_1 = 0

						else:
							print(single_clustering.scsuo_0.barrel_index_for_plot, single_clustering.scsuo_1.barrel_index_for_plot, single_clustering.scsuo_0.stimulated_whisker_index_for_plot)

			elif ((single_clustering.scsuo_0.barrel_index_for_plot == -1) | (single_clustering.scsuo_1.barrel_index_for_plot == -1)):

				atleast_one_in_septum_count += 1

			if (exc_inh_type_index_0 == exc_inh_type_index_1):
				exc_inh_pair_type_index = exc_inh_type_index_0
			elif (exc_inh_type_index_0 == 0):
				exc_inh_pair_type_index = 2
			elif (exc_inh_type_index_0 == 1):
				exc_inh_pair_type_index = 3


			key_0 = str(column_index_for_plot_0) + layer_string_0 + exc_inh_type_0
			key_1 = str(column_index_for_plot_1) + layer_string_1 + exc_inh_type_1
			pair_key_list = [key_0, key_1]
			pair_key_list.sort()

			pair_key = pair_key_list[0] + "_" + pair_key_list[1]
			if (pair_key not in dictionary_of_counts):
				dictionary_of_counts[pair_key] = 1
			else:
				dictionary_of_counts[pair_key] += 1

	return dictionary_of_counts, atleast_one_in_septum_count




"""
MISC PLOTTING FUNCTIONS
Misc collection of plotting functions
"""

def draw_basic_similarity_matrix(matrix, 
								file_name, 
								figsize=None):

	if (figsize == None):
		plt.figure()
	else:
		plt.figure(figsize=figsize)
	plt.imshow(matrix, cmap='hot')
	plt.savefig(file_name)
	plt.close()


def plot_cluster_reliability_plots(pca_overall_reliabilities, 
									pca_conj_reliabilities, 
									directory, 
									suf):

	pca_conj_reliabilities = np.asarray(pca_conj_reliabilities)
	pca_overall_reliabilities = np.asarray(pca_overall_reliabilities)

	plt.figure()
	plt.hist(pca_conj_reliabilities)
	plt.savefig(directory  + suf + "ConjReliabilityHist.pdf")
	plt.close()

	plt.figure()
	plt.hist(pca_overall_reliabilities)
	plt.savefig(directory + suf + "OverallReliabilityHist.pdf")
	plt.close()

	plt.figure()
	plt.scatter(pca_overall_reliabilities, pca_conj_reliabilities)
	plt.savefig(directory + suf + "OverallVsConjReliability.pdf")
	plt.close()
	
	normal_histo_plot([pca_overall_reliabilities, pca_conj_reliabilities], 
						directory + suf + "OverallAndConjHist", 
						bins=20, 
						histo_range=[0.0, 1.0], 
						x_axis_label="Reliability", 
						y_axis_label="Frequency", 
						custom_x_tick_locators=[1.0, 0.2], 
						custom_y_tick_locators=[20, 20], 
						alpha=0.78)


def draw_spike_count_histograms(dictionary_of_info_by_neuron_group, 
								file_name_and_directory):
	
	plt.figure(figsize=(3,16.5)) 

	sns.set()
	sns.set_style("ticks")


	r = []
	names = []
	raw_data = {'OneBars': [], 'TwoBars': [],'ThreeBars': [],'FourBars': [], 'FiveBars': []}
	

	for group_key_index in range(len(dictionary_of_info_by_neuron_group.keys())):
		r.append(group_key_index)

		group_key = list(dictionary_of_info_by_neuron_group.keys())[group_key_index]
		neurophys_neuron_group_dict = dictionary_of_info_by_neuron_group[group_key]

		names.append(group_key[:-2])
		# names.append(group_key)

		n,_ = np.histogram(neurophys_neuron_group_dict['non_zero_spike_counts'], bins=np.arange(1, 7) - 0.5)


		raw_data['OneBars'].append(float(n[0]))
		raw_data['TwoBars'].append(float(n[1]))
		raw_data['ThreeBars'].append(float(n[2]))
		raw_data['FourBars'].append(float(n[3]))
		raw_data['FiveBars'].append(float(n[4]))

	df = pd.DataFrame(raw_data)
	 
	totals = [i+j+k+l for i,j,k,l,m in zip(df['OneBars'], df['TwoBars'], df['ThreeBars'], df['FourBars'], df['FiveBars'])]
	OneBars = [i / j * 1.0 if j > 0 else 0.0 for i,j in zip(df['OneBars'], totals)]
	TwoBars  = [i / j * 1.0 if j > 0 else 0.0 for i,j in zip(df['TwoBars'], totals)]
	ThreeBars = [i / j * 1.0 if j > 0 else 0.0 for i,j in zip(df['ThreeBars'], totals)]
	FourBars = [i / j * 1.0 if j > 0 else 0.0 for i,j in zip(df['FourBars'], totals)]
	FiveBars = [i / j * 1.0 if j > 0 else 0.0 for i,j in zip(df['FiveBars'], totals)]
	 
	barWidth = 0.75
	plt.bar(r, OneBars, color='b', edgecolor='none', width=barWidth, alpha=0.88, linewidth=0.0)
	plt.bar(r, TwoBars, bottom=OneBars, color='g', edgecolor='none', width=barWidth, alpha=0.88, linewidth=0.0)
	plt.bar(r, ThreeBars, bottom=[i+j for i,j in zip(OneBars, TwoBars)], color='r', edgecolor='none', width=barWidth, alpha=0.88, linewidth=0.0)
	plt.bar(r, FourBars, bottom=[i+j+k for i,j,k in zip(OneBars, TwoBars, ThreeBars)], color='y', edgecolor='none', width=barWidth, alpha=0.88, linewidth=0.0)
	plt.bar(r, FiveBars, bottom=[i+j+k+l for i,j,k,l in zip(OneBars, TwoBars, ThreeBars, FourBars)], color='c', edgecolor='none', width=barWidth, alpha=0.88, linewidth=0.0)
	 
	plt.xticks(r, names, rotation=90)

	plt.gca().set_ylim([0.0, 1.0])


	for tick in plt.gca().xaxis.get_major_ticks():
		tick.label.set_fontsize(16) 
	for tick in plt.gca().yaxis.get_major_ticks():
		tick.label.set_fontsize(24) 

	plt.gca().yaxis.set_major_locator(MultipleLocator(1.0))
	plt.gca().yaxis.set_minor_locator(MultipleLocator(0.1))

	plt.gca().set_ylabel("Normalised count", fontsize=24)


	plt.savefig(file_name_and_directory, bbox_inches='tight')
	plt.close()




def draw_stimulus_frequency_histo(frequencies, 
								final_file_name_and_directory):
	
	plt.figure(figsize=(7, 5))
	sns.set()
	sns.set_style("ticks")
	freq_bins =   [0.0, 0.999, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5]
	freq_labels = ['0-1', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
	n, _ = np.histogram(frequencies, bins=freq_bins)
	plt.bar(freq_labels, n)
	plt.gca().set_xlabel('Stimulation frequency', fontsize=24)
	plt.gca().set_ylabel('Frequency', fontsize=24)

	plt.gca().yaxis.set_major_locator(MultipleLocator(20))
	plt.gca().yaxis.set_minor_locator(MultipleLocator(10))

	for tick in plt.gca().xaxis.get_major_ticks():
		tick.label.set_fontsize(24) 
	for tick in plt.gca().yaxis.get_major_ticks():
		tick.label.set_fontsize(24) 

	plt.savefig(final_file_name_and_directory, bbox_inches='tight')
	plt.close()



def plot_basic_raster_and_density_plot(first_spike_pairs_original, 
									analysis_recipe, 
									experiment_condition_units_string):

	buf = 0.0

	min_x = np.min([first_spike_pairs_original[:, 0].min() - buf, first_spike_pairs_original[:, 1].min() - buf])
	max_x = np.max([first_spike_pairs_original[:, 0].max() - buf, first_spike_pairs_original[:, 1].max() - buf])
	
	plt.figure()
	plt.scatter(first_spike_pairs_original[:, 0] - buf, first_spike_pairs_original[:, 1] - buf)
	plt.gca().set_xlim([min_x, max_x])
	plt.gca().set_ylim([min_x, max_x])

	plt.gca().set_xlabel("Time (ms)")
	plt.gca().set_ylabel("Time (ms)")

	plt.gca().spines['right'].set_visible(False)
	plt.gca().spines['top'].set_visible(False)

	plt.plot([min_x, max_x], [min_x, max_x], 'k-', alpha=0.75, zorder=1)
	plt.savefig(analysis_recipe['directory_holder'].basic_scatter_plot_directory + experiment_condition_units_string + "_SCATTER.pdf", bbox_inches='tight')
	plt.close()


	# k = kde.gaussian_kde(self.first_spike_pairs_original.T)
	# x_diff = max_x - min_x
	# xbins = x_diff * 2
	# xi, yi = np.mgrid[min_x:max_x:xbins*1j, min_x:max_x:xbins*1j]
	# zi = k(np.vstack([xi.flatten(), yi.flatten()]))

	# plt.figure()
	# plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.BuGn_r)
	# plt.savefig(self.directory_holder.basic_scatter_plot_directory + self.experiment_condition_units_string + "_DENSITY")
	# plt.close()

def plot_fa_comparison(compoarison_tests_tests_passed_not_normal_and_normal, 
							fa_tests_tests_passed_not_normal_and_normal,
							filename,
							y_equals_x_max,
							optional_y_max=-1.0,
							legend_labels=[],
							scatter_point_color_groups=['b', 'r', 'g']):

	basic_x_y_plot(compoarison_tests_tests_passed_not_normal_and_normal, 
										fa_tests_tests_passed_not_normal_and_normal, 
										filename, 
										s=4, 
										draw_y_equals_x=True, 
										y_equals_x_max=y_equals_x_max, 
										x_axis_label='$\sigma$', 
										y_axis_label='$\sigma$', 
										scatter_point_color_groups=scatter_point_color_groups, 
										custom_x_tick_locators=[6, 1],
										optional_y_max=optional_y_max,
										legend_labels=legend_labels)



"""
CUSTOM PAIRWISE RESPONSE DISTRIBUTION FIGURE
"""

import os
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1" 

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FixedLocator)

from shapely.geometry import Polygon
import descartes
import numpy as np
import glob

from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from statsmodels.graphics.gofplots import qqplot


def make_pairwise_distribution_analysis_plots(analysis_recipe, 
											first_spike_pairs_original, 
											experiment_condition_units_string, 
											single_clusterings_by_shuffle_option, 
											condition_index, 
											experiment_code, 
											tw_files_by_stwcboo_by_shuffle_option):


	if (analysis_recipe['do_plot_super_time_warp_cluster_busters']):

		# self.tw_files_by_stwcboo = []
		for stwcboo_index, stwcboo in enumerate(analysis_recipe['super_time_warp_cluster_buster_option_objects']):
			if (stwcboo['create_pdf'] | stwcboo['create_video']):
				plot_super_time_warp_cluster_busters(stwcboo_index, 
													stwcboo,
													analysis_recipe, 
													single_clusterings_by_shuffle_option, 
													condition_index, 
													experiment_code, 
													tw_files_by_stwcboo_by_shuffle_option, 
													experiment_condition_units_string)


	if (analysis_recipe['do_plot_basic_spike_pair_rasters']):
		plot_basic_raster_and_density_plot(first_spike_pairs_original, analysis_recipe, experiment_condition_units_string)



def plot_super_time_warp_cluster_busters(stwcboo_index, 
										stwcboo, 
										analysis_recipe, 
										single_clusterings_by_shuffle_option, 
										condition_index, 
										experiment_code, 
										tw_files_by_stwcboo_by_shuffle_option, 
										experiment_condition_units_string):	

	for shuffle_option in analysis_recipe['option_helper'].shuffle_options:
		shuffle_option_string = shuffle_option[0] + str(shuffle_option[1])

		list_of_files = []

		shuffle_option_pair_cluster_root_dir = analysis_recipe['directory_holder'].pair_cluster_plots_directory + shuffle_option_string + "/"
		sw.mkdir(shuffle_option_pair_cluster_root_dir)
		unique_pair_directory = shuffle_option_pair_cluster_root_dir + experiment_condition_units_string + '/'
		
		for single_clustering in single_clusterings_by_shuffle_option[shuffle_option_string]:

			if (single_clustering.do_use_clusters_in_analysis):

				do_make = analysis_recipe['do_plot_super_time_warp_cluster_buster_videos_SINGLE_PAIR']
				if (not do_make):
					for single_cluster_holder in single_clustering.secondary_single_cluster_holders:
						if (single_cluster_holder.final_stage_2_cluster):
							do_make = True

				if (do_make):

					unique_pair_and_options_file_full_path = [unique_pair_directory, stwcboo['type_string'] + "_" + single_clustering.option_string_with_experiment_condition_units_string]
					option_directory = shuffle_option_pair_cluster_root_dir + single_clustering.option_string  + "/"

					plot_super_time_warp_cluster_buster(stwcboo,
															single_clustering, 
															unique_pair_and_options_file_full_path, 
															option_directory, 
															single_clustering.experiment_condition_units_string_no_dots, 
															analysis_recipe['do_plot_super_time_warp_cluster_buster_videos_SINGLE_PAIR'], 
															condition_index, 
															experiment_code)

					if (stwcboo['create_video']):
						tw_files_by_stwcboo_by_shuffle_option[shuffle_option_string][stwcboo_index].append(sorted(glob.glob(unique_pair_and_options_file_full_path[0] + unique_pair_and_options_file_full_path[1] + "*"))[0])
						list_of_files.append(sorted(glob.glob(unique_pair_and_options_file_full_path[0] + unique_pair_and_options_file_full_path[1] + "*"))[0])

			
		if (analysis_recipe['do_plot_super_time_warp_cluster_buster_videos_SINGLE_PAIR'] & (len(list_of_files) > 0)):
			create_gif_from_list_of_image_file_names(list_of_files, unique_pair_directory + experiment_condition_units_string + '_' + shuffle_option_string + '.mp4')
			create_gif_from_list_of_image_file_names(list_of_files, analysis_recipe['directory_holder'].pair_cluster_plots_directory + experiment_condition_units_string + '_' + shuffle_option_string + '.mp4')




def draw_fa_line_and_variance_ellipse(ax, FA_mean_0_mean, FA_component_0_mean, FA_mean_1_mean, FA_component_1_mean, FA_noise_variance_0_mean, FA_noise_variance_1_mean, variance_scaling_factor):

	# Factor analysis line
	ax.plot([FA_mean_0_mean, FA_mean_0_mean + variance_scaling_factor*FA_component_0_mean], [FA_mean_1_mean, FA_mean_1_mean + variance_scaling_factor*FA_component_1_mean], 'm-', lw=1.5)
	ax.plot([FA_mean_0_mean, FA_mean_0_mean - variance_scaling_factor*FA_component_0_mean], [FA_mean_1_mean, FA_mean_1_mean - variance_scaling_factor*FA_component_1_mean], 'm-', lw=1.5)

	patch = sw.create_ellipse([FA_mean_0_mean, FA_mean_1_mean], variance_scaling_factor*np.sqrt(FA_noise_variance_0_mean), variance_scaling_factor*np.sqrt(FA_noise_variance_1_mean), 0.0)
	ax.add_patch(descartes.PolygonPatch(patch, fc=[0.0, 0.0, 0.0, 0.0], ec='m', linewidth=1.5, linestyle='-'))


def draw_single_analysis_type_column(axes, column_index, text_string, single_cluster_holder, indices_to_use, spike_pairs_to_use, analysis_dict_key, collective_scatter_lower_lims, collective_scatter_upper_lims, new_colour):

	scatter_point_size = 20

	axes[0][column_index].set_title(text_string, size='small')
	means = [single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['Mean_0'], single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['Mean_1']]
	axes[0][column_index].scatter(indices_to_use, spike_pairs_to_use[: , 0] - means[0], c='b', s=17)
	axes[0][column_index].scatter(indices_to_use, spike_pairs_to_use[: , 1] - means[1], c='r', s=17)
	collective_scatter_lower_lims[0].append(axes[0][column_index].get_ylim()[0])
	collective_scatter_upper_lims[0].append(axes[0][column_index].get_ylim()[1])

	axes[1][column_index].scatter(spike_pairs_to_use[: , 0] - means[0], spike_pairs_to_use[: , 1] - means[1], c=indices_to_use, s=scatter_point_size)
	
	number_of_points = 10
	x = np.linspace(0,45,number_of_points)
	y = np.linspace(0,45,number_of_points)
	axes[1][column_index].plot(x, y, '--k', linewidth=2)
	if (single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['tests_passed']):

		elipsoid_face_colour = [0.0, 0.0, 0.0, 0.0]
		axes[1][column_index].add_patch(descartes.PolygonPatch(single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['BS_PCA_ellipses_mean_zeroed'][1], fc=elipsoid_face_colour, ec=new_colour, linewidth=1.5, linestyle='-'))

		FA_mean_0_mean = 0.0
		FA_mean_1_mean = 0.0
		FA_component_0_mean = single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['FA_component_0_mean']
		FA_component_1_mean = single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['FA_component_1_mean']
		FA_noise_variance_0_mean = single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['FA_noise_variance_0_mean']
		FA_noise_variance_1_mean = single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['FA_noise_variance_1_mean']
		draw_fa_line_and_variance_ellipse(axes[1][column_index], FA_mean_0_mean, FA_component_0_mean, FA_mean_1_mean, FA_component_1_mean, FA_noise_variance_0_mean, FA_noise_variance_1_mean, 4)
	lims = [np.min([spike_pairs_to_use[: , 0] - means[0], spike_pairs_to_use[: , 1] - means[1]]) - 0.5, np.max([spike_pairs_to_use[: , 0] - means[0], spike_pairs_to_use[: , 1]  - means[1]]) + 0.5]
	axes[1][column_index].set_xlim(lims)
	axes[1][column_index].set_ylim(lims)
	# collective_scatter_lower_lims.append(lims[0])
	# collective_scatter_upper_lims.append(lims[1])
	axes[1][column_index].set_aspect('equal')
	collective_scatter_lower_lims[1].append(axes[1][column_index].get_ylim()[0])
	collective_scatter_upper_lims[1].append(axes[1][column_index].get_ylim()[1])

	lag_range = range(1, 11)
	plot_acf(spike_pairs_to_use[: , 0], ax=axes[2, column_index], c='b', lags=lag_range)
	plot_acf(spike_pairs_to_use[: , 1], ax=axes[2, column_index], c='r', lags=lag_range)
	collective_scatter_lower_lims[2].append(axes[2][column_index].get_ylim()[0])
	collective_scatter_upper_lims[2].append(axes[2][column_index].get_ylim()[1])
	axes[2, column_index].set_xticks(lag_range)

	lag_range = range(2, 11)
	plot_pacf(spike_pairs_to_use[: , 0], ax=axes[3, column_index], method='ols', c='b', lags=lag_range, alpha=0.025)
	plot_pacf(spike_pairs_to_use[: , 1], ax=axes[3, column_index], method='ols', c='r', lags=lag_range, alpha=0.025)
	collective_scatter_lower_lims[3].append(axes[3][column_index].get_ylim()[0])
	collective_scatter_upper_lims[3].append(axes[3][column_index].get_ylim()[1])
	axes[3, column_index].set_xticks(lag_range)


	qqplot(spike_pairs_to_use[:, 0], ax=axes[4, column_index], line='s', c='b', fit=True)
	qqplot(spike_pairs_to_use[:, 1], ax=axes[4, column_index], line='s', c='r', fit=True)
	collective_scatter_lower_lims[4].append(axes[4][column_index].get_ylim()[0])
	collective_scatter_upper_lims[4].append(axes[4][column_index].get_ylim()[1])


	return column_index



def plot_super_time_warp_cluster_buster(super_time_warp_cluster_buster_option_object, single_clustering, unique_pair_and_options_file_full_path, outlier_contamination_and_clustering_param_0_directory, experiment_condition_units_string_no_dots, do_plot_super_time_warp_cluster_buster_videos_SINGLE_PAIR, condition_index, experiment_code):

	stwcboo = super_time_warp_cluster_buster_option_object
	single_clustering = single_clustering

	sns.set()
	sns.set_style("ticks")
	
	if (stwcboo['do_plot_extra']):

		fig, axes = plt.subplots(5, 6, figsize=(24, 20)) 

		axes[3, 0].axis('off')
		axes[3, 1].axis('off')
		axes[3, 2].axis('off')
		axes[4, 0].axis('off')
		axes[4, 1].axis('off')
		axes[4, 2].axis('off')

		main_plot_position = (.0, .36, .52, .52) # [left, bottom, width, height]
		main_plot_ax = plt.subplot2grid((5,6), (0,0), colspan=3, rowspan=3)

	else:

		fig, axes = plt.subplots(2, 1, figsize=(10, 17)) 

		axes[1].axis('off')

		main_plot_ax = axes[0]
		main_plot_position = (.1, .3, .85, .65)


	time_warp_axis_and_text(main_plot_ax, single_clustering, stwcboo, condition_index, experiment_code, main_plot_position) #[left, bottom, width, height]


	if (stwcboo['do_plot_extra']):
		for secondary_cluster_index, single_cluster_holder in enumerate(single_clustering.secondary_single_cluster_holders):

			if (single_cluster_holder.final_stage_2_cluster):

				new_colour = single_clustering.rgbas_for_non_singleton_original_clusters[single_cluster_holder.non_singleton_index]

				start_column_index = 3

				collective_scatter_lower_lims = [[], [], [], [], []]
				collective_scatter_upper_lims = [[], [], [], [], []]
				text_string = ""
				draw_column_index = draw_single_analysis_type_column(axes, start_column_index, text_string, single_cluster_holder, single_cluster_holder.stage_2_cluster.cluster_indices, single_cluster_holder.stage_2_cluster.cluster_spike_pairs, 'Original', collective_scatter_lower_lims, collective_scatter_upper_lims, new_colour)

				if ((single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['tests_passed']) | (not single_cluster_holder.stage_2_cluster.use_selective_differences)):

					for row_index in range(0, 5):
						for column_index in range(4, 6):
							axes[row_index, column_index].axis('off')

				elif (single_cluster_holder.stage_2_cluster.use_selective_differences):
					indices_starting_index = single_cluster_holder.stage_2_cluster.number_of_samples_in_ellipse - single_cluster_holder.stage_2_cluster.number_of_samples_for_selective_differences
					draw_column_index = draw_single_analysis_type_column(axes, draw_column_index + 1, text_string, single_cluster_holder, single_cluster_holder.stage_2_cluster.cluster_indices[indices_starting_index:], single_cluster_holder.stage_2_cluster.selective_differences_scaled_down, 'SelectivelyDifferenced', collective_scatter_lower_lims, collective_scatter_upper_lims, new_colour)
					draw_column_index = draw_single_analysis_type_column(axes, draw_column_index + 1, text_string, single_cluster_holder, single_cluster_holder.stage_2_cluster.cluster_indices[indices_starting_index:], single_cluster_holder.stage_2_cluster.selective_differences_arima_residuals, 'SelectivelyDifferencedBoxJenkins', collective_scatter_lower_lims, collective_scatter_upper_lims, new_colour)
					
					

					for column_index in range(4, 6):
 						axes[4, column_index].set_ylabel('')


					for row_index in range(0, 5):

 						lower = np.min(collective_scatter_lower_lims[row_index])
 						upper = np.max(collective_scatter_upper_lims[row_index])
 						
 						for column_index in range(3, 6):

 							if (row_index == 1):
 								axes[row_index, column_index].set_xlim([lower, upper])

 							axes[row_index, column_index].set_ylim([lower, upper])

				axes[0, 3].set_ylabel('Time (ms)')
				axes[1, 3].set_ylabel('Time (ms)')
				axes[2, 3].set_ylabel('p-value')
				axes[3, 3 ].set_ylabel('p-value')


				for row_index in range(0, 5):
					for column_index in range(3, 6):
						axes[row_index, column_index].spines['right'].set_visible(False)
						axes[row_index, column_index].spines['top'].set_visible(False)

				for column_index in range(3, 6):
					axes[2, column_index].set_title('')
					axes[3, column_index].set_title('')

					axes[0, column_index].set_xlabel('Trial')
					axes[1, column_index].set_xlabel('Time (ms)')
					axes[2, column_index].set_xlabel('Lag')
					axes[3, column_index].set_xlabel('Lag')

				break




	if (stwcboo['create_pdf']):
		file_type_ending = '.png'
		if (stwcboo['is_mainz']):
			file_type_ending = '.pdf'
		sw.mkdir(outlier_contamination_and_clustering_param_0_directory)
		plt.savefig(outlier_contamination_and_clustering_param_0_directory + stwcboo['type_string'] + '_STWCB_' + experiment_condition_units_string_no_dots + "_" + single_clustering.option_string + file_type_ending, bbox_inches='tight')
	

	if(stwcboo['create_video']):
		sw.mkdir(unique_pair_and_options_file_full_path[0])
		unique_pair_and_options_file_full_path = unique_pair_and_options_file_full_path[0] + unique_pair_and_options_file_full_path[1]
		plt.savefig(unique_pair_and_options_file_full_path + '.png', bbox_inches='tight')

	plt.close()


def return_standard_text_string_for_analysis_dict(single_cluster_holder, analysis_dict_key):

	STs_acf_pvalues_0_string = ", PACF p: "
	STs_acf_pvalues_1_string = ", PACF p: "
	for i in range(5):
		STs_acf_pvalues_0_string += str(sw.round_to_2(single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['STs_pacf_pvalues_0'][i])) + ", "
		STs_acf_pvalues_1_string += str(sw.round_to_2(single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['STs_pacf_pvalues_1'][i])) + ", "
		
	text_string = "B LR p: " + str(sw.round_to_1(single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['TI_Vs_STs_LR_0_pvalue'])) + ", ADF p:" + str(sw.round_to_2(single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['ADFuller_STs_0_pvalue'])) + ", KPSS p: " + str(sw.round_to_2(single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['KPSS_STs_0_pvalue'])) + STs_acf_pvalues_0_string
	text_string = text_string + "\nR LR  p: " + str(sw.round_to_1(single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['TI_Vs_STs_LR_1_pvalue'])) + ", ADF p:" + str(sw.round_to_2(single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['ADFuller_STs_1_pvalue'])) + ", KPSS p: " + str(sw.round_to_2(single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['KPSS_STs_1_pvalue'])) + STs_acf_pvalues_1_string
	text_string = text_string + "\nMulti-out (R^2): " + str(sw.round_to_2(single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['TI_Vs_STs_LR_multiple_output_rsquared']))
	text_string = text_string + "\nAngle: " + str(sw.round_to_3(single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['BS_PCA_mean_angle_up_to_45'])) + " -" + str(sw.round_to_3(single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['PCA_BS_empirical_CI_lower'])) + " +" + str(sw.round_to_3(single_cluster_holder.stage_2_cluster.analysis_dicts[analysis_dict_key]['PCA_BS_empirical_CI_upper']))
	
	return text_string



def time_warp_axis_and_text(ax, single_clustering, stwcboo, condition_index, experiment_code, main_plot_position):

	lims = calculate_lims(single_clustering.spike_pair_distribution_spikes)
	if (not stwcboo['is_mainz']):
		lims = [2000.0, 2030.0]

	if (stwcboo['optional_fixed_lims'] != []):
		lims = stwcboo['optional_fixed_lims']
	
	# lims = [7.0, 28.0]
	ax.set_xlim(lims); ax.set_ylim(lims)
	ticks = np.arange(lims[0], lims[1], 2.0)
	ax.set_yticks(ticks)
	ax.set_xticks(ticks)
	ax.set_xlabel('Time (ms)', fontsize=24); ax.set_ylabel('Time (ms)', fontsize=24)
	ax.set_position(main_plot_position)
	ax.plot(lims, lims, '--k', linewidth=2, dashes=(8, 2))
	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(24) 
	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(24) 
	ax.xaxis.set_major_locator(FixedLocator(lims))
	ax.yaxis.set_major_locator(FixedLocator(lims))
	ax.xaxis.set_minor_locator(MultipleLocator(1))
	ax.yaxis.set_minor_locator(MultipleLocator(1))

	
	if (stwcboo['do_plot_3sd_conjunctive_ring']):
		patch = sw.create_ellipse([single_clustering.conjunctive_first_spike_0s_mean, single_clustering.conjunctive_first_spike_1s_mean], 3*single_clustering.conjunctive_first_spike_0s_standard_deviation, 3*single_clustering.conjunctive_first_spike_1s_standard_deviation, 0.0)
		ax.add_patch(descartes.PolygonPatch(patch, fc=[0.0, 0.0, 0.0, 0.0], ec='orange', linewidth=1.5, linestyle='-'))
	if (stwcboo['do_plot_4sd_conjunctive_ring']):
		patch = sw.create_ellipse([single_clustering.conjunctive_first_spike_0s_mean, single_clustering.conjunctive_first_spike_1s_mean], 4*single_clustering.conjunctive_first_spike_0s_standard_deviation, 4*single_clustering.conjunctive_first_spike_1s_standard_deviation, 0.0)
		ax.add_patch(descartes.PolygonPatch(patch, fc=[0.0, 0.0, 0.0, 0.0], ec='orange', linewidth=1.5, linestyle='-'))

		

	###### DRAW PRIMARY CLUSTER POINTS ######
	non_singleton_index = 0
	for original_cluster_index in range(single_clustering.number_of_original_clusters):
		first_spike_pairs_for_original_cluster = single_clustering.first_spike_pairs_for_each_original_cluster[original_cluster_index]
		number_of_samples_in_original_cluster = single_clustering.number_of_samples_for_each_original_cluster[original_cluster_index]

		scatter_colour = 'k'
		if (number_of_samples_in_original_cluster > 1):
			scatter_colour = single_clustering.rgbas_for_non_singleton_original_clusters[non_singleton_index]
			single_cluster_holder = single_clustering.non_singleton_single_cluster_holders[non_singleton_index]
			
			if (stwcboo['do_plot_3sd_flat_bounding_ellipses_dashed']):	
				ax.add_patch(descartes.PolygonPatch(single_cluster_holder.flat_intersection_ellipse, fc=[0.0, 0.0, 0.0, 0.0], ec=scatter_colour, linewidth=1.5, linestyle='--'))

			# Unangled ellipse
			if (single_cluster_holder.final_stage_1_cluster & (single_cluster_holder in single_clustering.non_intersecting_non_singleton_single_cluster_holders)):
			# if (single_cluster_holder.stage_1_cluster != None):
				if (stwcboo['do_plot_3sd_flat_bounding_ellipses']):	
					ax.add_patch(descartes.PolygonPatch(single_cluster_holder.stage_1_cluster.unangled_ellipses[0], fc=[0.0, 0.0, 0.0, 0.0], ec=scatter_colour, linewidth=1.5, linestyle='-'))
				if (stwcboo['do_plot_4sd_flat_bounding_ellipses']):	
					ax.add_patch(descartes.PolygonPatch(single_cluster_holder.stage_1_cluster.unangled_ellipses[1], fc=[0.0, 0.0, 0.0, 0.0], ec=scatter_colour, linewidth=1.5, linestyle='-'))

			non_singleton_index = non_singleton_index + 1

		ax.scatter(first_spike_pairs_for_original_cluster[:, 0], first_spike_pairs_for_original_cluster[:, 1], c=scatter_colour, marker='.', s=50)

		

	###### DRAW SECONDARY CLUSTER STUFF ######
	for secondary_cluster_index, single_cluster_holder in enumerate(single_clustering.secondary_single_cluster_holders):

		new_colour = single_clustering.rgbas_for_non_singleton_original_clusters[single_cluster_holder.non_singleton_index]
		elipsoid_face_colour = [0.0,0.0,0.0,0.0]

		if (single_cluster_holder.final_stage_2_cluster):

			if (single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['tests_passed']):

				if (stwcboo['do_plot_fa_line_and_3sd_ellipse'] | stwcboo['do_plot_fa_line_and_4sd_ellipse']):
					FA_mean_0_mean = single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['FA_mean_0_mean']
					FA_mean_1_mean = single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['FA_mean_1_mean']
					FA_component_0_mean = single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['FA_component_0_mean']
					FA_component_1_mean = single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['FA_component_1_mean']
					FA_noise_variance_0_mean = single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['FA_noise_variance_0_mean']
					FA_noise_variance_1_mean = single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['FA_noise_variance_1_mean']

					if (stwcboo['do_plot_fa_line_and_3sd_ellipse']):
						draw_fa_line_and_variance_ellipse(ax, FA_mean_0_mean, FA_component_0_mean, FA_mean_1_mean, FA_component_1_mean, FA_noise_variance_0_mean, FA_noise_variance_1_mean, 3)
					if (stwcboo['do_plot_fa_line_and_4sd_ellipse']):
						draw_fa_line_and_variance_ellipse(ax, FA_mean_0_mean, FA_component_0_mean, FA_mean_1_mean, FA_component_1_mean, FA_noise_variance_0_mean, FA_noise_variance_1_mean, 4)

			# Angled pca ellipse
			if (stwcboo['do_plot_3sd_angled_bounding_ellipses']):
				ax.add_patch(descartes.PolygonPatch(single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['BS_PCA_ellipses'][0], fc=elipsoid_face_colour, ec=new_colour, linewidth=1.5, linestyle='-'))
			if (stwcboo['do_plot_3sd_angled_bounding_ellipses_dashed']):
				ax.add_patch(descartes.PolygonPatch(single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['BS_PCA_ellipses'][0], fc=elipsoid_face_colour, ec=new_colour, linewidth=1.5, linestyle='--'))
			if (stwcboo['do_plot_4sd_angled_bounding_ellipses']):
				ax.add_patch(descartes.PolygonPatch(single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['BS_PCA_ellipses'][1], fc=elipsoid_face_colour, ec=new_colour, linewidth=1.5, linestyle='-'))



	####### FIGURE TEXT #######
	spike_pair_string = "Stim Freq: " + str(single_clustering.stimulation_frequency_of_condition) + 'Hz (' + str(experiment_code) + ", Condition: " + str(condition_index) + ') , EPS: ' + str(sw.round_to_3(single_clustering.clustering_param_0)) + '\n'
	spike_pair_string += single_clustering.neurophys_layer_strings_of_unit_pair[0] + ' ' + single_clustering.exc_inh_types_of_unit_pair[0] + ', ' + single_clustering.neurophys_layer_strings_of_unit_pair[1] + ' ' + single_clustering.exc_inh_types_of_unit_pair[1]
	spike_pair_string += ', PC: ' + str(single_clustering.principle_condition_bool_for_each_unit[0]) + ', PC: ' + str(single_clustering.principle_condition_bool_for_each_unit[1]) + ", "
	spike_pair_string += str(single_clustering.number_of_conjunctive_trials) + '/' + str(single_clustering.number_of_trials_for_condition) + ' conjunctive responses'
	 	
	fig_text_count = 0
	for secondary_cluster_index, single_cluster_holder in enumerate(single_clustering.secondary_single_cluster_holders):

		if (single_cluster_holder.final_stage_2_cluster):

			fig_text_string = 4 * (len(single_clustering.secondary_single_cluster_holders) - fig_text_count) * '\n' + '\n' # May want to use in future, if overlap!

			new_colour = single_clustering.rgbas_for_non_singleton_original_clusters[single_cluster_holder.non_singleton_index]

			precision = 6
			cluster_string = spike_pair_string
			cluster_string += "\nSamples - Stage 1: " + str(single_cluster_holder.stage_1_cluster.number_of_samples_in_ellipse) + ", Stage 2: " + str(single_cluster_holder.stage_2_cluster.number_of_samples_in_ellipse) + ", "
			cluster_string += "Criteria: " + str(single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['tests_passed']) + " , Criteria and Normal: " + str(single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['tests_passed_and_normal']) + '\n'
			cluster_string += "Θ\u2084\u2085: " + "{:.{}f}".format(single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['BS_PCA_mean_angle_up_to_45'], 2) + "+-(" + "{:.{}f}".format(single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['PCA_BS_empirical_CI_upper_original'], 2) + ", " + "{:.{}f}".format(single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['PCA_BS_empirical_CI_lower_original'], 2) + ")"
			# fig_text_string += "\n{:.{}f}".format(single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['BS_PCA_angle_sd'], precision)
			cluster_string += ", Different from 45 p=" + "{:.{}f}".format(single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['PCA_BS_empirical_pvalue_different_from_45'], precision) + ", from 0 p=" + "{:.{}f}".format(single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['PCA_BS_empirical_pvalue_different_from_0'], precision)
			cluster_string += "\nLinear Regression (Stage 1) - p: " + "{:.{}f}".format(single_cluster_holder.stage_1_cluster.analysis_dicts['FlatClusterStats']['FlatCluster_LR_pvalue'], precision) + ", r: " + "{:.{}f}".format(single_cluster_holder.stage_1_cluster.analysis_dicts['FlatClusterStats']['FlatCluster_LR_rvalue'], precision)
			cluster_string += "\n$σ_{{diff}}$: RDR 1: {:.{}f}".format(single_clustering.conjunctive_first_spike_differences_standard_deviation, 1) + ", RDR 2: {:.{}f}".format(single_cluster_holder.stage_2_cluster.cluster_spike_pairs_differences_sd, 1) + ", RDR 3: {:.{}f}".format(single_cluster_holder.stage_2_cluster.analysis_dicts['Original']['FA_1sd_estimated_difference_area_BS_mean'], 1) 
			# cluster_string += "\n[key + '_factor_correlation_matrix_eigen_value_0']"


			if (stwcboo['do_plot_extra'] & single_cluster_holder.stage_2_cluster.use_selective_differences):
				cluster_string +="\nBlue: ARIMA(" + str(single_cluster_holder.stage_2_cluster.AR_p_0) + "," + str(single_cluster_holder.stage_2_cluster.index_to_use_0) + "," + str(single_cluster_holder.stage_2_cluster.MA_q_0) + ")" + ", Red: ARIMA(" + str(single_cluster_holder.stage_2_cluster.AR_p_1) + "," + str(single_cluster_holder.stage_2_cluster.index_to_use_1) + "," + str(single_cluster_holder.stage_2_cluster.MA_q_1) + ")"
				cluster_string += "\n\nDifferenced: " + return_standard_text_string_for_analysis_dict(single_cluster_holder, 'SelectivelyDifferenced')
				cluster_string += "\n\nARIMA: " + return_standard_text_string_for_analysis_dict(single_cluster_holder, 'SelectivelyDifferencedBoxJenkins')


			fig_text_count = fig_text_count + 1
			plt.figtext(0.005, 0.31, cluster_string, fontsize='x-large', color=new_colour, va='top', fontweight='semibold')


 			# ax.plot(np.asarray(first_spike_pairs_for_cluster[:, 0]), linear_regression_for_cluster_of_first_spike_pairs.intercept + linear_regression_for_cluster_of_first_spike_pairs.slope * np.asarray(first_spike_pairs_for_cluster)[:, 0], c=regression_colour, label='fitted line', linewidth=1)

	

	plt.figtext(0.005, 0.31, spike_pair_string, fontsize='x-large', color='k', va='top', fontweight='semibold')	

	ax.set_aspect('equal', adjustable='box')



def calculate_lims(first_spike_pairs):

	min_x = np.min(first_spike_pairs[:, 0])
	max_x = np.max(first_spike_pairs[:, 0])

	min_y = np.min(first_spike_pairs[:, 1])
	max_y = np.max(first_spike_pairs[:, 1])

	# max_both = np.min([np.ceil(np.max([max_x, max_y]) / 2.0) * 2 + 1.0, 32.0]) #+ 0.01
	# min_both = np.max([np.floor( np.min([min_x, min_y]) / 2.0 ) * 2, 5.0])

	max_both = np.min([np.ceil(np.max([max_x, max_y])) + 1.0, 55.0]) + 1.0
	min_both = np.max([np.floor(np.min([min_x, min_y])) - 1.0, 5.0]) - 1.0

	# min_both = 7.0

	lims = [min_both, max_both]

	return lims







