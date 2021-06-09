"""
DIRECTORY HELPERS
Helper functions for creating directories
"""

import os
import errno

def mkdir(path):
	if (path != ''):
		try:
			os.mkdir(path)
		except OSError as exception:
			if exception.errno != errno.EEXIST:
				raise

def makedirs(path):
	if (path != ''):
		try:
			os.makedirs(path)
		except OSError as exception:
			if exception.errno != errno.EEXIST:
				raise



"""
JSON HELPERS
Helper functions for reading and writing JSON files
"""

import json

def load_json_recipe(fn):
	with open(fn, 'r') as fid:
	    recipe = json.load(fid)
	return recipe


def write_dict_to_json_file(dict, file):

	with open(file, 'w') as json_file:
	  json.dump(dict, json_file)



def save_str_to_file(str_to_save, file):

	with open(file, "w") as text_file:
   		print(str_to_save, file=text_file)



"""
PICKLE HELPERS
Helper functions for writing pickles to file
"""

import pickle

def pickle_clusters(analysis, single_cluster_holders_by_exp_cond_and_shuffle_option):

	with open(analysis['directory_holder'].collated_root_output_directory + 'single_cluster_holders_by_exp_cond_and_shuffle_option.pickle', 'wb') as handle:
		pickle.dump(single_cluster_holders_by_exp_cond_and_shuffle_option, handle, protocol=pickle.HIGHEST_PROTOCOL)


def pickle_directory_holders(analysis):
	
	with open(analysis['directory_holder'].collated_root_output_directory + 'directory_holder.pickle', 'wb') as handle:
		pickle.dump(analysis['directory_holder'], handle, protocol=pickle.HIGHEST_PROTOCOL)


def pickle_mcah(analysis, mcah_for_shuffle_type, shuffle_option_string):

	with open(analysis['directory_holder'].collated_root_output_directory + shuffle_option_string + '-mcah_for_shuffle_type.pickle', 'wb') as handle:
		pickle.dump(mcah_for_shuffle_type, handle, protocol=pickle.HIGHEST_PROTOCOL)