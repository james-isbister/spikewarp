import itertools

"""
Various constants and helpers. Particularly for neuron groups.
"""


principle_column_options = [True, False]
number_of_principle_column_types = len(principle_column_options)
exc_inh_types = ['EXC', 'INH']
exc_inh_type_colours = ['r', 'b']
number_of_exc_inh_types = len(exc_inh_types)

# ORIGINAL NEURO PHYS VARIABLES
neurophys_layer_original_data_indices = [2, 4, 5, 6]
neurophys_layer_original_data_names = ['L23', 'L4', 'L5A', 'L5B']
number_of_neurophys_layers = len(neurophys_layer_original_data_names)

# HELPERS
def return_neurophys_layer_string_given_neurophys_layer_index(neurophys_layer_index):
	index = neurophys_layer_original_data_indices.index(neurophys_layer_index)
	return neurophys_layer_original_data_names[index]


all_pairs_of_neurophys_layer_original_data_names = [[neurophys_layer_original_data_names[p1], neurophys_layer_original_data_names[p2]] for p1 in range(number_of_neurophys_layers) for p2 in range(p1,number_of_neurophys_layers)]
all_pairs_of_exc_inh_types = [['EXC', 'EXC'], ['INH', 'INH'], ['EXC', 'INH']]
all_pairs_of_principle_column_options = list(itertools.product(principle_column_options, principle_column_options))

combinations_of_all_pairs_of_neurophys_corresponding_layer_name_and_ei_types = [list(i) for i in itertools.product(all_pairs_of_neurophys_layer_original_data_names, all_pairs_of_exc_inh_types)]
number_of_combinations_of_all_pairs_of_neurophys_corresponding_layer_name_and_ei_types = len(combinations_of_all_pairs_of_neurophys_corresponding_layer_name_and_ei_types)

neurophys_layer_group_strings = ['L23', 'L4', 'L5A', 'L5B', 'L5']




def index_of_neurophys_layer_string_pair_ei_type_in_combination_list(layer_string_pair, exc_inh_type_pair): 

	for combination_index in range(number_of_combinations_of_all_pairs_of_neurophys_corresponding_layer_name_and_ei_types):

		layer_pair_and_ei_combination = combinations_of_all_pairs_of_neurophys_corresponding_layer_name_and_ei_types[combination_index]

		if (((layer_string_pair[0] == layer_pair_and_ei_combination[0][0]) & (layer_string_pair[1] == layer_pair_and_ei_combination[0][1])) | ((layer_string_pair[0] == layer_pair_and_ei_combination[0][1]) & (layer_string_pair[1] == layer_pair_and_ei_combination[0][0]))):

			if (((exc_inh_type_pair[0] == layer_pair_and_ei_combination[1][0]) & (exc_inh_type_pair[1] == layer_pair_and_ei_combination[1][1])) | ((exc_inh_type_pair[0] == layer_pair_and_ei_combination[1][1]) & (exc_inh_type_pair[1] == layer_pair_and_ei_combination[1][0]))):

				return combination_index


def index_of_principle_column_options_pair(principle_column_option_pair):

	for combination_index in range(4):

		combination = all_pairs_of_principle_column_options[combination_index]

		if ((combination[0] == principle_column_option_pair[0]) & (combination[1] == principle_column_option_pair[1])):

			return combination_index

def strings_for_combinations_of_all_pairs_of_neurophys_corresponding_layer_name_and_ei_types():

	string_list = []
	for index in range(number_of_combinations_of_all_pairs_of_neurophys_corresponding_layer_name_and_ei_types):

		string_list.append(str(combinations_of_all_pairs_of_neurophys_corresponding_layer_name_and_ei_types[index]))

	return string_list





