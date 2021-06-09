
"""
Helpers for adding cluster data to meta analysis dictionaries
"""


list_of_first_stage_data_names = [

'analysis_dicts',
'analysis_dict_member_keys',
	
]


list_of_second_stage_data_names = [

	'analysis_dicts',
	'analysis_dict_member_keys',

	'single_clustering',

	'PCA_ellipse_conj_reliability',
	'PCA_ellipse_overall_reliability',
	
	'angled_cluster_spike_pairs_differences_sd',

	'Unclustered_Conj_LR_pvalue', 'Unclustered_Conj_LR_rsquaredvalue',
	'Conjunctive_FS_SDs_N0', 'Conjunctive_FS_SDs_N1', 'Conjunctive_original_1sd_area', 'Conjunctive_original_FS_diff_SD',
	'NotOnlyConjunctive_FS_SDs_N0', 'NotOnlyConjunctive_FS_SDs_N1', 'NotOnlyConjunctive_original_1sd_area', 

	'stimulation_frequency_of_condition', 'neurophys_layer_strings_of_unit_pair', 'exc_inh_types_of_unit_pair', 'principle_condition_bool_for_each_unit'

]


def process(data_dict, ellipse_analysis_dicts, make_analysis_dict_keys):

	for (analysis_dict_key, analysis_dict) in ellipse_analysis_dicts.items():

		for (variable_key, variable_value) in analysis_dict.items():

			new_key = analysis_dict_key + "_" + variable_key

			if (make_analysis_dict_keys):
				if ('analysis_dict_member_keys' not in data_dict.keys()):
					data_dict['analysis_dict_member_keys'].append(variable_key)
				
			if (new_key not in data_dict.keys()):
				data_dict[new_key] = []

			data_dict[new_key].append(variable_value)



def add_single_cluster_holder_data_to_dict(data_dict, data_names, single_cluster_holder, single_clustering, cluster_type):

	for data_name in data_names:

		if (data_name == "analysis_dicts"):

			make_analysis_dict_keys = False
			if (data_dict['analysis_dict_member_keys'] == []):
				make_analysis_dict_keys = True

			if (cluster_type == "flat"):
				process(data_dict, single_cluster_holder.stage_1_cluster.analysis_dicts, make_analysis_dict_keys)
			elif (cluster_type == "angled"):
				process(data_dict, single_cluster_holder.stage_2_cluster.analysis_dicts, make_analysis_dict_keys)


		if (data_name == "single_clustering"):
			data_dict[data_name].append(single_clustering)

		if (data_name == "Unclustered_Conj_LR_pvalue"):
			data_dict[data_name].append(single_clustering.LR_all_conjunctive.pvalue)
		if (data_name == "Unclustered_Conj_LR_rsquaredvalue"):
			data_dict[data_name].append(single_clustering.LR_all_conjunctive.rvalue*single_clustering.LR_all_conjunctive.rvalue)

		if (data_name == "angled_cluster_spike_pairs_differences_sd"):
			data_dict[data_name].append(single_cluster_holder.stage_2_cluster.cluster_spike_pairs_differences_sd)


		if (data_name == "Conjunctive_FS_SDs_N0"):
			data_dict[data_name].append(single_clustering.conjunctive_first_spike_0s_standard_deviation)
		if (data_name == "Conjunctive_FS_SDs_N1"):
			data_dict[data_name].append(single_clustering.conjunctive_first_spike_1s_standard_deviation)
		if (data_name == "Conjunctive_original_1sd_area"):
			data_dict[data_name].append(single_clustering.conjunctive_first_spike_0s_standard_deviation * single_clustering.conjunctive_first_spike_1s_standard_deviation)
		if (data_name == "Conjunctive_original_FS_diff_SD"):
			data_dict[data_name].append(single_clustering.conjunctive_first_spike_differences_standard_deviation)
		if (data_name == "NotOnlyConjunctive_FS_SDs_N0"):
			data_dict[data_name].append(single_clustering.N0_not_only_conjunctive_original_FS_SD)
		if (data_name == "NotOnlyConjunctive_FS_SDs_N1"):
			data_dict[data_name].append(single_clustering.N1_not_only_conjunctive_original_FS_SD)
		if (data_name == "NotOnlyConjunctive_original_1sd_area"):
			data_dict[data_name].append(single_clustering.N0_not_only_conjunctive_original_FS_SD * single_clustering.N1_not_only_conjunctive_original_FS_SD)


		if (data_name == "stimulation_frequency_of_condition"):
			data_dict[data_name].append(single_clustering.stimulation_frequency_of_condition)
		if (data_name == "neurophys_layer_strings_of_unit_pair"):
			data_dict[data_name].append(single_clustering.neurophys_layer_strings_of_unit_pair)
		if (data_name == "exc_inh_types_of_unit_pair"):
			data_dict[data_name].append(single_clustering.exc_inh_types_of_unit_pair)
		if (data_name == "principle_condition_bool_for_each_unit"):
			data_dict[data_name].append(single_clustering.principle_condition_bool_for_each_unit)



		if (data_name == "PCA_ellipse_conj_reliability"):
			data_dict[data_name].append(single_cluster_holder.stage_2_cluster.ellipse_conj_reliability)
		if (data_name == "PCA_ellipse_overall_reliability"):
			data_dict[data_name].append(single_cluster_holder.stage_2_cluster.ellipse_overall_reliability)


