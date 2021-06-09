import spikewarp as sw

def unpack_parralel_results(list_of_results, meta_dict_by_clustering_type, exp_cond_units_holders_by_clustering_type, exp_cond_pop_state_procs_by_clustering_type, exp_cond_schs_for_conds_by_clustering_type, clustering_types, state_types, dh_by_clustering_type):
	
	# Function for unpacking results of parallel analysis of single experimental conditions

	for result in list_of_results:

		for key1, val1 in result.meta_dict.items():
			if isinstance(val1, list):
				meta_dict_by_clustering_type[result.clustering_type][key1].extend(result.meta_dict[key1])
			if isinstance(val1, dict):
				for key2 in result.meta_dict[key1].keys():
					meta_dict_by_clustering_type[result.clustering_type][key1][key2] += result.meta_dict[key1][key2]


		if (exp_cond_units_holders_by_clustering_type != {}):
			exp_cond_units_holders_by_clustering_type[result.clustering_type].update(result.exp_cond_units_holders)
			exp_cond_pop_state_procs_by_clustering_type[result.clustering_type].update(result.exp_cond_pop_state_procs)
			exp_cond_schs_for_conds_by_clustering_type[result.clustering_type].update(result.exp_cond_schs_for_cond)


	for clustering_type in clustering_types:	 
		sw.cluster_pair_plots(meta_dict_by_clustering_type[clustering_type], dh_by_clustering_type[clustering_type])
		for state_type in state_types:
			if hasattr(sw, 'population_plots'):
				sw.population_plots(meta_dict_by_clustering_type[clustering_type], state_type, dh_by_clustering_type[clustering_type])


class ClusteringTypeIterator(object):

	# Object for parallel analysis of single experimental condition

	def __init__(self, clustering_type, plot_opt, state_types, schs_for_cond, dh, exp_cond_str, premade_exp_cond_units_holder, cluster_pair_types):

		meta_dict = sw.create_meta_dict(cluster_pair_types, state_types)
		exp_cond_units_holders = {}
		exp_cond_pop_state_procs = {}
		exp_cond_schs_for_cond = {}

		exp_cond_schs_for_cond[exp_cond_str] = schs_for_cond

		if (len(schs_for_cond) > 0):	

			if (premade_exp_cond_units_holder == None):
				non_stationary_to_pop, unique_pair_strings, exp_cond_units_holder = sw.process_schs(schs_for_cond, exp_cond_str)
				meta_dict['all_clusters']['unique_pair_strings'].extend(unique_pair_strings)
				meta_dict['including_non_stationary']['angle_conf_intervals'].extend([abs(sch.stage_2_cluster.analysis_dicts['Original']['PCA_BS_empirical_CI_lower']) + abs(sch.stage_2_cluster.analysis_dicts['Original']['PCA_BS_empirical_CI_upper']) for sch in schs_for_cond])
				meta_dict['including_non_stationary']['unique_pair_strings'].extend([sch.unique_pair_str for sch in schs_for_cond])

				state_analysis = True
				if (state_analysis):
					non_stationary_to_pop.reverse()
					for sch_index in non_stationary_to_pop:
						schs_for_cond.pop(sch_index)

			else:
				exp_cond_units_holder = premade_exp_cond_units_holder
				
			if (len(schs_for_cond) > 0):
				exp_cond_units_holders[exp_cond_str] = exp_cond_units_holder

				sw.process_sch_combinations(schs_for_cond, exp_cond_units_holder, exp_cond_str, meta_dict, dh, plot_opt)
				
				pop_state_proc = None
				if hasattr(sw, 'population_plots'):
					pop_state_proc = sw.PopulationStateProcessor(exp_cond_units_holder, schs_for_cond, exp_cond_str, dh, meta_dict, plot_opt)
				exp_cond_pop_state_procs[exp_cond_str] = pop_state_proc							

		self.exp_cond_schs_for_cond = exp_cond_schs_for_cond
		self.meta_dict = meta_dict
		self.clustering_type = clustering_type
		self.exp_cond_units_holders = exp_cond_units_holders
		self.exp_cond_pop_state_procs = exp_cond_pop_state_procs