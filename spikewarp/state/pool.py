import copy

import spikewarp as sw

def create_conj_pool(exp_cond_units_holders_by_clustering_type, exp_cond_pop_state_procs_by_clustering_type, exp_cond_schs_for_conds_by_clustering_type, group_0_key, group_1_key, fig_dir, cluster_pair_types, state_types):

	"""
	Function for creating conjugative pool of neurons and clusters from two clustering type groups
	"""

	exp_cond_units_holders_0 = exp_cond_units_holders_by_clustering_type[group_0_key]
	exp_cond_units_holders_1 = exp_cond_units_holders_by_clustering_type[group_1_key]
	exp_cond_schs_for_conds_0s = exp_cond_schs_for_conds_by_clustering_type[group_0_key]
	exp_cond_schs_for_conds_1s = exp_cond_schs_for_conds_by_clustering_type[group_1_key]

	exp_cond_d = sw.check_intersect_and_union(list(exp_cond_units_holders_0.keys()), list(exp_cond_units_holders_1.keys()), optional_figure_path=fig_dir + 'EXP_COND_Venn.pdf', group_a_label=group_0_key, group_b_label=group_1_key)

	dh_by_conj = {'Conj':sw.create_dir_holder(fig_dir, "Conj")}
	meta_dict_conj = {'Conj':sw.create_meta_dict(cluster_pair_types, state_types)}


	exp_cond_unit_holders_conj_all = []
	schs_for_conds_conj_all = []
	exp_cond_key_conj_all = []
	for exp_cond_key in exp_cond_d['intersection']:
		
		exp_cond_unit_holder_0 = exp_cond_units_holders_0[exp_cond_key]	
		exp_cond_unit_holder_1 = exp_cond_units_holders_1[exp_cond_key]	

		exp_cond_schs_for_conds_0 = exp_cond_schs_for_conds_0s[exp_cond_key]
		exp_cond_schs_for_conds_1 = exp_cond_schs_for_conds_1s[exp_cond_key]


		exp_cond_unit_holder_conj = None
		if (exp_cond_unit_holder_0 != None):
			d = sw.check_intersect_and_union(exp_cond_unit_holder_0.d.keys(), exp_cond_unit_holder_1.d.keys())

			### REMOVE UNCLUSTERED FROM POOL ###
			exp_cond_unit_holder_conj = copy.copy(exp_cond_unit_holder_0)
			for unit_key, unit_dict in  exp_cond_unit_holder_1.d.items():
				if unit_key not in exp_cond_unit_holder_conj.d.keys():
					exp_cond_unit_holder_conj.d[unit_key] = unit_dict
		
		schs_for_cond_conj = copy.copy(exp_cond_schs_for_conds_0)

		unique_pair_strs_0 = [sch.unique_pair_str for sch in exp_cond_schs_for_conds_0]
		for sch in  exp_cond_schs_for_conds_1:
			if (sch.unique_pair_str not in unique_pair_strs_0):
				schs_for_cond_conj.append(sch)

		
		if (exp_cond_unit_holder_conj != None):
			exp_cond_unit_holder_conj.create_spikes_by_unit_array(schs_for_cond_conj)

		exp_cond_unit_holders_conj_all += [exp_cond_unit_holder_conj]
		schs_for_conds_conj_all += [schs_for_cond_conj]
		exp_cond_key_conj_all += [exp_cond_key]
		
	exp_cond_unit_holders_conj_all += [exp_cond_units_holders_0[exp_cond_key] for exp_cond_key in exp_cond_d['a_only_set']]
	exp_cond_unit_holders_conj_all += [exp_cond_units_holders_1[exp_cond_key] for exp_cond_key in exp_cond_d['b_only_set']]
	schs_for_conds_conj_all += [exp_cond_schs_for_conds_0s[exp_cond_key] for exp_cond_key in exp_cond_d['a_only_set']]
	schs_for_conds_conj_all += [exp_cond_schs_for_conds_1s[exp_cond_key] for exp_cond_key in exp_cond_d['b_only_set']]
	exp_cond_key_conj_all += exp_cond_d['a_only_set']
	exp_cond_key_conj_all += exp_cond_d['b_only_set']


	return exp_cond_unit_holders_conj_all, schs_for_conds_conj_all, exp_cond_key_conj_all, meta_dict_conj, dh_by_conj