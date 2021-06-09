import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def check_intersect_and_union_for_cluster_pair_type(clustering_type_a, clustering_type_b, meta_dict_by_clustering_type, cluster_pair_type, data_key, optional_figure_path=''):
	"""
	Function which tests for intersection and union of response distributions for two clustering types (i.e. DBSCAN and Unclustered)
	"""

	a = meta_dict_by_clustering_type[clustering_type_a][cluster_pair_type][data_key]
	b = meta_dict_by_clustering_type[clustering_type_b][cluster_pair_type][data_key]
	
	d = check_intersect_and_union(a, b, optional_figure_path=optional_figure_path, group_a_label=clustering_type_a, group_b_label=clustering_type_b)

	return d

def check_intersect_and_union(a, b, optional_figure_path='', group_a_label='', group_b_label=''):

	"""
	Function which calculates union and intersection of two lists
	"""

	d = {}

	d['a_set'] = set(a)
	d['b_set'] = set(b)
	d['intersection'] = d['a_set'].intersection(d['b_set'])
	d['union'] = d['a_set'].union(d['b_set'])
	d['a_only_set'] = d['a_set'] - d['intersection']
	d['b_only_set'] = d['b_set'] - d['intersection']


	d['a_set'] = list(d['a_set'])
	d['b_set'] = list(d['b_set'])
	d['intersection'] = list(d['intersection'])
	d['union'] = list(d['union'])
	d['a_only_set'] = list(d['a_only_set'])
	d['b_only_set'] = list(d['b_only_set'])

	d['a_set_num'] = len(d['a_set'])
	d['b_set_num'] = len(d['b_set'])
	d['intersection_num'] = len(d['intersection'])
	d['union_num'] = len(d['union'])
	d['a_only_set_num'] = len(d['a_only_set'])
	d['b_only_set_num'] = len(d['b_only_set'])

	if (optional_figure_path != ''):
		venn2(subsets=(d['a_only_set_num'], d['b_only_set_num'], d['intersection_num']), set_labels = (group_a_label, group_b_label), set_colors=('purple', 'skyblue'), alpha = 0.7)
		plt.savefig(optional_figure_path, bbox_inches='tight')
		plt.close()

	return d