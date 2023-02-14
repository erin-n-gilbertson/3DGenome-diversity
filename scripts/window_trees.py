import pandas as pd
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist

def main():

	#Old
	# tree_data = {}
	# comparisons = pd.read_csv('dataframes/HFF_comparisons.txt', sep='\t', header=0)
	# windows_list = comparisons['window'].unique().tolist()
	
	# for window in windows_list:
	# 	window_df = comparisons[comparisons['window'] == window]
	# 	idx = sorted(set(window_df['ind1']).union(window_df['ind2']))
	# 	array = window_df.pivot(index='ind1', columns='ind2', values='divergence').reindex(index=idx, columns=idx).fillna(0, downcast='infer').pipe(lambda x: x+x.values.T).to_numpy()
	# 	condensed = pdist(array)

	# 	dendrogram = sch.dendrogram(sch.linkage(condensed, method = 'complete'))

	# 	cluster_sample_IDs = dendrogram['leaves']
	# 	cluster_IDs = dendrogram['leaves_color_list']
	# 	clusters_dict = dict(zip(cluster_sample_IDs, cluster_IDs))
		
	# 	cluster_IDs_list = []
	# 	for key, value in sorted(clusters_dict.items()):
	# 		cluster_IDs_list.append(value)

	# 	top_tree_y = dendrogram['dcoord'][-1]
	# 	window_stats = cluster_IDs_list + top_tree_y
	# 	tree_data[window] = window_stats

	# dataframe = pd.DataFrame.from_dict(tree_data, orient='index')
	# dataframe.to_csv('dataframes/window_trees.txt', sep = '\t', header = False, index = True)

	# New
	tree_data = {}
	comparisons = pd.read_csv('dataframes/HFF_comparisons.txt', sep='\t', header=0)
	windows_list = comparisons['window'].unique().tolist()

	for window in windows_list:
			window_df = comparisons[comparisons['window'] == window]

			if window.startswith('chrX_'):
					length = 36
			else:
					length = 56

			idx = sorted(set(window_df['ind1']).union(window_df['ind2']))
			array = window_df.pivot(index='ind1', columns='ind2', values='divergence').reindex(index=idx, columns=idx).fillna(0, downcast='infer').pipe(lambda x: x+x.values.T).to_numpy()
			condensed = array[np.triu_indices(length, k = 1)]
			dendrogram = sch.dendrogram(sch.linkage(condensed, method = args.linkage))

			cluster_sample_IDs = dendrogram['leaves']
			cluster_IDs = dendrogram['leaves_color_list']
			clusters_dict = dict(zip(cluster_sample_IDs, cluster_IDs))

			cluster_IDs_list = []
			for key, value in sorted(clusters_dict.items()):
					cluster_IDs_list.append(value)

			top_tree_y = dendrogram['dcoord'][-1]
			window_stats = cluster_IDs_list + top_tree_y
			tree_data[window] = window_stats

	dataframe = pd.DataFrame.from_dict(tree_data, orient='index')
	dataframe.to_csv(args.out, sep = '\t', header = False, index = True)


if __name__ == '__main__':
	main()



