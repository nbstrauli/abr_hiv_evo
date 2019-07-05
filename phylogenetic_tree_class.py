import sys
from ete3 import Tree, TreeStyle, NodeStyle, CircleFace, faces
import math
from colour import Color

class phylo_tree(object):
	"""
	This class is for analyzing and manipulating pylogenetic trees. It mostly uses the python package 'ete3'.

	Input:
	Reads in a Newick formatted text file.

	Attributes:
	filepath - Gives the path to the input newick file.
	tree - This is an ete3 Tree object, where the provided newick file in 'filepath' is loaded into it. Each of the leaf nodes in the tree have a name, and there are attributes of a given node encoded in the name. The attributes are delimited by a '|'. Each attribute is coded as: "|attribute_name=attribute_value|"
	count_attribute_name - This is the name (as a string) of the attribute in each element name that gives the count information. If this is None, then it is assumed that each element has a count of one.
	freq_attribute_name - This is the name (as a string) of the attribute in each element name that gives the frequency information. If this is None, and 'count_attribute_name' is defined, then the freq is calculated as the count of a node/total counts. If this is None and 'count_attribute_name' is None, then it is assumed that each element has an equal frequency in the data.
	time_points - A set that gives all the time points found in the data. This attribute is only defined if the method 'add_time_info' has been run.
	outgroup_label - This gives the label for the node that is intended to be the outgroup. If it is not defined then the tree is assumed to be 'unrooted'. If it is defined, then this label (as a string) should be in the 'ID' portion of the outgroup node. So the naming outgroup node should look like this: '1229.0_[outgroup_label]|something=stuff|blah=blahblah|etc...'. So the outgroup_label should be the last part of the node ID, delimited by '_'.
	total_count - This gives the sum of the counts in each node. If each node has a count of one (default), then this will equal the number of nodes.
	cluster_attribute_name - If defined (default, None), this gives the name of the attribute that gives the cluster ID information for each of the leaves.
	extra_leaf_features - This is a dictionary of names for extra features that were added to each of the leaves of the tree. The def of each of the feature names is a set type of each of the unique values (as a string) for that feature. The default for this attribute is an empty dic. These extra features are typically added using the method 'add_attribute_to_leaves'.
	"""

	def __init__(self, filepath, count_attribute_name=None, freq_attribute_name=None, outgroup_label=None, cluster_attribute_name=None):
		self.filepath = filepath
		self.count_attribute_name = count_attribute_name
		self.freq_attribute_name = freq_attribute_name
		self.time_points = None
		self.tree = Tree(filepath)
		self.outgroup_label = outgroup_label
		self.cluster_attribute_name = cluster_attribute_name
		self.total_count = 0.
		self.extra_leaf_features = {}
		#get total counts before doing anything else
		for leaf in self.tree:
			if count_attribute_name:
				leaf_name = leaf.name.split('|')
				for i in leaf_name[1:]:
					attribute = i.split('=')
					attribute_name = attribute[0]
					if attribute_name == count_attribute_name:
						count = int(attribute[1])
			else:
				count = 1
			self.total_count += count
		for leaf in self.tree:
			leaf_name = leaf.name.split('|')
			if outgroup_label:
				if outgroup_label in leaf_name[0].split('_'):
					self.tree.set_outgroup(leaf.name)
			for i in leaf_name[1:]:
				attribute = i.split('=')
				attribute_name = attribute[0]
				if attribute_name == count_attribute_name:
					count = int(attribute[1])
				elif attribute_name == freq_attribute_name:
					freq = float(attribute[1])
				elif attribute_name == cluster_attribute_name:
					leaf.add_features(cluster_id=attribute[1])
			if not count_attribute_name:
				count = 1
			if not freq_attribute_name:
				freq = count / self.total_count
			leaf.add_features(count=count)
			leaf.add_features(freq=freq)
		return

	def add_time_info(self, time_series_info='start_of_id'):
		"""
		This method adds time information to each of the nodes in the tree.
		time_series_info - This tells where the information of time-point can be found in each of the element IDs. Acceptable values are:
			'start_of_id' - This means the time info is at the very beginning of each element ID, and is separated by a '_'. For example: '45.3_blahblahblah' would have a time point of 45.3.
		"""
		time_points = set()
		for leaf in self.tree:
			if time_series_info == 'start_of_id':
				tpoint = float(leaf.name.split("_")[0])
			time_points.update([tpoint])
			leaf.add_features(time_point=tpoint)
		self.time_points = time_points
		return

	def add_attribute_to_leaves(self, attribute_name):
		"""
		This method will cycle through all the leaves in the tree, and add an attribute to each of them. The attribute needs to be an attribute that is encoded in the input newick file.
		attribute_name - This gives the name of the attribute (that is in the newick file) to add to each of he leaves.
		"""
		attribute_value_set = set()
		for leaf in self.tree:
			attributes = leaf.name.split('|')
			for attribute in attributes:
				attribute_list = attribute.split('=')
				if attribute_name == attribute_list[0]:
					leaf.add_feature(attribute_name, attribute_list[1])
					attribute_value_set.update([attribute_list[1]])
					break
		self.extra_leaf_features[attribute_name] = attribute_value_set
		return

	def plot_tree(self, output_filepath, time_series_info=None, tree_style='horizontal_right', leaf_size_map_to=None, show_leaf_names=False, color_branches_by=None, line_width=1, start_color='red', end_color='purple', ladderize=True):
		"""
		This method plots the phylogenetic tree as a dendogram. Uses the ete3 package to do this.
		output_filepath - This is the path to the output image file.
		time_series_info - This tells where the information of time-point can be found (if at all) in each of the element IDs. Acceptable values are:
			None - (default) This means there is no time information in the tree
			'start_of_id' - This means the time info is at the very beginning of each element ID, and is separated by a '_'. For example: '45.3_blahblahblah' would have a time point of 45.3.
		tree_style - This give the style of tree plotting, i.e. circular, horizontal, etc. Acceptable values are:
			'half_circle' - Tree is plotted as a half circle, where branches are radiating outward and upward.
			'horizontal_right' - Tree is plotted as a normal dendogram where branching occurs from left to right.
		leaf_size_map_to - If defined (default, None), this will inform what attribute of the node the leaf size maps to. Acceptable values are:
			None - No leaf size info. All leaves the same size
			'count' - leaf size proportional to the count attribute
			'freq' - leaf size proportional to the freq attribute
		show_leaf_names - If True (default, False), then will plot the names of each of the leafs of the tree.
		color_branches_by - If defined (default, None) this will give how the leaf branches will be colored, if at all. Acceptable values are:
			None - Default. No branch coloring
			'time_point' - This means that the branches will be colored according to the 'time_point' attribute of each leaf node. There must be a 'time_point' attribute in the leaf nodes for this to work. So, one should run the 'add_time_info' method before doing this. Alternatively, one can define the 'time_series_info' parameter for this method, and this will be taken care of.
			any other string - This will give the name of any other attribute of the leaf nodes for which to map the value to the leaf branch color.
		then this will instruct to color the different time-points with different colors. This is ignored if 'time_series_info' is False
		line_width - controls the line width. Default, 1
		start_color - This gives the staring color. This should be a string the spells out the name of a color. Most simple colors should be fine. Uses the 'Color' module from 'colour'.
		end_color - This gives the ending color. This should be a string the spells out the name of a color. Most simple colors should be fine. Uses the 'Color' module from 'colour'. The colors used for each unique attribute that gives the colors will span the spectrom from 'start_color' to 'end_color'.
		ladderize - If True (default), this will ladderize the tree. That is it will sort the partitions of the internal nodes based upon number of decendant nodes in the child nodes.
		"""
		if not self.time_points and time_series_info:
			self.add_time_info(time_series_info=time_series_info)
		#set time-point colors, if desired
		if color_branches_by:
			start_color = Color(start_color)
			end_color = Color(end_color)
			if color_branches_by == 'time_point' and time_series_info:
				colors = list(start_color.range_to(end_color, len(self.time_points)))
				hex_colors = [i.hex_l for i in colors]
				tpoint_to_color_dic = {}
				for index, i in enumerate(sorted(self.time_points)):
					tpoint_to_color_dic[i] = hex_colors[index]
			else:
				#check if attribute already exists in the tree data. if not, add it
				for leaf in self.tree:
					if not color_branches_by in leaf.features:
						self.add_attribute_to_leaves(attribute_name=color_branches_by)
					break
				colors = list(start_color.range_to(end_color, len(self.extra_leaf_features[color_branches_by])))
				hex_colors = [i.hex_l for i in colors]
				attribute_to_color_dic = {}
				for index, i in enumerate(sorted(self.extra_leaf_features[color_branches_by])):
					attribute_to_color_dic[i] = hex_colors[index]
		#set node styles
		most_dist_leaf, size_to_tree_size_scaler = self.tree.get_farthest_leaf() #need to scale sizes by the length (divergence) of the tree
		for node in self.tree.traverse():
			node_style = NodeStyle()
			#do stuff to leaf nodes
			if node.is_leaf():
				if color_branches_by:
					if color_branches_by == 'time_point':
						color = tpoint_to_color_dic[node.time_point]
					else:
						color = attribute_to_color_dic[getattr(node, color_branches_by)]
				else:
					color = Color('black')
					color = color.hex_l
				node_style['hz_line_color'] = color
				node_style['vt_line_color'] = color
				if leaf_size_map_to:
					if leaf_size_map_to == 'count':
						radius = size_to_tree_size_scaler * 100 * math.log(node.count)
					elif leaf_size_map_to == 'freq':
						radius = size_to_tree_size_scaler * 100 * node.freq
					c = CircleFace(radius=radius, color=color, style='circle')
					c.opacity = 0.3
					node.add_face(c, 0, position='branch-right')
					node_style['size'] = 0
			node_style['hz_line_width'] = size_to_tree_size_scaler * 10 * line_width
			node_style['vt_line_width'] = size_to_tree_size_scaler * 10 * line_width
			node.set_style(node_style)
		#set tree style
		tree_steeze = TreeStyle()
		if tree_style == 'half_circle':
			tree_steeze.mode = 'c'
			tree_steeze.arc_start = -180
			tree_steeze.arc_span = 180
		elif tree_style == 'horizontal_right':
			pass
		if show_leaf_names:
			tree_steeze.show_leaf_name = True
		else:
			tree_steeze.show_leaf_name = False
		self.tree.ladderize()
		self.tree.render(output_filepath, w=700, h=700, units='mm', tree_style=tree_steeze)
		return

	def assign_outgroup(self, outgroup_node_name):
		"""
		This method takes a given node (within the existing tree structure) and assigns it as the root node to the tree.
		"""
		self.tree.set_outgroup(outgroup_node_name)
		return

	def find_outlier_clade(self, outlier_attribute_name, get_ancestor_of_outier_clade=False):
		"""
		This method will find the clade that contains all the outlier leaf nodes. It does this by finding the most recent common ancestor node to all the outlier seqs.
		outlier_attribute_name - This gives the name of the node attribute that tells if the leaf is an outlier or not. The value to this attribute should be boolean. So, 'True' if outlier, and 'False', if not. If the 'outlier_attribute_name' does not already exist as an attribute in the nodes, then the method 'add_attribute_to_leaves' is attempted to be used to add it.
		get_ancestor_of_outier_clade - If True (default, False), this will use the direct ancestor of the 'outlier_clade_root' to define the outlier clade. Use this if the most recent common ancestor of all the labeled outliers doesn't quite capture all the seqs that we want to label.

		OUPUT:
		Returns 'seq_ids' which is a list of the sequence IDs (or id component of each leaf-node) that is part of the outlier clade.
		"""
		for leaf in self.tree:
			if not outlier_attribute_name in leaf.features:
				self.add_attribute_to_leaves(attribute_name=outlier_attribute_name)
			break
		outlier_nodes = []
		for leaf in self.tree:
			if getattr(leaf, outlier_attribute_name) == 'True':
				outlier_nodes.append(leaf)
		#if there are no outlier leaves, return empty list
		if not outlier_nodes:
			self.extra_leaf_features[outlier_attribute_name].update(['maybe'])
			return []
		#if there is only one outlier node, then it is it's own clade. Nothing else to do
		elif len(outlier_nodes) == 1 and not get_ancestor_of_outier_clade:
			self.extra_leaf_features[outlier_attribute_name].update(['maybe'])
			seq_id = outlier_nodes[0].name.split('|')[0]
			return [seq_id]
		outlier_clade_root = self.tree.get_common_ancestor(outlier_nodes)
		if get_ancestor_of_outier_clade:
			outlier_clade_root = outlier_clade_root.get_ancestors()[0]
		seq_ids = []
		for leaf in outlier_clade_root.get_leaves():
			seq_id = leaf.name.split('|')[0]
			seq_ids.append(seq_id)
			if getattr(leaf, outlier_attribute_name) == 'False':
				setattr(leaf, outlier_attribute_name, 'maybe')
		self.extra_leaf_features[outlier_attribute_name].update(['maybe'])
		return seq_ids


# tree = phylo_tree(filepath='/Users/nstrauli/data/abr_hiv_coevo/phylogenetic_trees/hiv/fasttree/haplotypes/10.tree', count_attribute_name='total_count')
# tree.assign_outgroup(outgroup_node_name='1229.0_154|total_count=5|total_freq=0.000154564283285')
# tree.plot_tree(output_filepath='/Users/nstrauli/Desktop/test/test_tree.pdf', time_series_info='start_of_id', tree_style='half_circle', color_branches_by_time_point=True)
