from msa_class import msa
import numpy
import itertools
import math

class compare_MSAs(object):
	"""
	This is a class that will compare two MSA's to one another. 
	
	Input:
	Should be two filepaths that give the path to two different MSA's in fasta format. Ideally, both MSA's should contain time-point information, and both MSA's should have the same number of sequences. If not, not all the methods will work. Note that this method will read into memory the entire sequence information of both MSAs, so use caution if they are really big.

	Attributes:
	filepath_1 - Path the to fasta file that contains the first MSA
	filepath_2 - Path the to fasta file that contains the second MSA
	count_attribute_name - This is the name (as a string) of the attribute in the header that gives the count information for each sequence. If None, then it is assumed that each seq has a count of 1
	msa_1 - This is an object of class 'msa' for the first MSA (i.e. the MSA that corresponds to 'filepath_1'). It has all the methods contained in the msa_class (as well as the sequence_sample_class that it inherits from). See this class for explanations of all the attributes and methods associated with this class.
	msa_2 - This is the same as msa_1 but for the second input MSA.
	seq_[1|2]_array - This is a numpy array object for the sequences in the first or second MSA. Each sequence entry is a row, each column is a position.
	num_sites_[1|2] - Int. This gives the number of sites (i.e. columns) in MSA 1 or 2.
	remove_tpoints_that_arent_identical - This will cycle through the seq entries in each of the alignments, and will remove the entries that don't have an identical time-points in the cooresponding entry of the other alignment. Importantly, this assumes that i) the MSA's have time information, and ii) the MSA's entries are order in time.
	reduce_to_change_or_same - If True (default, False), this will cause the alignments to be altered upon initialization of the class such that a position will be coded as a '1' if that position is a changed relative to the same position in the row above it, or a '0' if there was no change. This altering of the alignment data only occurs in the 'seq_[1|2]_array' attributes. The rest of the data is unchanged.
	"""

	def __init__(self, filepath_1, filepath_2, count_attribute_name=None, remove_tpoints_that_arent_identical=False, reduce_to_change_or_same=False):
		self.filepath_1 = filepath_1
		self.filepath_2 = filepath_2
		self.count_attribute_name = count_attribute_name
		self.remove_tpoints_that_arent_identical = remove_tpoints_that_arent_identical
		self.msa_1 = msa(filepath=filepath_1, count_attribute_name=count_attribute_name, first_seq_ref=False)
		self.msa_2 = msa(filepath=filepath_2, count_attribute_name=count_attribute_name, first_seq_ref=False)
		self.reduce_to_change_or_same = reduce_to_change_or_same
		
		#remove the time-points that don't occur in both MSA's, if desired
		if remove_tpoints_that_arent_identical:
			tpoints_set_1 = set()
			for seq_entry in self.msa_1.data:
				tpoint = seq_entry['id'].split('_')[0]
				tpoints_set_1.update([tpoint])
			tpoints_set_2 = set()
			for seq_entry in self.msa_2.data:
				tpoint = seq_entry['id'].split('_')[0]
				tpoints_set_2.update([tpoint])
			set_1_unique = tpoints_set_1 - tpoints_set_2
			indeces_to_remove = []
			for index in xrange(len(self.msa_1.data)):
				tpoint = self.msa_1.data[index]['id'].split('_')[0]
				if tpoint in set_1_unique:
					indeces_to_remove.append(index)
			self.msa_1.remove_seq_entries(seq_indicators=indeces_to_remove, seq_ids=False)
			set_2_unique = tpoints_set_2 - tpoints_set_1
			indeces_to_remove = []
			for index in xrange(len(self.msa_2.data)):
				tpoint = self.msa_2.data[index]['id'].split('_')[0]
				if tpoint in set_2_unique:
					indeces_to_remove.append(index)
			self.msa_2.remove_seq_entries(seq_indicators=indeces_to_remove, seq_ids=False)
			if len(self.msa_1.data) < 2 or len(self.msa_2.data) < 2:
				print 'After removing unidentical time-point entries in the MSAs, there is not enough data. Aborting'
				return

		self.seq_1_array = []
		self.seq_2_array = []
		for index, msa_object in enumerate([self.msa_1, self.msa_2]):
			for seq_entry in msa_object.data:
				seq = [i for i in seq_entry['seq']]
				if index == 0:
					self.seq_1_array.append(seq)
				else:
					self.seq_2_array.append(seq)
		self.seq_1_array = numpy.array(self.seq_1_array)
		self.seq_2_array = numpy.array(self.seq_2_array)
		self.num_sites_1 = len(self.seq_1_array[0,])
		self.num_sites_2 = len(self.seq_2_array[0,])
		if reduce_to_change_or_same:
			new_seq_1_array = [['0' for i in self.seq_1_array[0]]]
			for row_index in xrange(1, len(self.seq_1_array)):
				new_seq_1_array.append([])
				for col_index in xrange(self.num_sites_1):
					if self.seq_1_array[row_index, col_index] != self.seq_1_array[row_index-1, col_index]:
						new_seq_1_array[-1].append('1')
					else:
						new_seq_1_array[-1].append('0')
			new_seq_2_array = [['0' for i in self.seq_2_array[0]]]
			for row_index in xrange(1, len(self.seq_2_array)):
				new_seq_2_array.append([])
				for col_index in xrange(self.num_sites_2):
					if self.seq_2_array[row_index, col_index] != self.seq_2_array[row_index-1, col_index]:
						new_seq_2_array[-1].append('1')
					else:
						new_seq_2_array[-1].append('0')
			self.seq_1_array = numpy.array(new_seq_1_array)
			self.seq_2_array = numpy.array(new_seq_2_array)
		return

	@staticmethod
	def calc_MI(vector1, vector2, log_base=2):
		"""
		This is a static method only meant to be used by other other methods. This will calculate the MI statistic for the two provided vectors. Each 'vector' is a list of characters. Each vector is meant to represent independent observations of two categorical random variables. 
		"""
		vec_1_states_prob_dic = {}
		vec_2_states_prob_dic = {}
		joint_prob_dic = {}
		num_obs = 0.
		for observation_1, observation_2 in itertools.izip(vector1, vector2):
			num_obs += 1
			try:
				vec_1_states_prob_dic[observation_1] += 1
			except KeyError:
				vec_1_states_prob_dic[observation_1] = 1
			try:
				vec_2_states_prob_dic[observation_2] += 1
			except KeyError:
				vec_2_states_prob_dic[observation_2] = 1
			try:
				joint_prob_dic['%s_%s' % (observation_1, observation_2)] += 1
			except KeyError:
				joint_prob_dic['%s_%s' % (observation_1, observation_2)] = 1
		for state in vec_1_states_prob_dic:
			vec_1_states_prob_dic[state] /= num_obs
		for state in vec_2_states_prob_dic:
			vec_2_states_prob_dic[state] /= num_obs
		for joint_state in joint_prob_dic:
			joint_prob_dic[joint_state] /= num_obs
		MI = 0.
		for state1 in vec_1_states_prob_dic:
			for state2 in vec_2_states_prob_dic:
				try:
					ratio = joint_prob_dic['%s_%s' % (state1, state2)] / (vec_1_states_prob_dic[state1]*vec_2_states_prob_dic[state2])
					MI += joint_prob_dic['%s_%s' % (state1, state2)] * math.log(ratio, log_base)
				except KeyError:
					continue
		return MI

	def calc_pairwise_mutual_information(self, output_array_filepath=None, output_distribution_filepath=None, ignore_sites_less_than_X_variants=None, ignore_sites_less_than_X_changes=None):
		"""
		This method will compare each column in one MSA to each column in the other MSA, and calculate the mutual information statistic for each comparison. 
		output_array_filepath - This is the path to the output file that will be a table containing the M.I. values for each column comparison. The names of the rows index the msa_1 sites, and the names of the columns index the msa_2 sites. If this is None, then this file is not written.
		output_distribution_filepath - This is the path to the file that will contain one long list (one column) of mutual information values. Each line is one MI value for one of the site comparisons.
		ignore_non_variable_sites - if True (default, False), this will cause the method to skip sites (columns) in either of the alignments that aren't variable, meaning that all seqs at this site are the same.
		ignore_sites_less_than_X_variants - Nonetype or Int. If defined (default, None), this will cause the method to ignore any sites (i.e. columns) that have a total number of variants that is less than this value.
		ignore_sites_less_than_X_changes - Nonetype or Int. If defined (default, None), this will cause the method to ignore any sites (i.e. columns) that have a total number of changes that is less than this value. A change is defined as a site that is different than the same site in the row directly preceding (above) it. This superceeds the 'ignore_sites_less_than_X_variants' parameter.

		Output:
		Returns a numpy array that is formatted just as 'output_array_filepath' is formatted (i.e. the rows index the msa_1 sites, and the columns index the msa_2 sites)
		"""
		#check to make sure that the MSAs have an equal number of entries
		if len(self.seq_1_array[:,0]) != len(self.seq_2_array[:,0]):
			print 'In order to calculate mutual information, the MSAs need to have the same number of sequence entries.'
			return
		#check to make sure that not using 'ignore_sites_less_than_X_variants' if this is change no change data
		if self.reduce_to_change_or_same and ignore_sites_less_than_X_variants:
			print "The 'ignore_sites_less_than_X_variants' is set but can't know the number of variants when the data has been reduced to change or no change. Aborting."
			return
		MI_values = []
		#find non-variable sites in MSA 2, if applicable
		if ignore_sites_less_than_X_changes:
			msa_2_sites_to_skip = set()
			for index, column_2 in enumerate(self.seq_2_array.T):
				if self.reduce_to_change_or_same:
					num_changes = 0
					for i in column_2:
						if i == '1':
							num_changes += 1
				else:
					num_changes = 0
					for row_index in xrange(1, len(column_2)):
						if column_2[row_index] != column_2[row_index-1]:
							num_changes += 1
						if num_changes >= ignore_sites_less_than_X_changes:
							break
				if num_changes < ignore_sites_less_than_X_changes:
					msa_2_sites_to_skip.update([index])
		elif ignore_sites_less_than_X_variants:
			msa_2_sites_to_skip = set()
			for index, column_2 in enumerate(self.seq_2_array.T):
				if len(set(column_2)) < ignore_sites_less_than_X_variants:
					msa_2_sites_to_skip.update([index])
		used_columns_1 = []
		used_columns_2 = []
		first_round = True
		for column_1_index, column_1 in enumerate(self.seq_1_array.T):
			#check if a variable site
			if ignore_sites_less_than_X_changes:
				if self.reduce_to_change_or_same:
					num_changes = 0
					for i in column_1:
						if i == '1':
							num_changes += 1
				else:
					num_changes = 0
					for row_index in xrange(1, len(column_1)):
						if column_1[row_index] != column_1[row_index-1]:
							num_changes += 1
						if num_changes >= ignore_sites_less_than_X_changes:
							break
				if num_changes < ignore_sites_less_than_X_changes:
					continue
			elif ignore_sites_less_than_X_variants:
				if len(set(column_1)) < ignore_sites_less_than_X_variants:
					continue
			used_columns_1.append(column_1_index)
			MI_values.append([])
			for column_2_index, column_2 in enumerate(self.seq_2_array.T):
				#check if a variable site
				if ignore_sites_less_than_X_changes or ignore_sites_less_than_X_variants:
					if column_2_index in msa_2_sites_to_skip:
						continue
				if first_round:
					used_columns_2.append(column_2_index)
				MI = self.calc_MI(column_1, column_2)
				MI_values[-1].append(MI)
			first_round = False
		MI_values = numpy.array(MI_values)
		if len(MI_values.flatten()) == 0:
			if ignore_sites_less_than_X_variants:
				threshold = ignore_sites_less_than_X_variants
			elif ignore_sites_less_than_X_changes:
				threshold = ignore_sites_less_than_X_changes
			print 'After removing sites with low variation (less than %s), there is no more data. Aborting' % threshold
			return None, None, None
		if output_array_filepath:
			fileout = open(output_array_filepath, 'w')
			fileout.write('\t%s\n' % '\t'.join([str(i+1) for i in used_columns_2]))
			for index, i in enumerate(MI_values):
				fileout.write('%s\t%s\n' % (used_columns_1[index]+1, '\t'.join([str(j) for j in i])))
			fileout.close()
		MI_values_list = []
		sites_1 = []
		sites_2 = []
		if output_distribution_filepath:
			fileout = open(output_distribution_filepath, 'w')
			fileout.write('MSA_1_position\tMSA_2_position\tmutual_information_values\n')
		for index_1, i in enumerate(MI_values):
			for index_2, j in enumerate(i):
				sites_1.append(used_columns_1[index_1]+1)
				sites_2.append(used_columns_2[index_2]+1)
				MI_values_list.append(j)
				if output_distribution_filepath:
					fileout.write('%s\t%s\t%s\n' % (used_columns_1[index_1]+1, used_columns_2[index_2]+1, j))
		if output_distribution_filepath:
			fileout.close()
		return sites_1, sites_2, MI_values_list

	def permute_MSA_rows(self):
		"""
		This method will permute the rows of both of the self.seq_1_array and self.seq_2_array objects. These are numpy arrays. None of the other attributes are affected, which means that the order of the sequence entries in self.msa_[1|2].data object will no longer correspond to the order of rows in the self.seq_[1|2]_array objects. Does one permutation for each MSA. The data in self.seq_[1|2]_array will be permanently changed.
		"""
		numpy.random.shuffle(self.seq_1_array)
		numpy.random.shuffle(self.seq_2_array)
		return

	def permute_rows_and_calc_MI(self, output_distribution_filepath, num_permute=10, ignore_sites_less_than_X_variants=None, ignore_sites_less_than_X_changes=None):
		"""
		This method will permute the rows of both of the seq arrays and then calculate the pairwise mutual information values between each of the sites after each round of permutations.
		output_distribution_filepath - This gives the path to the file that will contain all the MI values for each round of permutation.
		num_permute - Int. The number of rounds of permutations. Default, 10.
		ignore_sites_less_than_X_variants - Nonetype or Int. If defined (default, None), this will cause the method to ignore any sites (i.e. columns) that have a total number of variants that is less than this value.
		ignore_sites_less_than_X_changes - Nonetype or Int. If defined (default, None), this will cause the method to ignore any sites (i.e. columns) that have a total number of changes that is less than this value. A change is defined as a site that id different than the same site in the row directly preceding (above) it. This superceeds the 'ignore_sites_less_than_X_variants' parameter.
		"""
		fileout = open(output_distribution_filepath, "w")
		fileout.write('MSA_1_position\tMSA_2_position\tmutual_information_values\n')
		all_MI_values = []
		for permutation in xrange(num_permute):
			self.permute_MSA_rows()
			sites_1, sites_2, MI_values = self.calc_pairwise_mutual_information(output_array_filepath=None, output_distribution_filepath=None, ignore_sites_less_than_X_variants=ignore_sites_less_than_X_variants, ignore_sites_less_than_X_changes=ignore_sites_less_than_X_changes)
			all_MI_values += MI_values
			for index in xrange(len(MI_values)):
				fileout.write('%s\t%s\t%s\n' % (sites_1[index], sites_2[index], MI_values[index]))
		fileout.close()
		return all_MI_values

# msas = compare_MSAs(filepath_1='/Users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/abr/cdr3_time_series_alignments_of_representative_seqs_of_lineages/max_edit_dist_within_samps_6_across_samps_35/1/IGHV4-34_IGHJ6/554.fasta', filepath_2='/Users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/hiv/time_series_alignments_of_representative_seqs_of_wholePopulation/edit_dist_between_rep_seqs/1.fasta', count_attribute_name='DUPCOUNT')
# msas.calc_pairwise_mutual_information(output_array_filepath='/Users/nstrauli/Desktop/test_mi.txt', output_distribution_filepath='/Users/nstrauli/Desktop/test_mi_dstrb.txt', ignore_non_variable_sites=True)
# msas.permute_rows_and_calc_MI(output_distribution_filepath='/Users/nstrauli/Desktop/test_mi_null_dstrb.txt', num_permute=10, ignore_non_variable_sites=True)
