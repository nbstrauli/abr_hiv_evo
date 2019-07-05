import os
from sequence_sample_class import sequence_sample

class sequence_clusters(sequence_sample):
	"""
	This class deals with the manipulation of data types that consist of clusters of sequences within a sequence sample. The input data should be in fasta format, except formatted in a special way so that each fasta sequence entry is a representative sequence for a single unique cluster of seqs, and the informatin for each of the seqs within that cluster is contained in the header line. This class inherites all the methods and attributes of the "sequence_smaple" class.

	Attributes:
	filepath - path to the fasta or fastq formatted file. If fastq, however, the 3rd and 4th lines are ignored for each sequence entry.
	data - a list of dics, where each element has all the information for a given read (except the quality scores). Each read's info is the unique seq ID, the count of the sequence in the data, and the genetic sequence itself. It also contains any other information that is stored in the header. It stores this information as a dictionary called 'other' that is itself an element of the sequence's entry. This dictionary named 'other' has keys that are named after whatever other names are in the header line, and the keys values are whatever the value is listed as for that key. For example if there is "...|mut_freq=.012|..." in the header, than this information would be stored as 'other'={'mut_freq':'.012'}. Make sure that all headers in the sample have the same set of attributes listed in them.
		For the sequence_cluster_class, this attribute has some added features. In data[some_index]['other'] there are some sequence attributes that are going to be in the data. These are:
		'indiv_seqs' - This is a list of the individual seqs that make up the cluster.
		'indiv_seq_ids' - This is a list of the individual IDs for each of the seqs that make up the cluster.
		'indiv_seq_counts' - This is a list of the individual counts for each of the seqs that make up the cluster.
		'indiv_seq_freqs' - This is a list of the individual frequencies for each of the seqs that make up the cluster.
	counts - a list of counts for each of the sequences.
	total - total number of seqs represented in the file (not just unique seqs).
	diversity - diversity calculated as the mean pairwise genetic distance in the seqs in the fastq file. Must run the calc_diversity method before this is has a value (other than None type)
	downsamp_indicator - Boolean. If True, then it means the data has been down sampled. If False, then no downsampling has occured.
	count_attribute_name - This is the name (as a string) of the attribute in the header that gives the count information for all the reads.
	filetype - This is either 'fasta' or 'fastq'. It is determined by the suffix of the input file name.
	freq_attribute_name - This gives the name of the seq attribute that gives the relative frequency of each of the clusters in the dataset
    indiv_seqs_attribute_name - This gives the name of the seq attribute that lists individual sequences of the seqs for each cluster in the dataset.
    indiv_seq_ids_attribute_name - This gives the name of the seq attribute that lists individual sequence IDs of the seqs for each cluster in the dataset
    indiv_seq_counts_attribute_name - This gives the name of the seq attribute that lists individual sequence counts of the seqs for each cluster in the dataset
    indiv_seq_freqs_attribute_name - This gives the name of the seq attribute that lists individual sequence frequencies of the seqs for each cluster in the dataset
    min_freq - None (default), or Float. This gives the minimum frequency allowed for a sequence cluster to be included in the data. If None, then all clusters are included. This parameter is only considered if 'freq_attribute_name' is defined.
	"""

	def __init__(self, filepath, count_attribute_name=None, freq_attribute_name=None, indiv_seqs_attribute_name='indiv_seqs', indiv_seq_ids_attribute_name='indiv_seq_ids', indiv_seq_counts_attribute_name='indiv_seq_counts', indiv_seq_freqs_attribute_name='indiv_seq_freqs', min_freq=None):
		self.filepath = filepath
		suffix = os.path.basename(self.filepath).split('.')[-1]
		if suffix == 'fasta' or suffix == 'fa':
			self.filetype = 'fasta'
		elif suffix == 'fastq' or suffix == 'fq':
			self.filetype = 'fastq'
		else:
			print 'Unrecognized file type. Must be either fasta of fastq. Check file suffix.'
			return
		if not count_attribute_name and freq_attribute_name:
			print "'freq_attribute_name' cannont be defined if 'count_attribute_name' is not defined. Aborting."
			return
		self.count_attribute_name = count_attribute_name
		self.freq_attribute_name = freq_attribute_name
		self.indiv_seqs_attribute_name = indiv_seqs_attribute_name
		self.indiv_seq_ids_attribute_name = indiv_seq_ids_attribute_name
		self.indiv_seq_counts_attribute_name = indiv_seq_counts_attribute_name
		self.indiv_seq_freqs_attribute_name = indiv_seq_freqs_attribute_name
		self.min_freq = min_freq
		self.data = []
		self.counts = []
		self.pi = None
		self.downsamp_indicator = False
		filein = open(self.filepath, "r")
		count = 0
		missing_counts = False
		missing_freqs = False
		for i in filein:
			count += 1
			#get info in the header
			if count == 1:
				seq_entry_dic = {'id':'NA', 'count':1, 'seq':'NA', 'other':{}}
				line = i[1:-1].split('|')
				#the id is assumed to be the first entry in the header
				id = line[0]
				seq_entry_dic['id'] = id
				#get seq attribute info from header
				found_count = False
				found_freq = False
				freq_too_low = False
				for j in line[1:]:
					attribute = j.split('=')
					name = attribute[0]
					value = attribute[1]
					if name == self.count_attribute_name:
						seq_entry_dic['count'] = int(value)
						self.counts.append(int(value))
						found_count = True
					elif name == self.freq_attribute_name:
						freq = float(value)
						seq_entry_dic['other'][self.freq_attribute_name] = freq
						found_freq = True
						if min_freq:
							if freq < min_freq:
								freq_too_low = True
								break
					elif name == self.indiv_seqs_attribute_name:
						value = value.split(',')
						seq_entry_dic['other'][self.indiv_seqs_attribute_name] = value
					elif name == self.indiv_seq_ids_attribute_name:
						value = value.split(',')
						seq_entry_dic['other'][self.indiv_seq_ids_attribute_name] = value
					elif name == self.indiv_seq_counts_attribute_name:
						value = [int(k) for k in value.split(',')]
						seq_entry_dic['other'][self.indiv_seq_counts_attribute_name] = value
					elif name == self.indiv_seq_freqs_attribute_name:
						value = [float(k) for k in value.split(',')]
						seq_entry_dic['other'][self.indiv_seq_freqs_attribute_name] = value
					else:
				  		seq_entry_dic['other'][name] = value
				if not found_count and self.count_attribute_name:
					missing_counts = True
				if not found_freq and self.freq_attribute_name:
					missing_freqs = True
				if not self.count_attribute_name:
					seq_entry_dic['count'] = len(seq_entry_dic['other'][self.indiv_seqs_attribute_name])
					self.counts.append(seq_entry_dic['count'])
					seq_entry_dic['other'][self.indiv_seq_counts_attribute_name] = [1 for j in seq_entry_dic['other'][self.indiv_seqs_attribute_name]]
				if not self.freq_attribute_name:
					seq_entry_dic['other']['total_freq'] = float(seq_entry_dic['count'])
					seq_entry_dic['other']['indiv_seq_freqs'] = [float(j) for j in seq_entry_dic['other'][self.indiv_seq_counts_attribute_name]]
			elif count == 2:
				if self.filetype == 'fasta':
					count = 0
				if freq_too_low:
					if found_count:
						del self.counts[-1]
					continue
				seq = i[:-1]
				seq_entry_dic['seq'] = seq
				self.data.append(seq_entry_dic)
			elif count == 4:
				count = 0
		filein.close()
		if missing_counts:
			print "At least one sequence entry exists in the sample that does not have a count attribute, even though a count attribute was specified (i.e. not a None type). You might want to check this."
		if missing_freqs:
			print "At least one sequence entry exists in the sample that does not have a frequency attribute, even though a freq attribute was specified (i.e. not a None type). You might want to check this."
		self.total = float(sum(self.counts))
		if not self.freq_attribute_name:#normalize freqs by total if haven't already been calculated and stored in the file
			self.freq_attribute_name = 'total_freq'
			self.indiv_seq_freqs_attribute_name = 'indiv_seq_freqs'
			for index, i in enumerate(self.data):
				self.data[index]['total_freq'] = i['other']['total_freq'] / self.total
				self.data[index]['other']['indiv_seq_freqs'] = [j/self.total for j in i['other']['indiv_seq_freqs']]
		return

#seq_clusts = sequence_clusters(filepath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/sequence_clusters/hiv/clustered_by_edit_distance/most_numerous_seq_foreach_cluster/for_outlier_detection/max_edit_dist_4/1/3207.fasta', count_attribute_name='DUPCOUNT', freq_attribute_name='total_freq')
#print seq_clusts.counts
