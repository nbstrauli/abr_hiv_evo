import os
import sys
from sequence_sample_class import sequence_sample
from subprocess import Popen, PIPE
import subprocess
import random
import re
import itertools
from Bio import Seq
import re

class msa(sequence_sample):
	"""
	This class deals with analyzing and manipulating multiple sequence alignments (MSA's) in fasta format. Because the format is in fasta, this class inherits from the 'sequence_sample' class.

	Input:
	Should provide the filepath for the fasta formatted alignment. Also the attribute name that gives the count information for each of the entries

	Attributes:
	filepath - path to the fasta or fastq formatted file. If fastq, however, the 3rd and 4th lines are ignored for each sequence entry.
	data - a list of dics, where each element has all the information for a given read (except the quality scores). Each read's info is the unique seq ID, the count of the sequence in the data, and the genetic sequence itself. It also contains any other information that is stored in the header. It stores this information as a dictionary called 'other' that is itself an element of the sequence's entry. This dictionary named 'other' has keys that are named after whatever other names are in the header line, and the keys values are whatever the value is listed as for that key. For example if there is "...|mut_freq=.012|..." in the header, than this information would be stored as 'other'={'mut_freq':'.012'}. Make sure that all headers in the sample have the same set of attributes listed in them.
	counts - a list of counts for each of the sequences.
	total - total number of seqs represented in the file (not just unique seqs).
	diversity - diversity calculated as the mean pairwise genetic distance in the seqs in the fastq file. Must run the calc_diversity method before this is has a value (other than None type)
	downsamp_indicator - Boolean. If True, then it means the data has been down sampled. If False, then no downsampling has occured.
	count_attribute_name - This is the name (as a string) of the attribute in the header that gives the count information for all the reads.
	filetype - This is either 'fasta' or 'fastq'. It is determined by the suffix of the input file name. For the msa class, this will only be 'fasta'.
	variants_dic - this is a dic that is not defined until the method 'identitfy_variants' is run. Equals None type until then. The dic is indexed by variant positions, and the def of each index gives the identity and frequency of each of the alleles at that site.
	min_freq_threshold - This is undefined until the method 'identitfy_variants' is run. It gives the minimum frequency that an allele can have inorder to be included in the variants_dic.
	ref_seq - If the 'first_seq_ref' parameter is True, then this gives the reference sequence in the alignment. This will be None otherwise
	ref_seq_name - If the 'first_seq_ref' parameter is True, then this gives the name of the reference sequence in the alignment. This will be None otherwise.
	ref_seq_coding_frame - If the 'first_seq_ref' parameter is True, then this should be an integer, contained in [1,3]. It gives where one should start in the ref seq (indexed at 1) when translating to amino acids. This will be None if 'first_seq_ref' is not defined.
	"""

	def __init__(self, filepath, count_attribute_name=None, first_seq_ref=False):
		"""
		first_seq_ref - If this is True (default, False), then the first sequence in the fasta alignment file will be treated as a reference sequence (which is formatted differently). The header of this reference sequence must look like this: ">[name]|blah=blah_blah|coding_frame_start=[start]", where 'name' = the name of the reference sequence, and 'start' gives the coding frame of the sequence. This can be 1, 2, or 3. If it is 1, for example, then that means that the amino acids are translated starting at the very first nucleotide in the reference sequence.
		"""
		self.filepath = filepath
		self.variants_dic = None
		self.min_freq_threshold = None
		suffix = os.path.basename(self.filepath).split('.')[-1]
		if suffix == 'fasta' or suffix == 'fa':
			self.filetype = 'fasta'
		elif suffix == 'fastq' or suffix == 'fq':
			self.filetype = 'fastq'
		else:
			print 'Unrecognized file type. Must be either fasta of fastq. Check file suffix.'
			return
		self.count_attribute_name = count_attribute_name
		self.data = []
		self.counts = []
		self.pi = None
		self.downsamp_indicator = False
		filein = open(self.filepath, "r")
		count = 0
		missing_counts = False

		#This deals with getting the reference seq info. It is unique to the MSA class.
		if first_seq_ref:
			ref_header = filein.readline()[:-1].split('|')
			self.ref_seq_name = ref_header[0]
			for i in ref_header[1:]:
				attribute = i.split('=')
				name = attribute[0]
				value = attribute[1]
				if name == 'coding_frame_start':
					self.ref_seq_coding_frame = int(value)
			self.ref_seq = filein.readline()[:-1]
		else:
			self.ref_seq_name = None
			self.ref_seq_coding_frame = None
			self.ref_seq = None

		for i in filein:
			count += 1
				#get info in the header
			if count == 1:
				seq_entry_dic = {'id':'NA', 'count':1, 'seq':'NA', 'other':{}}
				line = i[1:-1].split('|')
				#the id is assumed to be the first entry in the header
				id = line[0]
				seq_entry_dic['id'] = id
				if len(id.split('_')) > 1:
					seq_entry_dic['other']['timepoint'] = id.split('_')[0]
				#look for the count attribute in the header
				if self.count_attribute_name != None:
					found_count = False
					for j in line[1:]:
						attribute = j.split('=')
						name = attribute[0]
						value = attribute[1]
						if name == self.count_attribute_name:
							seq_entry_dic['count'] = int(value)
							self.counts.append(int(value))
							found_count = True
						else:
							seq_entry_dic['other'][name] = value
					if found_count == False:
						missing_counts = True
						self.counts.append(1)
				else:
					self.counts.append(1)
			elif count == 2:
				seq = i[:-1]
				seq_entry_dic['seq'] = seq
				self.data.append(seq_entry_dic)
				if self.filetype == 'fasta':
					count = 0
			elif count == 4:
				count = 0
		filein.close()
		if missing_counts:
			print "At least one sequence entry exists in the sample that does not have a count attribute, even though a count attribute was specified (i.e. not a None type). You might want to check this."
		self.total = float(sum(self.counts))
		return

	def identitfy_variants(self, min_freq_threshold=0.0, time_series_data=False):
	 	"""
	 	This method goes through each column in the MSA and finds the positions that have significant variation. Significant means that the frequency is above 'min_freq_threshold'.
		min_freq_threshold - This gives the minimum value for a variant to be recorded. The default is 0, so that all variants, however rare, will be recorded.
		time_series_data - If True (default, False), then this tells the method that the data in the MSA is a time-series. This means that the data from across a time-series has been compiled and aligned all together. If this is the case then the time-point information for each of the reads needs to be encoded in a specific way. At the beggining of each seq header there needs to be: ">[X]_[seq_id]|etc=etc...", where X = the time-point (as a float), and seq_id = the sequence id.
	 	"""
	 	self.min_freq_threshold = min_freq_threshold
	 	#get time-point info
	 	tpoints = set()
	 	tpoint_totals = {}
	 	if time_series_data:
	 		for i in self.data:
	 			tpoint = float(i['id'].split('_')[0])
	 			tpoints.update([tpoint])
	 			try:
	 				tpoint_totals[tpoint] += i['count']
	 			except KeyError:
	 				tpoint_totals[tpoint] = float(i['count'])
	 		first_tpoint = min(tpoints)
	 	else:
	 		tpoints = set([1])
	 		tpoint_totals[1] = self.total
	 		first_tpoint = 1
	 	#this will be a 3 dimensional variable. 1st dim: a list, each element corresponds to a column in the alignment. 2nd dim: a dic, indexed by time-point, so each element corresponds to a time-point. 3rd dim: a dic, indexed by allele, so each element corresponds to each allele found at that position, at that time-point. Definition is the freq for that allele, at that time-point, in that position.
	 	variants = []
	 	num_cols = len(self.data[0]['seq'])
		#initialize the variants variable
		for i in xrange(num_cols):
			variants.append({})
			for j in tpoints:
				variants[-1][j] = {}
		#for each seq
		for i in self.data:
			tpoint = float(i['id'].split('_')[0])
			#for each column
			for j in xrange(num_cols):
				try:
					variants[j][tpoint][i['seq'][j]] += i['count'] / tpoint_totals[tpoint]
				except KeyError:
					variants[j][tpoint][i['seq'][j]] = i['count'] / tpoint_totals[tpoint]
		#now filter out the positions that have a minor allele freq less than 'min_freq_threshold'
		variants_dic = {}
		#cycle through each column
		for i in xrange(num_cols):
			found_variant = False
			#get identity of major allele
			starting_alleles = []
			for j in variants[i][first_tpoint]:
				starting_alleles.append([variants[i][first_tpoint], j])
			major_allele = sorted(starting_alleles)[-1][1]
			#for each time-point
			for j in tpoints:
				#look for minor alleles that rise above the min freq threshold
				for k in variants[i][j]:
					if variants[i][j][k] >= min_freq_threshold and k != major_allele:
						variants_dic[i] = variants[i]
						found_variant = True
						break
				if found_variant:
					break
		self.variants_dic = variants_dic
		return variants_dic

	def get_haplotypes(self, min_freq_threshold=0.0, time_series_data=False, freq_traj_output_filepath=None, fasta_output_filepath=None, write_seq_id_first=None):
		"""
		This method gets the haplotypes for each of the seqs in the data. The variants that it includes in each haplotype are determined by the method 'identitfy_variants'. If this method has not been run then it is run first.
		min_freq_threshold - This gives the minimum frequency threshold for including and allele in the haplotypes. This variable is only considered if 'identitfy_variants' has not been run yet. If it has than this is ignored.
		time_series_data - If True (default, False), then this tells the method that the data in the MSA is a time-series. This means that the data from across a time-series has been compiled and aligned all together. If this is the case then the time-point information for each of the reads needs to be encoded in a specific way. At the beggining of each seq header there needs to be: ">[X]_[seq_id]|etc=etc...", where X = the time-point (as a float), and seq_id = the sequence id.
		freq_traj_output_filepath - If defined (defined, None), then this will i) instruct the method to write the results to disk, and ii) give the filepath for which the output will be written.
		fasta_output_filepath - If defined (default, None), then this is the path to the fasta formatted file that will have the haplotypes found in the data, with some info about them, including their count, and what seq IDs have a given unique haplotype.
		write_seq_id_first - If defined (default, None), this will make the script write the haplotype that contains the provided sequence ID as the first entry in the fasta file. This should be a string that is identical to the sequence ID for which the haplotype that contains this seq ID will be written first. This is only considered if 'fasta_output_filepath' is defined.
		"""
		if self.min_freq_threshold == None:
			self.identitfy_variants(min_freq_threshold=min_freq_threshold, time_series_data=time_series_data)
		#if time-series data, get time-points
		tpoints = set()
		if time_series_data:
			for i in self.data:
				tpoint = float(i['id'].split('_')[0])
				tpoints.update([tpoint])
			tpoint_to_index_dic = {}
			for index, tpoint in enumerate(sorted(tpoints)):
				tpoint_to_index_dic[tpoint] = index
		#else, then treat as one time-point
		else:
			tpoints.update([1])
			tpoint_to_index_dic = {1:0}
		variant_positions = [i for i in sorted(self.variants_dic)]
		haplotype_freq_dic = {}
		tpoint_totals = [0. for i in tpoints]
		haplotype_seq_info_dic = {}
		#get haplotypes for each seq
		for i in self.data:
			haplotype = ''
			for j in variant_positions:
				haplotype += i['seq'][j]
			if time_series_data:
				tpoint = float(i['id'].split('_')[0])
			else:
				tpoint = 1
			tpoint_totals[tpoint_to_index_dic[tpoint]] += i['count']
			try:
				haplotype_freq_dic[haplotype][tpoint_to_index_dic[tpoint]] += i['count']
			except KeyError:
				haplotype_freq_dic[haplotype] = [0 for j in tpoints]
				haplotype_freq_dic[haplotype][tpoint_to_index_dic[tpoint]] = i['count']
			try:
				haplotype_seq_info_dic['%s_%s' % (haplotype, tpoint)]['ids'].append(i['id'])
				haplotype_seq_info_dic['%s_%s' % (haplotype, tpoint)]['freqs'].append(i['count'])
				haplotype_seq_info_dic['%s_%s' % (haplotype, tpoint)]['counts'].append(i['count'])
			except KeyError:
				haplotype_seq_info_dic['%s_%s' % (haplotype, tpoint)] = {'ids':[i['id']], 'freqs':[i['count']], 'counts':[i['count']]}
			if write_seq_id_first:
				if i['id'] == write_seq_id_first:
					first_haplotype_key = '%s_%s' % (haplotype, tpoint)
			
		#now normalize each entry by the total counts for that time-point
		for i in haplotype_freq_dic:
			for index in xrange(len(tpoints)):
				haplotype_freq_dic[i][index] /= tpoint_totals[index]
		for i in haplotype_seq_info_dic:
			tpoint = float(i.split('_')[-1])
			tpoint_total = tpoint_totals[tpoint_to_index_dic[tpoint]]
			for j in xrange(len(haplotype_seq_info_dic[i]['freqs'])):
				haplotype_seq_info_dic[i]['freqs'][j] /= tpoint_total
			haplotype_seq_info_dic[i]['total_count'] = sum(haplotype_seq_info_dic[i]['counts'])
			haplotype_seq_info_dic[i]['total_freq'] = sum(haplotype_seq_info_dic[i]['freqs'])

		if freq_traj_output_filepath:
			fileout = open(freq_traj_output_filepath, "w")
			header = '\t%s\n' % '\t'.join([str(i) for i in sorted(tpoints)])
			fileout.write(header)
			haplotype_freq_list = [[sum(haplotype_freq_dic[i]), i, haplotype_freq_dic[i]] for i in haplotype_freq_dic]
			for i in sorted(haplotype_freq_list):
				fileout.write('%s\t%s\n' % (i[1], '\t'.join([str(j) for j in i[2]])))
			fileout.close()

		if fasta_output_filepath:
			fileout = open(fasta_output_filepath, "w")
			if write_seq_id_first:
				tpoint = first_haplotype_key.split('_')[1]
				header = '>%s_%s|seq_ids=%s|counts=%s|freqs=%s|total_count=%s|total_freq=%s\n' % (tpoint, '0', ','.join(haplotype_seq_info_dic[first_haplotype_key]['ids']), ','.join([str(j) for j in haplotype_seq_info_dic[first_haplotype_key]['counts']]), ','.join([str(j) for j in haplotype_seq_info_dic[first_haplotype_key]['freqs']]), haplotype_seq_info_dic[first_haplotype_key]['total_count'], haplotype_seq_info_dic[first_haplotype_key]['total_freq'])
				fileout.write('%s%s\n' % (header, first_haplotype_key.split('_')[0]))
			index = 0
			for i in haplotype_seq_info_dic:
				if write_seq_id_first:
					if i == first_haplotype_key:
						continue
				index += 1
				tpoint = i.split('_')[1]
				header = '>%s_%s|seq_ids=%s|counts=%s|freqs=%s|total_count=%s|total_freq=%s\n' % (tpoint, index, ','.join(haplotype_seq_info_dic[i]['ids']), ','.join([str(j) for j in haplotype_seq_info_dic[i]['counts']]), ','.join([str(j) for j in haplotype_seq_info_dic[i]['freqs']]), haplotype_seq_info_dic[i]['total_count'], haplotype_seq_info_dic[i]['total_freq'])
				fileout.write('%s%s\n' % (header, i.split('_')[0]))
			fileout.close()
		return haplotype_freq_dic, haplotype_seq_info_dic

	def make_phylo_trees(self, output_filepath, method='fasttree', path_to_fasttree=None, nuc_or_aa='nucleotide', temp_dirpath=None, include_freq_info=False, freq_attribute_name=None, first_seq_is_outgroup=False, outgroup_label='outgroup', include_cluster_info=False, cluster_attribute_name=None, include_attributes=None):
		"""
		This method uses the MSA as input to make a phylogenetic tree.
		output_filepath - This gives the path to the output tree file.
		method - This tells what program to use to generate the tree. Acceptable values are:
			'fasttree' - This mean that the fastTree excecutable will be used to generate the tree.
		path_to_fasttree - This gives the path to the FastTree excecutable. If this is None (default), then it is assumed that FastTree is in the PATH.
		nuc_or_aa - This tells whether the MSA is made up of nucleotides or amino acids. Acceptable values are:
			'nucleotide' - made up of nucleotides
			'amino_acid' - made up of amino acids
		include_freq_info - If True (default, False), this will instruct the script to include the relative frequency information of each sequence in the data with the names of the seqs in the output newick file. The parameter 'freq_attribute_name' needs to be defined if this is True.
		freq_attribute_name - If defined (default, None), this gives the name (in the header or each seq) of the attribute that gives the relative frequency information. This is only considered if 'include_freq_info' is defined.
		first_seq_is_outgroup - If True (default, False), then the signifies that the first sequence in the MSA should be the outgroup for the phylogeny. NOTE, this will permanently change the ID of the fist seq to the value of 'outgroup_label'.
		outgroup_label - This should be a string (default is 'outgroup') that gives the identifier that will be in the newick tree formatted file in the node that signifies that the node is the outgroup. This label will be in place of the sequence ID. If there is time-point information, then this is retained in the outgroup node, but the 'ID' portion of the name will be replaced with this value.
		include_cluster_info - If True (default, False), this will include the cluster ID information in the leaf names in the output newick file.
		cluster_attribute_name - String. This gives the name of the attribute for each of the seqs that provides the cluster ID information. Ex: 'cluster_id'.
		include_attributes - If defined (default, None), this is a list of strings, where each string gives the name of an attribute in self.data[some_index]['other'] to be included in the name of the leaf node in the newick file. Make sure this is not too long because fasttree does not like long seq headers.
		"""
		if not temp_dirpath:
			temp_dirpath = os.path.dirname(output_filepath) + '/'
		elif temp_dirpath[-1] != '/':
			temp_dirpath += '/'
		if not os.path.exists(temp_dirpath):
			os.makedirs(temp_dirpath)
		random_suffix = str(random.random())
		if include_freq_info and not freq_attribute_name:
			print "'include_freq_info' parameter is defined but no attribute name is given for 'freq_attribute_name'. Aborting."
			return

		if first_seq_is_outgroup:
			if self.data[0]['id'].split('_') == 1:
				self.data[0]['id'] = outgroup_label
			else:
				self.data[0]['id'] = "%s_%s" % (self.data[0]['id'].split('_')[0], outgroup_label)

		if method == 'fasttree':
			if not path_to_fasttree:
				path_to_fasttree = 'FastTree'
			#fasttree seems to not like long headers, so need to shorten them to only necessary info: seq ID, count, and freq
			temp_fasta_input = temp_dirpath + 'temp_fasta_input_%s_%s.fasta' % (os.path.basename(self.filepath)[:-6], random_suffix)
			fileout = open(temp_fasta_input, "w")
			for i in self.data:
				header = '>%s|%s=%s' % (i['id'], self.count_attribute_name, i['count'])
				if include_freq_info:
					header += '|%s=%s' % (freq_attribute_name, i['other'][freq_attribute_name])
				if include_cluster_info:
					header += '|%s=%s' % (cluster_attribute_name, i['other'][cluster_attribute_name])
				if include_attributes:
					for j in include_attributes:
						header += '|%s=%s' % (j, i['other'][j])

				fileout.write('%s\n%s\n' % (header, i['seq']))
			fileout.close()
			p = Popen(['bash', 'call_fasttree.bash', path_to_fasttree, temp_fasta_input, output_filepath, nuc_or_aa], stdout=PIPE, stderr=PIPE)
			out, err = p.communicate()
			print out
			print err
			subprocess.call(['rm', temp_fasta_input])

		return

	def get_time_point_value_by_index(self, time_point_index):
		"""
		This method will return the value of a given time-point when provided the index (or rank). For example, if one wants to know the value of the 3rd time-point in the MSA, then they would make time_point_index=2, and it would return the value of this time-point (for example, if days).
		"""
		tpoints = set()
		for i in self.data:
			seq_id_info = i['id'].split('_')
			if len(seq_id_info) < 2:
				print 'There does not appear to be any time-point information if the sequence IDs. Aborting', i['id']
			tpoint = float(seq_id_info[0])
			tpoints.update([tpoint])
		return sorted(tpoints)[time_point_index]

	def get_consensus_seq(self, only_for_timepoint=False, output_filepath=None, method='cons', temp_dirpath=None, translate_by_ref_seq=False):
		"""
		This method will get the consensus seq for the sequences (or a subset of seqs) in the alignment.
		only_for_timepoint - This will tell which of the time-points the sequences that will be used to make the consensus seq will be restricted to. This should be a float. So, if it is 3.0 then only the seq from time-point 3.0 will be used to get the consensus seq. If this is False (default), then all sequences in the msa are used.
		output_filepath - This gives the path to the file for which the consensus sequence will be written
		method - This tells which tool will be used to get the consensus seq. Acceptable values are:
			'cons' - (default) This means that the 'cons' program from the EMBOSS package will be used. cons needs to be in the PATH
			'most_abundant' - This means that the most abundant seq in the alignment (or subset of the alignment) will be used as the 'consensus'. If there is a tie, then the seq that comes first (closer to the top) will arbitrarily be chosen.
		temp_dirpath - Gives the path to a directory for which temp data can be written. If None (default), then uses the master directory to 'output_filepath', if there is not output_filepath given then uses the current working directory.
		translate_by_ref_seq - If defined (default is False) this should be the path to a fasta formatted reference sequence for the seqs in the alignment. This reference seq needs to be annotated such that the coding frame is known. Meaning, it is known where to begin translating the ref seq. The annotation for the coding frame information must be in the header of the ref seq and should look like:'...|coding_frame_start=[X]|...', where 'X' gives the coding frame start position. The consensus seq is then aligned to this ref seq so that the coding frame can be mapped from the ref seq to the consensus seq. The alignment is done using 'needle' from the EMBOSS package. This needs to be in the PATH. The resulting translation is recorded in the 'output_filepath', so this should be defined too. If 'translate_by_ref_seq' equals False, then no translation occurs.

		Output:
		Returns the consensus seq as a string.
		"""
		if temp_dirpath:
			if temp_dirpath[-1] != '/':
				temp_dirpath += '/'
			if not os.path.exists(temp_dirpath):
				os.makedirs(temp_dirpath)
		else:
			if output_filepath:
				temp_dirpath = os.path.dirname(output_filepath) + '/'
			else:
				temp_dirpath = os.getcwd() + "/"
		if not self.count_attribute_name and method == 'most_abundant':
			print 'There is no reprted count information in the file, but the method for getting the consensus seq is "most_abundant". Aborting'
			return
		random_suffix = str(random.random())

		#if we are to restrict the sequence set to one of the time-points
		if only_for_timepoint:
			input_filepath = '%stemp_file_%s.fasta' % (temp_dirpath, random_suffix)
			fileout = open(input_filepath, "w")
			for i in self.data:
				tpoint = float(i['id'].split('_')[0])
				if tpoint == only_for_timepoint:
					other_headers = ['%s=%s' % (j, i['other'][j]) for j in i['other']]
					header = '>%s|%s=%s|%s\n' % (i['id'], self.count_attribute_name, i['count'], '|'.join(other_headers))
					fileout.write(header + i['seq'] + '\n')
			fileout.close()
		else:
			input_filepath = self.filepath()

		#if using 'cons' method
		if method == 'cons':
			output_filepath_temp = "%stemp_outfile_%s.fasta" % (temp_dirpath, random_suffix)
			p = Popen(['cons', '-sequence', input_filepath, '-outseq', output_filepath_temp], stdout=PIPE, stderr=PIPE)
			out, err = p.communicate()
			print out
			print err
			#delete temp input file if necessary
			if only_for_timepoint:
				subprocess.call(['rm', input_filepath])
			#read the consensus seq into memory
			filein = open(output_filepath_temp, "r")
			cons_seq = ''
			filein.readline()
			for i in filein:
				cons_seq += i[:-1]
			filein.close()
			subprocess.call(['rm', output_filepath_temp])

		#if usng the most abundant seq method
		elif method == 'most_abundant':
			filein = open(input_filepath, "r")
			counts = []
			seqs = []
			for i in filein:
				if i[0] == '>':
					header = i[:-1].split('|')
					for j in header:
						attribute = j.split('=')
						if attribute[0] == self.count_attribute_name:
							counts.append(int(attribute[1]))
							break
				else:
					seqs.append(re.sub('-', '', i[:-1]))
			filein.close()
			#this will find the most abundant seq
			cons_seq = seqs[counts.index(max(counts))]

		#if need to translate the consensus seq by a reference
		if translate_by_ref_seq:
			filein = open(translate_by_ref_seq, "r")
			header = filein.readline()[1:-1].split('|')
			for i in header:
				attribute = i.split('=')
				if attribute[0] == 'coding_frame_start':
					coding_frame_start = int(attribute[1])
			ref_seq = filein.readline()[:-1]
			filein.close()
			alignment_filepath = "%stemp_alignment_%s" % (temp_dirpath, random_suffix)
			p = Popen(['needle', '-asequence', 'asis:'+ref_seq, '-bsequence', 'asis:'+cons_seq, '-gapopen', '10.0', '-gapextend', '0.5', '-outfile', alignment_filepath, '-brief', 'Y'], stdout=PIPE, stderr=PIPE)
			out, err = p.communicate()
			print out
			print err
			filein = open(alignment_filepath, "r")
			ref_seq_alignment = ''
			cons_seq_alignment = ''
			first_seq = True
			for i in filein:
				if i[:4] == 'asis':
					seq = i[:-1].split()[2]
					if first_seq:
						ref_seq_alignment += seq
						first_seq = False
					else:
						cons_seq_alignment += seq
						first_seq = True
			filein.close()
			subprocess.call(['rm', alignment_filepath])
			ref_seq_pos = 0
			cons_seq_pos = 0
			for i, j in itertools.izip(ref_seq_alignment, cons_seq_alignment):
				if i != '-':
					ref_seq_pos += 1
				if j != '-':
					cons_seq_pos += 1
				if ref_seq_pos == coding_frame_start:
					cons_seq_coding_frame_start = cons_seq_pos
					break
			cons_aa_seq = Seq.translate(cons_seq[cons_seq_coding_frame_start-1:])
		else:
			cons_aa_seq = None

		if output_filepath:
			if cons_aa_seq:
				header = '>%s_consensus_seq|amino_acid_seq=%s|coding_frame_start=%s\n' % (os.path.basename(self.filepath)[:-6], cons_aa_seq, cons_seq_coding_frame_start)
			else:
				header = '>%s_consensus_seq\n' % os.path.basename(self.filepath)[:-6]
			fileout = open(output_filepath, "w")
			fileout.write("%s%s\n" % (header, cons_seq))
			fileout.close()

		return cons_seq

	def get_mean_dist_to_seq(self, seq, time_series_data=False, output_filepath=None, method='needle', temp_dirpath=None):
		"""
		This method calculates the mean distance to a given sequence for each seq (or subset of seqs) in the MSA.
		time_series_data - If True (default, False), then this tells the method that the data in the MSA is a time-series. This means that the data from across a time-series has been compiled and aligned all together. If this is the case then the time-point information for each of the reads needs to be encoded in a specific way. At the beggining of each seq header there needs to be: ">[X]_[seq_id]|etc=etc...", where X = the time-point (as a float), and seq_id = the sequence id. If time series data is True, then this method will calculate the mean genetic distance for each time-point individually.
		output_filepath - This gives the path to the file that the output will be written to. If this is None (default), no output is written, and is only returned via memory. If time_series_data=True, then one file will be written for each timepoint, in the same directory. This will then be treated as the base name for each of the output files. For example, if output_filepath equals 'blah/blah/divergence', then the output file for time-point 1046 will be 'blah/blah/divergence_1046'.
		method - This gives the method that will be used to make the pairwise alignments from which genetic distance will be calculated. Acceptable values are:
			'needle' - This means that the 'needle' program from the EMBOSS tool suite will be used. This program should be in the PATH.
		temp_dirpath - Gives the path to a directory for which temp data can be written. If None (default), then uses the master directory to 'output_filepath', if there is not output_filepath given then uses the current working directory.
		"""
		if temp_dirpath:
			if temp_dirpath[-1] != '/':
				temp_dirpath += '/'
			if not os.path.exists(temp_dirpath):
				os.makedirs(temp_dirpath)
		else:
			if output_filepath:
				temp_dirpath = os.path.dirname(output_filepath) + '/'
			else:
				temp_dirpath = os.getcwd() + "/"
		random_suffix = str(random.random())
		temp_aseq_filepath = '%s_temp_aseq_%s.fasta' % (temp_dirpath, random_suffix)
		fileout = open(temp_aseq_filepath, "w")
		fileout.write(">aseq\n%s\n" % seq)
		fileout.close()
		if time_series_data:
			#first get time-points and write to temp seq files
			temp_fileout_dic = {}
			for i in self.data:
				tpoint = float(i['id'].split('_')[0])
				header = '>count=%s_%s\n' % (i['count'], i['id'].split('_')[1])
				try:
					temp_fileout_dic[tpoint].write(header + i['seq'] + '\n')
				except KeyError:
					temp_output_filepath = '%stemp_%s_%s.fasta' % (temp_dirpath, tpoint, random_suffix)
					temp_fileout_dic[tpoint] = open(temp_output_filepath, "w")
					temp_fileout_dic[tpoint].write(header + i['seq'] + '\n')
			for i in temp_fileout_dic:
				temp_fileout_dic[i].close()
			#now calc dist between seqs in tpoints and query seq
			mean_percent_dist = {}
			for i in temp_fileout_dic:
				p = Popen(['needle', '-asequence', temp_aseq_filepath, '-bsequence', "%stemp_%s_%s.fasta" % (temp_dirpath, i, random_suffix), '-gapopen', '10.0', '-gapextend', '0.5', '-brief', 'Y', '-stdout', '-auto'], stdout=PIPE, stderr=PIPE)
				out, err = p.communicate()
				if output_filepath:
					fileout = open('%s_%s' % (output_filepath, i), "w")
					fileout.write('sequence_id\tsequence_count\tdivergence\n')
				#parse needle output and get seq counts too
				total_seqs = 0
				percent_dists = []
				for j in out.split('\n'):
					if j[:5] == '# 2: ':
						line = j.split('_')
						seq_id = line[-1]
						b_seq_count = int(line[-2].split('count=')[1])
						total_seqs += b_seq_count
					if j[:11] == '# Identity:':
						match = re.search('\(\s*([0-9]+\.[0-9])', j)
						percent_dist = 100 - float(match.group(1))
						percent_dist_count_scaled = percent_dist * b_seq_count
						percent_dists.append(percent_dist_count_scaled)
						if output_filepath:
							fileout.write('%s\t%s\t%s\n' % (seq_id, b_seq_count, percent_dist))
						#make sure we're not recycling seq_id and count values
						seq_id = None
						b_seq_count = None
				mean = sum(percent_dists) / total_seqs
				mean_percent_dist[i] = mean
				if output_filepath:
					fileout.write('mean\tNA\t%s\n' % mean)
					fileout.close()
				#delete temp bseq file
				subprocess.call(['rm', "%stemp_%s_%s.fasta" % (temp_dirpath, i, random_suffix)])
		else:
			#write temp sequence file. Need to do this so that the headers for each seq are formatted such that they can be parsed in the aligner output
			temp_output_filepath = '%stemp_%s.fasta' % (temp_dirpath, random_suffix)
			fileout = open(temp_output_filepath, "w")
			for i in self.data:
				tpoint = float(i['id'].split('_')[0])
				header = '>count=%s_%s\n' % (i['count'], i['id'].split('_')[1])
				fileout.write(header + i['seq'] + '\n')
			fileout.close()
			#now calc dist between seqs and query seq
			p = Popen(['needle', '-asequence', temp_aseq_filepath, '-bsequence', "%stemp_%s.fasta" % (temp_dirpath, random_suffix), '-gapopen', '10.0', '-gapextend', '0.5', '-brief', 'Y', '-stdout', '-auto'], stdout=PIPE, stderr=PIPE)
			out, err = p.communicate()
			if output_filepath:
				fileout = open(output_filepath, "w")
				fileout.write('sequence_id\tsequence_count\tdivergence\n')
			#parse needle output and get seq counts too
			total_seqs = 0
			percent_dists = []
			for j in out.split('\n'):
				if j[:5] == '# 2: ':
					line = j.split('_')
					seq_id = line[-1]
					b_seq_count = int(line[-2].split('count=')[1])
					total_seqs += b_seq_count
				if j[:11] == '# Identity:':
					match = re.search('\(\s*([0-9]+\.[0-9])', j)
					percent_dist = 100 - float(match.group(1))
					percent_dist_count_scaled = percent_dist * b_seq_count
					percent_dists.append(percent_dist_count_scaled)
					if output_filepath:
						fileout.write('%s\t%s\t%s\n' % (seq_id, b_seq_count, percent_dist))
					#make sure we're not recycling seq_id and count values
					seq_id = None
					b_seq_count = None
			mean_percent_dist = sum(percent_dists) / total_seqs
			if output_filepath:
				fileout.write('mean\tNA\t%s\n' % mean_percent_dist)
				fileout.close()
			#delete temp bseq file
			subprocess.call(['rm', "%stemp_%s.fasta" % (temp_dirpath, random_suffix)])
		#delete temp aseq file
		subprocess.call(['rm', temp_aseq_filepath])
		return mean_percent_dist

	def cluster_seqs_over_time(self, method='counting_diffs', diffs_cutoff=1, overwrite_data_object=False, cluster_representative_seq='most_abundant', freq_attribute_name=None, keep_first_seq_first=False):
		"""
		This method clusters the sequences in the multiple sequence alignment. Unique clusters are given sequenctial integer IDs, and the ID of the cluster that each sequence belongs to is stored in self.data[X]['other']['cluster_id'], where X is the index of the sequence.
		method - This gives the method of clustering that will be used. Acceptable values are:
			'counting_diffs' - This means that the number of differences between each pairwise sequence comparison are counted, and if they are below some threshold then they are deemed to be in the same cluster.
		diffs_cutoff - Int. This gives the cutoff for the max number of differences that are permissable for any two sequences to be in the same cluster.
		overwrite_data_object - If True (default, False), this will cause the information in self.data to be updated to reflect the information for each of the clusters, rather than the individual sequences. All other relevant attributes in the MSA object will be updated as well (such as counts).
		cluster_representative_seq - This gives the method for which the representative sequence for each cluster is determined. Acceptable values are:
			'most_abundant' - This means that the sequence with the highest frequency in the cluster is selected to represent the cluster.
		freq_attribute_name - If defined (default, None), this gives the name of the attribute that provides the frequency of each of the seqs in the data. If this is None type then frequency information (if present) is ignored in the data.
		keep_first_seq_first - If True (default, False), this will make sure that (after clustering) the cluster that owns the sequence that was first in the original data, will continue to be represented first in the data. This is because often times the first seq in a file has some significance, so this position is retained for the cluster that owns this seq. This paramter is only considered if 'overwrite_data_object' is True.
		"""
		new_cluster_id = 1
		list_of_indices = range(len(self.data))
		for seq1_index in list_of_indices[:-1]:

			#the algorithm below is flawed

			if 'cluster_id' in self.data[seq1_index]['other']:
				cluster_id = self.data[seq1_index]['other']['cluster_id']
			else:
				cluster_id = new_cluster_id
				new_cluster_id += 1
				self.data[seq1_index]['other']['cluster_id'] = cluster_id
			for seq2_index in list_of_indices[seq1_index+1:]:
				diffs = 0
				for i, j in itertools.izip(self.data[seq1_index]['seq'], self.data[seq2_index]['seq']):
					if i != j:
						diffs += 1
						if diffs > diffs_cutoff:
							break
				if diffs <= diffs_cutoff:
					self.data[seq2_index]['other']['cluster_id'] = cluster_id

		tpoint_freq_dic = {}
		for i in self.data:
			tpoint = i['id'].split('_')[0]
			try:
				tpoint_freq_dic[tpoint] += float(i['other']['total_freq'])
			except KeyError:
				tpoint_freq_dic[tpoint] = float(i['other']['total_freq'])
		for i in tpoint_freq_dic:
			print 'time-point:', i, '; freq sum:', tpoint_freq_dic[i]
		return

		#Now overwrite the data in memory, if desired.
		if overwrite_data_object:
			if keep_first_seq_first:
				if 'timepoint' in self.data[0]['other']:
					tpoint = self.data[0]['other']['timepoint']
					first_seq_key = '%s_%s' % (tpoint, self.data[0]['other']['cluster_id'])
				else:
					first_seq_key = self.data[0]['other']['cluster_id']
			#gather new data
			new_data = {}
			self.diversity = None#diversity statistic no longer applies (if calculated)
			for i in self.data:
				if 'timepoint' in i['other']:
					tpoint = i['id'].split('_')[0]
					cluster_seq_id = '%s_%s' % (tpoint, i['other']['cluster_id'])
				else:
					cluster_seq_id = i['other']['cluster_id']
				try:
					new_data[cluster_seq_id]['other']['indiv_seq_counts'].append(i['count'])
					new_data[cluster_seq_id]['other']['indiv_seqs'].append(i['seq'])
					new_data[cluster_seq_id]['count'] += i['count']
					if freq_attribute_name:
						new_data[cluster_seq_id]['other']['indiv_seq_freqs'].append(float(i['other'][freq_attribute_name]))
						new_data[cluster_seq_id]['other'][freq_attribute_name] += float(i['other'][freq_attribute_name])
				except KeyError:
					new_data[cluster_seq_id] = {'id':cluster_seq_id, 'count':i['count'], 'other':{'cluster_id':i['other']['cluster_id'], 'indiv_seqs':[i['seq']], 'indiv_seq_counts':[i['count']]}}
					if freq_attribute_name:
						new_data[cluster_seq_id]['other']['indiv_seq_freqs'] = [float(i['other'][freq_attribute_name])]
						new_data[cluster_seq_id]['other'][freq_attribute_name] = float(i['other'][freq_attribute_name])
					if 'timepoint' in i['other']:
						new_data[cluster_seq_id]['other']['timepoint'] = i['other']['timepoint']
			#now update data
			self.data = []
			self.counts = []
			for i in new_data:
				#deal with special case of the first sequence cluster, if desired
				if keep_first_seq_first:
					if i == first_seq_key:
						self.data.insert(0, new_data[i])
						self.counts.insert(0, new_data[i]['count'])
						index = 0
					else:
						self.data.append(new_data[i])
						self.counts.append(new_data[i]['count'])
						index = -1
				else:
					self.data.append(new_data[i])
					self.counts.append(new_data[i]['count'])
					index = -1
				#get the representative seq for the cluster
				if cluster_representative_seq == 'most_abundant':
					if freq_attribute_name:
						max_seq_index = new_data[i]['other']['indiv_seq_freqs'].index(max(new_data[i]['other']['indiv_seq_freqs']))
					else:
						max_seq_index = new_data[i]['other']['indiv_seq_counts'].index(max(new_data[i]['other']['indiv_seq_counts']))
					self.data[index]['seq'] = new_data[i]['other']['indiv_seqs'][max_seq_index]
				self.data[index]['other']['indiv_seq_counts'] = ','.join([str(j) for j in self.data[index]['other']['indiv_seq_counts']])
				self.data[index]['other']['indiv_seqs'] = ','.join(self.data[index]['other']['indiv_seqs'])
				if freq_attribute_name:
					self.data[index]['other']['indiv_seq_freqs'] = ','.join([str(j) for j in self.data[index]['other']['indiv_seq_freqs']])
		return

	def make_frequency_trajectories(self, output_filepath=None, attribute_that_groups_seqs='seq', freq_attribute_name=None, sort_by='sum'):
		"""
		This method will gather frequency trajectories of the data in the MSA over time. The seqs in the data can be grouped together based on any attribute in the data. For example, if on wants the freq trajectories of each unique cluster in the data, then they would group the seqs by 'cluster_id' (or whatever attribute name gives the cluster information). This method only works if there is time-point info in the alignment. It will return 'freq_traj_dic' which is a dictionary object, where each key is the unique identifier of the attribute that groups the seq, and the definition of each key is a list that gives the frequency values at each time-point in the data (i.e. a freq trajectory).
		output_filepath - This is the path the output file for whcich the freq trajs can be written. If this is None, then no output is written.
		attribute_that_groups_seqs - This gives the name of the attribute that will group the seqs (and thus the freqs that correspond to the seqs). Can be any attribute in the self.data object, or in the self.data['other'] object. Default id 'seq' which means that seq entries will be grouped based on there actual sequence identity (i.e. seqs have to be identical to be grouped)
		freq_attribute_name - If defined (default, None), then this gives the name of the attribute that provides the frequency information for each of the seqs. If this is None, then the freq values are calculated internally using the counts for each of the seqs and then normalizing by the total counts in a time-point.
		sort_by - This gives the sorting mechanism for sorting the freq trajs when writing to output. This parameter is only considered if 'output_filepath' is defined. Acceptable values are:
			'sum' - This means the entries will be sorted by the sum of their frequency trajectory.
		"""
		if not 'timepoint' in self.data[0]['other']:
			print 'There does not appear to be time-point information in the data, so cannot make trajectories. Aborting.'
			return
		tpoints = set()
		for i in self.data:
			tpoints.update([float(i['other']['timepoint'])])
		tpoint_to_index_dic = {}
		for index, tpoint in enumerate(sorted(tpoints)):
			tpoint_to_index_dic[tpoint] = index
		freq_traj_dic = {}
		if not freq_attribute_name:
			tpoint_totals_list = [0 for i in tpoint_to_index_dic]
		for i in self.data:
			tpoint = float(i['other']['timepoint'])
			try:
				key = i[attribute_that_groups_seqs]
			except KeyError:
				key = i['other'][attribute_that_groups_seqs]
			if freq_attribute_name:
				freq = float(i['other'][freq_attribute_name])
			else:
				freq = float(i['count'])
				tpoint_totals_list[tpoint_to_index_dic[tpoint]] += freq
			try:
				freq_traj_dic[key][tpoint_to_index_dic[tpoint]] += freq
			except KeyError:
				freq_traj_dic[key] = [0 for j in tpoint_to_index_dic]
				freq_traj_dic[key][tpoint_to_index_dic[tpoint]] += freq
		if not freq_attribute_name:
			for i in freq_traj_dic:
				for tpoint_index in xrange(len(tpoint_totals_list)):
					freq_traj_dic[i][tpoint_index] /= tpoint_totals_list[tpoint_index]
		if output_filepath:
			fileout = open(output_filepath, "w")
			fileout.write('\t%s\n' % '\t'.join([str(i) for i in sorted(tpoint_to_index_dic)]))
			freq_traj_list = []
			if sort_by == 'sum':
				freq_traj_list = [[sum(freq_traj_dic[i])] + [i] + freq_traj_dic[i] for i in freq_traj_dic]
			for i in sorted(freq_traj_list):
				fileout.write("%s\t%s\n" % (i[1], '\t'.join([str(j) for j in i[2:]])))
			fileout.close()
		return


# alignment = msa(filepath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/hiv/10.fasta', count_attribute_name='DUPCOUNT')
# con_seq = alignment.get_consensus_seq(only_for_timepoint=1229.0, output_filepath='/Users/nstrauli/Desktop/butt.fasta', method='cons')
# mean_percent_dist = alignment.get_mean_dist_to_seq(seq=con_seq, time_series_data=True, output_filepath='/Users/nstrauli/Desktop/butt', method='needle', temp_dirpath=None)
# print mean_percent_dist
