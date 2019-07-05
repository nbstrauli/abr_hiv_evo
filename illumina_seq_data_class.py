import os
import itertools
import re

class illumina_seq_data(object):
	"""
	A class for analyzing, and manipulating the raw fastq files that come strait out of the illumina sequencing machines. Importantly, this class does not read the entire file into memory, as these fastq files can be quite large. Instead, it walks through each entry, one at a time.

	Attributes:
	read1_filepath - A string. Path to the fastq formatted file for read1 (if paired end read data, otherwise just the read data) of the illumina output.
	read2_filepath - A string. Same as 'read1_filepath' but for the 2nd read. Default is that this is None, which implies that this is not paired end read data.
	"""

	def __init__(self, read1_fastq_filepath, read2_fastq_filepath=None):
		self.read1_filepath = read1_fastq_filepath
		self.read2_filepath = read2_fastq_filepath
		return

	def de_multiplex(self, barcodes, barcode_start_pos=1, sample_names=None, read_with_barcode='read1', output_dirpath=None, write_as_fasta=False, write_unknown_barcode_seqs=True, trim_barcode=True, barcode_counts_output_filepath=None):
		"""
		This method will take a list of barcode sequences and attempt to match each entry in the Illumina fastq files to one of those barcodes. If desired, it will write the sequence entries that match a given barcode to seperate output files.
		barcodes - a list of strings, where each string is a barcode (as in nucleic acid sequence). All barcodes should be the same length.
		barcode_start_pos - This is an integer that gives the position in each sequence where the barcode begins. Indexed at 1, so if the very first position in the each sequence is where the barcode begins, then this should be 1 (default).
		sample_names - If defined (default is that it is not), then this is a list of sample names that each of the barcodes represents. This should be a list of strings, of the same length as 'barcodes'. The order of sample names should be the same as 'barcodes', such that the 1st sample name corresponds to the 1st barcode, etc. If this is undefined (default), then the sample names are set to be equal to barcodes.
		read_with_barcode - This identifies the read that contains the barcode. This parameter is only considered if 'read2_fastq_filepath' in the object is defined. Exceptable values are:
			'read1' - This means the barcode is located on the 1st read.
			'read2' - This means the barcode is located on the 2nd read.
		output_dirpath - If defined (default is that it's not), this gives the path to the directory for which the output will be written. So, if defined, this instructs the script to write the output to disk. The structure of the output is a subdirectory for each sample (i.e. barcode) and two files within each sample-subdir, labeled 'read1.fast[qa], read2.fast[qa]'. 
			- This can also be a list of output directories, where each directory corresponds to the output directory for each of the sample names (in the same order). If 'write_unknown_barcode_seqs' is True, then the final element in the list of directory paths should correspond to the output directory path for which the unknown seqs should be written.
		write_as_fasta - If True (default False). This will instruct the method to write the output in fasta format. Otherwise the Q score information is retained and written in fastq format. This variable is only considered if output_dirpath is defined.
		write_unknown_barcode_seqs - If defined (default, True), then this means that any sequence (or pair of seqs) in the input fastq file(s) that don't have an identifyable barcode will be written to disk. This parameter is only considered if output_dirpath is defined. Exceptable values are:
			-True. If True then it will simply write the unknown seqs to a filepath that is named after the the read1 input fastq filepath.
			-any other string. It will use this string as the name of the fastq file that it will write the unknown seqs to. It will add a '.fast[qa]' to the end of this name, so don't include that suffix as input.
			-None/False. The unknown seqs will not be written to disk.
		trim_barcode - If True (default), then the barcode sequence, and sequence before (upstream, or 5') to the barcode will be trimmed off each seq where a barcode could be identified.
		barcode_counts_output_filepath - If defined (default, None), then this is the path to the file that the method will write the counts for each of the barcodes found. This is a tab delimited file. It will also write the counts for all the unknown barcodes found.
		"""
		#if no sample names, then make them equal to barcodes
		if not sample_names:
			sample_names = barcodes[:]
		#make sure barcodes and sample names have same length
		if len(barcodes) != len(sample_names):
			print 'barcodes and sample_names parameters must have same length. barcodes:', barcodes, 'sample_names:', sample_names
			return	
		#make sure that there is a directory path for the unknown seqs, if that is desired
		if type(output_dirpath) is list and write_unknown_barcode_seqs == True:
			if len(output_dirpath)-1 != len(sample_names):
				print '"write_unknown_barcode_seqs" parameter is set to True, but there does not appear to be an output directory set for the unknown seqs. "output_dirpath":', output_dirpath
				return

		#if writing the output to disk then set up the output directory and open all necessary output file handles. These file handles are saved in a dictionary
		if output_dirpath:
			if not type(output_dirpath) is list:
				if output_dirpath[-1] != '/':
					output_dirpath += '/'
				if not os.path.exists(output_dirpath):
					os.makedirs(output_dirpath)
			else:
				for i in xrange(len(output_dirpath)):
					if output_dirpath[i][-1] != '/':
						output_dirpath[i] += '/'
					if not os.path.exists(output_dirpath[i]):
						os.makedirs(output_dirpath[i])
			#make a dictionary that maps sample names to an output file handle
			fileout_dic = {}
			#set output file suffix
			if write_as_fasta:
				suffix = '.fasta'
			else:
				suffix = '.fastq'
			#cycle through sample names to open all output file handles
			if self.read2_filepath:
				for count, i in enumerate(sample_names):
					if type(output_dirpath) is list:
						output_dir = output_dirpath[count]
					else:
						output_dir = output_dirpath
					#make subdirectories, if paired read data
					if not os.path.exists(output_dir + i):
						os.makedirs(output_dir + i)
					output_filepath_read1 = output_dir + i + '/read1' + suffix
					output_filepath_read2 = output_dir + i + '/read2' + suffix
					fileout_dic[i] = [open(output_filepath_read1, "w"), open(output_filepath_read2, "w")]
				if write_unknown_barcode_seqs:
					if type(output_dirpath) is list:
						output_dir = output_dirpath[-1]
					else:
						output_dir = output_dirpath
					if write_unknown_barcode_seqs == True:
						if not os.path.exists(output_dir + re.sub('/', '_', os.path.dirname(self.read1_filepath)) + '_unknown'):
							os.makedirs(output_dir + re.sub('/', '_', os.path.dirname(self.read1_filepath)) + '_unknown')
						output_filepath_read1 = output_dir + re.sub('/', '_', os.path.dirname(self.read1_filepath)) + '_unknown/read1' + suffix
						output_filepath_read2 = output_dir + re.sub('/', '_', os.path.dirname(self.read1_filepath)) + '_unknown/read2' + suffix
					else:
						if not os.path.exists(output_dir + write_unknown_barcode_seqs):
							os.makedirs(output_dir + write_unknown_barcode_seqs)
						output_filepath_read1 = output_dir + write_unknown_barcode_seqs + '/read1' + suffix
						output_filepath_read2 = output_dir + write_unknown_barcode_seqs + '/read2' + suffix
					fileout_dic['unknown'] = [open(output_filepath_read1, "w"), open(output_filepath_read2, "w")]
			else:
				for count, i in enumerate(sample_names):
					if type(output_dirpath) is list:
						output_dir = output_dirpath[count]
					else:
						output_dir = output_dirpath
					output_filepath_read1 = output_dir + i + suffix
					fileout_dic[i] = [open(output_filepath_read1, "w")]
				if write_unknown_barcode_seqs:
					if type(output_dirpath) is list:
						output_dir = output_dirpath[-1]
					else:
						output_dir = output_dirpath
					if write_unknown_barcode_seqs == True:
						output_filepath_read1 = output_dir + re.sub('/', '_', os.path.dirname(self.read1_filepath) + '_unknown' + suffix)
					else:
						output_filepath_read1 = output_dir + write_unknown_barcode_seqs + suffix
					fileout_dic['unknown'] = [open(output_filepath_read1, "w")]
				
		filein_read1 = open(self.read1_filepath, "r")
		if self.read2_filepath:
			filein_read2 = open(self.read2_filepath, "r")
		#prime barcode dic
		barcode_count_dic = {}
		for i in barcodes:
			barcode_count_dic[i] = 0
		unknown_barcode_count_dic = {}
		#get length of each barcode
		barcode_len = len(barcodes[0])
		#map barcodes to sample names
		barcode_to_sample_name_dic = {}
		for i in xrange(len(barcodes)):
			barcode_to_sample_name_dic[barcodes[i]] = sample_names[i]

		#cycle through each line in input file(s)
		while True:
			read1_entry = []
			if self.read2_filepath:
				read2_entry = []
			for i in xrange(4):
				read1_entry.append(filein_read1.readline())
				if self.read2_filepath:
					read2_entry.append(filein_read2.readline())
			found_barcode = False
			if not read1_entry[0]:
				break
			else:
				#get barcode
				if read_with_barcode == 'read1':
					barcode = read1_entry[1][barcode_start_pos-1:barcode_start_pos+(barcode_len-1)]
				elif read_with_barcode == 'read2':
					barcode = read2_entry[1][barcode_start_pos-1:barcode_start_pos+(barcode_len-1)]
				try:
					barcode_count_dic[barcode] += 1
					sample_name = barcode_to_sample_name_dic[barcode]
				except KeyError:
					try:
						unknown_barcode_count_dic[barcode] += 1
					except KeyError:
						unknown_barcode_count_dic[barcode] = 1
					sample_name = 'unknown'
				#write output, if desired
				if output_dirpath:
					if sample_name == 'unknown' and not write_unknown_barcode_seqs:
						continue
					if not write_as_fasta:
						for i in xrange(4):
							if self.read2_filepath:
								if trim_barcode and (i == 1 or i == 3):
									fileout_dic[sample_name][0].write(read1_entry[i][barcode_start_pos+barcode_len-1:])
									fileout_dic[sample_name][1].write(read2_entry[i][barcode_start_pos+barcode_len-1:])
								else:
									fileout_dic[sample_name][0].write(read1_entry[i])
									fileout_dic[sample_name][1].write(read2_entry[i])
							else:
								if trim_barcode and (i == 1 or i == 3):
									fileout_dic[sample_name][0].write(read1_entry[i][barcode_start_pos+barcode_len-1:])
								else:
									fileout_dic[sample_name][0].write(read1_entry[i])
					else:
						fileout_dic[sample_name][0].write('>' + read1_entry[0][1:])
						if self.read2_filepath:
							fileout_dic[sample_name][1].write('>' + read2_entry[0][1:])
						if trim_barcode:
							fileout_dic[sample_name][0].write(read1_entry[1][barcode_start_pos+barcode_len-1:])
							if self.read2_filepath:
								fileout_dic[sample_name][1].write(read1_entry[1][barcode_start_pos+barcode_len-1:])
						else:
							fileout_dic[sample_name][0].write(read1_entry[1])
							if self.read2_filepath:
								fileout_dic[sample_name][1].write(read1_entry[1])

		filein_read1.close()
		if self.read2_filepath:
			filein_read2.close()
		for i in fileout_dic:
			for j in fileout_dic[i]:
				j.close()

		#write barcode counts to disk as well, if desired
		if barcode_counts_output_filepath:
			fileout = open(barcode_counts_output_filepath, "w")
			fileout.write('barcode\tcount\tknown_or_unknown_in_sample\n')
			barcode_count_list = [[barcode_count_dic[i], i] for i in barcode_count_dic]
			for i in sorted(barcode_count_list, reverse=True):
				fileout.write('%s\t%s\tknown\n' % (i[1], i[0]))
			unknown_barcode_count_list = [[unknown_barcode_count_dic[i], i] for i in unknown_barcode_count_dic]
			for i in sorted(unknown_barcode_count_list, reverse=True):
				fileout.write('%s\t%s\tunknown\n' % (i[1], i[0]))
			fileout.close()

		#print 'Number of reads for each query barcode:'
		#print barcode_count_dic
		#print 'Number of reads for each unknown barcode found:'
		#print unknown_barcode_count_dic
		return

#i = illumina_seq_data('/Users/nstrauli/Desktop/read1.fastq', '/Users/nstrauli/Desktop/read2.fastq')
#i.de_multiplex(barcodes=['ATGG', 'GTGG', 'CTGG'], barcode_start_pos=4, sample_names=['butt', 'weener', 'farts'], read_with_barcode='read1', output_dirpath=['/Users/nstrauli/Desktop/test_out/pat1', '/Users/nstrauli/Desktop/test_out/pat1', '/Users/nstrauli/Desktop/test_out/pat2', '/Users/nstrauli/Desktop/test_out/unknown'], write_as_fasta=False, write_unknown_barcode_seqs='huh?', trim_barcode=False, barcode_counts_output_filepath='/Users/nstrauli/Desktop/barcode_counts')
