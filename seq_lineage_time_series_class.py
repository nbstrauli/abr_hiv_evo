#!/usr/bin/python
#$ -S /usr/bin/python
#$ -e error
#$ -o out_2
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=200G
#$ -l h_rt=336:00:00

import sys
sys.path.insert(0, './')
import os
from master_directory_class import master_directory
from sequence_time_series_class import sequence_time_series
from sequence_sample_set_class import sequence_sample_set
from subprocess import Popen, PIPE
import subprocess
from sequence_sample_class import sequence_sample
import random
import time

class seq_lineage_time_series(sequence_sample_set):
	"""
	This class is for datasets of sequence data (in fasta or fastq format) that contain the seq data information of various lineages of sequences, over time. The data should be structured as a two level directory tree, where the directories in the first level represent each unique lineage, and the seq files in each of these 'lineage' directories give the seq data for in each of the time-points, for a given lineage. There should be no other files or information in this input directory tree (other than optional 'README' files). This class inherits from the 'sequence_sample_set' class.

	Attributes:
	dirpath - The path to the directory that contains all the sequence samples. There could be other sub-directories within this master-directory that contain other sequence samples. All ending in a fasta or fastq suffix will be included in the set. If this parameter equals None then 'list_of_filepaths' is considered
	sample_filepaths - A list of the absolute filepaths for each of the samples in the set. In no particular order.
	count_attribute_name - this gives the name that corresponds to the count of each sequence entity.
	file_suffix - This gives the suffix that is used in the input seq file data. This is determined by the first filepath that is listed in self.sample_filepaths
	lineages - This is a dictionary, where each key gives the name of a lineage, and the def for each key is a list of ints, where each element of this list gives the time-point (i.e. filename) of a seq data file. The timepoint is each of these lists are strings, but they are listed in numerical order.
	timepoints - This is a list of floats in ascending order that gives all the time-points found in the data. Not all the lineages will have sequences in all the time-points, so this attribute will give all the time-points in the dataset.
	lineages_grouped_into_dirs - Boolean. If True, this indicates that the structure of the input directory is such that each lineage directory is grouped into a higher directory that contians groups of lineages. So, the input directory will contain subdirectories or groups of lineages, that each contian sub-subdirectories that contain the time series data of lineages. If true, this will cause the 'lineages' attribute to be slightly different where each key to the dic will contain the name of the lineage group, followed by the name of the lineage, seperated by a '/'. The definition of each lineage remains the same. Ex: if a lineage is in group 'A' and has an ID of '1' then is key in self.lineages will be 'A/1' (like a mini path).
	just_one_lineage - Boolean. If True (default, False), this indicates that instead of the input 'dirpath' being a directory filled with the time-series seq data of multiple lineages, it instead is a directory that contains the time-series seq data of one lineage. This could be used for example if one wants to treat the time-series seq data of an entire population as one, large, lineage.
	lineage_groups - If 'lineages_grouped_into_dirs' is True, then this is a dictionary, where each key is the name of a group of lineages (defined by the name of the directory), and each definition is a set-type of the lineages that belong to that group. If 'lineages_grouped_into_dirs' is False, then this is None.
	min_num_tpoints - If defined (default, None), this gives the min number of time-points that a given lineage can have in order to be included in the data.
	"""

	def __init__(self, dirpath, count_attribute_name=None, list_of_filepaths=None, avoid_dirpath=None, avoid_timepoints_with_dot=False, lineages_grouped_into_dirs=False, min_num_tpoints=None):
		"""
		list_of_fileapths - This is a list of filepaths to the files that create the sequence set. They do not necessarily need to be in the same directory. This parameter is only considered if 'dirpath' is equal to None. If they are both equal to None, then throws an error.
		avoid_dirpath - This means that the subdirs and files in this directory path will be omitted from the set. If equals None (default) then nothing is avoided. This parameter is only considered if 'dirpath' is defined.
		avoid_timepoints_with_dot - If True (defalut, False), This will instruct the initialization to avoid including time-points that have a dot (i.e. a '.') in their values. The reason for this is that sometimes we will use a '.5' at the end of a time-point value to signify a duplicate, and we generally want to avoid duplicates here.
		"""
		if dirpath == None and list_of_filepaths == None:
			print "Both 'dirpath' and 'list_of_filepaths' variables are equal to None type. Aborting."
			return
		self.count_attribute_name = count_attribute_name
		self.lineages_grouped_into_dirs = lineages_grouped_into_dirs
		if self.lineages_grouped_into_dirs:
			self.lineage_groups = {}
		else:
			self.lineage_groups = None
		if dirpath and dirpath[-1] != '/':
			dirpath += '/'
		if avoid_dirpath:
			if avoid_dirpath[-1] != '/':
				avoid_dirpath += '/'
		self.dirpath = dirpath
		self.lineages = {}
		if dirpath:
			d = master_directory(self.dirpath)
			d.get_filepaths_containing_filetypes(['fasta', 'fa', 'fastq', 'fq'], avoid_dirpath=avoid_dirpath)
			self.sample_filepaths = d.files_containing_filetype
			#make these filepaths into a list
			self.sample_filepaths = [i for i in self.sample_filepaths]
			for lineage in os.listdir(dirpath):
				if lineage[0] == '.' or lineage[:6] == 'README' or dirpath+lineage+'/' == avoid_dirpath:
					continue

				if lineages_grouped_into_dirs:
					self.lineage_groups[lineage] = set()
					for sub_lineage in os.listdir(dirpath + lineage):
						if sub_lineage[0] == '.' or sub_lineage[:6] == 'README':
							continue
						self.lineage_groups[lineage].update([sub_lineage])
						self.lineages['%s/%s' % (lineage, sub_lineage)] = []
						for tpoint in os.listdir('%s%s/%s' % (dirpath, lineage, sub_lineage)):
							if avoid_timepoints_with_dot:
								if len(tpoint.split('.')) > 2:
									continue
							self.lineages['%s/%s' % (lineage, sub_lineage)].append(tpoint)
						if min_num_tpoints:
							if len(self.lineages['%s/%s' % (lineage, sub_lineage)]) < min_num_tpoints:
								del self.lineages['%s/%s' % (lineage, sub_lineage)]
				else:
					self.lineages[lineage] = []
					for tpoint in os.listdir('%s%s' % (dirpath, lineage)):
						if tpoint[0] == '.' or tpoint[:6] == 'README':
							continue
						if avoid_timepoints_with_dot:
							if len(tpoint.split('.')) > 2:
								continue
						self.lineages[lineage].append(tpoint)
					if min_num_tpoints:
						if len(self.lineages[lineage]) < min_num_tpoints:
							del self.lineages[lineage]
		else:
			self.sample_filepaths = list_of_filepaths
			for i in self.sample_filepaths:
				dirpath_list = dirpath.split('/')
				if lineages_grouped_into_dirs:
					lineage = '%s/%s' % (dirpath_list[-3], dirpath_list[-2])
					try:
						self.lineage_groups[dirpath_list[-3]].update([dirpath_list[-2]])
					except KeyError:
						self.lineage_groups[dirpath_list[-3]] = set([dirpath_list[-2]])
				else:
					lineage = dirpath_list[-2]
				tpoint = dirpath_list[-1]
				if avoid_timepoints_with_dot:
					if len(tpoint.split('.')) > 2:
						continue
				try:
					self.lineages[lineage].append(tpoint)
				except KeyError:
					self.lineages[lineage] = [tpoint]
			if min_num_tpoints:
				lineages_to_remove = []
				for lineage in self.lineages:
					if len(self.lineages[lineage]) < min_num_tpoints:
						lineages_to_remove.append(lineage)
				for lineage in lineages_to_remove:
					del self.lineages[lineage]
		if len(self.sample_filepaths) == 0:
			print 'No data in %s. Aborting.' % os.path.basename(self.dirpath[:-1])
			return
		self.file_suffix = self.sample_filepaths[0].split('.')[-1]
		#order each of the tpoint lists chronologically
		self.timepoints = set()
		for lineage in self.lineages:
			tpoint_floats = []
			for tpoint in self.lineages[lineage]:
				tpoint_float = float('.'.join(tpoint.split('.')[:-1]))
				tpoint_floats.append([tpoint_float, tpoint])
				self.timepoints.update([tpoint_float])
			self.lineages[lineage] = [i[1] for i in sorted(tpoint_floats)]
		self.timepoints = [i for i in sorted(self.timepoints)]

		return

	def calc_diversity_pi(self, full_seq_sample_dirpath, output_dirpath, num_alignments_per_job=100000, method='needle', path_to_needle='needle', path_to_vsearch=None, temp_dirpath=None, num_parallel_cores=12, try_again=False, one_job_per='lineage', min_freq_cutoff=0.0):
		"""
		This method will calculate diversity (pi) for each of the seq time-points, for each of the lineages.
		full_seq_sample_dirpath - This gives the directory path for the full sequence samples from which all of the lineages came from. We need this to get the total count for the entire sample.
		path_to_needle - This is the path to the needle program by EMBOSS that impliments the needleman-wunsch global alignment algorithm. This is what we use to get the genetic distance between two seqs.
        method - This gives the method that will be used to get the genetic distance from the pairs of seqs. The exceptable values for this are:
            'needle' - This means the needleman-wunsch global alignment algorithm will be used in a program called 'needle' in the EMBOSS package
            'pairwise2' - This means the 'pairwise2' Biopython package will be used, which also implements the needleman-wunsch global alignment algorithm.
            'vsearch' - This means the vsearch program will be used. This program is capable of doing an all pairwise alignment by itself, so the flow of this method is significantly changed if this is selected.
        path_to_vsearch - This is the path to the vsearch excecutable. If this is None (default) then the path is assumed to be 'vsearch' (i.e. in $PATH). Ignored if method!='vsearch'
        temp_dirpath - This gives the path to the temporary directory that will contain the fasta files to be submitted to the aligner (if the chosen aligner needs this), and the sub_sum_of_distances values outputed by 'pi_calculator_compCluster'. The temp fasta files can get quite large, so make sure that there is enough space where ever this path leads. If temp_dirpath=None (default) then a temp dir in the current working directory is made for this.
        num_parallel_cores - This gives the number of parallel jobs that vsearch can use. vsearch's default is to use as many cores as it can, and this slows down the cluster for others. So, we need to set this to a fixed amount, and then let the cluster know when submitting array job. Default is 12.
        try_again - Boolean. If True (not default), then the script will keep trying to calculate pi until it works. It seems that vsearch is a bit buggy and can write uninterpretable lines of output. This setting (when set to True) will simply keep trying until it works. Only relevant if method='vsearch'.
        one_job_per - This tells the method how to split up the array jobs. Acceptable values are:
        	'timepoint' - This makes it so one cycles through lineages, and for each lineages, submits and array job where each job is one time-point. The script waits for each lineage to complete before starting the next.
        	'lineage' - Default. This makes it so there is one job per lineage, and within each job the script will cycle through each time-point. This is most likely faster.
        min_freq_cutoff - Float. This gives the frequency threshold that each lineage must pass (at least once) in order to be recorded in the output.
		"""
		if full_seq_sample_dirpath[-1] != '/':
			full_seq_sample_dirpath += '/'
		if output_dirpath[-1] != '/':
			output_dirpath += '/'
		if not os.path.exists(output_dirpath):
			os.makedirs(output_dirpath)
		#get full sample total counts foreach tpoint
		total_counts_dic = {}
		for i in os.listdir(full_seq_sample_dirpath):
			if i[0] == '.' or i[:6] == 'README':
				continue
			tpoint_float = float('.'.join(i.split('.')[:-1]))
			full_seq_sample_filepath = '%s%s' % (full_seq_sample_dirpath, i)
			full_seq_sample = sequence_sample(filepath=full_seq_sample_filepath, count_attribute_name=self.count_attribute_name)
			total_counts_dic[tpoint_float] = full_seq_sample.total
		total_counts_string = ''
		total_counts_list = []
		for tpoint in sorted(total_counts_dic):
			total_counts_string += '%s:%s,' % (tpoint, total_counts_dic[tpoint])
			total_counts_list.append(total_counts_dic[tpoint])
		total_counts_string = total_counts_string[:-1]
		if one_job_per == 'lineage':
			num_jobs = 0
			random_suffix = str(random.random())
			for i in self.lineages:
				num_jobs += 1
				temp_file = open('%s%s_%s' % (temp_dirpath, random_suffix, num_jobs), "w")
				input_dirpath = '%s%s/' % (self.dirpath, i)
				temp_file.write('%s\t%s\n' % (input_dirpath, i))
				temp_file.close()
			num_jobs = len(self.lineages)
			timepoints_string = ','.join([str(i) for i in self.timepoints])
			p = Popen(['qsub', '-t', '1-%s' % num_jobs, 'seq_lineage_time_series_class.py', 'calc_diversity_pi_compcluster', self.dirpath, self.count_attribute_name, str(num_alignments_per_job), method, path_to_needle, path_to_vsearch, temp_dirpath, str(num_parallel_cores), str(try_again), output_dirpath, timepoints_string, total_counts_string, str(min_freq_cutoff), random_suffix], stdout=PIPE, stderr=PIPE)
			out, err = p.communicate()
			print out
			print err
		elif one_job_per == 'timepoint':
			for lineage in self.lineages:
				lineage_dirpath = '%s%s/' % (self.dirpath, lineage)
				seq_samp_time_series = sequence_time_series(dirpath=lineage_dirpath, count_attribute_name=self.count_attribute_name)
				sample_counts = seq_samp_time_series.get_time_series_sample_counts()
				freqs_too_low = True
				for index in xrange(len(sample_counts)):
					freq = sample_counts[index] / total_counts_list[index]
					if freq >= min_freq_cutoff:
						freqs_too_low = False
						break
				if freqs_too_low:
					print 'lineage from %s is composed of seq samples that are too low in frequency, so skipping.' % lineage_dirpath
					continue
				diversity_pi_values = seq_samp_time_series.down_sample_and_calc_pi_compCluster(downsamp_to='Nope', num_downsamp_trials=1, num_alignments_per_job=num_alignments_per_job, method=method, path_to_needle=path_to_needle, path_to_vsearch=path_to_vsearch, temp_dirpath=temp_dirpath, num_parallel_cores=num_parallel_cores, try_again=try_again, no_array_job=False)
				pi_vals_dic = {}
				for i in diversity_pi_values:
					tpoint = float('.'.join(os.path.basename(i).split('.')[:-1]))
					pi_vals_dic[tpoint] = diversity_pi_values[i][0]
				output_filepath = '%s%s.txt' % (output_dirpath, lineage)
				fileout = open(output_filepath, "w")
				fileout.write('timepoint\tdiversity\n')
				for i in self.timepoints:
					if i in pi_vals_dic:
						pi = pi_vals_dic[i]
					else:
						pi = "NA"
					fileout.write('%s\t%s\n' % (i, pi))
				fileout.close()
		return

	def get_mean_attribute(self, full_seq_sample_dirpath, attribute_name, output_dirpath, temp_dirpath, weight_by_counts=True, treat_NAs='ignore', min_freq_cutoff=0.0):
		"""
		This method will calculate the mean value of a given attribute for the seq data in each of the lineages. It will then collect this data as trajectories of the attribute value over time. Writes this data to output.
		attribute_name - String. This gives the name of the attribute as it is recorded in the seq data
		output_dirpath - This gives the path to the directory that each of the trajectories for each of the lineages will be written.
		weight_by_counts - If True (default), then this weights each of the sequence entries values by the 'count' of that sequence in the data (i.e. weighted average).
		treat_NAs - This informs the method how to treat 'NA' values. Acceptable values are:
            'ignore' - Default. This means that seq entries that have a value of 'NA' are simply ignored. Neither their value nor their count is included in the numerator nor the denominator of the mean calculation, respectively.
            'zero' - This means that a seq entry with a value of 'NA' will be treated as if this value were 0.
        min_freq_cutoff - Float. This gives the frequency threshold that each lineage must pass (at least once) in order to be recorded in the output.
		"""
		if full_seq_sample_dirpath[-1] != '/':
			full_seq_sample_dirpath += '/'
		if output_dirpath[-1] != '/':
			output_dirpath += '/'
		if not os.path.exists(output_dirpath):
			os.makedirs(output_dirpath)
		if temp_dirpath[-1] != '/':
			temp_dirpath += '/'
		if not os.path.exists(temp_dirpath):
			os.makedirs(temp_dirpath)
		num_jobs = 0
		random_suffix = str(random.random())
		for i in self.lineages:
			num_jobs += 1
			temp_filepath = '%s%s_%s' % (temp_dirpath, random_suffix, num_jobs)
			temp_file = open(temp_filepath, "w")
			input_dirpath = '%s%s/' % (self.dirpath, i)
			temp_file.write('%s\t%s\n' % (input_dirpath, i))
			temp_file.close()
		p = Popen(['qsub', '-t', '1-%s' % (num_jobs), 'seq_lineage_time_series_class.py', 'get_mean_attribute_compcluster', attribute_name, output_dirpath, self.dirpath, self.count_attribute_name, str(weight_by_counts), ','.join([str(i) for i in self.timepoints]), treat_NAs, full_seq_sample_dirpath, str(min_freq_cutoff), random_suffix, temp_dirpath], stdout=PIPE, stderr=PIPE)
		out, err = p.communicate()
		print out
		print err
		return

	def get_freq_trajectories(self, full_seq_sample_dirpath, output_dirpath=None, min_freq_cutoff=0.0):
		"""
		This method gets the total count of each of the lineages and then divides this by the total count of the sample to arrive at the relative frequency of each of the lineages. Writes the trajectory of this frequency over time.
		full_seq_sample_dirpath - This gives the directory path for the full sequence samples from which all of the lineages came from. We need this to get the total count for the entire sample.
		output_dirpath - This gives the path to the directory that the frequency trajectories for each of the lineages will be written. If this is None, than no output written.
		"""
		if full_seq_sample_dirpath[-1] != '/':
			full_seq_sample_dirpath += '/'
		if output_dirpath:
			if output_dirpath[-1] != '/':
				output_dirpath += '/'
			if not os.path.exists(output_dirpath):
				os.makedirs(output_dirpath)
		#get full sample total counts foreach tpoint
		total_counts_dic = {}
		for i in os.listdir(full_seq_sample_dirpath):
			if i[0] == '.' or i[:6] == 'README':
				continue
			tpoint_float = float('.'.join(i.split('.')[:-1]))
			full_seq_sample_filepath = '%s%s' % (full_seq_sample_dirpath, i)
			full_seq_sample = sequence_sample(filepath=full_seq_sample_filepath, count_attribute_name=self.count_attribute_name)
			total_counts_dic[tpoint_float] = full_seq_sample.total
		lineage_freq_trajs = {}
		for lineage in self.lineages:
			print '\tlineage:', lineage
			lineage_dirpath = '%s%s/' % (self.dirpath, lineage)
			lineage_freq_traj = {}
			freqs_too_low = True
			for tpoint in self.lineages[lineage]:
				sample_filepath = '%s%s' % (lineage_dirpath, tpoint)
				seq_sample = sequence_sample(filepath=sample_filepath, count_attribute_name=self.count_attribute_name)
				tpoint_float = float('.'.join(tpoint.split('.')[:-1]))
				freq = seq_sample.total / total_counts_dic[tpoint_float]
				if freq >= min_freq_cutoff:
					freqs_too_low = False
				lineage_freq_traj[tpoint_float] = seq_sample.total / total_counts_dic[tpoint_float]
			if freqs_too_low:
				print 'lineage from %s is composed of seq samples that are too low in frequency, so skipping.' % lineage_dirpath
				continue
			lineage_freq_trajs[lineage] = lineage_freq_traj
			if output_dirpath:
				output_filepath = '%s%s.txt' % (output_dirpath, lineage)
				if not os.path.exists(os.path.dirname(output_filepath)):
					os.makedirs(os.path.dirname(output_filepath))
				fileout = open(output_filepath, "w")
				fileout.write('timepoint\tfrequency\n')
				for i in self.timepoints:
					if i in lineage_freq_traj:
						freq = lineage_freq_traj[i]
					else:
						freq = 0.
					fileout.write('%s\t%s\n' % (i, freq))
				fileout.close()
		return lineage_freq_trajs

	def make_time_series_alingments(self, output_dirpath, path_to_mafft=None, temp_dirpath=None, align_attribute_value_instead=None, wait_till_cluster_done=False):
		"""
		This method will make alignments of the representative seqs for each lineage at each time-point. In other words, it will find the most numerous sequnece for each of the time-points of a lineage, and then align these representative seqs. It then writes these alignments in fasta format to an output directory: one output file for each lineage.
		output_dirpath - This gives the path to the directory that will contain all the MSAs
		path_to_mafft - This gives the path to the executable that will execute the MAFFT multiple aligner. If this is None, then it is assumed that MAFFT is in $PATH
		align_attribute_value_instead - If defined (default, None), this will align (and write) the value of the provided attribute. If defined this must be the name of an attribute in the seq data, and the value of this attribute should be a sequence of some sort. This is essentially a hack to allow the alignment of CDR3 seq data.
		wait_till_cluster_done - Boolean. If True (default, False), this will cause the method to not exit until the job on the cluster has completed.
		"""
		if output_dirpath[-1] != '/':
			output_dirpath += '/'
		if not os.path.exists(output_dirpath):
			os.makedirs(output_dirpath)
		if temp_dirpath:
			if temp_dirpath[-1] != '/':
				temp_dirpath += '/'
			if not os.path.exists(temp_dirpath):
				os.makedirs(temp_dirpath)
		else:
			temp_dirpath = os.getcwd() + '/'
		random_suffix = str(random.random())

		num_jobs = 0
		for i in self.lineages:
			num_jobs += 1
			temp_filepath = '%s%s_%s' % (temp_dirpath, random_suffix, num_jobs)
			temp_file = open(temp_filepath, "w")
			input_dirpath = '%s%s/' % (self.dirpath, i)
			temp_file.write('%s\t%s\n' % (input_dirpath, i))
			temp_file.close()
		p = Popen(['qsub', '-t', '1-%s' % num_jobs, 'seq_lineage_time_series_class.py', 'make_time_series_alingments_compcluster', self.dirpath, random_suffix, output_dirpath, str(self.count_attribute_name), temp_dirpath, str(path_to_mafft), str(align_attribute_value_instead)], stdout=PIPE, stderr=PIPE)
		out, err = p.communicate()
		print out
		print err
		if wait_till_cluster_done:
			submit_job_id = set([out.split()[2].split('.')[0]])
			#sleep till job is done
			good_to_go = False
			while good_to_go == False:
				time.sleep(10)
				p = Popen(['qstat'], stdout=PIPE, stderr=PIPE)
				jobs = p.communicate()[0].split('\n')
				job_ids = set()
				if len(jobs) >= 1:
					for i in jobs[2:-1]:
						job_id = i.split()[0]
						job_ids.update([job_id])
						if len(job_ids.intersection(submit_job_id)) == 0:
							good_to_go = True
				else:
					good_to_go = True
		return

	def get_lineage_stats(self, return_lineage_length_dic=False):
		"""
		This method retrieves various statistics about the lineage structure.
		return_lineage_length_dic - Boolean. If True (default, False), this will return a dic whose index is the name of a lineage, and definition is the length of that lineage. Returns this in addition to 'lineage_stats'.

		RETURNS:
		lineage_stats. This is a dic where each index gives that value of some statistic. This statistics returned here are:
		'mean_length' - This is the average length (i.e. number of time-points) of the lineages
		'length_range' - Two element list. Gives the higest and lowest values that went into 'mean_length'
		'mean_proportion_new_lineages' - Each time-point has a certain proportion of the lineages in it that are novel (i.e. started at that time-point). This gives the mean proportion of novel lineages across all time-points. 
		'proportion_new_lineages_range' - Two element list. Gives the range of 'mean_proportion_new_lineages'
		'mean_total_lineages' - Gives the mean number of lineages across the time-points.
		'total_lineages_range' - Two element list. Gives the range of 'mean_total_lineages'.
		"""
		lineage_stats = {'mean_length':[], 'length_range':[], 'mean_proportion_new_lineages':[], 'proportion_new_lineages_range':[], 'mean_total_lineages':[], 'total_lineages_range':[]}
		first_tpoint = self.timepoints[0]
		tpoint_to_new_lin_dic = {}
		for tpoint in self.timepoints[1:]:
			tpoint_to_new_lin_dic[tpoint] = 0.
		tpoint_to_num_lin_dic = {}
		for tpoint in self.timepoints:
			tpoint_to_num_lin_dic[tpoint] = 0.
		lineage_length_dic = {}
		for lineage in self.lineages:
			tpoints = sorted([float(tpoint[:-6]) for tpoint in os.listdir(self.dirpath+lineage) if tpoint[0] != '.'])
			length = len(tpoints)
			lineage_stats['mean_length'].append(length)
			lineage_length_dic[lineage] = length
			if tpoints[0] != first_tpoint:
				tpoint_to_new_lin_dic[tpoints[0]] += 1
			for tpoint in tpoints:
				tpoint_to_num_lin_dic[tpoint] += 1
		for tpoint in self.timepoints:
			if tpoint != first_tpoint:
				new_lin_prop = tpoint_to_new_lin_dic[tpoint] / tpoint_to_num_lin_dic[tpoint]
				lineage_stats['mean_proportion_new_lineages'].append(new_lin_prop)
			lineage_stats['mean_total_lineages'].append(tpoint_to_num_lin_dic[tpoint])
		lineage_stats['length_range'] = [min(lineage_stats['mean_length']), max(lineage_stats['mean_length'])]
		lineage_lengths = lineage_stats['mean_length'][:]
		lineage_stats['mean_length'] = float(sum(lineage_stats['mean_length'])) / len(lineage_stats['mean_length'])
		lineage_stats['proportion_new_lineages_range'] = [min(lineage_stats['mean_proportion_new_lineages']), max(lineage_stats['mean_proportion_new_lineages'])]
		lineage_stats['mean_proportion_new_lineages'] = sum(lineage_stats['mean_proportion_new_lineages']) / len(lineage_stats['mean_proportion_new_lineages'])
		lineage_stats['total_lineages_range'] = [min(lineage_stats['mean_total_lineages']), max(lineage_stats['mean_total_lineages'])]
		lineage_stats['mean_total_lineages'] = sum(lineage_stats['mean_total_lineages']) / len(lineage_stats['mean_total_lineages'])
		if return_lineage_length_dic:
			return lineage_stats, lineage_length_dic
		else:
			return lineage_stats


#############################################
# Below are scripts that are used for running
# jobs on the computational cluster. They
# should only be used by the methods above.
#############################################

def calc_diversity_pi_compcluster(lineage_master_dirpath, count_attribute_name, num_alignments_per_job, method, path_to_needle, path_to_vsearch, temp_dirpath, num_parallel_cores, try_again, output_dirpath, all_timepoints, total_counts_string, min_freq_cutoff, random_suffix):
	num_alignments_per_job = int(num_alignments_per_job)
	num_parallel_cores = int(num_parallel_cores)
	all_timepoints = [float(i) for i in all_timepoints.split(',')]
	total_counts_dic = {}
	for i in total_counts_string.split(','):
		tpoint = float(i.split(':')[0])
		total_counts_dic[tpoint] = float(i.split(':')[1])
	if try_again == 'True':
		try_again = True
	else:
		try_again = False
	min_freq_cutoff = float(min_freq_cutoff)
	sge_task_id = int(os.environ['SGE_TASK_ID'])
	index = 1
	temp_filepath = '%s%s_%s' % (temp_dirpath, random_suffix, sge_task_id)
	temp_file = open('%s%s_%s' % (temp_dirpath, random_suffix, sge_task_id), "r")
	filepath_info = temp_file.readline()[:-1].split('\t')
	lineage_dirpath = filepath_info[0]
	lineage = filepath_info[1]
	temp_file.close()
	subprocess.call(['rm', temp_filepath])
	print lineage
	sys.stdout.flush()

	# if lineage != 'IGHV4-34_IGHJ4':
	# 	return

	seq_samp_time_series = sequence_time_series(dirpath=lineage_dirpath, count_attribute_name=count_attribute_name)
	sample_counts = seq_samp_time_series.get_time_series_sample_counts()
	freqs_too_low = True
	for index in xrange(len(sample_counts)):
		freq = sample_counts[index] / total_counts_dic[seq_samp_time_series.timepoints[index]]
		if freq >= min_freq_cutoff:
			freqs_too_low = False
			break
	if freqs_too_low:
		print 'lineage from %s is composed of seq samples that are too low in frequency, so aborting.' % lineage_dirpath
		return
	diversity_pi_values = seq_samp_time_series.down_sample_and_calc_pi_compCluster(downsamp_to='Nope', num_downsamp_trials=1, num_alignments_per_job=num_alignments_per_job, method=method, path_to_needle=path_to_needle, path_to_vsearch=path_to_vsearch, temp_dirpath='/scratch/', num_parallel_cores=num_parallel_cores, try_again=try_again, no_array_job=True)
	pi_vals_dic = {}
	for i in diversity_pi_values:
		tpoint = float('.'.join(os.path.basename(i).split('.')[:-1]))
		pi_vals_dic[tpoint] = diversity_pi_values[i][0]
	output_filepath = '%s%s.txt' % (output_dirpath, lineage)
	if not os.path.exists(os.path.dirname(output_filepath)):
		os.makedirs(os.path.dirname(output_filepath))
	fileout = open(output_filepath, "w")
	fileout.write('timepoint\tdiversity\n')
	for i in all_timepoints:
		if i in pi_vals_dic:
			pi = pi_vals_dic[i]
		else:
			pi = "NA"
		fileout.write('%s\t%s\n' % (i, pi))
	fileout.close()
	return

def get_mean_attribute_compcluster(attribute_name, output_dirpath, lineage_master_dirpath, count_attribute_name, weight_by_counts, all_timepoints, treat_NAs, full_seq_sample_dirpath, min_freq_cutoff, random_suffix, temp_dirpath):
	all_timepoints = [float(i) for i in all_timepoints.split(',')]
	if weight_by_counts == 'True':
		weight_by_counts = True
	else:
		weight_by_counts = False
	min_freq_cutoff = float(min_freq_cutoff)
	sge_task_id = int(os.environ['SGE_TASK_ID'])
	temp_filepath = '%s%s_%s' % (temp_dirpath, random_suffix, sge_task_id)
	temp_file = open(temp_filepath, "r")
	filepath_info = temp_file.readline()[:-1].split('\t')
	lineage_dirpath = filepath_info[0]
	lineage = filepath_info[1]
	temp_file.close()
	subprocess.call(['rm', temp_filepath])
	#get full sample time-series counts
	seq_time_series_sample = sequence_time_series(dirpath=full_seq_sample_dirpath, count_attribute_name=count_attribute_name)
	total_counts_dic = seq_time_series_sample.get_time_series_sample_counts(return_dic=True)
	#check that lineage has high enough frequency, abort if not
	seq_time_series_sample = sequence_time_series(dirpath=lineage_dirpath, count_attribute_name=count_attribute_name)
	lineage_counts_dic = seq_time_series_sample.get_time_series_sample_counts(return_dic=True)
	freqs_too_low = True
	for tpoint_float in lineage_counts_dic:
		freq = lineage_counts_dic[tpoint_float] / total_counts_dic[tpoint_float]
		if freq >= min_freq_cutoff:
			freqs_too_low = False
			break
	if freqs_too_low:
		print 'lineage from %s is composed of seq samples that are too low in frequency, so aborting.' % lineage_dirpath
		return
	tpoint_value_dic = {}
	for tpoint_sample in os.listdir(lineage_dirpath):
		if tpoint_sample[0] == '.' or tpoint_sample[:6] == 'README':
			continue
		tpoint_float = float('.'.join(tpoint_sample.split('.')[:-1]))
		sample_filepath = '%s/%s' % (lineage_dirpath, tpoint_sample)
		seq_sample = sequence_sample(filepath=sample_filepath, count_attribute_name=count_attribute_name)
		mean_attribute_value = seq_sample.get_mean_attribute_value(query_attribute_name=attribute_name, weight_by_counts=True, treat_NAs=treat_NAs)
		tpoint_value_dic[tpoint_float] = mean_attribute_value
	output_filepath = '%s%s.txt' % (output_dirpath, lineage)
	fileout = open(output_filepath, "w")
	fileout.write('timepoint\t%s\n' % attribute_name)
	for tpoint in all_timepoints:
		if tpoint in tpoint_value_dic:
			attribute_value = tpoint_value_dic[tpoint]
		else:
			attribute_value = 'NA'
		fileout.write('%s\t%s\n' % (tpoint, attribute_value))
	fileout.close()
	return

def make_time_series_alingments_compcluster(input_dirpath, random_suffix, output_dirpath, count_attribute_name, temp_dirpath, path_to_mafft, align_attribute_value_instead):
	if align_attribute_value_instead == 'None':
		align_attribute_value_instead = None
	if count_attribute_name == 'None':
		count_attribute_name = None
	if path_to_mafft == 'None':
		path_to_mafft = None
	sge_task_id = int(os.environ['SGE_TASK_ID'])

	temp_filepath = '%s%s_%s' % (temp_dirpath, random_suffix, sge_task_id)
	temp_file = open(temp_filepath, "r")
	filepath_info = temp_file.readline()[:-1].split('\t')
	lineage_input_dirpath = filepath_info[0]
	lineage = filepath_info[1]
	temp_file.close()
	subprocess.call(['rm', temp_filepath])
	temp_representative_seqs_dirpath = '/scratch/%s_%s/' % (random_suffix, sge_task_id)
	os.makedirs(temp_representative_seqs_dirpath)
	output_msa_filepath = '%s%s.fasta' % (output_dirpath, lineage)

	# if sge_task_id != 1:
	# 	return
	print lineage_input_dirpath
	print lineage
	print output_msa_filepath

	if not os.path.exists(os.path.dirname(output_msa_filepath)):
		os.makedirs(os.path.dirname(output_msa_filepath))
	seq_time_series = sequence_time_series(dirpath=lineage_input_dirpath, count_attribute_name=count_attribute_name)
	seq_time_series.write_most_abundant_seq_foreach_tpoint(output_dirpath=temp_representative_seqs_dirpath, write_attribute_value_instead=align_attribute_value_instead)
	seq_time_series = sequence_time_series(dirpath=temp_representative_seqs_dirpath, count_attribute_name=count_attribute_name)
	seq_time_series.make_MSA(output_filepath=output_msa_filepath, temp_dirpath='/scratch/', method='mafft', path_to_mafft=path_to_mafft, add_ref_seq=None, ref_seq_name='ref_seq', write_seq_id_first=None)
	subprocess.call(['rm', '-r', temp_representative_seqs_dirpath])
	return

if __name__ == '__main__':
	if sys.argv[1] == 'test':
		seq_lins = seq_lineage_time_series(dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/abr_seq_unique_VJ/changeo/1', count_attribute_name='DUPCOUNT', list_of_filepaths=None, avoid_dirpath=None)
		print seq_lins.lineages
	elif sys.argv[1] == 'calc_diversity_pi_compcluster':
		calc_diversity_pi_compcluster(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13], sys.argv[14], sys.argv[15])
	elif sys.argv[1] == 'get_mean_attribute_compcluster':
		get_mean_attribute_compcluster(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12])
	elif sys.argv[1] == 'make_time_series_alingments_compcluster':
		make_time_series_alingments_compcluster(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8])
