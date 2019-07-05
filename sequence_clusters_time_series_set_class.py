#!/usr/bin/python
#$ -S /usr/bin/python
#$ -e error
#$ -o out_1
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=16G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=200G
#$ -l h_rt=336:00:00

import sys
sys.path.insert(0, './')
import os
from sequence_sample_set_class import sequence_sample_set
from sequence_clusters_time_series_class import sequence_clusters_time_series
from master_directory_class import master_directory
import random
from subprocess import Popen, PIPE
import subprocess
import time

class sequence_clusters_time_series_set(sequence_sample_set):
	"""
	This is a class for managing groups of directories, where each directory contains a time-series set of clustered sequence samples, such that they can be loaded as a 'sequence_clusters_time_series_class'.

	Attributes:
	dirpath - The path to the directory that contains all the sequence samples. There could be other sub-directories within this master-directory that contain other sequence samples. All ending in a fasta or fastq suffix will be included in the set. If this parameter equals None then 'list_of_filepaths' is considered
	sample_filepaths - A list of the absolute filepaths for each of the samples in the set. In no particular order.
	count_attribute_name - this gives the name that corresponds to the count of each sequence entity.
	file_suffix - This gives the suffix that is used in the input seq file data. This is determined by the first filepath that is listed in self.sample_filepaths
	freq_attribute_name - This gives the name of the seq attribute that gives the relative frequency of each of the clusters in the dataset
	indiv_seqs_attribute_name - This gives the name of the seq attribute that lists individual sequences of the seqs for each cluster in the dataset.
	indiv_seq_ids_attribute_name - This gives the name of the seq attribute that lists individual sequence IDs of the seqs for each cluster in the dataset
	indiv_seq_counts_attribute_name - This gives the name of the seq attribute that lists individual sequence counts of the seqs for each cluster in the dataset
	indiv_seq_freqs_attribute_name - This gives the name of the seq attribute that lists individual sequence frequencies of the seqs for each cluster in the dataset
	"""
	def __init__(self, dirpath, count_attribute_name=None, freq_attribute_name=None, indiv_seqs_attribute_name='indiv_seqs', indiv_seq_ids_attribute_name='indiv_seq_ids', indiv_seq_counts_attribute_name='indiv_seq_counts', indiv_seq_freqs_attribute_name='indiv_seq_freqs'):
		if dirpath[-1] != '/':
			dirpath += '/'
		self.dirpath = dirpath
		self.count_attribute_name = count_attribute_name
		self.freq_attribute_name = freq_attribute_name
		self.indiv_seqs_attribute_name = indiv_seqs_attribute_name
		self.indiv_seq_ids_attribute_name = indiv_seq_ids_attribute_name
		self.indiv_seq_counts_attribute_name = indiv_seq_counts_attribute_name
		self.indiv_seq_freqs_attribute_name = indiv_seq_freqs_attribute_name
		d = master_directory(self.dirpath)
		d.get_filepaths_containing_filetypes(['fasta', 'fa', 'fastq', 'fq'])
		self.sample_filepaths = d.files_containing_filetype
		#make these filepaths into a list
		self.sample_filepaths = [i for i in self.sample_filepaths]
		return

	def create_lineages(self, output_fasta_lineage_dirpath=None, temp_dirpath=None, hold_for_job_id=None, distance_metric='closest_seq_pair', path_to_needle=None, distance_units='edit_distance', compare_tpoint_to_all_previous=False, genetic_dist_dirpath=None, write_full_original_seqs=None, muller_plot_output_dirpath=None, start_end_colors=None, max_distance=None, min_cluster_freq=None, lineage_min_freq_cutoff=None, make_master_network_plot=None, indiv_sample_network_dirpath=None, mut_count_attribute_name=None, alignment_method='needle', sim_null_lineages_fasta_dirpath=None, sim_null_lineages_muller_plot_dirpath=None, multiply_num_sim_lineages_by=1, only_simulate=False, prob_of_new_lineage_dirpath=None):
		"""
		This method uses the sequence_cluster_time_series_class to create cluster the clustered sequences within a large array of samples. Specifically, it creates lineages for a large group of sequence directories. The data MUST be structured as a master directory which contains a bunch of subdirectories, where each subdirectory contains a time-series of clustered sequences (i.e. each sub-directory must be compatible with the 'sequence_cluster_time_series_class').
		temp_dirpath - Gives the path to the dir that will hold temporary data
		output_fasta_lineage_dirpath - Gives the path to the directory for which the fasta formatted output lineages will be written.
		hold_for_job_id - If defined (default, None), this should be a string, and will hold starting the array job until a job with the provided ID is finished.
		sim_null_lineages_fasta_dirpath - If defined (default, None), this means that null lineages will be simulated and then written in fasta format (similar to those written in 'output_fasta_dirpath'). The lineages are simulated by randomly choosing an ancestor for a given cluster, including the possibility of no ancestor, and building this out moving down the time-series.
		sim_null_lineages_muller_plot_dirpath - If defined (default, False), this will be the path to the directory that will contain the information to make Muller plots of the simulated null lineages. This is only considered if 'sim_null_lineages_fasta_dirpath' is defined.
		multiply_num_sim_lineages_by - Int. This gives the fold increase for the number of clusters per time-point in the null lineages. This effectively increases the size of the null lineages data by X fold
		only_simulate - If True (default, False), this will cause the method to only create and record simulated lineages (see 'sim_null_lineages_fasta_dirpath'). If this is True, then no 'true' lineages will be created and none of that information will be recorded, so the parameters like 'muller_plot_output_dirpath' and 'output_fasta_dirpath' will be ignored. If this is true, then the parameter 'prob_of_new_lineage_filepath' MUST be defined. This is because these probabilities are learned from the observed data.
        prob_of_new_lineage_dirpath - If defined (default, None), then this is the path to the directory that contains the information for what the probability for a new lineage occurring in each of the time-points (the 1st time-point, it is of course, 1), for each of the seq sets. If 'only_simulate' is True, then it is assumed that this file already exists, and if it is False, then this file is made from the data. The data files in this directory MUST end in '.txt'
		"""
		if not output_fasta_lineage_dirpath and not sim_null_lineages_fasta_dirpath:
			print "Either 'sim_null_lineages_fasta_dirpath' or 'output_fasta_lineage_dirpath' must be defined. Aborting."
			return

		if temp_dirpath:
			if temp_dirpath[-1] != '/':
				temp_dirpath += '/'
			if not os.path.exists(temp_dirpath):
				os.makedirs(temp_dirpath)
		else:
			temp_dirpath = os.getcwd() + "/"
		if output_fasta_lineage_dirpath:
			if output_fasta_lineage_dirpath[-1] != '/':
				output_fasta_lineage_dirpath += '/'
			if not os.path.exists(output_fasta_lineage_dirpath):
				os.makedirs(output_fasta_lineage_dirpath)
		if genetic_dist_dirpath:
			if genetic_dist_dirpath[-1] != '/':
				genetic_dist_dirpath += '/'
			if not os.path.exists(genetic_dist_dirpath):
				os.makedirs(genetic_dist_dirpath)
		if write_full_original_seqs:
			if write_full_original_seqs[-1] != '/':
				write_full_original_seqs += '/'
		if muller_plot_output_dirpath:
			if muller_plot_output_dirpath[-1] != '/':
				muller_plot_output_dirpath += '/'
			if not os.path.exists(muller_plot_output_dirpath):
				os.makedirs(muller_plot_output_dirpath)
		if make_master_network_plot:
			if make_master_network_plot[-1] != '/':
				make_master_network_plot += '/'
			if not os.path.exists(make_master_network_plot):
				os.makedirs(make_master_network_plot)
		if indiv_sample_network_dirpath:
			if indiv_sample_network_dirpath[-1] != '/':
				indiv_sample_network_dirpath += '/'
		if start_end_colors:
			start_end_colors = ','.join(start_end_colors)
		if sim_null_lineages_fasta_dirpath:
			if sim_null_lineages_fasta_dirpath[-1] != '/':
				sim_null_lineages_fasta_dirpath += '/'
			if not os.path.exists(sim_null_lineages_fasta_dirpath):
				os.makedirs(sim_null_lineages_fasta_dirpath)
		if prob_of_new_lineage_dirpath:
			if prob_of_new_lineage_dirpath[-1] != '/':
				prob_of_new_lineage_dirpath += '/'
			if not os.path.exists(prob_of_new_lineage_dirpath) and not only_simulate:
				os.makedirs(prob_of_new_lineage_dirpath)
		random_suffix = str(random.random())
		temp_filepath_info_filepath = '%sfilepath_info_%s' % (temp_dirpath, random_suffix)
		fileout = open(temp_filepath_info_filepath, "w")
		num_jobs = 0
		for i in os.listdir(self.dirpath):
			if i[0] == '.' or i[:6] == 'README':
				continue
			num_jobs += 1
			input_dirpath = '%s%s/' % (self.dirpath, i)
			if output_fasta_lineage_dirpath:
				output_dirpath = '%s%s/' % (output_fasta_lineage_dirpath, i)
			else:
				output_dirpath = '%s%s/' % (sim_null_lineages_fasta_dirpath, i)
			fileout.write('%s\t%s\n' % (input_dirpath, output_dirpath))
		fileout.close()
		job_name = 'clust_across_%s' % random_suffix
		if hold_for_job_id:
			p = Popen(['qsub', '-hold_jid', hold_for_job_id, '-N', job_name, '-t', '1-%s' % num_jobs, 'sequence_clusters_time_series_set_class.py', 'create_lineages_compcluster', temp_filepath_info_filepath, self.count_attribute_name, str(self.freq_attribute_name), self.indiv_seqs_attribute_name, self.indiv_seq_ids_attribute_name, self.indiv_seq_counts_attribute_name, self.indiv_seq_freqs_attribute_name, distance_metric, '/scratch/', str(path_to_needle), distance_units, str(compare_tpoint_to_all_previous), str(genetic_dist_dirpath), str(write_full_original_seqs), str(muller_plot_output_dirpath), str(start_end_colors), str(max_distance), str(min_cluster_freq), str(lineage_min_freq_cutoff), str(make_master_network_plot), str(indiv_sample_network_dirpath), str(mut_count_attribute_name), alignment_method, str(sim_null_lineages_fasta_dirpath), str(sim_null_lineages_muller_plot_dirpath), str(multiply_num_sim_lineages_by), str(only_simulate), str(prob_of_new_lineage_dirpath)], stdout=PIPE, stderr=PIPE)
		else:
			p = Popen(['qsub', '-N', job_name, '-t', '1-%s' % num_jobs, 'sequence_clusters_time_series_set_class.py', 'create_lineages_compcluster', temp_filepath_info_filepath, self.count_attribute_name, str(self.freq_attribute_name), self.indiv_seqs_attribute_name, self.indiv_seq_ids_attribute_name, self.indiv_seq_counts_attribute_name, self.indiv_seq_freqs_attribute_name, distance_metric, '/scratch/', str(path_to_needle), distance_units, str(compare_tpoint_to_all_previous), str(genetic_dist_dirpath), str(write_full_original_seqs), str(muller_plot_output_dirpath), str(start_end_colors), str(max_distance), str(min_cluster_freq), str(lineage_min_freq_cutoff), str(make_master_network_plot), str(indiv_sample_network_dirpath), str(mut_count_attribute_name), alignment_method, str(sim_null_lineages_fasta_dirpath), str(sim_null_lineages_muller_plot_dirpath), str(multiply_num_sim_lineages_by), str(only_simulate), str(prob_of_new_lineage_dirpath)], stdout=PIPE, stderr=PIPE)
		out, err = p.communicate()
		return job_name

#############################################
# Below are scripts that are used for running
# jobs on the computational cluster. They
# should only be used by the methods above.
#############################################

def create_lineages_compcluster(temp_filepath_info_filepath, count_attribute_name, freq_attribute_name, indiv_seqs_attribute_name, indiv_seq_ids_attribute_name, indiv_seq_counts_attribute_name, indiv_seq_freqs_attribute_name, distance_metric, temp_dirpath, path_to_needle, distance_units, compare_tpoint_to_all_previous, genetic_dist_dirpath, write_full_original_seqs, muller_plot_output_dirpath, start_end_colors, max_distance, min_cluster_freq, lineage_min_freq_cutoff, make_master_network_plot, indiv_sample_network_dirpath, mut_count_attribute_name, alignment_method, sim_null_lineages_fasta_dirpath, sim_null_lineages_muller_plot_dirpath, multiply_num_sim_lineages_by, only_simulate, prob_of_new_lineage_dirpath):
	sge_task_id = int(os.environ['SGE_TASK_ID'])
	#try to open the file 3 times. The server can be quite slow, and can cause an IO error sometimes, so we try this.
	try:
		filein = open(temp_filepath_info_filepath, "r")
	except IOError:
		time.sleep(10)
		try:
		   filein = open(temp_filepath_info_filepath, "r")
		except IOError:
			time.sleep(10)
			filein = open(temp_filepath_info_filepath, "r")
	count = 0
	for i in filein:
		count += 1
		if count == sge_task_id:
			line = i[:-1].split('\t')
			input_dirpath = line[0]
			output_dirpath = line[1]
			sub_dir_name = os.path.basename(input_dirpath[:-1])
			break
	filein.close()

	# if sge_task_id != 1:
	# 	return
	# if output_dirpath != '/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/lineage_fasta_files/abr/edit_dist_between_rep_seqs/max_edit_dist_within_samps_6_across_samps_30/2/IGHV4-31_IGHJ6/':
	# 	return
	print input_dirpath
	print output_dirpath
	sys.stdout.flush()

	if freq_attribute_name == 'None':
		freq_attribute_name = None
	if path_to_needle == 'None':
		path_to_needle = None
	if compare_tpoint_to_all_previous == 'True':
		compare_tpoint_to_all_previous = True
	elif compare_tpoint_to_all_previous == 'False':
		compare_tpoint_to_all_previous = False
	if genetic_dist_dirpath == 'None':
		genetic_dist_filepath = None
	else:
		genetic_dist_filepath = '%s%s.txt' % (genetic_dist_dirpath, sub_dir_name)
	if write_full_original_seqs == 'None':
		write_full_original_seqs = None
	else:
		write_full_original_seqs = '%s%s/' % (write_full_original_seqs, sub_dir_name)
	if muller_plot_output_dirpath == 'None':
		muller_plot_output_dirpath = None
	else:
		muller_plot_output_dirpath = '%s%s/' % (muller_plot_output_dirpath, sub_dir_name)
	if start_end_colors == 'None':
		start_end_colors = None
	else:
		start_end_colors = start_end_colors.split(',')
	if max_distance == 'None':
		max_distance = None
	elif max_distance != 'parsimony':
		max_distance = float(max_distance)
	if min_cluster_freq == 'None':
		min_cluster_freq = None
	else:
		min_cluster_freq = float(min_cluster_freq)
	if lineage_min_freq_cutoff == 'None':
		lineage_min_freq_cutoff = None
	else:
		lineage_min_freq_cutoff = float(lineage_min_freq_cutoff)
	if make_master_network_plot == 'None':
		make_master_network_plot = None
	else:
		make_master_network_plot = '%s%s/' % (make_master_network_plot, sub_dir_name)
	if indiv_sample_network_dirpath == 'None':
		indiv_sample_network_dirpath = None
	else:
		indiv_sample_network_dirpath = '%s%s/' % (indiv_sample_network_dirpath, sub_dir_name)
	if mut_count_attribute_name == 'None':
		mut_count_attribute_name = None
	if sim_null_lineages_fasta_dirpath == 'None':
		sim_null_lineages_fasta_dirpath = None
	else:
		sim_null_lineages_fasta_dirpath = '%s%s/' % (sim_null_lineages_fasta_dirpath, sub_dir_name)
	if sim_null_lineages_muller_plot_dirpath == 'None':
		sim_null_lineages_muller_plot_dirpath = None
	else:
		sim_null_lineages_muller_plot_dirpath = '%s%s/' % (sim_null_lineages_muller_plot_dirpath, sub_dir_name)
	multiply_num_sim_lineages_by = int(multiply_num_sim_lineages_by)
	if only_simulate == 'True':
		only_simulate = True
	elif only_simulate == 'False':
		only_simulate = False
	if prob_of_new_lineage_dirpath == 'None':
		prob_of_new_lineage_filepath = None
	else:
		prob_of_new_lineage_filepath = '%s%s.txt' % (prob_of_new_lineage_dirpath, sub_dir_name)

	seq_clust_tser = sequence_clusters_time_series(dirpath=input_dirpath, count_attribute_name=count_attribute_name, freq_attribute_name=freq_attribute_name, indiv_seqs_attribute_name=indiv_seqs_attribute_name, indiv_seq_ids_attribute_name=indiv_seq_ids_attribute_name, indiv_seq_counts_attribute_name=indiv_seq_counts_attribute_name, indiv_seq_freqs_attribute_name=indiv_seq_freqs_attribute_name)
	seq_clust_tser.create_lineages(distance_metric=distance_metric, temp_dirpath=temp_dirpath, path_to_needle=path_to_needle, distance_units=distance_units, compare_tpoint_to_all_previous=compare_tpoint_to_all_previous, genetic_dist_filepath=genetic_dist_filepath, output_fasta_dirpath=output_dirpath, write_full_original_seqs=write_full_original_seqs, muller_plot_output_dirpath=muller_plot_output_dirpath, start_end_colors=start_end_colors, max_distance=max_distance, min_cluster_freq=min_cluster_freq, use_comp_cluster=False, lineage_min_freq_cutoff=lineage_min_freq_cutoff, make_master_network_plot=make_master_network_plot, indiv_sample_network_dirpath=indiv_sample_network_dirpath, mut_count_attribute_name=mut_count_attribute_name, alignment_method=alignment_method, sim_null_lineages_fasta_dirpath=sim_null_lineages_fasta_dirpath, sim_null_lineages_muller_plot_dirpath=sim_null_lineages_muller_plot_dirpath, multiply_num_sim_lineages_by=multiply_num_sim_lineages_by, only_simulate=only_simulate, prob_of_new_lineage_filepath=prob_of_new_lineage_filepath)

	return

if __name__ == '__main__':
	if sys.argv[1] == 'create_lineages_compcluster':
		create_lineages_compcluster(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13], sys.argv[14], sys.argv[15], sys.argv[16], sys.argv[17], sys.argv[18], sys.argv[19], sys.argv[20], sys.argv[21], sys.argv[22], sys.argv[23], sys.argv[24], sys.argv[25], sys.argv[26], sys.argv[27], sys.argv[28], sys.argv[29])
	else:
		print 'Wha!?'
