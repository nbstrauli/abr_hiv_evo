#!/usr/bin/python
#$ -S /usr/bin/python
#$ -e error
#$ -o out_3
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=200G
#$ -l h_rt=336:00:00

import sys
sys.path.insert(0, './')
from compare_MSAs_class import compare_MSAs
import os
from scipy import stats
import random
import subprocess
from subprocess import Popen, PIPE

class compare_sets_of_MSAs(object):
	"""This class is for comparing groups of multiple sequence alignments to one another. For example, if one wants to compare the MSA of each lineage of a population, to the MSAs of each lineage of another population, one could use this class.

	Input:
	Should be two directory paths that give the paths to the 1st and 2nd directories (respectively) that contain all the fasta formatted files of the sets of alignments.

	Attributes:
	dirpath_[1|2] - This is the path to the 1st or 2nd directory that contain all the MSAs for the first or second set.
	count_attribute_name - This is the name (as a string) of the attribute in the header that gives the count information for each sequence. If None, then it is assumed that each seq has a count of 1.
	MSAs_[1|2] - This is a list of each of the filenames of the MSAs contained in the 1st or 2nd directory. Includes the file suffix (ex. '.fasta')
	MSA_[1|2]_is_filepath - if True (default, False), this means that the 'directory' path for the 1st or 2nd 'sets' of MSAs is in fact a filepath, which means that there is only a single MSA in this set. This will mean the the self.dirpath_[1|2] will be the directory that contains the filepath (provided by the input of dirpath_[1|2]), and will NOT equal dirpath_[1|2].
	MSA_[1|2]_is_grouped_into_subdirs - If True (default, False), this means that the input directory for 1st or 2nd set of MSA's is grouped into subdirectories, which means that the input directory will have a 2 layer directory tree. This also changes the naming scheme for MSAs_[1|2], which will be of the form: 'subdir/msa_filename', instead of just 'msa_filename'.
	reduce_to_change_or_same - If True (default, False), this will cause the alignments to be altered upon initialization of the class such that a position will be coded as a '1' if that position is a changed relative to the same position in the row above it, or a '0' if there was no change. This altering of the alignment data only occurs in the 'seq_[1|2]_array' attributes. The rest of the data is unchanged.
	"""

	def __init__(self, dirpath_1, dirpath_2, count_attribute_name=None, MSA_1_is_filepath=False, MSA_2_is_filepath=False, MSA_1_is_grouped_into_subdirs=False, MSA_2_is_grouped_into_subdirs=False, reduce_to_change_or_same=False):
		self.count_attribute_name = count_attribute_name
		self.MSA_1_is_filepath = MSA_1_is_filepath
		self.MSA_2_is_filepath = MSA_2_is_filepath
		self.MSA_1_is_grouped_into_subdirs = MSA_1_is_grouped_into_subdirs
		self.MSA_2_is_grouped_into_subdirs = MSA_2_is_grouped_into_subdirs
		self.reduce_to_change_or_same = reduce_to_change_or_same
		if MSA_1_is_filepath:
			self.dirpath_1 = os.path.dirname(dirpath_1) + '/'
			self.MSAs_1 = [os.path.basename(dirpath_1)]
		else:
			if dirpath_1[-1] != '/':
				dirpath_1 += '/'
			self.dirpath_1 = dirpath_1
			self.MSAs_1 = []
			for i in os.listdir(dirpath_1):
				if i[0] == '.' or i[:6] == 'README':
					continue
				if MSA_1_is_grouped_into_subdirs:
					for j in os.listdir('%s%s' % (dirpath_1, i)):
						if j[0] == '.' or j[:6] == 'README':
							continue
						self.MSAs_1.append('%s/%s' % (i, j))
				else:
					self.MSAs_1.append(i)
		if MSA_2_is_filepath:
			self.dirpath_2 = os.path.dirname(dirpath_2) + '/'
			self.MSAs_2 = [os.path.basename(dirpath_2)]
		else:
			if dirpath_2[-1] != '/':
				dirpath_2 += '/'
			self.dirpath_2 = dirpath_2
			self.MSAs_2 = []
			for i in os.listdir(dirpath_2):
				if i[0] == '.' or i[:6] == 'README':
					continue
				if MSA_2_is_grouped_into_subdirs:
					for j in os.listdir('%s%s' % (dirpath_2, i)):
						if j[0] == '.' or j[:6] == 'README':
							continue
						self.MSAs_2.append('%s/%s' % (i, j))
				else:
					self.MSAs_2.append(i)
		return

	def calc_MI_foreach_MSA_pair(self, output_array_dirpath=None, output_distribution_dirpath=None, ignore_sites_less_than_X_variants=None, ignore_sites_less_than_X_changes=None, min_num_of_seq_entries=2, perform_permutation_null=False, num_permute=10, permutation_null_comparison_results_dirpath=None, use_comp_cluster=False, temp_dirpath=None):
		"""
		This method will use the 'calc_pairwise_mutual_information' method from the 'compare_MSAs_class' to calculate the mutual information values for each pair of sites, between each pair of MSAs.
		output_array_dirpath - This is the path to the output directory that will be a directory containing tables containing the M.I. values for each column comparison for a given MSA pair. The names of the rows index the msa_1 sites, and the names of the columns index the msa_2 sites. If this is None, then this file is not written.
		output_distribution_dirpath - This is the path to the directory that will contain the files that have one long list (one column) of mutual information values. Each line is one MI value for one of the site comparisons.
		ignore_sites_less_than_X_variants - Nonetype or Int. If defined (default, None), this will cause the method to ignore any sites (i.e. columns) that have a total number of variants that is less than this value.
		ignore_sites_less_than_X_changes - Nonetype or Int. If defined (default, None), this will cause the method to ignore any sites (i.e. columns) that have a total number of changes that is less than this value. A change is defined as a site that id different than the same site in the row directly preceding (above) it. This superceeds the 'ignore_sites_less_than_X_variants' parameter.
		min_num_of_seq_entries - Int. This gives the minimum number of sequence entries in any of the MSAs for them to be considered in this method.
		perform_permutation_null - If True (default, False), this will instruct the method to also calculate the mutual information values after permuting the rows of each of the MSA pairs.
		num_permute - Int. The number of rounds of permutations. Default, 10. Only considered if perform_permutation_null is True
		permutation_null_comparison_results_dirpath - If defined (default, None), this will give the path to the directory that will contain the file that gives the results of using the Mann-Whitney U test comparing the observed MI values to the permuted null MI values. This is only considered if 'perform_permutation_null' is True.
		use_comp_cluster - If True (default, False), this indicates that the computational cluster (utilizing SGE) shall be used. One job for each MSA pair.
		temp_dirpath - path to directory where temporary data will be written.
		"""
		if output_array_dirpath:
			if output_array_dirpath[-1] != '/':
				output_array_dirpath += '/'
			if not os.path.exists(output_array_dirpath):
				os.makedirs(output_array_dirpath)
		if output_distribution_dirpath:
			if output_distribution_dirpath[-1] != '/':
				output_distribution_dirpath += '/'
			if not os.path.exists(output_distribution_dirpath):
				os.makedirs(output_distribution_dirpath)
		if temp_dirpath:
			if temp_dirpath[-1] != '/':
				temp_dirpath += '/'
			if not os.path.exists(temp_dirpath):
				os.makedirs(temp_dirpath)
		else:
			temp_dirpath = os.getcwd() + '/'
		if permutation_null_comparison_results_dirpath:
			if permutation_null_comparison_results_dirpath[-1] != '/':
				permutation_null_comparison_results_dirpath += '/'
			if not os.path.exists(permutation_null_comparison_results_dirpath):
				os.makedirs(permutation_null_comparison_results_dirpath)
		results = []
		random_suffix = str(random.random())
		num_jobs = 0

		if use_comp_cluster:
			for MSA_file_1 in self.MSAs_1:
				MSA_filepath_1 = self.dirpath_1 + MSA_file_1
				MSA_name_1 = '.'.join(MSA_file_1.split('.')[:-1])
				for MSA_file_2 in self.MSAs_2:
					num_jobs += 1
					MSA_filepath_2 = self.dirpath_2 + MSA_file_2
					MSA_name_2 = '.'.join(MSA_file_2.split('.')[:-1])
					temp_filepath = '%s%s_%s' % (temp_dirpath, num_jobs, random_suffix)
					temp_file = open(temp_filepath, "w")
					temp_file.write('%s\t%s\t%s\t%s\n' % (MSA_filepath_1, MSA_filepath_2, MSA_name_1, MSA_name_2))
					temp_file.close()
					if output_array_dirpath:
						output_array_filepath = '%s%s;%s.txt' % (output_array_dirpath, MSA_name_1, MSA_name_2)
						if not os.path.exists(os.path.dirname(output_array_filepath)):
								os.makedirs(os.path.dirname(output_array_filepath))
					if output_distribution_dirpath:
						output_distribution_filepath = '%s%s;%s.txt' % (output_distribution_dirpath, MSA_name_1, MSA_name_2)
						if not os.path.exists(os.path.dirname(output_distribution_filepath)):
							os.makedirs(os.path.dirname(output_distribution_filepath))
			job_name = 'calc_MI_%s' % random_suffix
			p = Popen(['qsub', '-N', job_name, '-t', '1-%s' % num_jobs, 'compare_sets_of_MSAs_class.py', 'calc_MI_foreach_MSA_pair_comp_cluster', temp_dirpath, random_suffix, str(output_array_dirpath), str(output_distribution_dirpath), str(self.count_attribute_name), str(min_num_of_seq_entries), str(ignore_sites_less_than_X_variants), str(perform_permutation_null), str(num_permute), str(permutation_null_comparison_results_dirpath), str(self.reduce_to_change_or_same), str(ignore_sites_less_than_X_changes)], stdout=PIPE, stderr=PIPE)
			out, err = p.communicate()
			print out
			print err
			if permutation_null_comparison_results_dirpath:
				p = Popen(['qsub', '-hold_jid', job_name, 'compare_sets_of_MSAs_class.py', 'gather_test_results_comp_cluster', permutation_null_comparison_results_dirpath, random_suffix, str(num_jobs)], stdout=PIPE, stderr=PIPE)
				out, err = p.communicate()
				print out
				print err

		else:
			for MSA_file_1 in self.MSAs_1:
				MSA_filepath_1 = self.dirpath_1 + MSA_file_1
				MSA_name_1 = '.'.join(MSA_file_1.split('.')[:-1])
				for MSA_file_2 in self.MSAs_2:
					num_jobs += 1
					MSA_filepath_2 = self.dirpath_2 + MSA_file_2
					MSA_name_2 = '.'.join(MSA_file_2.split('.')[:-1])
					if output_array_dirpath:
						output_array_filepath = '%s%s_%s.txt' % (output_array_dirpath, MSA_name_1, MSA_name_2)
					else:
						output_array_filepath = None
					if output_distribution_dirpath:
						output_distribution_filepath = '%s%s_%s.txt' % (output_distribution_dirpath, MSA_name_1, MSA_name_2)
					else:
						output_distribution_filepath = None
					MSA_comp = compare_MSAs(filepath_1=MSA_filepath_1, filepath_2=MSA_filepath_2, count_attribute_name=self.count_attribute_name, remove_tpoints_that_arent_identical=True, reduce_to_change_or_same=self.reduce_to_change_or_same)
					if len(MSA_comp.msa_1.data) < min_num_of_seq_entries:
						# print 'Not enough seq entries: aborting'
						continue
					sites_1, sites_2, MI_values = MSA_comp.calc_pairwise_mutual_information(output_array_filepath=output_array_filepath, output_distribution_filepath=output_distribution_filepath, ignore_sites_less_than_X_variants=ignore_sites_less_than_X_variants, ignore_sites_less_than_X_changes=ignore_sites_less_than_X_changes)
					if perform_permutation_null and MI_values:
						if not output_distribution_dirpath:
							print 'If performing a null permutation, then the parameter "output_distribution_dirpath" must be defined. Aborting'
							return
						output_distribution_filepath = '%s%s_%s_null_permute.txt' % (output_distribution_dirpath, MSA_name_1, MSA_name_2)
						null_MI_values = MSA_comp.permute_rows_and_calc_MI(output_distribution_filepath=output_distribution_filepath, num_permute=num_permute, ignore_sites_less_than_X_variants=ignore_sites_less_than_X_variants, ignore_sites_less_than_X_changes=ignore_sites_less_than_X_changes)
						test_results = stats.mannwhitneyu(MI_values, null_MI_values, alternative='greater')
						results.append([test_results[1], test_results[0], MSA_name_1, MSA_name_2])
			if perform_permutation_null and permutation_null_comparison_results_dirpath:
				results_output_filepath = '%scombined_results.txt' % permutation_null_comparison_results_dirpath
				fileout = open(results_output_filepath, "w")
				fileout.write('MSA_1\tMSA_2\tmannwhitneyu_test_stat\tp_value\n')
				for result in sorted(results):
					fileout.write('%s\t%s\t%s\t%s\n' % (result[2], result[3], result[1], result[0]))
				fileout.close()
		return

#############################################
# Below are scripts that are used for running
# jobs on the computational cluster. They
# should only be used by the methods above.
#############################################

def calc_MI_foreach_MSA_pair_comp_cluster(temp_dirpath, random_suffix, output_array_dirpath, output_distribution_dirpath, count_attribute_name, min_num_of_seq_entries, ignore_sites_less_than_X_variants, perform_permutation_null, num_permute, permutation_null_comparison_results_dirpath, reduce_to_change_or_same, ignore_sites_less_than_X_changes):
	sge_task_id = int(os.environ['SGE_TASK_ID'])

	# if sge_task_id > 10:
	# 	subprocess.call(['rm', filepath_info_file])
	# 	return

	filepath_info_file = '%s%s_%s' % (temp_dirpath, sge_task_id, random_suffix)
	temp_file = open(filepath_info_file, "r")
	filepath_info = temp_file.readline()[:-1].split('\t')
	MSA_filepath_1 = filepath_info[0]
	MSA_filepath_2 = filepath_info[1]
	MSA_name_1 = filepath_info[2]
	MSA_name_2 = filepath_info[3]
	temp_file.close()
	subprocess.call(['rm', filepath_info_file])
	if output_array_dirpath != 'None':
		output_array_filepath = '%s%s;%s.txt' % (output_array_dirpath, MSA_name_1, MSA_name_2)
	else:
		output_array_filepath = None
	if output_distribution_dirpath != 'None':
		output_distribution_filepath = '%s%s;%s.txt' % (output_distribution_dirpath, MSA_name_1, MSA_name_2)
	else:
		output_distribution_filepath = None
	if count_attribute_name == 'None':
		count_attribute_name = None
	min_num_of_seq_entries = int(min_num_of_seq_entries)
	if ignore_sites_less_than_X_variants == 'None':
		ignore_sites_less_than_X_variants = None
	else:
		ignore_sites_less_than_X_variants = int(ignore_sites_less_than_X_variants)
	if perform_permutation_null == 'False':
		perform_permutation_null = False
	elif perform_permutation_null == 'True':
		perform_permutation_null = True
	num_permute = int(num_permute)
	if permutation_null_comparison_results_dirpath == 'None':
		permutation_null_comparison_results_dirpath = None
	if reduce_to_change_or_same == 'False':
		reduce_to_change_or_same = False
	elif reduce_to_change_or_same == 'True':
		reduce_to_change_or_same = True
	if ignore_sites_less_than_X_changes == 'None':
		ignore_sites_less_than_X_changes = None
	else:
		ignore_sites_less_than_X_changes = int(ignore_sites_less_than_X_changes)
	
	print MSA_filepath_1
	print MSA_filepath_2
	print permutation_null_comparison_results_dirpath

	MSA_comp = compare_MSAs(filepath_1=MSA_filepath_1, filepath_2=MSA_filepath_2, count_attribute_name=count_attribute_name, remove_tpoints_that_arent_identical=True, reduce_to_change_or_same=reduce_to_change_or_same)
	if len(MSA_comp.msa_1.data) < min_num_of_seq_entries:
		# print 'Not enough seq entries: aborting'
		return
	sites_1, sites_2, MI_values = MSA_comp.calc_pairwise_mutual_information(output_array_filepath=output_array_filepath, output_distribution_filepath=output_distribution_filepath, ignore_sites_less_than_X_variants=ignore_sites_less_than_X_variants, ignore_sites_less_than_X_changes=ignore_sites_less_than_X_changes)
	if perform_permutation_null and MI_values:
		if not output_distribution_dirpath:
			print 'If performing a null permutation, then the parameter "output_distribution_dirpath" must be defined. Aborting'
			return
		output_distribution_filepath = '%s%s;%s_null_permute.txt' % (output_distribution_dirpath, MSA_name_1, MSA_name_2)
		null_MI_values = MSA_comp.permute_rows_and_calc_MI(output_distribution_filepath=output_distribution_filepath, num_permute=num_permute, ignore_sites_less_than_X_variants=ignore_sites_less_than_X_variants)
		test_results = stats.mannwhitneyu(MI_values, null_MI_values) #alternative='greater' --> can't use the 'alternative' argument here b/c cluster has old version of scipy. so we do this hack below
		if sum(MI_values)/float(len(MI_values)) <= sum(null_MI_values)/float(len(null_MI_values)):
			p_val = min(test_results[1] + 0.5, 1.)
		else:
			p_val = test_results[1]
		if permutation_null_comparison_results_dirpath:
			fileout = open('%s%s_%s' % (permutation_null_comparison_results_dirpath, sge_task_id, random_suffix), "w")
			fileout.write('%s\t%s\t%s\t%s\n' % (p_val, test_results[0], MSA_name_1, MSA_name_2))
			fileout.close()
	return

def gather_test_results_comp_cluster(permutation_null_comparison_results_dirpath, random_suffix, num_jobs):
	num_jobs = int(num_jobs)
	results = []
	for i in xrange(num_jobs):
		temp_filepath = '%s%s_%s' % (permutation_null_comparison_results_dirpath, i+1, random_suffix)
		if not os.path.exists(temp_filepath):
			continue
		temp_file = open(temp_filepath, "r")
		result_line = temp_file.readline()[:-1].split('\t')
		temp_file.close()
		subprocess.call(['rm', temp_filepath])
		results.append([float(result_line[0]), float(result_line[1]), result_line[2], result_line[3]])
	output_filepath = '%scombined_results.txt' % permutation_null_comparison_results_dirpath
	fileout = open(output_filepath, "w")
	fileout.write('MSA_1\tMSA_2\tmannwhitneyu_test_stat\tp_value\n')
	for result in sorted(results):
		fileout.write('%s\t%s\t%s\t%s\n' % (result[2], result[3], result[1], result[0]))
	fileout.close()
	return


if __name__ == '__main__':
	if sys.argv[1] == 'calc_MI_foreach_MSA_pair_comp_cluster':
		calc_MI_foreach_MSA_pair_comp_cluster(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13])
	elif sys.argv[1] == 'gather_test_results_comp_cluster':
		gather_test_results_comp_cluster(sys.argv[2], sys.argv[3], sys.argv[4])
	else:
		print '\n\nWha!?\n\n'
