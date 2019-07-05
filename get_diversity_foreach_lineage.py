from seq_lineage_time_series_class import seq_lineage_time_series
import os
import subprocess

def get_lineage_diversities(input_dirpath, input_full_sample_dirpath, output_dirpath, abr_lineage=False):

	###### parameters ######
	count_attribute_name = 'DUPCOUNT'
	avoid_timepoints_with_dot = True
	num_alignments_per_job = 100000
	method = 'vsearch'
	path_to_needle = '/netapp/home/nstrauli/tools_c/EMBOSS-6.6.0/emboss/needle'
	path_to_vsearch = '/netapp/home/nstrauli/tools_c/vsearch-2.4.3-linux-x86_64/bin/vsearch'
	temp_dirpath = '/netapp/home/nstrauli/data/abr_hiv_coevo/temp_stuff'
	num_parallel_cores = 1
	try_again = True
	one_job_per = 'lineage'
	if abr_lineage:
		min_freq_cutoff = 0.0
	else:
		min_freq_cutoff = 0.01
	########################

	if input_dirpath[-1] != '/':
		input_dirpath += '/'
	if input_full_sample_dirpath[-1] != '/':
		input_full_sample_dirpath += '/'
	if output_dirpath[-1] != '/':
		output_dirpath += '/'
	if not os.path.exists(output_dirpath):
		os.makedirs(output_dirpath)
	for patient in os.listdir(input_dirpath):
		if patient[0] == '.' or patient[:6] == 'README':
			continue
		print patient

		# if patient != '7':
		# 	continue

		patient_full_sample_dirpath = '%s%s/' % (input_full_sample_dirpath, patient)
		if abr_lineage:
			for gene_pair in os.listdir('%s%s' % (input_dirpath, patient)):
				if gene_pair[0] == '.' or gene_pair[:6] == 'README':
					continue

				# if gene_pair != 'IGHV6-1_IGHJ5':
				# 	continue

				output_lineage_dirpath = '%s%s/%s/' % (output_dirpath, patient, gene_pair)
				if os.path.exists(output_lineage_dirpath):
					subprocess.call(['rm', '-r', output_lineage_dirpath])
				seq_lin_series = seq_lineage_time_series(dirpath='%s%s/%s' % (input_dirpath, patient, gene_pair), count_attribute_name=count_attribute_name, list_of_filepaths=None, avoid_dirpath=None, avoid_timepoints_with_dot=avoid_timepoints_with_dot)
				seq_lin_series.calc_diversity_pi(full_seq_sample_dirpath=patient_full_sample_dirpath, output_dirpath=output_lineage_dirpath, num_alignments_per_job=num_alignments_per_job, method=method, path_to_needle=path_to_needle, path_to_vsearch=path_to_vsearch, temp_dirpath=temp_dirpath, num_parallel_cores=num_parallel_cores, try_again=try_again, one_job_per=one_job_per, min_freq_cutoff=min_freq_cutoff)
		else:
			output_patient_dirpath = '%s%s/' % (output_dirpath, patient)
			if os.path.exists(output_patient_dirpath):
				subprocess.call(['rm', '-r', output_patient_dirpath])
			seq_lin_series = seq_lineage_time_series(dirpath='%s%s' % (input_dirpath, patient), count_attribute_name=count_attribute_name, list_of_filepaths=None, avoid_dirpath=None, avoid_timepoints_with_dot=avoid_timepoints_with_dot)
			seq_lin_series.calc_diversity_pi(full_seq_sample_dirpath=patient_full_sample_dirpath, output_dirpath=output_patient_dirpath, num_alignments_per_job=num_alignments_per_job, method=method, path_to_needle=path_to_needle, path_to_vsearch=path_to_vsearch, temp_dirpath=temp_dirpath, num_parallel_cores=num_parallel_cores, try_again=try_again, one_job_per=one_job_per, min_freq_cutoff=min_freq_cutoff)
	return

if __name__ == '__main__':

	#below are hard coded examples of how we used this code. One will need to replace these files with there own data if they want to use this code.
	
	#for HIV
	#get_lineage_diversities(input_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/lineage_fasta_files/hiv/edit_dist_between_rep_seqs', input_full_sample_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_fasta_files', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/diversity_lineages/hiv/pi/edit_dist_1')
	#for AbR, unique VJ gene pairs
	# get_lineage_diversities(input_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/abr_seq_unique_VJ/changeo', input_full_sample_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/changeo_annotated_fasta_files', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/diversity_lineages/abr/pi/unique_gene_pairs')
	#for AbR, lineages
	# get_lineage_diversities(input_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/lineage_fasta_files/abr/edit_dist_between_rep_seqs/max_edit_dist_within_samps_6_across_samps_35', input_full_sample_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/changeo_annotated_fasta_files', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/diversity_lineages/abr/pi/lineages_max_edit_dist_within_samps_6_across_samps_35', abr_lineage=True)
