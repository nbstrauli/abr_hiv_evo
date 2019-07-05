import os
import itertools
from seq_lineage_time_series_class import seq_lineage_time_series
import subprocess

def get_mean_trajectories(input_lineage_fasta_dirpath, input_full_sample_dirpath, attibute_names, output_trajectory_dirpaths, abr_lineage=False):

	###### parameters ######
	count_attribute_name = 'DUPCOUNT'
	avoid_timepoints_with_dot = True
	weight_by_counts = True
	treat_NAs = 'ignore'
	if abr_lineage:
		min_freq_cutoff = 0.0
	else:
		min_freq_cutoff = 0.01
	temp_dirpath = '/scrapp/'
	########################

	if input_lineage_fasta_dirpath[-1] != '/':
		input_lineage_fasta_dirpath += '/'
	if input_full_sample_dirpath[-1] != '/':
		input_full_sample_dirpath += '/'
	for patient in os.listdir(input_lineage_fasta_dirpath):
		if patient[0] == '.' or patient[:6] == 'README':
			continue
		print patient

		# if patient != '6':
		# 	continue

		patient_full_sample_dirpath = '%s%s/' % (input_full_sample_dirpath, patient)
		
		if abr_lineage:
			for gene_pair in os.listdir('%s%s' % (input_lineage_fasta_dirpath, patient)):
				if gene_pair[0] == '.' or gene_pair[:6] == 'README':
					continue

				# if gene_pair != 'IGHV6-1_IGHJ5':
				# 	continue

				gene_pair_input_dirpath = '%s%s/%s/' % (input_lineage_fasta_dirpath, patient, gene_pair)
				for attribute, output_dirpath in itertools.izip(attibute_names, output_trajectory_dirpaths):
					if output_dirpath[-1] != '/':
						output_dirpath += '/'
					output_gene_pair_dirpath = '%s%s/%s/' % (output_dirpath, patient, gene_pair)
					if os.path.exists(output_gene_pair_dirpath):
						subprocess.call(['rm', '-r', output_gene_pair_dirpath])
					os.makedirs(output_gene_pair_dirpath)
					seq_lineages = seq_lineage_time_series(dirpath=gene_pair_input_dirpath, count_attribute_name=count_attribute_name, list_of_filepaths=None, avoid_dirpath=None, avoid_timepoints_with_dot=avoid_timepoints_with_dot)
					seq_lineages.get_mean_attribute(full_seq_sample_dirpath=patient_full_sample_dirpath, attribute_name=attribute, output_dirpath=output_gene_pair_dirpath, temp_dirpath=temp_dirpath, weight_by_counts=weight_by_counts, treat_NAs=treat_NAs, min_freq_cutoff=min_freq_cutoff)

		else:
			patient_input_dirpath = '%s%s/' % (input_lineage_fasta_dirpath, patient)
			for attribute, output_dirpath in itertools.izip(attibute_names, output_trajectory_dirpaths):
				if output_dirpath[-1] != '/':
					output_dirpath += '/'
				output_patient_dirpath = '%s%s/' % (output_dirpath, patient)
				if os.path.exists(output_patient_dirpath):
					subprocess.call(['rm', '-r', output_patient_dirpath])
				os.makedirs(output_patient_dirpath)
				seq_lineages = seq_lineage_time_series(dirpath=patient_input_dirpath, count_attribute_name=count_attribute_name, list_of_filepaths=None, avoid_dirpath=None, avoid_timepoints_with_dot=avoid_timepoints_with_dot)
				seq_lineages.get_mean_attribute(full_seq_sample_dirpath=patient_full_sample_dirpath, attribute_name=attribute, output_dirpath=output_patient_dirpath, temp_dirpath=temp_dirpath, weight_by_counts=weight_by_counts, treat_NAs=treat_NAs, min_freq_cutoff=min_freq_cutoff)

	return

if __name__ == '__main__':

	#below are hard coded examples of how we used this code. One will need to replace these files with there own data if they want to use this code.
	
	#for HIV
	#get_mean_trajectories(input_lineage_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/lineage_fasta_files/hiv/edit_dist_between_rep_seqs', input_full_sample_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_fasta_files', attibute_names=['dN_dS', 'divergence_syn', 'divergence_non_syn'], output_trajectory_dirpaths=['/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/selection_lineages/hiv/edit_dist_1/dN_dS_by_counting_changes_relative_to_most_numerouse_1st_tpoint_seq', '/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/divergence_lineages/hiv/edit_dist_1/counting_relative_to_most_numerous_1st_tpoint_seq/syn', '/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/divergence_lineages/hiv/edit_dist_1/counting_relative_to_most_numerous_1st_tpoint_seq/non_syn'])
	#for AbR, unique VJ gene pairs
	# get_mean_trajectories(input_lineage_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/abr_seq_unique_VJ/changeo', input_full_sample_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/changeo_annotated_fasta_files', attibute_names=['CDR_baseline_selection', 'FWR_baseline_selection', 'total_mut_freq'], output_trajectory_dirpaths=['/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/selection_lineages/abr/unique_gene_pairs/baseline/mean_across_indiv_seqs/cdr', '/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/selection_lineages/abr/unique_gene_pairs/baseline/mean_across_indiv_seqs/fwr', '/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/divergence_lineages/abr/unique_gene_pairs/mean_divergence'])
	#for AbR, lineages
	within_dist = '6'
	across_dist = '35'
	# get_mean_trajectories(input_lineage_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/lineage_fasta_files/abr/edit_dist_between_rep_seqs/max_edit_dist_within_samps_%s_across_samps_%s' % (within_dist, across_dist), input_full_sample_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/changeo_annotated_fasta_files', attibute_names=['CDR_baseline_selection', 'FWR_baseline_selection', 'total_mut_freq'], output_trajectory_dirpaths=['/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/selection_lineages/abr/lineages_max_edit_dist_within_samps_%s_across_samps_%s/baseline/mean_across_indiv_seqs/cdr' % (within_dist, across_dist), '/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/selection_lineages/abr/lineages_max_edit_dist_within_samps_%s_across_samps_%s/baseline/mean_across_indiv_seqs/fwr' % (within_dist, across_dist), '/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/divergence_lineages/abr/lineages_max_edit_dist_within_samps_%s_across_samps_%s/mean_divergence' % (within_dist, across_dist)], abr_lineage=True)
