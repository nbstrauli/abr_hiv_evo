import os
from seq_lineage_time_series_class import seq_lineage_time_series
import subprocess

def get_freq_trajs(input_lineage_dirpath, input_full_sample_dirpath, output_freq_traj_dirpath, abr_lineage=False):

	###### parameters ######
	count_attribute_name = 'DUPCOUNT'
	avoid_timepoints_with_dot = True
	if abr_lineage:
		min_freq_cutoff = 0.0
	else:
		min_freq_cutoff = 0.01
	########################

	if input_lineage_dirpath[-1] != '/':
		input_lineage_dirpath += '/'
	if input_full_sample_dirpath[-1] != '/':
		input_full_sample_dirpath += '/'
	if output_freq_traj_dirpath[-1] != '/':
		output_freq_traj_dirpath +='/'
	for patient in os.listdir(input_lineage_dirpath):
		if patient[0] == '.' or patient[:6] == 'README':
			continue
		print ''
		print 'patient:', patient

		# if patient != '6':
		# 	continue

		patient_full_sample_dirpath = '%s%s/' % (input_full_sample_dirpath, patient)

		if abr_lineage:
			for gene_pair in os.listdir('%s%s' % (input_lineage_dirpath, patient)):
				if gene_pair[0] == '.' or gene_pair[:6] == 'README':
					continue

				# if gene_pair != 'IGHV6-1_IGHJ5':
				# 	continue

				gene_pair_input_dirpath = '%s%s/%s/' % (input_lineage_dirpath, patient, gene_pair)
				output_gene_pair_dirpath = '%s%s/%s/' % (output_freq_traj_dirpath, patient, gene_pair)
				if os.path.exists(output_gene_pair_dirpath):
					subprocess.call(['rm', '-r', output_gene_pair_dirpath])
				os.makedirs(output_gene_pair_dirpath)
				seq_lins = seq_lineage_time_series(dirpath=gene_pair_input_dirpath, count_attribute_name=count_attribute_name, list_of_filepaths=None, avoid_dirpath=None, avoid_timepoints_with_dot=avoid_timepoints_with_dot)
				seq_lins.get_freq_trajectories(full_seq_sample_dirpath=patient_full_sample_dirpath, output_dirpath=output_gene_pair_dirpath, min_freq_cutoff=min_freq_cutoff)

		else:
			patient_input_dirpath = '%s%s/' % (input_lineage_dirpath, patient)
			output_patient_dirpath = '%s%s/' % (output_freq_traj_dirpath, patient)
			if os.path.exists(output_patient_dirpath):
				subprocess.call(['rm', '-r', output_patient_dirpath])
			os.makedirs(output_patient_dirpath)
			seq_lins = seq_lineage_time_series(dirpath=patient_input_dirpath, count_attribute_name=count_attribute_name, list_of_filepaths=None, avoid_dirpath=None, avoid_timepoints_with_dot=avoid_timepoints_with_dot)
			seq_lins.get_freq_trajectories(full_seq_sample_dirpath=patient_full_sample_dirpath, output_dirpath=output_patient_dirpath, min_freq_cutoff=min_freq_cutoff)

	return

if __name__ == '__main__':

	#below are hard coded examples of how we used this code. One will need to replace these files with there own data if they want to use this code.
	
	#for HIV
	#get_freq_trajs(input_lineage_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/lineage_fasta_files/hiv/edit_dist_between_rep_seqs', input_full_sample_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_fasta_files', output_freq_traj_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/freq_trajectory_lineages/hiv/edit_dist_1')
	#for AbR, unique VJ gene pairs
	# get_freq_trajs(input_lineage_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/abr_seq_unique_VJ/changeo', input_full_sample_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/changeo_annotated_fasta_files', output_freq_traj_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/freq_trajectory_lineages/abr/unique_gene_pairs')
	#for AbR, lineages
	# get_freq_trajs(input_lineage_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/lineage_fasta_files/abr/edit_dist_between_rep_seqs/max_edit_dist_within_samps_6_across_samps_35', input_full_sample_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/changeo_annotated_fasta_files', output_freq_traj_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/freq_trajectory_lineages/abr/lineages_max_edit_dist_within_samps_6_across_samps_35', abr_lineage=True)
