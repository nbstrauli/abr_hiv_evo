import os
from seq_lineage_time_series_class import seq_lineage_time_series
import time
import subprocess

def run(input_lineage_dirpath, output_alignment_dirpath, output_AA_alignment_dirpath, output_alignment_dirpath_cdr3, output_AA_alignment_dirpath_cdr3, type_of_data):

	##### parameters #####
	# temp_dirpath = '/scrapp'
	temp_dirpath = '/netapp/home/nstrauli/data/abr_hiv_coevo/temp_stuff/'
	path_to_mafft = '/netapp/home/nstrauli/bin/mafft'
	time_to_sleep = 1000
	min_num_tpoints = 7
	######################

	if input_lineage_dirpath[-1] != '/':
		input_lineage_dirpath += '/'
	if output_alignment_dirpath[-1] != '/':
		output_alignment_dirpath += '/'
	if output_AA_alignment_dirpath[-1] != '/':
		output_AA_alignment_dirpath += '/'
	if output_alignment_dirpath_cdr3:
		if output_alignment_dirpath_cdr3[-1] != '/':
			output_alignment_dirpath_cdr3 += '/'
	if output_AA_alignment_dirpath_cdr3:
		if output_AA_alignment_dirpath_cdr3[-1] != '/':
			output_AA_alignment_dirpath_cdr3 += '/'
	if type_of_data == 'hiv_population':
		patients = ['nuh_bro']
	else:
		patients = os.listdir(input_lineage_dirpath)
	for patient in patients:
		if patient[0] == '.' or patient[:6] == 'README':
			continue
		print patient

		# if patient != '1':
		# 	continue

		if type_of_data == 'hiv_lineage' or type_of_data == 'hiv_population' or type_of_data == 'abr_partition':
			gene_pairs = ['nuh_bra']
		elif type_of_data == 'abr_lineage':
			# gene_pairs = os.listdir('%s%s/' % (input_lineage_dirpath, patient))
			gene_pairs = ['nuh_bra']
		for gene_pair in gene_pairs:
			if gene_pair[0] == '.' or gene_pair[:6] == 'README':
				continue

			# if gene_pair != 'IGHV3-15_IGHJ4':
			# 	continue

			if type_of_data == 'abr_lineage':
				input_dirpath = '%s%s/' % (input_lineage_dirpath, patient)
				output_dirpath = '%s%s/' % (output_alignment_dirpath, patient)
				output_dirpath_AA = '%s%s/' % (output_AA_alignment_dirpath, patient)
				output_dirpath_cdr3 = '%s%s/' % (output_alignment_dirpath_cdr3, patient)
				output_dirpath_cdr3_AA = '%s%s/' % (output_AA_alignment_dirpath_cdr3, patient)
				wait_till_cluster_done = False
				AA_seq_attribute_name = 'VDJ_AA_seq'
				cdr3_AA_seq_attribute_name = 'CDR3_IGBLAST_AA'
				lineages_grouped_into_dirs = True
			elif type_of_data == 'abr_partition':
				input_dirpath = '%s%s/' % (input_lineage_dirpath, patient)
				output_dirpath = '%s%s/' % (output_alignment_dirpath, patient)
				output_dirpath_AA = '%s%s/' % (output_AA_alignment_dirpath, patient)
				output_dirpath_cdr3 = '%s%s/' % (output_alignment_dirpath_cdr3, patient)
				output_dirpath_cdr3_AA = '%s%s/' % (output_AA_alignment_dirpath_cdr3, patient)
				wait_till_cluster_done = False
				AA_seq_attribute_name = 'VDJ_AA_seq'
				cdr3_AA_seq_attribute_name = 'CDR3_IGBLAST_AA'
				lineages_grouped_into_dirs = False
			elif type_of_data == 'hiv_lineage':
				input_dirpath = '%s%s/' % (input_lineage_dirpath, patient)
				output_dirpath = '%s%s/' % (output_alignment_dirpath, patient)
				output_dirpath_AA = '%s%s/' % (output_AA_alignment_dirpath, patient)
				wait_till_cluster_done = False
				AA_seq_attribute_name = 'amino_acid_seq'
				lineages_grouped_into_dirs = False
			elif type_of_data == 'hiv_population':
				input_dirpath = input_lineage_dirpath
				output_dirpath = output_alignment_dirpath
				output_dirpath_AA = output_AA_alignment_dirpath
				wait_till_cluster_done = False
				AA_seq_attribute_name = 'amino_acid_seq'
				lineages_grouped_into_dirs = False

			seq_lineage = seq_lineage_time_series(dirpath=input_dirpath, count_attribute_name='DUPCOUNT', list_of_filepaths=None, avoid_dirpath=None, avoid_timepoints_with_dot=None, lineages_grouped_into_dirs=lineages_grouped_into_dirs, min_num_tpoints=min_num_tpoints)
			#align nucleotide seqs
			if os.path.exists(output_dirpath):
				subprocess.call(['rm', '-r', output_dirpath])
			os.makedirs(output_dirpath)
			seq_lineage.make_time_series_alingments(output_dirpath=output_dirpath, path_to_mafft=path_to_mafft, temp_dirpath=temp_dirpath, align_attribute_value_instead=None, wait_till_cluster_done=wait_till_cluster_done)
			#align amino acid seqs
			if os.path.exists(output_dirpath_AA):
				subprocess.call(['rm', '-r', output_dirpath_AA])
			os.makedirs(output_dirpath_AA)
			seq_lineage.make_time_series_alingments(output_dirpath=output_dirpath_AA, path_to_mafft=path_to_mafft, temp_dirpath=temp_dirpath, align_attribute_value_instead=AA_seq_attribute_name, wait_till_cluster_done=wait_till_cluster_done)
			if type_of_data == 'abr_lineage' or type_of_data == 'abr_partition':
				#align nucleotide seqs
				if os.path.exists(output_dirpath_cdr3):
					subprocess.call(['rm', '-r', output_dirpath_cdr3])
				os.makedirs(output_dirpath_cdr3)
				seq_lineage.make_time_series_alingments(output_dirpath=output_dirpath_cdr3, path_to_mafft=path_to_mafft, temp_dirpath=temp_dirpath, align_attribute_value_instead='CDR3_IMGT', wait_till_cluster_done=wait_till_cluster_done)
				#align amino acid seqs
				if os.path.exists(output_dirpath_cdr3_AA):
					subprocess.call(['rm', '-r', output_dirpath_cdr3_AA])
				os.makedirs(output_dirpath_cdr3_AA)
				seq_lineage.make_time_series_alingments(output_dirpath=output_dirpath_cdr3_AA, path_to_mafft=path_to_mafft, temp_dirpath=temp_dirpath, align_attribute_value_instead=cdr3_AA_seq_attribute_name, wait_till_cluster_done=wait_till_cluster_done)

			# time.sleep(time_to_sleep)

	return

if __name__ == '__main__':

	#below are hard coded examples of how we used this code. One will need to replace these files with there own data if they want to use this code.
	
	#for AbR lineage data
	# run(input_lineage_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/lineage_fasta_files/abr/edit_dist_between_rep_seqs/max_edit_dist_within_samps_6_across_samps_30', output_alignment_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/abr/time_series_alignments_of_representative_seqs_of_lineages/max_edit_dist_within_samps_6_across_samps_30', output_AA_alignment_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/abr/time_series_alignments_of_representative_AA_seqs_of_lineages/max_edit_dist_within_samps_6_across_samps_30', output_alignment_dirpath_cdr3='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/abr/cdr3_time_series_alignments_of_representative_seqs_of_lineages/max_edit_dist_within_samps_6_across_samps_30', output_AA_alignment_dirpath_cdr3='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/abr/cdr3_time_series_alignments_of_representative_AA_seqs_of_lineages/max_edit_dist_within_samps_6_across_samps_30', type_of_data='abr_lineage')

	#for AbR lineage data, null simulation
	# run(input_lineage_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/lineage_fasta_files_sim_null/abr/edit_dist_between_rep_seqs/max_edit_dist_within_samps_6_across_samps_30', output_alignment_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments_sim_null_2/abr/time_series_alignments_of_representative_seqs_of_lineages/max_edit_dist_within_samps_6_across_samps_30', output_AA_alignment_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments_sim_null_2/abr/time_series_alignments_of_representative_AA_seqs_of_lineages/max_edit_dist_within_samps_6_across_samps_30', output_alignment_dirpath_cdr3='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments_sim_null_2/abr/cdr3_time_series_alignments_of_representative_seqs_of_lineages/max_edit_dist_within_samps_6_across_samps_30', output_AA_alignment_dirpath_cdr3='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments_sim_null_2/abr/cdr3_time_series_alignments_of_representative_AA_seqs_of_lineages/max_edit_dist_within_samps_6_across_samps_30', type_of_data='abr_lineage')

	#for AbR lineage select data
	# run(input_lineage_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/lineage_fasta_files/abr/edit_dist_between_rep_seqs_select/max_edit_dist_within_samps_6_across_samps_30', output_alignment_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/abr/time_series_alignments_of_representative_seqs_of_lineages_select/max_edit_dist_within_samps_6_across_samps_30', output_AA_alignment_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/abr/time_series_alignments_of_representative_AA_seqs_of_lineages_select/max_edit_dist_within_samps_6_across_samps_30', output_alignment_dirpath_cdr3='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/abr/cdr3_time_series_alignments_of_representative_seqs_of_lineages_select/max_edit_dist_within_samps_6_across_samps_30', output_AA_alignment_dirpath_cdr3='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/abr/cdr3_time_series_alignments_of_representative_AA_seqs_of_lineages_select/max_edit_dist_within_samps_6_across_samps_30', type_of_data='abr_lineage')

	#for AbR lineage data, null simulation, select data
	# run(input_lineage_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/lineage_fasta_files_sim_null/abr/edit_dist_between_rep_seqs_select/max_edit_dist_within_samps_6_across_samps_30', output_alignment_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments_sim_null/abr/time_series_alignments_of_representative_seqs_of_lineages_select/max_edit_dist_within_samps_6_across_samps_30', output_AA_alignment_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments_sim_null/abr/time_series_alignments_of_representative_AA_seqs_of_lineages_select/max_edit_dist_within_samps_6_across_samps_30', output_alignment_dirpath_cdr3='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments_sim_null/abr/cdr3_time_series_alignments_of_representative_seqs_of_lineages_select/max_edit_dist_within_samps_6_across_samps_30', output_AA_alignment_dirpath_cdr3='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments_sim_null/abr/cdr3_time_series_alignments_of_representative_AA_seqs_of_lineages_select/max_edit_dist_within_samps_6_across_samps_30', type_of_data='abr_lineage')

	#for AbR whole partition data
	# run(input_lineage_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/abr_seq_unique_VJ/changeo/', output_alignment_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/abr/time_series_alignments_of_representative_seqs_of_VJgenePair_partitions', output_AA_alignment_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/abr/time_series_alignments_of_representative_AA_seqs_of_VJgenePair_partitions', output_alignment_dirpath_cdr3='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/abr/cdr3_time_series_alignments_of_representative_seqs_of_VJgenePair_partitions', output_AA_alignment_dirpath_cdr3='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/abr/cdr3_time_series_alignments_of_representative_AA_seqs_of_VJgenePair_partitions', type_of_data='abr_partition')

	#for HIV lineage data
	# run(input_lineage_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/lineage_fasta_files/hiv/edit_dist_between_rep_seqs', output_alignment_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/hiv/time_series_alignments_of_representative_seqs_of_lineages/edit_dist_between_rep_seqs', output_AA_alignment_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/hiv/time_series_alignments_of_representative_AA_seqs_of_lineages/edit_dist_between_rep_seqs', output_alignment_dirpath_cdr3=None, output_AA_alignment_dirpath_cdr3=None, type_of_data='hiv_lineage')

	#for HIV whole population
	# run(input_lineage_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_annotated_fasta_files/', output_alignment_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/hiv/time_series_alignments_of_representative_seqs_of_wholePopulation/edit_dist_between_rep_seqs', output_AA_alignment_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/hiv/time_series_alignments_of_representative_AA_seqs_of_wholePopulation/edit_dist_between_rep_seqs', output_alignment_dirpath_cdr3=None, output_AA_alignment_dirpath_cdr3=None, type_of_data='hiv_population')
