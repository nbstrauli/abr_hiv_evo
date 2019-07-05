import os
from immune_time_series_class import immune_time_series
import subprocess

def cluster(input_seqs_dirpath, output_VJ_seqs_dirpath):
	"""This script uses 'immune_time_series_class.py' to go through the immune seqs and partition them into groups based upon their germline V and J gene assignment."""

	###### parameters ######
	count_attribute_name = 'DUPCOUNT'
	vgene_name = 'V_CALL'
	dgene_name = 'D_CALL'
	jgene_name = 'J_CALL'
	cdr3_name = 'CDR3_IMGT'
	ignor_allele_info = True
	use_cluster = True
	replace_slash_with = '.'
	###### parameters ######

	if input_seqs_dirpath[-1] != '/':
		input_seqs_dirpath += '/'
	if output_VJ_seqs_dirpath[-1] != '/':
		output_VJ_seqs_dirpath += '/'
	if not os.path.exists(output_VJ_seqs_dirpath):
		os.makedirs(output_VJ_seqs_dirpath)
	for patient in os.listdir(input_seqs_dirpath):
		if patient[0] == '.' or patient[:6] == 'README' or patient == 'unknown':
			continue
		
		print '####################'
		print patient
		print '####################'
		print ''

		input_seqs_patient_dirpath = input_seqs_dirpath + patient + '/'
		output_VJ_seqs_patient_dirpath = output_VJ_seqs_dirpath + patient + '/'
		if os.path.exists(output_VJ_seqs_patient_dirpath):
			subprocess.call(['rm', '-r', output_VJ_seqs_patient_dirpath])
		os.makedirs(output_VJ_seqs_patient_dirpath)

		imm_t_series = immune_time_series(dirpath=input_seqs_patient_dirpath, count_attribute_name=count_attribute_name, vgene_name=vgene_name, dgene_name=dgene_name, jgene_name=jgene_name, cdr3_name=cdr3_name, ignor_allele_info=ignor_allele_info, replace_slash_with=replace_slash_with)
		imm_t_series.divide_data_by_VJ_combo(output_dirpath=output_VJ_seqs_patient_dirpath, use_cluster=use_cluster)

	return

if __name__ == '__main__':

	#below are hard coded examples of how we used this code. One will need to replace these files with there own data if they want to use this code.
	
	# cluster(input_seqs_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/changeo_annotated_fasta_files/', output_VJ_seqs_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/abr_seq_unique_VJ/changeo')
	# cluster(input_seqs_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/changeo_annotated_fasta_files', output_VJ_seqs_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/abr_seq_unique_VJ/changeo')
