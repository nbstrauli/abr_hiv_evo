import os
import sys
from sequence_sample_set_class import sequence_sample_set

def run(seq_input_dirpath, immune_reads_output_dirpath):
	"""
	Uses the 'sequence_sample_set' class to get the immune reads from a sequence sample set
	"""
	###### parameters ######
	count_attribute_name = 'DUPCOUNT'
	path_to_changeo = '/netapp/home/nstrauli/tools_c/changeo-0.3.9'
	#path_to_changeo = '/Users/nstrauli/tools/changeo-0.3.9'
	path_to_igblast = '/netapp/home/nstrauli/tools_c/ncbi-igblast-1.8.0'
	#path_to_igblast = '/Users/nstrauli/tools/ncbi-igblast-1.8.0'
	temp_dirpath = '/netapp/home/nstrauli/data/abr_hiv_coevo/temp_stuff/'
	#temp_dirpath = '/Users/nstrauli/data/abr_hiv_coevo/temp_stuff/'
	add_germline = ['dmask', 'full']
	add_selection = True
	use_comp_cluster = True
	add_divergence = 'python'
	overwrite_output = True
	translate_VDJ = True
	remove_seqs_with_stop = True
	########################

	sample_set = sequence_sample_set(seq_input_dirpath, count_attribute_name=count_attribute_name)
	sample_set.get_immune_reads_with_changeo_foreach_sample(output_dirpath=immune_reads_output_dirpath, path_to_changeo=path_to_changeo, path_to_igblast=path_to_igblast, temp_dirpath=temp_dirpath, add_germline=add_germline, add_selection=add_selection, use_comp_cluster=use_comp_cluster, add_divergence=add_divergence, overwrite_output=overwrite_output, translate_VDJ=translate_VDJ, remove_seqs_with_stop=remove_seqs_with_stop)
	return

if __name__ == '__main__':

	#below are hard coded examples of how we used this code. One will need to replace these files with there own data if they want to use this code.
	
	# run(seq_input_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/QCed_data_clean_final/abr/', immune_reads_output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/changeo_annotated_fasta_files/')
	#run(seq_input_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/temp_stuff/abr/', immune_reads_output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/temp_stuff/changeo_annotated_fasta_files/')
	#run(seq_input_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/QCed_data_clean_final/abr/', immune_reads_output_dirpath='/users/nstrauli/data/abr_hiv_coevo/seq_data/changeo_annotated_fasta_files')
	#run(seq_input_dirpath='/Users/nstrauli/data/abr_hiv_coevo/temp_stuff/abr/', immune_reads_output_dirpath='/Users/nstrauli/data/abr_hiv_coevo/temp_stuff/changeo_annotated_fasta_files/')

