import os
import sys
from sequence_time_series_class import sequence_time_series
import random

def run(input_fasta_dirpath, consensus_output_dirpath, fasta_output_dirpath):

	####### parameters #######
	count_attribute_name = 'DUPCOUNT'
	alignment_method = 'mafft'
	consensus_method = 'most_abundant'
	temp_dirpath = '/netapp/home/nstrauli/data/abr_hiv_coevo/temp_stuff'
	translate_by_ref_seq = '/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/ref_seqs/hiv/HXB2_amplicon.fasta'
	first_seq_ref = True
	remove_seqs_with_stop_codons = True
	divergence_method = 'counting'
	coding_frame_start = 'relative_to_reference'
	method_counting = 'max_synonymous'
	path_to_needle = '/netapp/home/nstrauli/tools_c/EMBOSS-6.6.0/emboss/needle'
	####### parameters #######

	if input_fasta_dirpath[-1] != '/':
		input_fasta_dirpath += '/'
	if consensus_output_dirpath[-1] != '/':
		consensus_output_dirpath += '/'
	if not os.path.exists(consensus_output_dirpath):
		os.makedirs(consensus_output_dirpath)
	if fasta_output_dirpath[-1] != '/':
		fasta_output_dirpath += '/'
	if not os.path.exists(fasta_output_dirpath):
		os.makedirs(fasta_output_dirpath)
	random_suffix = str(random.random())
	for i in os.listdir(input_fasta_dirpath):
		if i[0] == '.' or i[:6] == 'README':
			continue
		patient = i
		print '########################'
		print "running for patient:", patient
		print '########################'
		#make consensus seq from first time point, by getting the seq that is most abundant in that time-point
		consensus_output_filepath = consensus_output_dirpath + i + '.fasta'
		time_series_seq_data = sequence_time_series(dirpath=input_fasta_dirpath+i, count_attribute_name=count_attribute_name)
		consensus_seq = time_series_seq_data.get_most_abundant_seq_from_tpoint(tpoint_index=0, output_filepath=consensus_output_filepath, translate_by_ref_seq=translate_by_ref_seq, temp_dirpath=temp_dirpath, append_to_output=False, path_to_needle=path_to_needle)
		#patient 6 seems to have had a superinfection in the 2nd tpoint, so we are getting two consensus seqs for that patient
		if patient == '6':
			consensus_seq = time_series_seq_data.get_most_abundant_seq_from_tpoint(tpoint_index=1, output_filepath=consensus_output_filepath, translate_by_ref_seq=translate_by_ref_seq, temp_dirpath=temp_dirpath, append_to_output=True, path_to_needle=path_to_needle)
		fasta_output_patient_dirpath = fasta_output_dirpath + i + '/'
		#translate and compare each seq to consensus to get divergence (syn and non-syn)
		time_series_seq_data.translate_and_get_dN_dS_divergence(output_dirpath=fasta_output_patient_dirpath, ref_seq=consensus_output_filepath, coding_frame_start=coding_frame_start, temp_dirpath=temp_dirpath, remove_seqs_with_stop_codons=remove_seqs_with_stop_codons, method=divergence_method, method_counting=method_counting, path_to_needle=path_to_needle)
	return

if __name__ == '__main__':

	#below are hard coded examples of how we used this code. One will need to replace these files with there own data if they want to use this code.
	
	# run(input_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/QCed_data_clean_final/hiv', consensus_output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/1st_time_point_consensus_seqs/hiv/most_abundant_method_patient6', fasta_output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_annotated_fasta_files_patient6')
