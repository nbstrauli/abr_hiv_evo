import os
from compare_sets_of_MSAs_class import compare_sets_of_MSAs
import time
import subprocess

def run(input_hiv_MSA_dirpath, input_abr_MSA_dirpath, output_dirpath, hiv_is_whole_pop, abr_is_whole_pop, abr_partition=False):

	##### parameteres #####
	min_num_of_seq_entries = 7
	num_permute = 10
	ignore_sites_less_than_X_changes = 2
	# temp_dirpath = '/scrapp'
	temp_dirpath = '/netapp/home/nstrauli/data/abr_hiv_coevo/temp_stuff'
	time_to_sleep = 100
	#######################

	if input_hiv_MSA_dirpath[-1] != '/':
		input_hiv_MSA_dirpath += '/'
	if input_abr_MSA_dirpath[-1] != '/':
		input_abr_MSA_dirpath += '/'
	if output_dirpath[-1] != '/':
		output_dirpath += '/'
	if hiv_is_whole_pop:
		MSA_1_is_filepath = True
	else:
		MSA_1_is_filepath = False
	if abr_is_whole_pop:
		MSA_2_is_filepath = True
	else:
		MSA_2_is_filepath = False
	for patient in os.listdir(input_abr_MSA_dirpath):
		if patient[0] == '.' or patient[:6] == 'README':
			continue
		print patient

		# if patient != '6' and patient != '7' and patient != '8' and patient != '9' and patient != '10':
			# continue

		hiv_MSAs_dirpath = '%s%s.fasta' % (input_hiv_MSA_dirpath, patient)
		abr_MSAs_dirpath = '%s%s/' % (input_abr_MSA_dirpath, patient)
		output_patient_dirpath = '%s%s/' % (output_dirpath, patient)
		if os.path.exists(output_patient_dirpath):
			subprocess.call(['rm', '-r', output_patient_dirpath])
		os.makedirs(output_patient_dirpath)
		output_patient_dirpath_change_or_same = '%s_change_or_same/%s/' % (output_dirpath[:-1], patient)
		if os.path.exists(output_patient_dirpath_change_or_same):
			subprocess.call(['rm', '-r', output_patient_dirpath_change_or_same])
		os.makedirs(output_patient_dirpath_change_or_same)
		if abr_partition:
			MSA_set_comp = compare_sets_of_MSAs(dirpath_1=hiv_MSAs_dirpath, dirpath_2=abr_MSAs_dirpath, count_attribute_name='DUPCOUNT', MSA_1_is_filepath=MSA_1_is_filepath, MSA_2_is_filepath=MSA_2_is_filepath, MSA_1_is_grouped_into_subdirs=False, MSA_2_is_grouped_into_subdirs=False)
			MSA_set_comp_change_or_same = compare_sets_of_MSAs(dirpath_1=hiv_MSAs_dirpath, dirpath_2=abr_MSAs_dirpath, count_attribute_name='DUPCOUNT', MSA_1_is_filepath=MSA_1_is_filepath, MSA_2_is_filepath=MSA_2_is_filepath, MSA_1_is_grouped_into_subdirs=False, MSA_2_is_grouped_into_subdirs=False, reduce_to_change_or_same=True)
		else:
			MSA_set_comp = compare_sets_of_MSAs(dirpath_1=hiv_MSAs_dirpath, dirpath_2=abr_MSAs_dirpath, count_attribute_name='DUPCOUNT', MSA_1_is_filepath=MSA_1_is_filepath, MSA_2_is_filepath=MSA_2_is_filepath, MSA_1_is_grouped_into_subdirs=False, MSA_2_is_grouped_into_subdirs=True)
			MSA_set_comp_change_or_same = compare_sets_of_MSAs(dirpath_1=hiv_MSAs_dirpath, dirpath_2=abr_MSAs_dirpath, count_attribute_name='DUPCOUNT', MSA_1_is_filepath=MSA_1_is_filepath, MSA_2_is_filepath=MSA_2_is_filepath, MSA_1_is_grouped_into_subdirs=False, MSA_2_is_grouped_into_subdirs=True, reduce_to_change_or_same=True)
		print '\tcalculating MI'
		MSA_set_comp.calc_MI_foreach_MSA_pair(output_array_dirpath=None, output_distribution_dirpath=output_patient_dirpath, ignore_sites_less_than_X_changes=ignore_sites_less_than_X_changes, min_num_of_seq_entries=min_num_of_seq_entries, perform_permutation_null=True, num_permute=num_permute, permutation_null_comparison_results_dirpath=output_patient_dirpath, use_comp_cluster=True, temp_dirpath=temp_dirpath)
		print '\tcalculating MI, for change or same data'
		MSA_set_comp_change_or_same.calc_MI_foreach_MSA_pair(output_array_dirpath=None, output_distribution_dirpath=output_patient_dirpath_change_or_same, ignore_sites_less_than_X_changes=ignore_sites_less_than_X_changes, min_num_of_seq_entries=min_num_of_seq_entries, perform_permutation_null=True, num_permute=num_permute, permutation_null_comparison_results_dirpath=output_patient_dirpath_change_or_same, use_comp_cluster=True, temp_dirpath=temp_dirpath)

		# time.sleep(time_to_sleep)
		
	return

if __name__ == '__main__':

	#below are hard coded examples of how we used this code. One will need to replace these files with there own data if they want to use this code.

	#for HIV AA population against AbR AA lineages
	# run(input_hiv_MSA_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/hiv/time_series_alignments_of_representative_AA_seqs_of_wholePopulation/edit_dist_between_rep_seqs', input_abr_MSA_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/abr/time_series_alignments_of_representative_AA_seqs_of_lineages/max_edit_dist_within_samps_6_across_samps_30', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/MI_values/hiv_AA_pop_against_abr_AA_lineages', hiv_is_whole_pop=True, abr_is_whole_pop=False)

	#for HIV AA population against AbR AA lineages, null simulation
	# run(input_hiv_MSA_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/hiv/time_series_alignments_of_representative_AA_seqs_of_wholePopulation/edit_dist_between_rep_seqs', input_abr_MSA_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments_sim_null/abr/time_series_alignments_of_representative_AA_seqs_of_lineages/max_edit_dist_within_samps_6_across_samps_30', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/MI_values_null_sim/hiv_AA_pop_against_abr_AA_lineages', hiv_is_whole_pop=True, abr_is_whole_pop=False)

	#for HIV AA population against AbR AA lineages, select
	# run(input_hiv_MSA_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/hiv/time_series_alignments_of_representative_AA_seqs_of_wholePopulation/edit_dist_between_rep_seqs', input_abr_MSA_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/abr/time_series_alignments_of_representative_AA_seqs_of_lineages_select/max_edit_dist_within_samps_6_across_samps_30', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/MI_values/hiv_AA_pop_against_abr_AA_lineages_select', hiv_is_whole_pop=True, abr_is_whole_pop=False)

	#for HIV AA population against AbR AA lineages, null simulation, select
	# run(input_hiv_MSA_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/hiv/time_series_alignments_of_representative_AA_seqs_of_wholePopulation/edit_dist_between_rep_seqs', input_abr_MSA_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments_sim_null/abr/time_series_alignments_of_representative_AA_seqs_of_lineages_select/max_edit_dist_within_samps_6_across_samps_30', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/MI_values_null_sim/hiv_AA_pop_against_abr_AA_lineages_select', hiv_is_whole_pop=True, abr_is_whole_pop=False)

	#for HIV nt population against AbR nt lineages
	# run(input_hiv_MSA_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/hiv/time_series_alignments_of_representative_seqs_of_wholePopulation/edit_dist_between_rep_seqs/', input_abr_MSA_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/abr/time_series_alignments_of_representative_seqs_of_lineages/max_edit_dist_within_samps_6_across_samps_30/', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/MI_values/hiv_nt_pop_against_abr_nt_lineages', hiv_is_whole_pop=True, abr_is_whole_pop=False)

	#for HIV AA population against AbR AA CDR3 lineages
	# run(input_hiv_MSA_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/hiv/time_series_alignments_of_representative_AA_seqs_of_wholePopulation/edit_dist_between_rep_seqs', input_abr_MSA_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/abr/cdr3_time_series_alignments_of_representative_AA_seqs_of_lineages/max_edit_dist_within_samps_6_across_samps_30', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/MI_values/hiv_AA_pop_against_abr_AA_CDR3_lineages', hiv_is_whole_pop=True, abr_is_whole_pop=False)

	#for HIV nt population against AbR nt CDR3 lineages
	# run(input_hiv_MSA_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/hiv/time_series_alignments_of_representative_seqs_of_wholePopulation/edit_dist_between_rep_seqs/', input_abr_MSA_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/abr/cdr3_time_series_alignments_of_representative_seqs_of_lineages/max_edit_dist_within_samps_6_across_samps_30', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/MI_values/hiv_nt_pop_against_abr_nt_CDR3_lineages', hiv_is_whole_pop=True, abr_is_whole_pop=False)

	#for HIV population against AbR partitions
	# run(input_hiv_MSA_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/hiv/time_series_alignments_of_representative_AA_seqs_of_wholePopulation/edit_dist_between_rep_seqs', input_abr_MSA_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/abr/time_series_alignments_of_representative_AA_seqs_of_VJgenePair_partitions', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/MI_values/hiv_pop_against_abr_partitions', hiv_is_whole_pop=True, abr_is_whole_pop=False, abr_partition=True)
