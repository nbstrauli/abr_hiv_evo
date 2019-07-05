import os
import sys
from sequence_sample_set_class import sequence_sample_set

def run(input_fasta_dirpath, output_sequence_cluster_dirpath):
	
	####### parameters #######
	count_attribute_name = 'DUPCOUNT'
	sample_clustering_method = 'by_edit_distance'
	max_edit_distance_within_sample = 4
	use_comp_cluster = True
	path_to_needle = '/netapp/home/nstrauli/tools_c/EMBOSS-6.6.0/emboss/needle'
	temp_dirpath = '/netapp/home/nstrauli/data/abr_hiv_coevo/temp_stuff'
	##########################

	if input_fasta_dirpath[-1] != '/':
		input_fasta_dirpath += '/'
	if output_sequence_cluster_dirpath[-1] != '/':
		output_sequence_cluster_dirpath += '/'
	if not os.path.exists(output_sequence_cluster_dirpath):
		os.makedirs(output_sequence_cluster_dirpath)

	#add the right base name to each output dirpath, depending on max edit dist parameters
	output_sequence_cluster_dirpath += 'max_edit_dist_%s/' % max_edit_distance_within_sample

	#cluster within samples
	seq_sample_set = sequence_sample_set(dirpath=input_fasta_dirpath, count_attribute_name=count_attribute_name)
	seq_sample_set.cluster_seqs_within_samples(output_dirpath=output_sequence_cluster_dirpath, method=sample_clustering_method, max_edit_distance=max_edit_distance_within_sample, use_comp_cluster=use_comp_cluster, wait_till_cluster_done=True, path_to_needle=path_to_needle, temp_dirpath=temp_dirpath)

	return

if __name__ == '__main__':

	#below are hard coded examples of how we used this code. One will need to replace these files with there own data if they want to use this code.
	
	#on cluster for lineage assignment
	# run(input_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_fasta_files', output_sequence_cluster_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/sequence_clusters/hiv/clustered_by_edit_distance/most_numerous_seq_foreach_cluster/for_lineage_assignment')
	#on cluster for outlier detection
	# run(input_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_annotated_fasta_files', output_sequence_cluster_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/sequence_clusters/hiv/clustered_by_edit_distance/most_numerous_seq_foreach_cluster/for_outlier_detection')
