import sys
sys.path.insert(0, './')
import os
from sequence_sample_set_class import sequence_sample_set
from sequence_sample_class import sequence_sample
import remove_outlier_clade_from_hiv_seqs

def run(input_clustered_samples_fasta_dirpath, input_indiv_seq_fasta_dirpath, output_network_filepath, outlier_filepath, output_msa_dirpath, output_newick_tree_dirpath, output_outliers_removed_fasta_dirpath, output_seqs_in_outlier_clade_dirpath):
	"""This script clusters all the HIV data across all the patient/time-points. It does this by first clustering within a sample, and then clustering across samples. The purpose is to see if some sequences cluster more closely with other than there assigned patient."""

	########### parameters ###########
	count_attribute_name = 'DUPCOUNT'
	sample_clustering_method = 'by_edit_distance'
	max_edit_distance_within_sample = 4
	max_edit_distance_across_samples = 4
	temp_dirpath = '/Users/nstrauli/data/abr_hiv_coevo/temp_stuff'
	use_comp_cluster = True
	# path_to_needle = '/netapp/home/nstrauli/tools_c/EMBOSS-6.6.0/emboss/needle'
	path_to_needle = '/Users/nstrauli/tools/EMBOSS-6.6.0/emboss/needle'
	node_attributes_to_add = ['element_0_of_seq_id', 'element_1_of_seq_id', 'total_freq', 'count']
	names_for_node_attributes_in_header = ['patient_ID', 'time_point', 'frequency', 'count']
	outlier_freq_cutoff = 0.001
	freq_attribute_name = 'total_freq'
	########### parameters ###########

	if input_clustered_samples_fasta_dirpath[-1] != '/':
		input_clustered_samples_fasta_dirpath += '/'
	if input_indiv_seq_fasta_dirpath[-1] != '/':
		input_indiv_seq_fasta_dirpath += '/'
	if output_network_filepath[-1] != '/':
		output_network_filepath += '/'
	if outlier_filepath[-1] != '/':
		outlier_filepath += '/'
	if output_msa_dirpath[-1] != '/':
		output_msa_dirpath += '/'
	if output_newick_tree_dirpath[-1] != '/':
		output_newick_tree_dirpath += '/'
	if output_outliers_removed_fasta_dirpath[-1] != '/':
		output_outliers_removed_fasta_dirpath += '/'
	if output_seqs_in_outlier_clade_dirpath[-1] != '/':
		output_seqs_in_outlier_clade_dirpath += '/'
	#add appropriate basename to each of the file/dir paths
	base_name = 'max_edit_dist_witin_samp_%s_across_samps_%s' % (max_edit_distance_within_sample, max_edit_distance_across_samples)
	input_clustered_samples_fasta_dirpath += 'max_edit_dist_%s/' % max_edit_distance_within_sample
	output_network_filepath += base_name + '.sif'
	output_node_attribute_filepath = output_network_filepath[:-4] + '_node_attributes.txt'
	outlier_filepath += 'outliers_%s.fasta' % base_name
	output_msa_dirpath += base_name + '/'
	output_newick_tree_dirpath += base_name + '/'
	output_seqs_in_outlier_clade_dirpath += base_name + '/'
	if not os.path.exists(input_clustered_samples_fasta_dirpath):
		os.makedirs(input_clustered_samples_fasta_dirpath)
	if not os.path.exists(output_msa_dirpath):
		os.makedirs(output_msa_dirpath)
	if not os.path.exists(output_newick_tree_dirpath):
		os.makedirs(output_newick_tree_dirpath)
	if not os.path.exists(output_outliers_removed_fasta_dirpath):
		os.makedirs(output_outliers_removed_fasta_dirpath)
	if not os.path.exists(output_seqs_in_outlier_clade_dirpath):
		os.makedirs(output_seqs_in_outlier_clade_dirpath)
	
	#cluster across samples
	seq_sample_set = sequence_sample_set(dirpath=input_clustered_samples_fasta_dirpath, count_attribute_name=count_attribute_name)
	seq_sample_set.cluster_seqs_across_samples(output_network_filepath=output_network_filepath, output_node_attribute_filepath=output_node_attribute_filepath, node_attributes_to_add=node_attributes_to_add, names_for_node_attributes_in_header=names_for_node_attributes_in_header, temp_dirpath=temp_dirpath, add_string_to_seq_ids='dirname_and_filename', add_to_start_or_end='start', max_edit_distance=max_edit_distance_across_samples, path_to_needle=path_to_needle, id_outliers=outlier_filepath, outlier_def='clusters_with_dif_directory', outlier_freq_cutoff=outlier_freq_cutoff, freq_attribute_name=freq_attribute_name)

	#label outliers seqs in sequence clusters data
	filein = open(outlier_filepath, "r")
	sample_filepath_outlier_dic_seq_clusts = {}
	seq_sample_set_seq_clusts = sequence_sample_set(dirpath=input_clustered_samples_fasta_dirpath, count_attribute_name=count_attribute_name)
	for i in seq_sample_set_seq_clusts.sample_filepaths:
		sample_filepath_outlier_dic_seq_clusts[i] = [] #initialize the dic
	for i in filein:
		if i[0] == '>':
			header = i[1:-1].split('|')
			seq_data = header[0].split('_')
			patient_ID = seq_data[0]
			tpoint = seq_data[1]
			seq_id = seq_data[2]
			#append cluster ID to cluster outlier dic
			sample_filepath_seq_clusts = '%s%s/%s.fasta' % (input_clustered_samples_fasta_dirpath, patient_ID, tpoint)
			sample_filepath_outlier_dic_seq_clusts[sample_filepath_seq_clusts].append(seq_id)
	filein.close()
	#lable seqs in seq clusters data
	for i in sample_filepath_outlier_dic_seq_clusts:
		sample = sequence_sample(filepath=i, count_attribute_name=count_attribute_name)
		sample.add_boolean_attribute_to_seqs(attribute_name='is_outlier', seq_ids_that_are_True=sample_filepath_outlier_dic_seq_clusts[i], overwrite_existing=True)
		sample.write_full_data_to_disk_fasta(output_filepath=i, append_to_file=False, seq_id_fileout_dic=None)

	#remove outlier seqs and seqs that are in 'outlier clade' using a phylogenetic approach
	remove_outlier_clade_from_hiv_seqs.run(input_seq_clusters_dirpath=input_clustered_samples_fasta_dirpath, output_msa_dirpath=output_msa_dirpath, output_newick_tree_dirpath=output_newick_tree_dirpath, input_indiv_seq_fasta_dirpath=input_indiv_seq_fasta_dirpath, output_outliers_removed_fasta_dirpath=output_outliers_removed_fasta_dirpath, output_outlier_seqs_dirpath=output_seqs_in_outlier_clade_dirpath)

	return

if __name__ == '__main__':

	#below are hard coded examples of how we used this code. One will need to replace these files with there own data if they want to use this code.
	
	#on local computer
	# run(input_clustered_samples_fasta_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/sequence_clusters/hiv/clustered_by_edit_distance/most_numerous_seq_foreach_cluster/for_outlier_detection', input_indiv_seq_fasta_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_annotated_fasta_files', output_network_filepath='/Users/nstrauli/data/abr_hiv_coevo/network_files/hiv/by_edit_distance_pooled_all_data', outlier_filepath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/outlier_seqs/hiv/by_clustering_with_other_samples', output_msa_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/hiv/alignments_of_seq_clusters/for_outlier_detection', output_newick_tree_dirpath='/Users/nstrauli/data/abr_hiv_coevo/phylogenetic_trees/hiv/fasttree/for_outlier_detection', output_outliers_removed_fasta_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_fasta_files', output_seqs_in_outlier_clade_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/outlier_seqs/hiv/by_being_in_clade_with_miss_IDed_barcode_seqs')
	#on cluster
	# run(input_clustered_samples_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/sequence_clusters/hiv/clustered_by_edit_distance/most_numerous_seq_foreach_cluster/for_outlier_detection', input_indiv_seq_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_annotated_fasta_files', output_network_filepath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/network_files/hiv/by_edit_distance_pooled_all_data', outlier_filepath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/outlier_seqs/hiv/by_clustering_with_other_samples', output_msa_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/hiv/alignments_of_seq_clusters/for_outlier_detection', output_newick_tree_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/phylogenetic_trees/hiv/fasttree/for_outlier_detection', output_outliers_removed_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_fasta_files', output_seqs_in_outlier_clade_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/outlier_seqs/hiv/by_being_in_clade_with_miss_IDed_barcode_seqs')
