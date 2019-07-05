import os
from sequence_time_series_class import sequence_time_series
from msa_class import msa
from phylogenetic_tree_class import phylo_tree
from sequence_sample_class import sequence_sample
import subprocess

def run(input_seq_clusters_dirpath, output_msa_dirpath, output_newick_tree_dirpath, input_indiv_seq_fasta_dirpath, output_outliers_removed_fasta_dirpath, output_outlier_seqs_dirpath):
	"""
	This script will take fasta files that are of clustered HIV data, and first make MSA's of them for each patient, and then make phylogeneic trees for each patient.
	"""

	########### parameters ###########
	count_attribute_name = 'DUPCOUNT'
	temp_dirpath = '/Users/nstrauli/data/abr_hiv_coevo/temp_stuff'
	alignment_method = 'mafft'
	ref_seq_dirpath = '/Users/nstrauli/data/abr_hiv_coevo/seq_data/1st_time_point_consensus_seqs/hiv/most_abundant_method/'
	phylo_method = 'fasttree'
	path_to_fasttree = '/Users/nstrauli/tools/fastTree/FastTree'
	nuc_or_aa = 'nucleotide'
	include_freq_info = True
	freq_attribute_name = 'total_freq'
	outgroup_label = 'outgroup'
	info_included_in_trees = ['is_outlier']
	#tree plotting parameters
	outlier_attribute_name = 'is_outlier'
	tree_plot_file_suffix = '_full_tree_half_circ_plot.pdf'
	time_series_info = 'start_of_id'
	tree_style = 'half_circle'
	leaf_size_map_to = 'freq'
	show_leaf_names = False
	color_branches_by = 'is_outlier'
	line_width = 3
	start_color = 'red'
	end_color = 'blue'
	ladderize = True
	########### parameters ###########

	if input_seq_clusters_dirpath[-1] != '/':
		input_seq_clusters_dirpath += '/'
	if output_msa_dirpath[-1] != '/':
		output_msa_dirpath += '/'
	if not os.path.exists(output_msa_dirpath):
		os.makedirs(output_msa_dirpath)
	if output_newick_tree_dirpath[-1] != '/':
		output_newick_tree_dirpath += '/'
	if not os.path.exists(output_newick_tree_dirpath):
		os.makedirs(output_newick_tree_dirpath)
	if input_indiv_seq_fasta_dirpath[-1] != '/':
		input_indiv_seq_fasta_dirpath += '/'
	if output_outliers_removed_fasta_dirpath[-1] != '/':
		output_outliers_removed_fasta_dirpath += '/'
	if not os.path.exists(output_outliers_removed_fasta_dirpath):
		os.makedirs(output_outliers_removed_fasta_dirpath)
	if output_outlier_seqs_dirpath[-1] != '/':
		output_outlier_seqs_dirpath += '/'
	if not os.path.exists(output_outlier_seqs_dirpath):
		os.makedirs(output_outlier_seqs_dirpath)
	for patient in os.listdir(input_seq_clusters_dirpath):
		if patient[0] == '.' or patient[:6] == 'README':
			continue

		print ''
		print "#####################"
		print "patient:", patient
		print "#####################"
		print ''

		patient_alignment_filepath = output_msa_dirpath + patient + '.fasta'
		patient_ref_seq_filepath = ref_seq_dirpath + patient + '.fasta'
		newick_output_filepath = '%s%s.tree' % (output_newick_tree_dirpath, patient)

		#get seq ID for the consensous seq for the patient
		filein = open(patient_ref_seq_filepath, "r")
		ref_seq_id = filein.readline()[1:].split('|')[0]
		filein.close()

		#make MSA
		sample_time_series = sequence_time_series(dirpath=input_seq_clusters_dirpath+patient, count_attribute_name=count_attribute_name)
		sample_time_series.make_MSA(output_filepath=patient_alignment_filepath, temp_dirpath=temp_dirpath, method=alignment_method, add_ref_seq=None, write_seq_id_first='%s_%s' % (sample_time_series.timepoints[0], ref_seq_id))

		#make tree in newick format
		time_series_msa = msa(filepath=patient_alignment_filepath, count_attribute_name=count_attribute_name)
		time_series_msa.make_phylo_trees(output_filepath=newick_output_filepath, method=phylo_method, path_to_fasttree=path_to_fasttree, nuc_or_aa=nuc_or_aa, temp_dirpath=temp_dirpath, include_freq_info=include_freq_info, freq_attribute_name=freq_attribute_name, first_seq_is_outgroup=True, outgroup_label=outgroup_label, include_attributes=info_included_in_trees)

		#plot tree, and get outlier clade
		tree = phylo_tree(filepath=newick_output_filepath, count_attribute_name=count_attribute_name, freq_attribute_name=freq_attribute_name, outgroup_label=outgroup_label)
		if patient == '5' or patient == '6':
			outlier_seq_ids = tree.find_outlier_clade(outlier_attribute_name=outlier_attribute_name, get_ancestor_of_outier_clade=True)
		else:
			outlier_seq_ids = tree.find_outlier_clade(outlier_attribute_name=outlier_attribute_name)
		tree.plot_tree(output_filepath=newick_output_filepath[:-5]+tree_plot_file_suffix, time_series_info=time_series_info, tree_style=tree_style, leaf_size_map_to=leaf_size_map_to, show_leaf_names=show_leaf_names, color_branches_by=color_branches_by, line_width=line_width, start_color=start_color, end_color=end_color, ladderize=ladderize)

		#remove seqs that are in outlier clade from data
		sequence_clusters = sequence_sample(filepath=patient_alignment_filepath, count_attribute_name=count_attribute_name)
		outlier_cluster_to_seq_ids_dic = sequence_clusters.get_seqID_to_attribute_dic(query_attribute_name='indiv_seq_ids', only_for_seq_ids=outlier_seq_ids)
		#init dic
		tpoint_to_outlier_seq_ids_dic = {}
		for i in sample_time_series.timepoints:
			tpoint_to_outlier_seq_ids_dic[str(int(i))] = []
		for i in outlier_cluster_to_seq_ids_dic:
			tpoint = i.split('_')[0].split('.')[0]
			seq_ids = outlier_cluster_to_seq_ids_dic[i].split(',')
			tpoint_to_outlier_seq_ids_dic[tpoint] += seq_ids
		input_indiv_seqs_patient_dirpath = input_indiv_seq_fasta_dirpath + patient + '/'
		output_fasta_dirpath = output_outliers_removed_fasta_dirpath + patient + '/'
		if not os.path.exists(output_fasta_dirpath):
			os.makedirs(output_fasta_dirpath)
		outlier_seqs_patient_filepath = output_outlier_seqs_dirpath + patient + '.fasta'
		if os.path.exists(outlier_seqs_patient_filepath):
			subprocess.call(['rm', outlier_seqs_patient_filepath])
		for tpoint in tpoint_to_outlier_seq_ids_dic:
			input_idiv_seqs_filepath = input_indiv_seqs_patient_dirpath + tpoint + '.fasta'
			sample = sequence_sample(filepath=input_idiv_seqs_filepath, count_attribute_name=count_attribute_name)
			if tpoint_to_outlier_seq_ids_dic[tpoint]:
				sample.write_subset_of_seqs_to_disk(seq_indices=None, output_filepath=outlier_seqs_patient_filepath, append_to_file=True, seq_ids=tpoint_to_outlier_seq_ids_dic[tpoint], add_string_to_each_id_being_written=tpoint+'_')
			sample.remove_seq_entries(seq_indicators=tpoint_to_outlier_seq_ids_dic[tpoint], seq_ids=True)
			output_fasta_filepath = output_fasta_dirpath + tpoint + '.fasta'
			sample.write_full_data_to_disk_fasta(output_filepath=output_fasta_filepath, append_to_file=False, seq_id_fileout_dic=None)

	return

if __name__ == '__main__':

	#below are hard coded examples of how we used this code. One will need to replace these files with there own data if they want to use this code.
	
	# run(input_seq_clusters_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/sequence_clusters/hiv/clustered_by_edit_distance/max_edit_dist_4', output_msa_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/hiv/alignments_of_seq_clusters/max_edit_dist_4', output_newick_tree_dirpath='/Users/nstrauli/data/abr_hiv_coevo/phylogenetic_trees/hiv/fasttree/clusters_edit_dist_4_within_tpoints', input_indiv_seq_fasta_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_annotated_fasta_files', output_outliers_removed_fasta_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_fasta_files', output_outlier_seqs_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/outlier_seqs/hiv/by_being_in_clade_with_miss_IDed_barcode_seqs')
