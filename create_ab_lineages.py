import os
from sequence_sample_class import sequence_sample
from sequence_sample_set_class import sequence_sample_set
from sequence_clusters_time_series_class import sequence_clusters_time_series
from subprocess import Popen, PIPE
import subprocess
import time
from sequence_clusters_time_series_set_class import sequence_clusters_time_series_set
from seq_lineage_time_series_class import seq_lineage_time_series

def run(input_VJ_pair_seq_dirpath, input_VJ_pair_freq_traj_dirpath, output_sequence_cluster_fasta_dirpath, output_lineage_fasta_dirpath, genetic_dist_dirpath, output_muller_plot_dirpath, genetic_dists_withinSamps_dirpath, network_dirpath, pooled_network_dirpath, max_edit_distance_within_sample, max_edit_distance_across_samples, patient_VJ_pair_of_interest_dic, sim_null_lineages_fasta_dirpath, sim_null_lineages_muller_plot_dirpath, prob_of_new_lineage_dirpath):

	####### parameters #######
	count_attribute_name = 'DUPCOUNT'
	sample_clustering_method = 'by_edit_distance_seqanpy'
	use_comp_cluster = True
	path_to_needle = '/netapp/home/nstrauli/tools_c/EMBOSS-6.6.0/emboss/needle'
	temp_dirpath = '/netapp/home/nstrauli/data/abr_hiv_coevo/temp_stuff'
	# temp_dirpath = '/scrapp'
	distance_metric = 'dist_between_rep_seqs'
	distance_units = 'edit_distance'
	node_attributes_to_add = ['count', 'freq', 'total_mut_count', 'total_mut_freq', 'CDR_baseline_selection', 'FWR_baseline_selection']
	names_for_node_attributes_in_header = ['count', 'frequency', 'number_mutations', 'frequency_of_mutations', 'selection_CDR', 'selection_FWR']
	min_cluster_freq = None
	# min_cluster_freq = 0.0001
	compare_tpoint_to_all_previous = False
	start_end_colors = ['red', 'purple']
	lineage_min_freq_cutoff = 0.0001
	cluster_within_job_name = None
	alignment_method = 'seqanpy'
	multiply_num_sim_lineages_by = 100
	cluster_across_job_name = None
	##########################

	if input_VJ_pair_seq_dirpath[-1] != '/':
		input_VJ_pair_seq_dirpath += '/'
	if input_VJ_pair_freq_traj_dirpath[-1] != '/':
		input_VJ_pair_freq_traj_dirpath += '/'
	if output_sequence_cluster_fasta_dirpath[-1] != '/':
		output_sequence_cluster_fasta_dirpath += '/'
	if output_lineage_fasta_dirpath[-1] != '/':
		output_lineage_fasta_dirpath += '/'
	if genetic_dist_dirpath[-1] != '/':
		genetic_dist_dirpath += '/'
	if output_muller_plot_dirpath[-1] != '/':
		output_muller_plot_dirpath += '/'
	if genetic_dists_withinSamps_dirpath[-1] != '/':
		genetic_dists_withinSamps_dirpath += '/'
	if network_dirpath[-1] != '/':
		network_dirpath += '/'
	if pooled_network_dirpath[-1] != '/':
		pooled_network_dirpath += '/'
	if sim_null_lineages_fasta_dirpath[-1] != '/':
		sim_null_lineages_fasta_dirpath += '/'
	if sim_null_lineages_muller_plot_dirpath[-1] != '/':
		sim_null_lineages_muller_plot_dirpath
	if prob_of_new_lineage_dirpath[-1] != '/':
		prob_of_new_lineage_dirpath += '/'

	print 'max edit dist within samp:', max_edit_distance_within_sample
	print 'max edit dist across samps:', max_edit_distance_across_samples

	for patient in patient_VJ_pair_of_interest_dic:
		print 'patient:', patient

		# if patient != '7':
		# 	continue

		for VJ_pair in patient_VJ_pair_of_interest_dic[patient]:

			# if VJ_pair != 'IGHV3-15_IGHJ4':
			# 	continue

			print '\tV,J gene pair:', VJ_pair

			input_patient_VJ_pair_seq_dirpath = '%s%s/%s/' % (input_VJ_pair_seq_dirpath, patient, VJ_pair)

			#get relative frequencies of VJ gene pair over time
			VJ_freq_traj_filepath = '%s%s/%s.txt' % (input_VJ_pair_freq_traj_dirpath, patient, VJ_pair)
			filein = open(VJ_freq_traj_filepath, "r")
			filein.readline()
			tpoint_freq_dic = {}
			for i in filein:
				line = i[:-1].split('\t')
				tpoint = line[0].split('.')[0]
				freq = float(line[1])
				sample_filepath = '%s%s.fasta' % (input_patient_VJ_pair_seq_dirpath, tpoint)
				tpoint_freq_dic[sample_filepath] = freq
			filein.close()

			#cluster within time-points
			patient_seq_cluster_output_dirpath = "%smax_edit_dist_%s/%s/%s/" % (output_sequence_cluster_fasta_dirpath, max_edit_distance_within_sample, patient, VJ_pair)
			# if os.path.exists(patient_seq_cluster_output_dirpath):
			# 	subprocess.call(['rm', '-r', patient_seq_cluster_output_dirpath])
			# os.makedirs(patient_seq_cluster_output_dirpath)
			patient_within_dist_dirpath = '%s%s/%s/' % (genetic_dists_withinSamps_dirpath, patient, VJ_pair)
			if not os.path.exists(patient_within_dist_dirpath):
				os.makedirs(patient_within_dist_dirpath)
			patient_network_dirpath = '%smax_edit_dist_%s/%s/%s/' % (network_dirpath, max_edit_distance_within_sample, patient, VJ_pair)
			if not os.path.exists(patient_network_dirpath):
				os.makedirs(patient_network_dirpath)
			# seq_sample_set = sequence_sample_set(dirpath=input_patient_VJ_pair_seq_dirpath, count_attribute_name=count_attribute_name, avoid_files_with='.')
			# cluster_within_job_name = seq_sample_set.cluster_seqs_within_samples(output_dirpath=patient_seq_cluster_output_dirpath, method=sample_clustering_method, max_edit_distance=max_edit_distance_within_sample, use_comp_cluster=use_comp_cluster, wait_till_cluster_done=False, path_to_needle=path_to_needle, temp_dirpath=temp_dirpath, scale_freqs_by=tpoint_freq_dic, dists_dstrb_output_dirpath=patient_within_dist_dirpath, output_network_dirpath=patient_network_dirpath, output_node_attribute_dirpath=patient_network_dirpath, node_attributes_to_add=node_attributes_to_add, names_for_node_attributes_in_header=names_for_node_attributes_in_header, add_freq_prior_to_clustering='freq', add_indiv_seq_attribute=['total_mut_count', 'total_mut_count_non_D_masked'])

			#cluster across time-points
			patient_lineage_output_dirpath = '%smax_edit_dist_within_samps_%s_across_samps_%s/%s/%s/' % (output_lineage_fasta_dirpath, max_edit_distance_within_sample, max_edit_distance_across_samples, patient, VJ_pair)
			if os.path.exists(patient_lineage_output_dirpath):
				subprocess.call(['rm', '-r', patient_lineage_output_dirpath])
			os.makedirs(patient_lineage_output_dirpath)
			patient_genetic_dist_filepath = '%smax_edit_dist_within_samps_%s_across_samps_%s/%s/%s.txt' % (genetic_dist_dirpath, max_edit_distance_within_sample, max_edit_distance_across_samples, patient, VJ_pair)
			if not os.path.exists(os.path.dirname(patient_genetic_dist_filepath)):
				os.makedirs(os.path.dirname(patient_genetic_dist_filepath))
			patient_output_muller_plot_dirpath = '%smax_edit_dist_within_samps_%s_across_samps_%s/%s/%s/' % (output_muller_plot_dirpath, max_edit_distance_within_sample, max_edit_distance_across_samples, patient, VJ_pair)
			if not os.path.exists(patient_output_muller_plot_dirpath):
				os.makedirs(patient_output_muller_plot_dirpath)
			patient_pooled_network_dirpath = '%smax_edit_dist_within_samps_%s_across_samps_%s/%s/%s/' % (pooled_network_dirpath, max_edit_distance_within_sample, max_edit_distance_across_samples, patient, VJ_pair)
			if not os.path.exists(patient_pooled_network_dirpath):
				os.makedirs(patient_pooled_network_dirpath)
			patient_prob_of_new_lineage_filepath = '%smax_edit_dist_within_samps_%s_across_samps_%s/%s/%s.txt' % (prob_of_new_lineage_dirpath, max_edit_distance_within_sample, max_edit_distance_across_samples, patient, VJ_pair)
			if not os.path.exists(os.path.dirname(patient_prob_of_new_lineage_filepath)):
				os.makedirs(os.path.dirname(patient_prob_of_new_lineage_filepath))
			patient_sim_null_lineages_fasta_dirpath = '%smax_edit_dist_within_samps_%s_across_samps_%s/%s/%s/' % (sim_null_lineages_fasta_dirpath, max_edit_distance_within_sample, max_edit_distance_across_samples, patient, VJ_pair)
			if not os.path.exists(patient_sim_null_lineages_fasta_dirpath):
				os.makedirs(patient_sim_null_lineages_fasta_dirpath)
			patient_sim_null_lineages_muller_plot_dirpath = '%smax_edit_dist_within_samps_%s_across_samps_%s/%s/%s' % (sim_null_lineages_muller_plot_dirpath, max_edit_distance_within_sample, max_edit_distance_across_samples, patient, VJ_pair)
			if not os.path.exists(patient_sim_null_lineages_muller_plot_dirpath):
				os.makedirs(patient_sim_null_lineages_muller_plot_dirpath)
			seq_clusters_time_series = sequence_clusters_time_series(dirpath=patient_seq_cluster_output_dirpath, count_attribute_name=count_attribute_name, freq_attribute_name='total_freq', indiv_seqs_attribute_name='indiv_seqs', indiv_seq_ids_attribute_name='indiv_seq_ids', indiv_seq_counts_attribute_name='indiv_seq_counts', indiv_seq_freqs_attribute_name='indiv_seq_freqs')
			#create 'real' lineages
			cluster_across_job_name = seq_clusters_time_series.create_lineages(distance_metric=distance_metric, temp_dirpath=temp_dirpath, path_to_needle=path_to_needle, distance_units=distance_units, compare_tpoint_to_all_previous=compare_tpoint_to_all_previous, genetic_dist_filepath=patient_genetic_dist_filepath, output_fasta_dirpath=patient_lineage_output_dirpath, write_full_original_seqs=input_patient_VJ_pair_seq_dirpath, muller_plot_output_dirpath=patient_output_muller_plot_dirpath, start_end_colors=start_end_colors, max_distance=max_edit_distance_across_samples, min_cluster_freq=min_cluster_freq, use_comp_cluster=use_comp_cluster, lineage_min_freq_cutoff=lineage_min_freq_cutoff, make_master_network_plot=patient_pooled_network_dirpath, indiv_sample_network_dirpath=patient_network_dirpath, mut_count_attribute_name='indiv_seq_total_mut_count_non_D_masked', alignment_method=alignment_method, sim_null_lineages_fasta_dirpath=None, sim_null_lineages_muller_plot_dirpath=None, multiply_num_sim_lineages_by=1, only_simulate=False, prob_of_new_lineage_filepath=patient_prob_of_new_lineage_filepath, hold_for_job_id=cluster_within_job_name)
			#create simulated null lineages
			# seq_clusters_time_series.create_lineages(distance_metric=distance_metric, temp_dirpath=temp_dirpath, path_to_needle=path_to_needle, distance_units=distance_units, compare_tpoint_to_all_previous=compare_tpoint_to_all_previous, genetic_dist_filepath=None, output_fasta_dirpath=None, write_full_original_seqs=input_patient_VJ_pair_seq_dirpath, muller_plot_output_dirpath=None, start_end_colors=start_end_colors, max_distance=None, min_cluster_freq=min_cluster_freq, use_comp_cluster=use_comp_cluster, lineage_min_freq_cutoff=lineage_min_freq_cutoff, make_master_network_plot=None, indiv_sample_network_dirpath=None, mut_count_attribute_name='indiv_seq_total_mut_count_non_D_masked', alignment_method=alignment_method, sim_null_lineages_fasta_dirpath=patient_sim_null_lineages_fasta_dirpath, sim_null_lineages_muller_plot_dirpath=patient_sim_null_lineages_muller_plot_dirpath, multiply_num_sim_lineages_by=multiply_num_sim_lineages_by, only_simulate=True, prob_of_new_lineage_filepath=patient_prob_of_new_lineage_filepath, hold_for_job_id=cluster_across_job_name)

			# time.sleep(100)
	return

def run_for_all_AbR_partitions(input_VJ_pair_seq_dirpath, input_full_sample_dirpath, output_sequence_cluster_fasta_dirpath, output_lineage_fasta_dirpath, genetic_dist_dirpath, output_muller_plot_dirpath, genetic_dists_withinSamps_dirpath, network_dirpath, pooled_network_dirpath, max_edit_distance_within_sample, max_edit_distance_across_samples, sim_null_lineages_fasta_dirpath, sim_null_lineages_muller_plot_dirpath, prob_of_new_lineage_dirpath):

	####### parameters #######
	count_attribute_name = 'DUPCOUNT'
	sample_clustering_method = 'by_edit_distance_seqanpy'
	use_comp_cluster = True
	path_to_needle = '/netapp/home/nstrauli/tools_c/EMBOSS-6.6.0/emboss/needle'
	temp_dirpath = '/netapp/home/nstrauli/data/abr_hiv_coevo/temp_stuff'
	# temp_dirpath = '/scrapp'
	distance_metric = 'dist_between_rep_seqs'
	distance_units = 'edit_distance'
	node_attributes_to_add = ['count', 'freq', 'total_mut_count', 'total_mut_freq', 'CDR_baseline_selection', 'FWR_baseline_selection']
	names_for_node_attributes_in_header = ['count', 'frequency', 'number_mutations', 'frequency_of_mutations', 'selection_CDR', 'selection_FWR']
	min_cluster_freq = None
	# min_cluster_freq = 0.0001
	compare_tpoint_to_all_previous = False
	start_end_colors = ['red', 'purple']
	lineage_min_freq_cutoff = 0.0001
	cluster_within_job_name = None
	alignment_method = 'seqanpy'
	multiply_num_sim_lineages_by = 100
	cluster_across_job_name = None
	##########################

	if input_VJ_pair_seq_dirpath[-1] != '/':
		input_VJ_pair_seq_dirpath += '/'
	if input_full_sample_dirpath[-1] != '/':
		input_full_sample_dirpath += '/'
	if output_sequence_cluster_fasta_dirpath[-1] != '/':
		output_sequence_cluster_fasta_dirpath += '/'
	if output_lineage_fasta_dirpath[-1] != '/':
		output_lineage_fasta_dirpath += '/'
	if genetic_dist_dirpath[-1] != '/':
		genetic_dist_dirpath += '/'
	if output_muller_plot_dirpath[-1] != '/':
		output_muller_plot_dirpath += '/'
	if genetic_dists_withinSamps_dirpath[-1] != '/':
		genetic_dists_withinSamps_dirpath += '/'
	if network_dirpath[-1] != '/':
		network_dirpath += '/'
	if pooled_network_dirpath[-1] != '/':
		pooled_network_dirpath += '/'
	if sim_null_lineages_fasta_dirpath[-1] != '/':
		sim_null_lineages_fasta_dirpath += '/'
	if sim_null_lineages_muller_plot_dirpath[-1] != '/':
		sim_null_lineages_muller_plot_dirpath += '/'
	if prob_of_new_lineage_dirpath[-1] != '/':
		prob_of_new_lineage_dirpath += '/'

	print 'max edit dist within samp:', max_edit_distance_within_sample
	print 'max edit dist across samps:', max_edit_distance_across_samples

	for patient in os.listdir(input_VJ_pair_seq_dirpath):
		if patient[0] == '.' or patient[:6] == 'README':
			continue

		print 'patient:', patient

		# if patient != '5':
		# 	continue

		#get relative frequencies of VJ gene pair over time
		input_patient_VJ_pair_seq_dirpath = '%s%s/' % (input_VJ_pair_seq_dirpath, patient)
		input_patient_full_sample_dirpath = '%s%s/' % (input_full_sample_dirpath, patient)
		seq_lins = seq_lineage_time_series(dirpath=input_patient_VJ_pair_seq_dirpath, count_attribute_name=count_attribute_name, list_of_filepaths=None, avoid_dirpath=None, avoid_timepoints_with_dot=True)
		lineage_freq_trajs = seq_lins.get_freq_trajectories(full_seq_sample_dirpath=input_patient_full_sample_dirpath, output_dirpath=None, min_freq_cutoff=0.)
		sample_to_freq_dic = {}
		for VJ_pair in lineage_freq_trajs:
			for tpoint_float in lineage_freq_trajs[VJ_pair]:
				tpoint = str(tpoint_float).split('.')[0]
				sample_filepath = '%s%s/%s.fasta' % (input_patient_VJ_pair_seq_dirpath, VJ_pair, tpoint)
				freq = lineage_freq_trajs[VJ_pair][tpoint_float]
				sample_to_freq_dic[sample_filepath] = freq

		#cluster within time-points
		patient_seq_cluster_output_dirpath = "%smax_edit_dist_%s/%s/" % (output_sequence_cluster_fasta_dirpath, max_edit_distance_within_sample, patient)
		# if os.path.exists(patient_seq_cluster_output_dirpath):
		# 	subprocess.call(['rm', '-r', patient_seq_cluster_output_dirpath])
		# os.makedirs(patient_seq_cluster_output_dirpath)
		patient_within_dist_dirpath = '%s%s/' % (genetic_dists_withinSamps_dirpath, patient)
		if not os.path.exists(patient_within_dist_dirpath):
			os.makedirs(patient_within_dist_dirpath)
		patient_network_dirpath = '%smax_edit_dist_%s/%s/' % (network_dirpath, max_edit_distance_within_sample, patient)
		if not os.path.exists(patient_network_dirpath):
			os.makedirs(patient_network_dirpath)
		# seq_sample_set = sequence_sample_set(dirpath=input_patient_VJ_pair_seq_dirpath, count_attribute_name=count_attribute_name, avoid_files_with='.')
		# cluster_within_job_name = seq_sample_set.cluster_seqs_within_samples(output_dirpath=patient_seq_cluster_output_dirpath, method=sample_clustering_method, max_edit_distance=max_edit_distance_within_sample, use_comp_cluster=use_comp_cluster, wait_till_cluster_done=False, path_to_needle=path_to_needle, temp_dirpath=temp_dirpath, scale_freqs_by=sample_to_freq_dic, dists_dstrb_output_dirpath=patient_within_dist_dirpath, output_network_dirpath=patient_network_dirpath, output_node_attribute_dirpath=patient_network_dirpath, node_attributes_to_add=node_attributes_to_add, names_for_node_attributes_in_header=names_for_node_attributes_in_header, add_freq_prior_to_clustering='freq', add_indiv_seq_attribute=['total_mut_count', 'total_mut_count_non_D_masked'])

		#cluster across time-points
		patient_lineage_output_dirpath = '%smax_edit_dist_within_samps_%s_across_samps_%s/%s/' % (output_lineage_fasta_dirpath, max_edit_distance_within_sample, max_edit_distance_across_samples, patient)
		if os.path.exists(patient_lineage_output_dirpath):
			subprocess.call(['rm', '-r', patient_lineage_output_dirpath])
		os.makedirs(patient_lineage_output_dirpath)
		patient_sim_null_lineages_fasta_dirpath = '%smax_edit_dist_within_samps_%s_across_samps_%s/%s/' % (sim_null_lineages_fasta_dirpath, max_edit_distance_within_sample, max_edit_distance_across_samples, patient)
		if os.path.exists(patient_sim_null_lineages_fasta_dirpath):
			subprocess.call(['rm', '-r', patient_sim_null_lineages_fasta_dirpath])
		os.makedirs(patient_sim_null_lineages_fasta_dirpath)
		patient_genetic_dist_dirpath = '%smax_edit_dist_within_samps_%s_across_samps_%s/%s/' % (genetic_dist_dirpath, max_edit_distance_within_sample, max_edit_distance_across_samples, patient)
		if not os.path.exists(os.path.dirname(patient_genetic_dist_dirpath)):
			os.makedirs(os.path.dirname(patient_genetic_dist_dirpath))
		patient_output_muller_plot_dirpath = '%smax_edit_dist_within_samps_%s_across_samps_%s/%s/' % (output_muller_plot_dirpath, max_edit_distance_within_sample, max_edit_distance_across_samples, patient)
		if not os.path.exists(patient_output_muller_plot_dirpath):
			os.makedirs(patient_output_muller_plot_dirpath)
		patient_pooled_network_dirpath = '%smax_edit_dist_within_samps_%s_across_samps_%s/%s/' % (pooled_network_dirpath, max_edit_distance_within_sample, max_edit_distance_across_samples, patient)
		if not os.path.exists(patient_pooled_network_dirpath):
			os.makedirs(patient_pooled_network_dirpath)
		patient_sim_null_lineages_muller_plot_dirpath = '%smax_edit_dist_within_samps_%s_across_samps_%s/%s/' % (sim_null_lineages_muller_plot_dirpath, max_edit_distance_within_sample, max_edit_distance_across_samples, patient)
		if not os.path.exists(patient_sim_null_lineages_muller_plot_dirpath):
			os.makedirs(patient_sim_null_lineages_muller_plot_dirpath)
		patient_prob_of_new_lineage_dirpath = '%smax_edit_dist_within_samps_%s_across_samps_%s/%s/' % (prob_of_new_lineage_dirpath, max_edit_distance_within_sample, max_edit_distance_across_samples, patient)
		if not os.path.exists(patient_prob_of_new_lineage_dirpath):
			os.makedirs(patient_prob_of_new_lineage_dirpath)
		seq_clusts_set = sequence_clusters_time_series_set(dirpath=patient_seq_cluster_output_dirpath, count_attribute_name=count_attribute_name, freq_attribute_name='total_freq', indiv_seqs_attribute_name='indiv_seqs', indiv_seq_ids_attribute_name='indiv_seq_ids', indiv_seq_counts_attribute_name='indiv_seq_counts', indiv_seq_freqs_attribute_name='indiv_seq_freqs')
		#create 'real' lineages
		cluster_across_job_name = seq_clusts_set.create_lineages(output_fasta_lineage_dirpath=patient_lineage_output_dirpath, temp_dirpath=temp_dirpath, hold_for_job_id=cluster_within_job_name, distance_metric=distance_metric, path_to_needle=path_to_needle, distance_units=distance_units, compare_tpoint_to_all_previous=compare_tpoint_to_all_previous, genetic_dist_dirpath=patient_genetic_dist_dirpath, write_full_original_seqs=input_patient_VJ_pair_seq_dirpath, muller_plot_output_dirpath=patient_output_muller_plot_dirpath, start_end_colors=start_end_colors, max_distance=max_edit_distance_across_samples, min_cluster_freq=min_cluster_freq, lineage_min_freq_cutoff=lineage_min_freq_cutoff, make_master_network_plot=patient_pooled_network_dirpath, indiv_sample_network_dirpath=patient_network_dirpath, mut_count_attribute_name='indiv_seq_total_mut_count_non_D_masked', alignment_method=alignment_method, sim_null_lineages_fasta_dirpath=None, sim_null_lineages_muller_plot_dirpath=None, multiply_num_sim_lineages_by=1, only_simulate=False, prob_of_new_lineage_dirpath=patient_prob_of_new_lineage_dirpath)
		#create simulated null lineages
		seq_clusts_set.create_lineages(output_fasta_lineage_dirpath=None, temp_dirpath=temp_dirpath, hold_for_job_id=cluster_across_job_name, compare_tpoint_to_all_previous=compare_tpoint_to_all_previous, genetic_dist_dirpath=None, write_full_original_seqs=input_patient_VJ_pair_seq_dirpath, muller_plot_output_dirpath=None, start_end_colors=start_end_colors, max_distance=None, min_cluster_freq=min_cluster_freq, lineage_min_freq_cutoff=lineage_min_freq_cutoff, make_master_network_plot=None, indiv_sample_network_dirpath=None, mut_count_attribute_name='indiv_seq_total_mut_count_non_D_masked', sim_null_lineages_fasta_dirpath=patient_sim_null_lineages_fasta_dirpath, sim_null_lineages_muller_plot_dirpath=patient_sim_null_lineages_muller_plot_dirpath, multiply_num_sim_lineages_by=multiply_num_sim_lineages_by, only_simulate=True, prob_of_new_lineage_dirpath=patient_prob_of_new_lineage_dirpath)

		# time.sleep(100)

	return


# max_edit_distance_within_sample_list = [5,6,7]
# max_edit_distance_across_samples_list = [25,30,35,'parsimony']

max_edit_distance_within_sample_list = [6]
max_edit_distance_across_samples_list = [30]

#This is for making lineages of the AbR partitions that seem to be associated with HIV
# HIV_associated_VJpairs_dic = {'1':['IGHV3-23_IGHJ6', 'IGHV4-34_IGHJ6'], '2':['IGHV3-48_IGHJ6', 'IGHV7-4-1_IGHJ4'], '3':['IGHV3-30_IGHJ3'], '4':['IGHV1-46_IGHJ6', 'IGHV4-59_IGHJ4'], '5':['IGHV4-31_IGHJ4', 'IGHV3-74_IGHJ4'], '6':['IGHV3-11_IGHJ4'], '7':['IGHV4-31_IGHJ5', 'IGHV4-59_IGHJ4'], '8':['IGHV6-1_IGHJ5', 'IGHV4-61_IGHJ5'], '9':['IGHV6-1_IGHJ4'], '10':[]}
HIV_associated_VJpairs_dic = {'3':['IGHV3-30_IGHJ3'], '7':['IGHV2-70_IGHJ6', 'IGHV3-15_IGHJ4', 'IGHV4-31_IGHJ5'], '8':['IGHV6-1_IGHJ4', 'IGHV6-1_IGHJ5']}
for edit_dist_within in max_edit_distance_within_sample_list:
	for edit_dist_across in max_edit_distance_across_samples_list:
		run(input_VJ_pair_seq_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/abr_seq_unique_VJ/changeo', input_VJ_pair_freq_traj_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/freq_trajectory_lineages/abr/unique_gene_pairs', output_sequence_cluster_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/sequence_clusters/abr/clustered_by_edit_distance/most_numerous_seq_foreach_cluster_select', output_lineage_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/lineage_fasta_files/abr/edit_dist_between_rep_seqs_select', genetic_dist_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/genetic_dist_between_clusters_across_samples/abr/edit_dist_between_rep_seqs_select', output_muller_plot_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/muller_plots/abr/edit_dist_between_rep_seqs_select', genetic_dists_withinSamps_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/genetic_dist_between_seqs_within_samples/abr/edit_dist_between_rep_seqs_select', network_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/network_files/abr/each_tpoint_cluster_by_edit_dist_select', pooled_network_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/network_files/abr/all_tpoints_pooled_clustered_by_edit_dist_select', max_edit_distance_within_sample=edit_dist_within, max_edit_distance_across_samples=edit_dist_across, patient_VJ_pair_of_interest_dic=HIV_associated_VJpairs_dic, sim_null_lineages_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/lineage_fasta_files_sim_null/abr/edit_dist_between_rep_seqs_select', sim_null_lineages_muller_plot_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/muller_plots_sim_null/abr/edit_dist_between_rep_seqs_select', prob_of_new_lineage_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/prob_of_new_lineage_foreach_tpoint_select/abr')

#This is for making lineages of the AbR partitions that seem to NOT be ascociated with HIV
HIV_NonAssociated_VJpairs_dic = {'1':['IGHV4-30-2_IGHJ6', 'IGHV4-31_IGHJ2'], '2':['IGHV1-8_IGHJ6', 'IGHV3-9_IGHJ3'], '3':['IGHV1-18_IGHJ5', 'IGHV3-7_IGHJ4'], '4':['IGHV6-1_IGHJ4', 'IGHV4-59_IGHJ5'], '5':['IGHV3-64_IGHJ3', 'IGHV4-39_IGHJ6'], '6':['IGHV1-8_IGHJ4', 'IGHV2-70_IGHJ4'], '7':['IGHV1-8_IGHJ5', 'IGHV1-46_IGHJ6'], '8':['IGHV2-70D_IGHJ5', 'IGHV3-11_IGHJ4'], '9':['IGHV7-4-1_IGHJ6', 'IGHV7-4-1_IGHJ5'], '10':[]}
# for edit_dist_within in max_edit_distance_within_sample_list:
# 	for edit_dist_across in max_edit_distance_across_samples_list:
# 		run(input_VJ_pair_seq_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/abr_seq_unique_VJ/changeo', input_VJ_pair_freq_traj_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/freq_trajectory_lineages/abr/unique_gene_pairs', output_sequence_cluster_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/sequence_clusters/abr/clustered_by_edit_distance/most_numerous_seq_foreach_cluster', output_lineage_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/lineage_fasta_files/abr/edit_dist_between_rep_seqs', genetic_dist_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/genetic_dist_between_clusters_across_samples/abr/edit_dist_between_rep_seqs', output_muller_plot_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/muller_plots/abr/edit_dist_between_rep_seqs', genetic_dists_withinSamps_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/genetic_dist_between_seqs_within_samples/abr/edit_dist_between_rep_seqs', network_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/network_files/abr/each_tpoint_cluster_by_edit_dist', pooled_network_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/network_files/abr/all_tpoints_pooled_clustered_by_edit_dist', max_edit_distance_within_sample=edit_dist_within, max_edit_distance_across_samples=edit_dist_across, patient_VJ_pair_of_interest_dic=HIV_NonAssociated_VJpairs_dic)

#This is for making lineages of all the other AbR partitions (ad hoc added later)
# all_other_VJpairs_dic = {'1':[], '2':[], '3':[], '4':[], '5':[], '6':[], '7':[], '8':[], '9':[], '10':[]}
# gene_pair_dirpath = '/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/freq_trajectory_lineages/abr/unique_gene_pairs/'
# for patient in os.listdir(gene_pair_dirpath):
# 	if patient[0] == '.' or patient[:6] == 'README' or patient == '10':
# 		continue
# 	for gene_pair in os.listdir('%s%s' % (gene_pair_dirpath, patient)):
# 		if gene_pair[0] == '.' or gene_pair[:6] == 'README':
# 			continue
# 		gene_pair = gene_pair.split('.')[0]
# 		if gene_pair in HIV_associated_VJpairs_dic[patient] or gene_pair in HIV_NonAssociated_VJpairs_dic[patient]:
# 			continue
# 		all_other_VJpairs_dic[patient].append(gene_pair)
# for edit_dist_within in max_edit_distance_within_sample_list:
# 	for edit_dist_across in max_edit_distance_across_samples_list:
# 		run(input_VJ_pair_seq_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/abr_seq_unique_VJ/changeo', input_VJ_pair_freq_traj_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/freq_trajectory_lineages/abr/unique_gene_pairs', output_sequence_cluster_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/sequence_clusters/abr/clustered_by_edit_distance/most_numerous_seq_foreach_cluster', output_lineage_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/lineage_fasta_files/abr/edit_dist_between_rep_seqs', genetic_dist_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/genetic_dist_between_clusters_across_samples/abr/edit_dist_between_rep_seqs', output_muller_plot_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/muller_plots/abr/edit_dist_between_rep_seqs', genetic_dists_withinSamps_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/genetic_dist_between_seqs_within_samples/abr/edit_dist_between_rep_seqs', network_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/network_files/abr/each_tpoint_cluster_by_edit_dist', pooled_network_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/network_files/abr/all_tpoints_pooled_clustered_by_edit_dist', max_edit_distance_within_sample=edit_dist_within, max_edit_distance_across_samples=edit_dist_across, patient_VJ_pair_of_interest_dic=all_other_VJpairs_dic)

# for edit_dist_within in max_edit_distance_within_sample_list:
# 	for edit_dist_across in max_edit_distance_across_samples_list:
# 		run_for_all_AbR_partitions(input_VJ_pair_seq_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/abr_seq_unique_VJ/changeo', input_full_sample_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/changeo_annotated_fasta_files', output_sequence_cluster_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/sequence_clusters/abr/clustered_by_edit_distance/most_numerous_seq_foreach_cluster', output_lineage_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/lineage_fasta_files/abr/edit_dist_between_rep_seqs', genetic_dist_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/genetic_dist_between_clusters_across_samples/abr/edit_dist_between_rep_seqs', output_muller_plot_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/muller_plots/abr/edit_dist_between_rep_seqs', genetic_dists_withinSamps_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/genetic_dist_between_seqs_within_samples/abr/edit_dist_between_rep_seqs', network_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/network_files/abr/each_tpoint_cluster_by_edit_dist', pooled_network_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/network_files/abr/all_tpoints_pooled_clustered_by_edit_dist', max_edit_distance_within_sample=edit_dist_within, max_edit_distance_across_samples=edit_dist_across, sim_null_lineages_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/lineage_fasta_files_sim_null_2/abr/edit_dist_between_rep_seqs', sim_null_lineages_muller_plot_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/muller_plots_sim_null/abr/edit_dist_between_rep_seqs', prob_of_new_lineage_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/prob_of_new_lineage_foreach_tpoint/abr')
