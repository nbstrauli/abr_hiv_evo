import os

def run(input_fasta_dirpath, divergence_dstrb_output_dirpath, divergence_mean_output_dirpath, dNdS_dstrb_output_dirpath, dNdS_mean_output_dirpath):
	if input_fasta_dirpath[-1] != '/':
		input_fasta_dirpath += '/'
	if divergence_dstrb_output_dirpath[-1] != '/':
		divergence_dstrb_output_dirpath += '/'
	if not os.path.exists(divergence_dstrb_output_dirpath):
		os.makedirs(divergence_dstrb_output_dirpath)
	if divergence_mean_output_dirpath[-1] != '/':
		divergence_mean_output_dirpath += '/'
	if not os.path.exists(divergence_mean_output_dirpath):
		os.makedirs(divergence_mean_output_dirpath)
	if dNdS_dstrb_output_dirpath[-1] != '/':
		dNdS_dstrb_output_dirpath += '/'
	if not os.path.exists(dNdS_dstrb_output_dirpath):
		os.makedirs(dNdS_dstrb_output_dirpath)
	if dNdS_mean_output_dirpath[-1] != '/':
		dNdS_mean_output_dirpath += '/'
	if not os.path.exists(dNdS_mean_output_dirpath):
		os.makedirs(dNdS_mean_output_dirpath)
	for i in os.listdir(input_fasta_dirpath):
		if i[0] == '.' or i[:6] == 'README':
			continue
		divergence_dstrb_patient_dirpath = divergence_dstrb_output_dirpath + i + '/'
		if not os.path.exists(divergence_dstrb_patient_dirpath):
			os.makedirs(divergence_dstrb_patient_dirpath)
		divergence_mean_output_filepath = divergence_mean_output_dirpath + i
		fileout_mean_divergence = open(divergence_mean_output_filepath, "w")
		fileout_mean_divergence.write('time_point\tmean_synonymous_divergence\tmean_non_synonymous_divergence\n')
		dNdS_dstrb_patient_dirpath = dNdS_dstrb_output_dirpath + i + '/'
		if not os.path.exists(dNdS_dstrb_patient_dirpath):
			os.makedirs(dNdS_dstrb_patient_dirpath)
		dNdS_mean_output_filepath = dNdS_mean_output_dirpath + i
		fileout_mean_dnds = open(dNdS_mean_output_filepath, "w")
		fileout_mean_dnds.write('time_point\tmean_dN_dS\n')
		#order files by time-point
		tpoint_filepaths = []
		for j in os.listdir(input_fasta_dirpath + i):
			if j[0] == '.' or j[:6] == 'README':
				continue
			tpoint_filepaths.append([float(j[:-6]), input_fasta_dirpath + i + '/' + j])
		input_fasta_filepaths = [j[1] for j in sorted(tpoint_filepaths)]
		#this will be a list of lists, where each sub list has [dS, dN, count, seq_id]
		seq_info = []
		for input_fasta_filepath in input_fasta_filepaths:
			tpoint = os.path.basename(input_fasta_filepath)[:-6]
			filein = open(input_fasta_filepath, "r")
			divergence_dstrb_output_filepath = divergence_dstrb_patient_dirpath + tpoint
			fileout_dstrb_divergence = open(divergence_dstrb_output_filepath, "w")
			fileout_dstrb_divergence.write('sequence_id\tsynonymous_divergence\tnon_synonymous_divergence\n')
			dNdS_dstrb_output_filepath = dNdS_dstrb_patient_dirpath + tpoint
			fileout_dstrb_dnds = open(dNdS_dstrb_output_filepath, "w")
			fileout_dstrb_dnds.write('sequence_id\tdN_dS\n')
			syn_divegence_sum = 0.
			non_syn_divergence_sum = 0.
			divergence_seq_count = 0.
			dnds_sum = 0.
			dnds_seq_count = 0.
			for k in filein:
				if k[0] == '>':
					header = k[1:-1].split('|')
					seq_id = header[0]
					for l in header:
						attribute = l.split('=')
						if attribute[0] == 'dN_dS':
							dnds = attribute[1]
						elif attribute[0] == 'divergence_syn':
							syn_divg = attribute[1]
						elif attribute[0] == 'divergence_non_syn':
							nonsyn_divg = attribute[1]
						elif attribute[0] == 'DUPCOUNT':
							count = float(attribute[1])
					for l in xrange(int(count)):
						fileout_dstrb_divergence.write('%s\t%s\t%s\n' % (seq_id, syn_divg, nonsyn_divg))
						fileout_dstrb_dnds.write('%s\t%s\n' % (seq_id, dnds))
					if dnds != 'NA':
						dnds_seq_count += count
						dnds_sum += float(dnds) * count
					if syn_divg != 'NA':
						divergence_seq_count += count
						syn_divegence_sum += float(syn_divg) * count
						non_syn_divergence_sum += float(nonsyn_divg) * count
						seq_info.append([syn_divg, nonsyn_divg, count, tpoint+'_'+seq_id])
				else:
					dnds = None
					syn_divg = None
					nonsyn_divg = None
					count = None
			filein.close()
			fileout_dstrb_divergence.close()
			fileout_dstrb_dnds.close()
			fileout_mean_divergence.write('%s\t%s\t%s\n' % (tpoint, syn_divegence_sum/divergence_seq_count, non_syn_divergence_sum/divergence_seq_count))
			fileout_mean_dnds.write('%s\t%s\n' % (tpoint, dnds_sum/dnds_seq_count))

		divergence_dstrb_allTpoints_filepath = divergence_dstrb_output_dirpath + i + '_pooled.txt'
		fileout_divergence_dstrb_pooled = open(divergence_dstrb_allTpoints_filepath, "w")
		fileout_divergence_dstrb_pooled.write('sequence_id\tsynonymous_divergence\tnon_synonymous_divergence\tsequence_count\n')
		for j in sorted(seq_info, reverse=True):
			fileout_divergence_dstrb_pooled.write('%s\t%s\t%s\t%s\n' % (j[3], j[0], j[1], j[2]))
		fileout_divergence_dstrb_pooled.close()

		fileout_mean_divergence.close()
		fileout_mean_dnds.close()
	return

if __name__ == '__main__':

	#below are hard coded examples of how we used this code. One will need to replace these files with there own data if they want to use this code.
	
	#on local comp
	# run(input_fasta_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_fasta_files', divergence_dstrb_output_dirpath='/Users/nstrauli/data/abr_hiv_coevo/divergence/hiv/counting_relative_to_most_numerous_1st_tpoint_seq/divergence_distrbs', divergence_mean_output_dirpath='/Users/nstrauli/data/abr_hiv_coevo/divergence/hiv/counting_relative_to_most_numerous_1st_tpoint_seq/mean_divergence_trajectories', dNdS_dstrb_output_dirpath='/Users/nstrauli/data/abr_hiv_coevo/selection/hiv/dN_dS_by_counting_changes_relative_to_most_numerouse_1st_tpoint_seq/selection_distrbs', dNdS_mean_output_dirpath='/Users/nstrauli/data/abr_hiv_coevo/selection/hiv/dN_dS_by_counting_changes_relative_to_most_numerouse_1st_tpoint_seq/mean_selection_trajectories')
	#on cluster
	# run(input_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_fasta_files', divergence_dstrb_output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/divergence/hiv/counting_relative_to_most_numerous_1st_tpoint_seq/divergence_distrbs', divergence_mean_output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/divergence/hiv/counting_relative_to_most_numerous_1st_tpoint_seq/mean_divergence_trajectories', dNdS_dstrb_output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/selection/hiv/dN_dS_by_counting_changes_relative_to_most_numerouse_1st_tpoint_seq/selection_distrbs', dNdS_mean_output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/selection/hiv/dN_dS_by_counting_changes_relative_to_most_numerouse_1st_tpoint_seq/mean_selection_trajectories')
	# run(input_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_annotated_fasta_files_patient6', divergence_dstrb_output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/divergence/hiv/counting_relative_to_most_numerous_1st_tpoint_seq_patient6/divergence_distrbs', divergence_mean_output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/divergence/hiv/counting_relative_to_most_numerous_1st_tpoint_seq_patient6/mean_divergence_trajectories', dNdS_dstrb_output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/selection/hiv/dN_dS_by_counting_changes_relative_to_most_numerouse_1st_tpoint_seq_patient6/selection_distrbs', dNdS_mean_output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/selection/hiv/dN_dS_by_counting_changes_relative_to_most_numerouse_1st_tpoint_seq_patient6/mean_selection_trajectories')
