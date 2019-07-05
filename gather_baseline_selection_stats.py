import os

def run(input_fasta_dirpath, output_selection_dstrb_dirpath, output_selection_mean_dirpath):
	if input_fasta_dirpath[-1] != '/':
		input_fasta_dirpath += '/'
	if output_selection_dstrb_dirpath[-1] != '/':
		output_selection_dstrb_dirpath += '/'
	if not os.path.exists(output_selection_dstrb_dirpath):
		os.makedirs(output_selection_dstrb_dirpath)
	if output_selection_mean_dirpath[-1] != '/':
		output_selection_mean_dirpath += '/'
	if not os.path.exists(output_selection_mean_dirpath):
		os.makedirs(output_selection_mean_dirpath)
	for i in os.listdir(input_fasta_dirpath):
		if i[0] == '.' or i[:6] == 'README' or i[-4:] == '.pdf' or i == 'unknown':
			continue

		print 'patient:', i

		output_selection_dstrb_patient_dirpath = output_selection_dstrb_dirpath + i + '/'
		if not os.path.exists(output_selection_dstrb_patient_dirpath):
			os.makedirs(output_selection_dstrb_patient_dirpath)
		output_selection_mean_filepath = output_selection_mean_dirpath + i
		fileout_mean = open(output_selection_mean_filepath, "w")
		fileout_mean.write('time_point\tmean_FWR_baseline_sigma\tmean_CDR_baseline_sigma\n')
		#order files by time-point
		tpoint_filepaths = []
		for j in os.listdir(input_fasta_dirpath + i):
			if j[0] == '.' or j[:6] == 'README' or j[-4:] == '.pdf':
				continue
			tpoint_filepaths.append([float(j[:-6]), input_fasta_dirpath + i + '/' + j])
		input_fasta_filepaths = [j[1] for j in sorted(tpoint_filepaths)]
		for input_fasta_filepath in input_fasta_filepaths:
			tpoint = os.path.basename(input_fasta_filepath)[:-6]

			print '\ttime-point:', tpoint

			filein = open(input_fasta_filepath, "r")
			output_selection_dstrb_filepath = output_selection_dstrb_patient_dirpath + tpoint
			fileout_dstrb = open(output_selection_dstrb_filepath, "w")
			fileout_dstrb.write('sequence_ID\tFWR_baseline_sigma\tFWR_baseline_pvalue\tCDR_baseline_sigma\tCDR_baseline_pvalue\n')
			cdr_sigma_sum = 0.
			cdr_sigma_count = 0.
			fwr_sigma_sum = 0.
			fwr_sigma_count = 0.
			for k in filein:
				if k[0] == '>':
					header = k[1:-1].split('|')
					seq_id = header[0]
					for l in header:
						attribute = l.split('=')
						if attribute[0] == 'FWR_baseline_selection':
							fwr_sigma = attribute[1]
						elif attribute[0] == 'FWR_baseline_selection_pval':
							fwr_pval = attribute[1]
						elif attribute[0] == 'CDR_baseline_selection':
							cdr_sigma = attribute[1]
						elif attribute[0] == 'CDR_baseline_selection_pval':
							cdr_pval = attribute[1]
						elif attribute[0] == 'DUPCOUNT':
							count = float(attribute[1])
					for l in xrange(int(count)):
						fileout_dstrb.write('%s\t%s\t%s\t%s\t%s\n' % (seq_id, fwr_sigma, fwr_pval, cdr_sigma, cdr_pval))
					if fwr_sigma != 'NA':
						fwr_sigma_sum += float(fwr_sigma) * count
						fwr_sigma_count += count
					if cdr_sigma != 'NA':
						cdr_sigma_sum += float(cdr_sigma) * count
						cdr_sigma_count += count
				else:
					seq_id = None
					fwr_sigma = None
					fwr_pval = None
					cdr_sigma = None
					cdr_pval = None
					count = None
			filein.close()
			fileout_dstrb.close()
			fileout_mean.write('%s\t%s\t%s\n' % (tpoint, fwr_sigma_sum/fwr_sigma_count, cdr_sigma_sum/cdr_sigma_count))
		fileout_mean.close()
	return

if __name__ == '__main__':

	#below are hard coded examples of how we used this code. One will need to replace these files with there own data if they want to use this code.
	
	#On Gibbon
	# run(input_fasta_dirpath='/Volumes/nicolas_ext_harddrive/data/abr_hiv_coevo/seq_data/changeo_annotated_fasta_files', output_selection_dstrb_dirpath='/Users/nstrauli/data/abr_hiv_coevo/selection/abr/changeo/baseline/selection_dstrbs', output_selection_mean_dirpath='/Users/nstrauli/data/abr_hiv_coevo/selection/abr/changeo/baseline/mean_selection_trajectories')
	#On netapp
	# run(input_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/changeo_annotated_fasta_files', output_selection_dstrb_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/selection/abr/changeo/baseline/selection_dstrbs', output_selection_mean_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/selection/abr/changeo/baseline/mean_selection_trajectories')

