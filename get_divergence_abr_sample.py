import os
import sys
from sequence_sample_class import sequence_sample

def get_mean_divergence(input_fasta_dirpath, output_dirpath, dataset='partis'):
	"""
	dataset - gives the type of preprocessing that was done on the raw seq data. The different processing have different formatted fasta headers. Exceptable values are:
		'partis' - means that the processing was PRESTO->partis->mixcr
		'changeo' - mean that the processing was PRESTO->IgBLAST->changeo
	"""
	if input_fasta_dirpath[-1] != '/':
		input_fasta_dirpath += '/'
	if output_dirpath[-1] != '/':
		output_dirpath += '/'
	if not os.path.exists(output_dirpath):
		os.makedirs(output_dirpath)
	for i in os.listdir(input_fasta_dirpath):
		if i[0] == '.' or i[:6] == 'README':
			continue

		print 'patient:', i

		fileout = open(output_dirpath + i, "w")
		fileout.write('time_point\tmean_divergence\n')
		for j in os.listdir(input_fasta_dirpath + i):
			if j[0] == '.' or j[:6] == 'README':
				continue

			print '\ttime-point:', j

			time_point = j[:-6]
			input_fasta_filepath = input_fasta_dirpath + i + '/' + j
			if dataset == 'partis':
				sample = sequence_sample(input_fasta_filepath, count_attribute_name='count')
				mean_divergence = sample.get_mean_attribute_value('mut_freq')
			elif dataset == 'changeo':
				sample = sequence_sample(input_fasta_filepath, count_attribute_name='DUPCOUNT')
				mean_divergence = sample.get_mean_attribute_value('total_mut_freq')
			fileout.write('%s\t%s\n' % (time_point, mean_divergence))
		fileout.close()
	return

def get_divergence_distrbs(input_fasta_dirpath, output_dirpath, dataset='partis'):
	"""
	dataset - gives the type of preprocessing that was done on the raw seq data. The different processing have different formatted fasta headers. Exceptable values are:
		'partis' - means that the processing was PRESTO->partis->mixcr
		'changeo' - mean that the processing was PRESTO->IgBLAST->changeo
	"""
	if input_fasta_dirpath[-1] != '/':
		input_fasta_dirpath += '/'
	if output_dirpath[-1] != '/':
		output_dirpath += '/'
	for i in os.listdir(input_fasta_dirpath):
		if i[0] == '.' or i[:6] == 'README':
			continue

		print 'patient:', i

		if not os.path.exists(output_dirpath + i):
			os.makedirs(output_dirpath + i)
		for j in os.listdir(input_fasta_dirpath + i):
			if j[0] == '.' or j[:6] == 'README':
				continue

			print '\ttime-point:', j

			time_point = j[:-6]
			input_fasta_filepath = input_fasta_dirpath + i + '/' + j
			output_filepath = output_dirpath + i + '/' + time_point
			# print output_filepath
			# if os.path.exists(output_filepath):
			# 	print "output_filepath exists, so skipping"
			# 	continue
			fileout = open(output_filepath, "w")
			fileout.write('sequence_ID\tdivergence\n')
			if dataset == 'partis':
				sample = sequence_sample(input_fasta_filepath, count_attribute_name='count')
				for k in sample.data:
					for l in xrange(k['count']):
						fileout.write('%s\t%s\n' % (k['id'], k['other']['mut_freq']))
			elif dataset == 'changeo':
				sample = sequence_sample(input_fasta_filepath, count_attribute_name='DUPCOUNT')
				for k in sample.data:
					for l in xrange(k['count']):
						fileout.write('%s\t%s\n' % (k['id'], k['other']['total_mut_freq']))
			fileout.close()
	return

if __name__ == '__main__':

	#below are hard coded examples of how we used this code. One will need to replace these files with there own data if they want to use this code.
	
	# get_mean_divergence(input_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/changeo_annotated_fasta_files', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/divergence/abr/changeo/mean_divergence_trajectories', dataset='changeo')
	# get_divergence_distrbs(input_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/changeo_annotated_fasta_files', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/divergence/abr/changeo/divergence_distrbs', dataset='changeo')
	# get_mean_divergence(input_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/partis_annotated_fasta_files', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/divergence/abr/partis/mean_divergence_trajectories', dataset='partis')
	# get_divergence_distrbs(input_fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/partis_annotated_fasta_files', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/divergence/abr/partis/divergence_distrbs', dataset='partis')
