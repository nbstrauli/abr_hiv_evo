#!/usr/bin/python
#$ -S /usr/bin/python
#$ -e error
#$ -o out
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=4G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=336:00:00

import sys
sys.path.insert(0, './')
from illumina_seq_data_class import illumina_seq_data
import os
import re
from subprocess import Popen, PIPE

def de_multiplex_all_seq_runs(input_illumina_data_dirpath, sample_stats_filepath, output_dirpath, barcode_count_output_dirpath, unzip=True):
	"""
	This script will cycle through the raw output from the Illumina sequencer machine, and use the 'illumina_seq_data' class to demultiplex the samples that are in each of the sequencing runs. It does this by using the barcodes seqs for each of the samples as input.
	input_illumina_data_dirpath - The path to the directory that contains the raw MiSeq output. This should be the directory that contains the FULL output, with all files that are associated with each run.
	unzip - If True (default), the script will assume that the fastq sequence files have not been unzipped, and will do this first.
	sample_stats_filepath - The path to the file that contains a whole bunch of misc stats for each of the samples in the sample set. Included in these stats are sequencing date (to identify the sequencing run), barcode (to get get barcodes for each seq run), and sample names (to make names for output files).
	output_dirpath - The path to the directory for which the de-multiplexed seq data will be written.
	barcode_count_output_dirpath - This gives the path to the output directory for which the counts for each of the barcodes found for each of the seq runs, will be written.
	"""

	#############parameters####################
	barcode_start_pos = 4
	trim_barcode = True
	###########################################

	if input_illumina_data_dirpath[-1] != '/':
		input_illumina_data_dirpath += '/'
	if output_dirpath[-1] != '/':
		output_dirpath += '/'
	if not os.path.exists(output_dirpath):
		os.makedirs(output_dirpath)
	if barcode_count_output_dirpath[-1] != '/':
		barcode_count_output_dirpath += '/'
	if not os.path.exists(barcode_count_output_dirpath):
		os.makedirs(barcode_count_output_dirpath)

	#if desired, first unzip the fastq files for each seq run
	if unzip:
		for i in os.listdir(input_illumina_data_dirpath):
			if i[0] == '.' or i[:6] == 'README':
				continue
			#this is the directory that contains the fastq files (amoung other things)
			base_dirpath = '%s%s/Data/Intensities/BaseCalls/' % (input_illumina_data_dirpath, i)
			#need to find the fastq files we want from this directory
			for j in os.listdir(base_dirpath):
				match = re.search('_S1_L001_R1', j)
				if match:
					read1_filepath = base_dirpath + j
				match = re.search('_S1_L001_R2', j)
				if match:
					read2_filepath = base_dirpath + j
			p = Popen(['gunzip', read1_filepath], stderr=PIPE, stdout=PIPE)
			print p.communicate()
			p = Popen(['gunzip', read2_filepath], stderr=PIPE, stdout=PIPE)
			print p.communicate()

	#make a dic that will get the barcodes and sample names that belong to each seq run
	filein = open(sample_stats_filepath, "r")
	filein.readline()
	seq_date_dic = {}
	for i in filein:
		line = i[:-1].split('\t')
		patient_id = line[0]
		time_point = line[1]
		sample_id = line[2]
		if sample_id == '1_17' or sample_id == '5_2':
			continue
		barcode = line[5]
		seq_date = line[7]
		try:
			seq_date_dic[seq_date]['barcodes'].append(barcode)
			seq_date_dic[seq_date]['sample_names'].append(patient_id + '_' + time_point)
			seq_date_dic[seq_date]['sample_ids'].append(sample_id)
			seq_date_dic[seq_date]['time_points'].append(time_point)
			seq_date_dic[seq_date]['output_dirpaths'].append(output_dirpath + patient_id + '/')
		except KeyError:
			seq_date_dic[seq_date] = {}
			seq_date_dic[seq_date]['barcodes'] = [barcode]
			seq_date_dic[seq_date]['sample_names'] = [patient_id + '_' + time_point]
			seq_date_dic[seq_date]['sample_ids'] = [sample_id]
			seq_date_dic[seq_date]['time_points'] = [time_point]
			seq_date_dic[seq_date]['output_dirpaths'] = [output_dirpath + patient_id + '/']
		#if this sample has a duplicate, get the data for the duplicate as well
		if sample_id == '10_10' or sample_id == '2_17' or sample_id == '3_8':
			time_point = time_point + '.5'
			sample_id = sample_id + '-2'
			barcode = line[6]
			seq_date = line[8]
			try:
				seq_date_dic[seq_date]['barcodes'].append(barcode)
				seq_date_dic[seq_date]['sample_names'].append(patient_id + '_' + time_point)
				seq_date_dic[seq_date]['sample_ids'].append(sample_id)
				seq_date_dic[seq_date]['time_points'].append(time_point)
				seq_date_dic[seq_date]['output_dirpaths'].append(output_dirpath + patient_id + '/')
			except KeyError:
				seq_date_dic[seq_date] = {}
				seq_date_dic[seq_date]['barcodes'] = [barcode]
				seq_date_dic[seq_date]['sample_names'] = [patient_id + '_' + time_point]
				seq_date_dic[seq_date]['sample_ids'] = [sample_id]
				seq_date_dic[seq_date]['time_points'] = [time_point]
				seq_date_dic[seq_date]['output_dirpaths'] = [output_dirpath + patient_id + '/']
	for i in seq_date_dic:
		print i
		print seq_date_dic[i]['sample_names']
		sys.stdout.flush()
	print len(seq_date_dic)

    #now de-multiplex each seq run's data
	for i in os.listdir(input_illumina_data_dirpath):
		if i[0] == '.' or i[:6] == 'README':
			continue
		print i
		barcode_count_filepath = barcode_count_output_dirpath + i
		filename = i.split('_')
		seq_date = '/'.join([filename[-3], filename[-2], filename[-1][-2:]])
		#this is the directory that contains the fastq files (amoung other things)
		base_dirpath = '%s%s/Data/Intensities/BaseCalls/' % (input_illumina_data_dirpath, i)
		#need to find the fastq files we want from this directory
		for j in os.listdir(base_dirpath):
			match = re.search('_S1_L001_R1', j)
			if match:
				read1_filepath = base_dirpath + j
			match = re.search('_S1_L001_R2', j)
			if match:
				read2_filepath = base_dirpath + j
		barcodes = seq_date_dic[seq_date]['barcodes']
		sample_names = seq_date_dic[seq_date]['time_points']
		#add sample name for unknown seqs
		unknown_seqs_sample_name = i
		output_dirpaths = seq_date_dic[seq_date]['output_dirpaths']
		#add the directory for the seqs for which barcodes cannont be found
		unknown_seqs_dirpath = output_dirpath + 'unknown/'
		output_dirpaths.append(unknown_seqs_dirpath)
		seq_data = illumina_seq_data(read1_filepath, read2_filepath)
		seq_data.de_multiplex(barcodes=barcodes, barcode_start_pos=barcode_start_pos, sample_names=sample_names, read_with_barcode='read1', output_dirpath=output_dirpaths, write_as_fasta=False, write_unknown_barcode_seqs=unknown_seqs_sample_name, trim_barcode=trim_barcode, barcode_counts_output_filepath=barcode_count_filepath)

	return

if __name__ == '__main__':

	#below is a hard coded example of how we used this code. One will need to replace these files with there own data if they want to use this code.

	de_multiplex_all_seq_runs(input_illumina_data_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/raw_data_from_MiSeq/abr/successful_runs', sample_stats_filepath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/sample_info/abr.txt', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/raw_data_organized/abr', barcode_count_output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/abr_barcode_analysis/barcode_counts_foreach_seq_run', unzip=False)
