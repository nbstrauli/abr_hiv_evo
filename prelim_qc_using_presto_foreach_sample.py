#!/usr/bin/python
#$ -S /usr/bin/python
#$ -e error
#$ -o out
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=16G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=200G
#$ -l h_rt=336:00:00

import sys
sys.path.insert(0, './')
import os
from subprocess import Popen, PIPE
from sequence_sample_class import sequence_sample
import get_divergence_and_dNdS_hiv_sample
from sequence_sample_set_class import sequence_sample_set
import find_outliers_by_clustering

def QC_all_samples(input_sequence_dirpath, fwd_primer_fasta_filepath, rev_primer_fasta_filepath, output_dirpath, clean_data_output_dirpath, abr_or_hiv='hiv', fasta_dirpath=None, consensus_output_dirpath=None):
    """This script will walk through each of the patient/timepoint samples and run the QC script 'prelim_qc_using_presto.bash' in order to do basic quality control on the data.
    'clean_data_output_dirpath' = this is the data that only the final fastq file (after all the QC) will be written.
    '[fwd|rev]_primer_fasta_filepath' = Paths to the forward and reverse fasta formatted files. If these = None, then it is assumed that no primer information exists and the porgram 'prelim_qc_using_presto_noPrimerMask.bash' will be used (this is what we do for the AbR data).
    'abr_or_hiv' = Determines what data type we are dealing with. If 'hiv', then the cluster is not used, and the FWD and REV primer seqs are known. If 'abr' then we use the cluster (and all the complications that go with that), and the primers seqs are unknown.
    'fasta_dirpath' = If defined (default, None) then this should be the path to a directory that will contain the sequence data (post QC) in fasta format.
    'consensus_output_dirpath' - If defined (default, None), this will be the path to the directory that will contain the concesus seq for each of the patients. This is currently only applicable to HIV data.
    """

    ######################
    ##### parameters #####
    ######################

    #This is name of the file that is produced at the 
    #final step of the QC process
    final_QCed_filename = 'QC_atleast-2.fastq'
    #this is the basename of the read 1 and read 2 fastq input files
    # read1_basename = 'read1'
    read1_basename = '1'
    # read2_basename = 'read2'
    read2_basename = '2'
    count_attribute_name = 'DUPCOUNT'
    # path_to_blast = '/Users/nstrauli/tools/ncbi-blast-2.7.1+/bin'
    path_to_blast = '/netapp/home/nstrauli/tools_c/ncbi-blast-2.2.28+/bin'
    #path_to_blast_db = '/Users/nstrauli/tools/ncbi-blast-2.7.1+/databases/hxb2_env/env'
    # path_to_blast_db = '/Users/nstrauli/tools/ncbi-blast-2.7.1+/databases/hiv_subtype_ref_seqs/env'
    path_to_blast_db = '/netapp/home/nstrauli/tools_c/ncbi-blast-2.2.28+/databases/hiv_subtype_ref_seqs/env'
    # path_to_blast_db_all_env = '/Users/nstrauli/tools/ncbi-blast-2.7.1+/databases/all_hiv_env_lanl/all_env'
    path_to_blast_db_all_env = '/netapp/home/nstrauli/tools_c/ncbi-blast-2.2.28+/databases/all_hiv_env_lanl/all_env'
    percent_ident_to_ref_cutoff = '70'
    max_percent_ident_to_any_env_cutoff = '99'
    # path_to_hxb2_ref_seq = '/Users/nstrauli/data/abr_hiv_coevo/seq_data/ref_seqs/hiv/HXB2_env_nucleotide_reference_seq.fasta'
    path_to_hxb2_ref_seq = '/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/ref_seqs/hiv/HXB2_env_nucleotide_reference_seq.fasta'
    #Below are the boundaries to the target reference sequence (HXB2) that will be searched using a pairwise alignment program. These boundaries were found by aligning our fwd and rev primers to the HXB2 ref seq and finding the positions that these primers map to the seq. Did this manually using the EMBL gui for their 'water' program.
    start_of_hxb2_to_search = '704'
    end_of_hxb2_to_search = '1120'
    #vsearch path for rocinante
    # path_to_vsearch = '/Users/nstrauli/tools/vsearch-2.4.0-macos-x86_64/bin/vsearch'
    #vsearch path for gibbon
    # path_to_vsearch = '/Users/nstrauli/tools/vsearch-2.4.3-macos-x86_64/bin/vsearch'
    #vsearch path for cluster
    path_to_vsearch = '/netapp/home/nstrauli/tools_c/vsearch-2.4.3-linux-x86_64/bin/vsearch'
    ######################
    ##### parameters #####
    ######################

    if input_sequence_dirpath[-1] != '/':
        input_sequence_dirpath += '/'
    if output_dirpath[-1] != '/':
        output_dirpath += '/'
    if not os.path.exists(output_dirpath):
        os.makedirs(output_dirpath)
    if clean_data_output_dirpath[-1] != '/':
        clean_data_output_dirpath += '/'
    if not os.path.exists(clean_data_output_dirpath):
        os.makedirs(clean_data_output_dirpath)
    if fasta_dirpath:
        if fasta_dirpath[-1] != '/':
            fasta_dirpath += '/'
        if not os.path.exists(fasta_dirpath):
            os.makedirs(fasta_dirpath)
    if consensus_output_dirpath:
        if consensus_output_dirpath[-1] != '/':
            consensus_output_dirpath += '/'
        if not os.path.exists(consensus_output_dirpath):
            os.makedirs(consensus_output_dirpath)
    #do QC
    for i in os.listdir(input_sequence_dirpath):
        if i[0] == '.' or i == 'README':
            continue
        output_dirpath_patient = output_dirpath + i + '/'
        if not os.path.exists(output_dirpath_patient):
            os.makedirs(output_dirpath_patient)
        clean_output_dirpath_patient = clean_data_output_dirpath + i + '/'
        if not os.path.exists(clean_output_dirpath_patient):
            os.makedirs(clean_output_dirpath_patient)
        if fasta_dirpath:
            fasta_dirpath_patient = fasta_dirpath + i + '/'
            if not os.path.exists(fasta_dirpath_patient):
                os.makedirs(fasta_dirpath_patient)
        for j in os.listdir(input_sequence_dirpath + i):
            if j[0] == '.' or j[:6] == 'README':
                continue

            print ''
            print '################################'
            print 'patient:', i
            print 'time-point:', j
            print '################################'
            print ''

            output_dirpath_tpoint = output_dirpath_patient + j + '/'
            if not os.path.exists(output_dirpath_tpoint):
                os.makedirs(output_dirpath_tpoint)
            if fasta_dirpath:
                fasta_filepath_tpoint = fasta_dirpath_patient + j + '.fasta'
            input_filepath_read1 = input_sequence_dirpath + i + '/' + j + '/' + read1_basename + '.fastq'
            input_filepath_read2 = input_sequence_dirpath + i + '/' + j + '/' + read2_basename + '.fastq'
            if abr_or_hiv == 'hiv':
                clean_output_filepath_tpoint = clean_output_dirpath_patient + j + '.fasta'
                p = Popen(['bash', 'prelim_qc_using_presto.bash', input_filepath_read1, input_filepath_read2, fwd_primer_fasta_filepath, rev_primer_fasta_filepath, output_dirpath_tpoint, 'QC', clean_output_filepath_tpoint, path_to_blast, path_to_blast_db, percent_ident_to_ref_cutoff, path_to_hxb2_ref_seq, start_of_hxb2_to_search, end_of_hxb2_to_search, path_to_vsearch, path_to_blast_db_all_env, max_percent_ident_to_any_env_cutoff, '/netapp/home/nstrauli/tools_c/presto-0.5.4/bin'], stderr=PIPE, stdout=PIPE)
                out_err = p.communicate()
                print out_err[0]
                print out_err[1]
            elif abr_or_hiv == 'abr':
                clean_output_filepath_tpoint = clean_output_dirpath_patient + j + '.fastq'
                p = Popen(['qsub', 'prelim_qc_using_presto_abr.bash', input_filepath_read1, input_filepath_read2, output_dirpath_tpoint, 'QC', clean_output_filepath_tpoint, '/netapp/home/nstrauli/tools_c/presto-0.5.4/bin'], stderr=PIPE, stdout=PIPE)
                out_err = p.communicate()
                print out_err[0]
                print out_err[1]
                if fasta_dirpath:
                    sample = sequence_sample(clean_output_filepath_tpoint, count_attribute_name='DUPCOUNT')
                    sample.write_full_data_to_disk_fasta(output_filepath=fasta_filepath_tpoint)

    if abr_or_hiv == 'hiv' and consensus_output_dirpath:
        #get divergence and selection values for each seq using in-house scripts, and write to final output file. This script has a lot of hardcoded parameters of its own. See the script to view these parameters.
        get_divergence_and_dNdS_hiv_sample.run(input_fasta_dirpath=clean_data_output_dirpath, consensus_output_dirpath=consensus_output_dirpath, fasta_output_dirpath=fasta_dirpath)
    
    return

if __name__ == '__main__':
    
    #below are hard coded examples of how we used this code. One will need to replace these files with there own data if they want to use this code. 
    
    #for HIV, cluster:
    # QC_all_samples(input_sequence_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/raw_data_organized/hiv', fwd_primer_fasta_filepath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/primers/hiv/fwd.fasta', rev_primer_fasta_filepath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/primers/hiv/rev.fasta', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/QCed_data/hiv', clean_data_output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/QCed_data_clean_final/hiv', abr_or_hiv='hiv', fasta_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_annotated_fasta_files', consensus_output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/1st_time_point_consensus_seqs/hiv/most_abundant_method')
    #for HIV, local:
    # QC_all_samples(input_sequence_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/raw_data_organized/hiv', fwd_primer_fasta_filepath='/Users/nstrauli/data/abr_hiv_coevo/primers/hiv/fwd.fasta', rev_primer_fasta_filepath='/Users/nstrauli/data/abr_hiv_coevo/primers/hiv/rev.fasta', output_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/QCed_data/hiv', clean_data_output_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/QCed_data_clean_final/hiv', abr_or_hiv='hiv', fasta_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_annotated_fasta_files', consensus_output_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/1st_time_point_consensus_seqs/hiv/most_abundant_method')
    #for AbR:
    # QC_all_samples(input_sequence_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/raw_data_organized/abr', fwd_primer_fasta_filepath=None, rev_primer_fasta_filepath=None, output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/QCed_data/abr', clean_data_output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/QCed_data_clean_final/abr', abr_or_hiv='abr', fasta_dirpath=None, consensus_output_dirpath=None)
