from sequence_sample_set_class import sequence_sample_set
import os
from master_directory_class import master_directory
import itertools
import sys

def get_pi(input_dirpath, output_dirpath, temp_dirpath, count_attribute_name='count', method='vsearch', try_again=True):
    """
    """
    
    ###### parameters ######
    #min sample size for partis dataset is 116394.0
    #min sample size for changeo dataset is 116408.0
    #min sample size for irepertoire dataset is 140899.0
    #min sample size for hiv dataset is 2274.0
    #downsamp_to = None
    #downsamp_to = 116394.0
    # downsamp_to = 116408.0
    #downsamp_to = 140899.0
    downsamp_to = 2274.0
    method = 'vsearch'
    try_again = True
    num_downsamp_trials = 10
    num_alignments_per_job = 1000000
    path_to_needle = '/netapp/home/nstrauli/tools_c/EMBOSS-6.6.0/emboss/needle'
    path_to_vsearch = '/netapp/home/nstrauli/tools_c/vsearch-2.4.3-linux-x86_64/bin/vsearch'
    #num_parallel_cores = 12
    num_parallel_cores = 1
    avoid_dirpath = '/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/partis_annotated_fasta_files/unknown'
    ########################

    if input_dirpath[-1] != '/':
        input_dirpath += '/'
    if output_dirpath[-1] != '/':
        output_dirpath += '/'
    if not os.path.exists(output_dirpath):
        os.makedirs(output_dirpath)
    if not os.path.exists(temp_dirpath):
        os.makedirs(temp_dirpath)
    d = master_directory(input_dirpath)
    d.mirror_directory_tree(output_dirpath, avoid_dirpath=avoid_dirpath)

    #if using needle, then do one patient at a time
    if method == 'needle':
        for i in os.listdir(input_dirpath):
            if i[0] == '.' or i[:6] == 'README':
                continue
            patient_input_dirpath = input_dirpath + i + '/'
            sample_set = sequence_sample_set(patient_input_dirpath, count_attribute_name, avoid_dirpath=avoid_dirpath)
            diversity_pi_values = sample_set.down_sample_and_calc_pi_compCluster(downsamp_to=downsamp_to, num_downsamp_trials=num_downsamp_trials, num_alignments_per_job=num_alignments_per_job, method=method, path_to_needle=path_to_needle, path_to_vsearch=path_to_vsearch, temp_dirpath=temp_dirpath, try_again=try_again)
            for j in diversity_pi_values:
                path = j.split('/')
                patient = path[-2]
                tpoint = path[-1][:-6]
                output_filepath = '%s%s/%s' % (output_dirpath, patient, tpoint)
                fileout = open(output_filepath, "w")
                fileout.write('trial\tpi\n')
                count = 0
                for k in diversity_pi_values[j]:
                    count += 1
                    fileout.write('%s\t%s\n' % (count, k))
                mean_pi = sum(diversity_pi_values[j]) / len(diversity_pi_values[j])
                fileout.write('mean\t%s\n' % mean_pi)
                fileout.close()

    #else if using vsearch then do all patients at once
    elif method == 'vsearch':
        sample_set = sequence_sample_set(input_dirpath, count_attribute_name, avoid_dirpath=avoid_dirpath)
        diversity_pi_values = sample_set.down_sample_and_calc_pi_compCluster(downsamp_to=downsamp_to, num_downsamp_trials=num_downsamp_trials, num_alignments_per_job=num_alignments_per_job, num_parallel_cores=num_parallel_cores, method=method, path_to_needle=path_to_needle, path_to_vsearch=path_to_vsearch, temp_dirpath=temp_dirpath, try_again=try_again)

        print diversity_pi_values

        for j in diversity_pi_values:
            path = j.split('/')
            patient = path[-2]
            tpoint = path[-1][:-6]
            output_filepath = '%s%s/%s' % (output_dirpath, patient, tpoint)
            fileout = open(output_filepath, "w")
            fileout.write('trial\tpi\n')
            count = 0
            for k in diversity_pi_values[j]:
                count += 1
                fileout.write('%s\t%s\n' % (count, k))
            mean_pi = sum(diversity_pi_values[j]) / len(diversity_pi_values[j])
            fileout.write('mean\t%s\n' % mean_pi)
            fileout.close()
    return

if __name__ == '__main__':

    #below are hard coded examples of how we used this code for 4 different datasets (3 are antibody data that were processed with differnet pipelines). One will need to replace these files with there own data if they want to use this code.

    # if sys.argv[1] == 'partis':
    #     get_pi(input_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/partis_annotated_fasta_files', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/diversity/abr/partis/pi', temp_dirpath='/netapp/home/nstrauli/data/abr_hiv_coevo/temp_stuff', count_attribute_name='count')
    # elif sys.argv[1] == 'changeo':
    #     get_pi(input_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/changeo_annotated_fasta_files', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/diversity/abr/changeo/pi', temp_dirpath='/netapp/home/nstrauli/data/abr_hiv_coevo/temp_stuff', count_attribute_name='DUPCOUNT')
    # elif sys.argv[1] == 'irep':
    #     get_pi(input_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/irep_fasta_files/abr', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/diversity/abr/irep/pi', temp_dirpath='/netapp/home/nstrauli/data/abr_hiv_coevo/temp_stuff', count_attribute_name='count')
    # elif sys.argv[1] == 'hiv':
    #     get_pi(input_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_fasta_files', output_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/diversity/hiv/pi', temp_dirpath='/netapp/home/nstrauli/data/abr_hiv_coevo/temp_stuff', count_attribute_name='DUPCOUNT')
    # else:
    #     print 'Wha?'
