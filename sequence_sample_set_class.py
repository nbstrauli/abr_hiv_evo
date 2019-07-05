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
from master_directory_class import master_directory
from sequence_sample_class import sequence_sample
import time
from subprocess import Popen, PIPE
import subprocess
import tempfile
import itertools
import random
import numpy

class sequence_sample_set(object):
    """
    This class is for performing functions on sets of sequence files. The sequence file set could all be from the same patient, study, or what have you. They just need to be conceptually linked in some way. The methods use the 'sequence_sample' class quite a bit to perform functions on each individual sample. Each file representing a sequence sample, should be either in fasta or fastq format (quality information is ignored however). The header for each sequence should start with an ID (could be anything, ideally unique) and then should be followed by a series of attributes. Each attribute should be delimited by a '|'. Each attribute should have a name followed by and '=', followed by the value of the attribute. An example of a well formatted header: ">seq_id|vgene=IGHV1|dgene=NA|jgene=IGHJ3|count=65\n". Here the attributes are 'vgene', 'dgene', 'jgene', and 'count'. The first entry in the header is assumed to be the sequence unique ID.

    Attributes:
    dirpath - The path to the directory that contains all the sequence samples. There could be other sub-directories within this master-directory that contain other sequence samples. All ending in a fasta or fastq suffix will be included in the set. If this parameter equals None then 'list_of_filepaths' is considered
    sample_filepaths - A list of the absolute filepaths for each of the samples in the set. In no particular order.
    count_attribute_name - this gives the name that corresponds to the count of each sequence entity.
    file_suffix - This gives the suffix that is used in the input seq file data. This is determined by the first filepath that is listed in self.sample_filepaths
    avoid_files_with - If defined (defualt, None), this means that any file's basename that has this string in it (not including the file suffix, i.e. ',fasta') will be omitted from the sample set. This parameter is included because often there will be duplicates in the data set that aren't meant to be included for some methods.
    """
    def __init__(self, dirpath, count_attribute_name=None, list_of_filepaths=None, avoid_dirpath=None, avoid_files_with=None):
        """
        list_of_fileapths - This is a list of filepaths to the files that create the sequence set. They do not necessarily need to be in the same directory. This parameter is only considered if 'dirpath' is equal to None. If they are both equal to None, then throws an error.
        avoid_dirpath - This means that the subdirs and files in this directory path will be omitted from the set. If equals None (default) then nothing is avoided. This parameter is only considered if 'dirpath' is defined.
        avoid_files_with - If defined (defualt, None), this means that any file's basename that has this string in it (not including the file suffix, i.e. '.fasta') will be omitted from the sample set. This parameter is included because often there will be duplicates in the data set that aren't meant to be included for some methods.
        """

        if dirpath == None and list_of_filepaths == None:
            print "Both 'dirpath' and 'list_of_filepaths' variables are equal to None type. Aborting."
            return
        self.count_attribute_name = count_attribute_name
        if dirpath and dirpath[-1] != '/':
            dirpath += '/'
        self.dirpath = dirpath
        self.avoid_files_with = avoid_files_with
        if dirpath:
            d = master_directory(self.dirpath)
            d.get_filepaths_containing_filetypes(['fasta', 'fa', 'fastq', 'fq'], avoid_dirpath=avoid_dirpath, avoid_files_with=avoid_files_with)
            self.sample_filepaths = d.files_containing_filetype
            #make these filepaths into a list
            self.sample_filepaths = [i for i in self.sample_filepaths]
        else:
            self.sample_filepaths = list_of_filepaths
        self.file_suffix = self.sample_filepaths[0].split('.')[-1]
        return

    def down_sample_and_calc_pi_compCluster(self, downsamp_to=None, num_downsamp_trials=10, num_alignments_per_job=100000, method='needle', path_to_needle='needle', path_to_vsearch=None, temp_dirpath=None, num_parallel_cores=12, try_again=False, no_array_job=False):
        """
        This method first finds the lowest sample size in the sample set, then downsamples each sample to that size, then finds the diversity of the samples after the downsampling event.
        downsamp_to - The number to downsample each of the samples to. If this is None type (default) the method will cycle through each of the samples to find the one that has the lowest sample size and use this value as that which all samples should be down sampled to. Make sure this value is lower or equal to the sample size of all samples in the set. If this equals 'Nope', then no down-sampling occurs, and thus no downsamp trials are necessary. So, if this is the case, 'num_downsamp_trials' will be set to 1.
        num_downsamp_trials - This gives the number of downsampling events that will be carried out for each sample.
        'path_to_needle' - This is the path to the needle program by EMBOSS that impliments the needleman-wunsch global alignment algorithm. This is what we use to get the genetic distance between two seqs.
        'method' - This gives the method that will be used to get the genetic distance from the pairs of seqs. The exceptable values for this are:
            'needle' - This means the needleman-wunsch global alignment algorithm will be used in a program called 'needle' in the EMBOSS package
            'pairwise2' - This means the 'pairwise2' Biopython package will be used, which also implements the needleman-wunsch global alignment algorithm.
            'vsearch' - This means the vsearch program will be used. This program is capable of doing an all pairwise alignment by itself, so the flow of this method is significantly changed if this is selected.
        path_to_vsearch - This is the path to the vsearch excecutable. If this is None (default) then the path is assumed to be 'vsearch' (i.e. in $PATH). Ignored if method!='vsearch'
        temp_dirpath - This gives the path to the temporary directory that will contain the fasta files to be submitted to the aligner (if the chosen aligner needs this), and the sub_sum_of_distances values outputed by 'pi_calculator_compCluster'. The temp fasta files can get quite large, so make sure that there is enough space where ever this path leads. If temp_dirpath=None (default) then a temp dir in the current working directory is made for this.
        num_parallel_cores - This gives the number of parallel jobs that vsearch can use. vsearch's default is to use as many cores as it can, and this slows down the cluster for others. So, we need to set this to a fixed amount, and then let the cluster know when submitting array job. Default is 12.
        try_again - Boolean. If True (not default), then the script will keep trying to calculate pi until it works. It seems that vsearch is a bit buggy and can write uninterpretable lines of output. This setting (when set to True) will simply keep trying until it works. Only relevant if method='vsearch'.
        no_array_job - If True (default, False), this will inform the method to not submit an array job to calculate pi. This function is currently only supported if method='vsearch'.

        Output:
        Returns 'diversity_pi_values', which is a dictionary, where each index is the filepath of a given sequence sample, and the definition of that index is a list pi values that were calculated for that sample. Each individual pi value calculated is the result of one of the downsampling trials for that sample.
        """
        #first get lowest sample size (as in lowest num of seqs)
        if downsamp_to == None:
            sample_sizes = []
            for i in self.sample_filepaths:
                sample = sequence_sample(i, self.count_attribute_name)
                sample_size = sample.total
                sample_sizes.append(sample_size)
                print i
                print 'sample size:', sample_size
                print '################################################'
            min_sample_size = min(sample_sizes)
            print 'Min sample size is:', min_sample_size
            print ''
        elif downsamp_to == 'Nope':
            num_downsamp_trials = 1
            min_sample_size = 'Nope'
        else:
            min_sample_size = float(downsamp_to)

        #now downsample and calculate pi

        #if using needle or pairwise2
        if method == 'needle' or method == 'pairwise2':
            diversity_pi_values = {}
            for i in self.sample_filepaths:
                print 'calculating pi for sample:'
                print i
                pis = []
                for j in xrange(num_downsamp_trials):
                    print '\ttrial:', j+1
                    sample = sequence_sample(i, self.count_attribute_name)
                    #if the sample is the min, then no
                    #need to downsample
                    if sample.total == min_sample_size or downsamp_to == 'Nope':
                        print '\t\tcalculating pi'
                        start_time = time.time()
                        pi = sample.calc_diversity_pi_compCluster(num_alignments_per_job, method, path_to_needle, temp_dirpath)
                        end_time = time.time()
                        print '\t\tpi:', pi
                        print '\t\ttime elapsed while calculating pi:', end_time-start_time
                        pis.append(pi)
                        break
                    #else, downsample
                    else:
                        print '\t\tdownsampling'
                        start_time = time.time()
                        sample.down_sample(min_sample_size)
                        end_time = time.time()
                        print '\t\ttime elapsed while downsampling:', end_time-start_time
                        print ''
                        print '\t\tcalculating pi'
                        start_time = time.time()
                        pi = sample.calc_diversity_pi_compCluster(num_alignments_per_job, method, path_to_needle, temp_dirpath)
                        end_time = time.time()
                        print '\t\tpi:', pi
                        print '\t\ttime elapsed while calculating pi:', end_time-start_time
                        pis.append(pi)
                mean_pi = sum(pis) / float(len(pis))
                diversity_pi_values[i] = pis
                print 'mean pi across downsampling trials:', mean_pi
                print '#################################################'
        
        #else if using vsearch
        elif method == 'vsearch':
            #If no temp dirpath given, make one in CWD
            if temp_dirpath == None:
                temp_dirpath = tempfile.mkdtemp(dir=os.getcwd()) + '/'
            #Else, make the temp dir within the provided directory
            else:
                temp_dirpath = tempfile.mkdtemp(dir=temp_dirpath) + '/'
            num_jobs = 0
            for i in self.sample_filepaths:
                for j in xrange(num_downsamp_trials):
                    num_jobs += 1
                    indiv_job_filepath = '%s%s_input_filepath' % (temp_dirpath, num_jobs)
                    fileout = open(indiv_job_filepath, "w")
                    fileout.write(i)
                    fileout.close()
            #if it is desired to not use array
            if no_array_job:
                for i in xrange(num_jobs):
                    p = Popen(['python', 'sequence_sample_set_class.py', 'vsearch_pi_calculator_compCluster', str(min_sample_size), self.count_attribute_name, path_to_vsearch, temp_dirpath, str(num_parallel_cores), str(i+1)], stdout=PIPE, stderr=PIPE)
                    out, err = p.communicate()
                    print out
                    print err
            #else, submit array job
            else:
                if num_parallel_cores == 1:
                    p = Popen(['qsub', '-t', '1-%s' % num_jobs, 'sequence_sample_set_class.py', 'vsearch_pi_calculator_compCluster', str(min_sample_size), self.count_attribute_name, path_to_vsearch, temp_dirpath, str(num_parallel_cores), 'None'], stdout=PIPE, stderr=PIPE)
                else:
                    p = Popen(['qsub', '-t', '1-%s' % num_jobs, '-pe', 'smp', str(num_parallel_cores), 'sequence_sample_set_class.py', 'vsearch_pi_calculator_compCluster', str(min_sample_size), self.count_attribute_name, path_to_vsearch, temp_dirpath, str(num_parallel_cores), 'None'], stdout=PIPE, stderr=PIPE)
                out, err = p.communicate()

                print out
                print err

                job_id = out.split()[2].split('.')[0]
                #sleep till job is done
                good_to_go = False
                while good_to_go == False:
                    p = Popen(['qstat'], stdout=PIPE, stderr=PIPE)
                    jobs = p.communicate()[0].split('\n')
                    jobs = [i.split()[0] for i in jobs[2:-1]]
                    if job_id in jobs:
                        time.sleep(10)
                    else:
                        good_to_go = True
            #when array jobs are done, retrieve pi values
            diversity_pi_values = {}
            failed_filepaths = []
            for i in os.listdir(temp_dirpath):
                if i.split('_')[-1] != 'pi':
                    continue
                filein = open('%s%s' % (temp_dirpath, i), "r")
                line = filein.readline().split('\t')
                filein.close()
                sample_filepath = line[0]
                try: #try to retrieve pi
                    pi = float(line[1])
                    try: #if pi calc was successfull add it to 'diversity_pi_values'
                        diversity_pi_values[sample_filepath].append(pi)
                    except KeyError:
                        diversity_pi_values[sample_filepath] = [pi]
                except ValueError: #except if pi calc was unsuccessful
                    failed_filepaths.append(sample_filepath)

            #keep trying to get the failed pi values, using recursion
            if try_again and len(failed_filepaths)>0:
                print "Some attempts to calculate pi were unsuccessful. Vsearch likely failed, so will try again."
                print "Number of failed samples:", len(failed_filepaths)
                for i in sorted(failed_filepaths):
                    print i
                sys.stdout.flush()
                samples = sequence_sample_set(dirpath=None, count_attribute_name=self.count_attribute_name, list_of_filepaths=failed_filepaths)
                diversity_pi_values_failed = samples.down_sample_and_calc_pi_compCluster(downsamp_to=min_sample_size, num_downsamp_trials=1, num_alignments_per_job=num_alignments_per_job, method=method, path_to_needle=path_to_needle, path_to_vsearch=path_to_vsearch, temp_dirpath=temp_dirpath, num_parallel_cores=num_parallel_cores, try_again=try_again)

                #upadate the diversity_pi_values dic
                for i in diversity_pi_values_failed:
                    try:
                        diversity_pi_values[i] += diversity_pi_values_failed[i]
                    except KeyError:
                        diversity_pi_values[i] = diversity_pi_values_failed[i]

            subprocess.call(['rm', '-r', temp_dirpath])
            #calc mean pi foreach sample, and also make sure
            #there are the correct number of pi values in each
            for i in diversity_pi_values:
                len_pis = len(diversity_pi_values[i])
                if len_pis != num_downsamp_trials:
                    print 'found sample that did not have pi values for each downsampling trial:', i
                mean_pi = sum(diversity_pi_values[i]) / float(len_pis)

        return diversity_pi_values

    def get_immune_reads_with_changeo_foreach_sample(self, output_dirpath, path_to_changeo='/Users/nstrauli/tools/changeo-0.3.9/', path_to_igblast='/Users/nstrauli/tools/ncbi-igblast-1.8.0/', temp_dirpath=None, imgt_ref_seqs=None, add_germline=None, add_selection=False, use_comp_cluster=False, add_divergence=False, overwrite_output=True, translate_VDJ=False, remove_seqs_with_stop=False):
        """
        This method uses the similarily named method 'get_immune_reads_with_changeo' in the sequence_sample class to get the immune reads for each of the samples in the dataset. Running this method will create a new dataset that will consist of immune reads.
        output_dirpath - This is the path to the directory that will contain the immune data. It will be fasta formatted.
        path_to_changeo - This is the path to the change-o excecutable. It should be a directory that contains 'bin', 'changeo/', etc.
        path_to_igblast - This gives the path to the directory that contains the excecutables for IgBLAST. That is, the IgBLAST excecutables should be in a directory called 'bin' in this directory. The other file such as 'database', 'internal_data', and 'optional_file' should be in this directory as well.
        temp_dirpath - This gives the directory for which temporary files will be written. If None (default) then the temp dirpath is the same as 'output_dirpath'.
        imgt_ref_seqs - This is the path to the files that contain the germline reference sequences for the genes that were used in the IgBLAST database. If this is None (default), then this will equal path_to_igblast + 'ref_seqs_from_imgt/'
        add_germline - This instructs if and what type of germline sequence to add for each sequence. If defined, this should be a list of strings. Acceptable values are:
            'dmask' - This means that the germline seq created will have the D-region masked (i.e. all N's). This is needed in order to calculate selection values for each seq using Baseline.
            'full' - This gives the complete germline (i.e. includes V gene).
            'vonly' - This will return only the germline of the V gene.
            None - This will make it so the germline is not returned (default).
        add_selection - This indicates whether or not the method should invoke Baseline to add selection values to each of the seqs in the data. If this is True (default, False), then the method will run Baseline using the 'shazam' R package to get selection values.
        use_comp_cluster - This indicates if the computational cluster should be used. If True (default, False), then the SGE system will be used to submit jobs to a computational cluster. One job will be submitted for each sample (regardless of the number of seqs in each file).
        add_divergence - If defined (default, False), this will instruct the method to also include a divergence from the inferred unmutated sequence in the immune reads info. Acceptable values are:
            'shazam' - Causes divergence to be calculated using the 'calcObservedMutations' function to be run in shazam. This is an R script and is quite slow, but will give the number of replacement and silent mutations in addition to the total.
            'python' - Causes divergence to be calculated using a simple in house python script that gets the total number of muts between the germline and observed seq. These seqs must be the same length (i.e. does not do an alignment), but they should be because they are imgt numbered from changeo. It also ignores positions that are '.', '-', or 'N' in either the observed or germline seqs.
        overwrite_output - If True (default), then this will write the results to 'output_filepath' regardless of whether or not it exists. If False, the this method will first check to see if the output_filepath exists, and if it does, then it will abort with an error message.
        translate_VDJ - Boolean. If True (default, Flase), this will cause the method to add a field to the header that will be called 'AA_VDJ_sequence' that will give the translated amino acid sequence of the nucleotride sequence contained in the 'SEQUENCE_VDJ' field.
        remove_seqs_with_stop - Boolean. If True (default, False), this will cause the method to not include any seqs that were labled as containing stop codons (by looking at the 'STOP' field given by changeo.)
        """
        if output_dirpath[-1] != '/':
            output_dirpath += '/'
        if not os.path.exists(output_dirpath):
            os.makedirs(output_dirpath)
        if temp_dirpath:
            if temp_dirpath[-1] != '/':
                temp_dirpath += '/'
            if not os.path.exists(temp_dirpath):
                os.makedirs(temp_dirpath)
        else:
            temp_dirpath = output_dirpath
        if self.dirpath:
            dir = master_directory(self.dirpath)
            dir.mirror_directory_tree_with_files(out_dirpath=output_dirpath, only_include_filetypes=['fasta', 'fa', 'fastq', 'fq'], include_file_suffix=False)
            output_filepaths = [i + '.fasta' for i in dir.mirrored_filepaths]
        else:
            output_filepaths = [output_dirpath+i for i in self.sample_filepaths]
        if use_comp_cluster:
            #write input and output filepaths in temp file
            temp_filepath = temp_dirpath + 'temp_input_file_' + str(random.random())
            fileout = open(temp_filepath, "w")
            for i, j in itertools.izip(sorted(self.sample_filepaths), sorted(output_filepaths)):
                fileout.write('%s\t%s\n' % (i, j))
            fileout.close()
            p = Popen(['qsub', '-t', '1-'+str(len(output_filepaths)), 'sequence_sample_set_class.py', 'get_immune_reads_with_changeo_compCluster', temp_filepath, self.count_attribute_name, path_to_changeo, path_to_igblast, '/scratch/', str(imgt_ref_seqs), ','.join(add_germline), str(add_selection), str(use_comp_cluster), str(add_divergence), str(overwrite_output), str(translate_VDJ), str(remove_seqs_with_stop)], stdout=PIPE, stderr=PIPE)
            submit_job_id = set([p.communicate()[0].split()[2].split('.')[0]])
            #sleep till job is done
            good_to_go = False
            while good_to_go == False:
                p = Popen(['qstat'], stdout=PIPE, stderr=PIPE)
                jobs = p.communicate()[0].split('\n')
                job_ids = set()
                if len(jobs) >= 1:
                    for i in jobs[2:-1]:
                        job_id = i.split()[0]
                        job_ids.update([job_id])
                    if len(job_ids.intersection(submit_job_id)) == 0:
                        good_to_go = True
                else:
                    good_to_go = True
                time.sleep(10)
            subprocess.call(['rm', temp_filepath])
        else:
            for i, j in itertools.izip(sorted(self.sample_filepaths), sorted(output_filepaths)):
                print '$$$$$$$$$$$$$$$$$$'
                print i
                print j
                sample = sequence_sample(i, count_attribute_name=self.count_attribute_name)
                sample.get_immune_reads_with_changeo(output_filepath=j, path_to_changeo=path_to_changeo, path_to_igblast=path_to_igblast, temp_dirpath=temp_dirpath, imgt_ref_seqs=imgt_ref_seqs, add_germline=add_germline, add_selection=add_selection, add_divergence=add_divergence, overwrite_output=overwrite_output, translate_VDJ=translate_VDJ, remove_seqs_with_stop=remove_seqs_with_stop)
        return

    def translate_each_seq(self, output_dirpath, coding_frame_start=1, ref_seq=None, temp_dirpath=None, remove_seqs_with_stop_codons=False):
        """
        This method will translate each of the nucleotide seqs in the data. A challenge with this is to know where this coding frame is in the sequence. The coding frame can start at either the 1st, 2nd, or 3rd nucleotide in a given seq. The resulting translated amino acid seqs will be stored in memory as a string in self.data['other']['amino_acid_seq'].
        output_dirpath - This gives the path to the directory that will contain the sequence data for all the samples, with the translated amino acid seqs included in the data.
        coding_frame_start - This gives the coding frame of each of the sequences. Acceptable values are:
            [1,2,3] - (read as 1, 2, or 3). This gives the location of the start of the coding frame in the sequence. For example if this is 1, then that means that the coding frame starts at the 1st nucleotide in the seq.
            'relative_to_reference' - This means that the start of the coding frame should be determined by mapping each sequence to a reference sequence, where the coding frame is known. This is done by using 'needle' from the emboss package to globally align each seq (pairwise) to the reference seq. Using this alignment, the coding frame start location from the reference is then mapped onto a given seq in the data. If this approach is used, then the 'ref_seq' parameter must be defined.
        ref_seq - If defined (default, None), then this will be the path to the fasta file that contains the reference sequence for each of the seqs in the data. This will be used as a reference for translation.
        temp_dirpath - Gives the path to a directory for which temp data can be written. If None (default), then uses the current working directory.
        remove_seqs_with_stop_codons - If True (default, False) this will remove seq entries in the data that, when translated, have stop codons.
        """
        if output_dirpath[-1] != '/':
            output_dirpath += '/'
        if not os.path.exists(output_dirpath):
            os.makedirs(output_dirpath)
        file_structure = master_directory(self.dirpath)
        file_structure.mirror_directory_tree_with_files(out_dirpath=output_dirpath, only_include_filetypes=['fasta', 'fa', 'fastq', 'fq'], include_file_suffix=False)
        output_filepaths = file_structure.mirrored_filepaths
        for file_out, file_in in itertools.izip(sorted(output_filepaths), sorted(self.sample_filepaths)):
            print file_out
            print file_in
            print '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
            sample = sequence_sample(filepath=file_in, count_attribute_name=self.count_attribute_name)
            sample.translate_each_seq(coding_frame_start=coding_frame_start, ref_seq=ref_seq, temp_dirpath=temp_dirpath, remove_seqs_with_stop_codons=remove_seqs_with_stop_codons)
            output_filepath = file_out + '.fasta'
            sample.write_full_data_to_disk_fasta(output_filepath=output_filepath)
        return

    def get_dN_dS_and_divergence(self, ref_seq, output_dirpath, method='counting', temp_dirpath=None, method_counting='equal_weights'):
        """
        This method will get both the non-synonymous and synonymous divergence, as well as dN/dS, relative to some reference sequence, for all seqs in the data set. It does this for each of the seqs in the data. It will store the resulting values in memory as self.data['other']['divergence_non_syn'] and self.data['other']['divergence_syn']
        ref_seq - This gives the path to the fasta file that has the reference sequence that each of the sequences in the data will be compared to to get the divergence values. This fasta file should contain one sequence with no line breaks (two lines total in the whole file). It should have an attribute in the header of the sq that looks like this: '>blah_blah|amino_acid_seq=[X]|blah_blah_blah...', where [X] is a string of amino acid characters.
        output_dirpath - This gives the path to the directory that will contain the sequence data for all the samples, with the divergence info included with them.
        method - This gives the method for which divergence is calculated. Acceptable values are:
            'counting' - This means that non-syn. and syn. divergence is calculated by first codon aligning the query seq to the ref seq and then counting all the non-syn and syn changes between the two seqs, and then normalizing by all the possible non-syn and syn changes that could have taken place between the two seqs. A tutorial that describes how this was implemented can be found here: http://bioinformatics.cvr.ac.uk/blog/calculating-dnds-for-ngs-datasets/ .
        temp_dirpath - Gives the path to a directory for which temp data can be written. If None (default), then uses the current working directory.
        method='counting'. Acceptable values are:
            'equal_weights' - Default. This means that the different possible mutational paths between two codons are weigted equally. I.e. non-syn and syn mutations are treated as equally likely.
            'max_synonymous' - This means that the mutational path with the maximum amount of syn muts will be chosen for each reference/query codon pair.
        """
        if output_dirpath[-1] != '/':
            output_dirpath += '/'
        if not os.path.exists(output_dirpath):
            os.makedirs(output_dirpath)
        file_structure = master_directory(self.dirpath)
        file_structure.mirror_directory_tree_with_files(out_dirpath=output_dirpath, only_include_filetypes=['fasta', 'fa', 'fastq', 'fq'], include_file_suffix=False)
        output_filepaths = file_structure.mirrored_filepaths
        for file_out, file_in in itertools.izip(sorted(output_filepaths), sorted(self.sample_filepaths)):
            print file_out
            print file_in
            print '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
            sample = sequence_sample(filepath=file_in, count_attribute_name=self.count_attribute_name)
            sample.get_dN_dS_and_divergence(ref_seq=ref_seq, method=method, temp_dirpath=temp_dirpath, method_counting=method_counting)
            output_filepath = file_out + '.fasta'
            sample.write_full_data_to_disk_fasta(output_filepath=output_filepath)
        return

    def translate_and_get_dN_dS_divergence(self, output_dirpath, ref_seq, coding_frame_start=1, temp_dirpath=None, remove_seqs_with_stop_codons=False, method='counting', method_counting='equal_weights', path_to_needle='needle'):
        """
        This method runs 'translate_each_seq' and 'get_divergence_nonSynAndSyn' for each of the sequence samples in the dataset. See these two methods for descriptions of the parameters.
        """
        if output_dirpath[-1] != '/':
            output_dirpath += '/'
        if not os.path.exists(output_dirpath):
            os.makedirs(output_dirpath)
        file_structure = master_directory(self.dirpath)
        file_structure.mirror_directory_tree_with_files(out_dirpath=output_dirpath, only_include_filetypes=['fasta', 'fa', 'fastq', 'fq'], include_file_suffix=False)
        output_filepaths = file_structure.mirrored_filepaths
        for file_out, file_in in itertools.izip(sorted(output_filepaths), sorted(self.sample_filepaths)):
            print file_out
            print file_in
            print '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'

#            if file_in != '/Users/nstrauli/data/abr_hiv_coevo/seq_data/QCed_data_clean_final/hiv/2/4656.fasta':
#                continue

            sample = sequence_sample(filepath=file_in, count_attribute_name=self.count_attribute_name)
            sample.translate_each_seq(coding_frame_start=coding_frame_start, ref_seq=ref_seq, temp_dirpath=temp_dirpath, remove_seqs_with_stop_codons=remove_seqs_with_stop_codons, path_to_needle=path_to_needle)
            sample.get_dN_dS_and_divergence(ref_seq=ref_seq, method=method, temp_dirpath=temp_dirpath, method_counting=method_counting, path_to_needle=path_to_needle)
            output_filepath = file_out + '.fasta'
            sample.write_full_data_to_disk_fasta(output_filepath=output_filepath)
        return

    def identify_outlier_seqs(self, align_method='vsearch', path_to_vsearch=None, temp_dirpath=None, num_parallel_cores=1, mean_pairwise_distance_output_filepath=None, pairwise_distance_output_filepath=None, outlier_id_method='MAD', MAD_threshold=3.5, outlier_fasta_output_dirpath=None, remove_outliers_from_data=False):
        """
        This method will find outlier sequences (by genetic distance) in the sample set. It does this by doing an all-by-all pairwise alignment, so excersize caution when using with large datasets.
        align_method - This gives the method that will be used to align and calculate the edit distance between each pair of seqs. Acceptable values are:
            'vsearch' - This means that the vsearch program will be used to run the all-by-all alignment.
        path_to_vsearch - If method id 'vsearch' then this gives the path to the vsearch excecutable, otherwise this parameter is not considered. If this is None (default) then it is assumed that vsearch excecutables are in the PATH
        temp_dirpath - If defined (default, None), this gives the directory that will store the temporary data to run this method. If this is None, then this will be set to the current working directory.
        num_parallel_cores - This gives the number of cores to be used, if the aligner method allows for parrallel processing (ex: 'vsearch'). If not, this is ignored.
        mean_pairwise_distance_output_filepath - If defined (default, None), this gives the filepath to the output file that will contain a list of the mean pairwise distance for each sequence against all the other seqs in the dataset.
        pairwise_distance_output_filepath - If defined (default, None), this gives the filepath to the output file that will contain a long list of the pairwise genetic distances in the data. They will be sorted (decending) by distance.
        outlier_id_method - This gives the statistical method for identifying the outliers. Acceptable values are:
            'MAD' - This means that the median-absolute-deviation (MAD) approach will be used. A reference for this method is "Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and Handle Outliers", The ASQC Basic References in Quality Control: Statistical Techniques, Edward F. Mykytka, Ph.D., Editor." Can also refer to "https://stackoverflow.com/questions/22354094/pythonic-way-of-detecting-outliers-in-one-dimensional-observation-data".
        MAD_threshold - This gives the modified z-score threshold for the MAD test. Default if 3.5. This is only considered if 'outlier_id_method' is 'MAD'.
        outlier_fasta_output_dirpath - If defined (default, None), this gives the path to the directory for which the outlier sequence will be written in fasta format. The outliers for each of the sequence files will be written together and will have the same filename as the seq file from which they came.
        remove_outliers_from_data - If True (default, False) this will instruct the method to actually remove the outlier sequences from the data. BEWARE, this will permanently overwrite your sequence files.
        """
        if temp_dirpath == None:
            temp_dirpath = os.getcwd() + '/'
        else:
            if temp_dirpath[-1] != '/':
                temp_dirpath += '/'
            if not os.path.exists(temp_dirpath):
                os.makedirs(temp_dirpath)
        random_suffix = str(random.random())
        temp_fasta_file = temp_dirpath + 'temp_input_' + random_suffix + '.fasta'
        for i in self.sample_filepaths:
            filename = os.path.basename(i)[:-6]
            sample = sequence_sample(filepath=i, count_attribute_name=self.count_attribute_name)
            sample.add_string_to_each_id(string_to_add=filename+'_', add_to_start_or_end='start')
            sample.write_to_disk(output_filepath=temp_fasta_file, append_to_file=True)
        temp_output_filepath = temp_dirpath + 'temp_output_' + random_suffix
        if align_method == 'vsearch':
            p = Popen([path_to_vsearch, '--allpairs_global', temp_fasta_file, '--userout', temp_output_filepath, '--acceptall', '--userfields', 'query+target+id', '--threads', str(num_parallel_cores)], stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()
            subprocess.call(['rm', temp_fasta_file])
            filein = open(temp_output_filepath, "r")
            distances_dic = {}
            if pairwise_distance_output_filepath:
                fileout_distances = open(pairwise_distance_output_filepath, "w")
                fileout_distances.write('pairwise_distances\n')
            for i in filein:
                line = i.split('\t')
                header1 = line[0].split('|')
                header2 = line[1].split('|')
                count1 = int(header1[1].split('=')[1])
                count2 = int(header2[1].split('=')[1])
                seq_id1 = header1[0]
                seq_id2 = header2[0]
                distance = 100 - float(line[-1])
                num_comparisons = count1 * count2
                try:
                    distances_dic[seq_id1][0] += distance*num_comparisons
                    distances_dic[seq_id1][1] += num_comparisons
                except KeyError:
                    distances_dic[seq_id1] = [distance*num_comparisons, num_comparisons]
                try:
                    distances_dic[seq_id2][0] += distance*num_comparisons
                    distances_dic[seq_id2][1] += num_comparisons
                except KeyError:
                    distances_dic[seq_id2] = [distance*num_comparisons, num_comparisons]
                if pairwise_distance_output_filepath:
                    for j in [distance] * num_comparisons:
                        fileout_distances.write('%s\n' % j)
            filein.close()
            subprocess.call(['rm', temp_output_filepath])
            if pairwise_distance_output_filepath:
                fileout_distances.close()
        distances_ids_list = sorted([(distances_dic[i][0] / distances_dic[i][1], i) for i in distances_dic], reverse=True)
        distances_list, ids_list = zip(*distances_ids_list)
        if mean_pairwise_distance_output_filepath:
            fileout_mean_distances = open(mean_pairwise_distance_output_filepath, "w")
            fileout_mean_distances.write('sequence_ID\tmean_pairwise_genetic_distance\n')
            for ID, distance in zip(ids_list, distances_list):
                fileout_mean_distances.write('%s\t%s\n' % (ID, distance))
            fileout_mean_distances.close()
        #identify outliers
        if outlier_id_method == 'MAD':
            median = numpy.median(distances_list)
            diffs = [numpy.abs(i - median) for i in distances_list]
            med_abs_deviation = numpy.median(diffs)
            modified_z_score = [0.6745*i/med_abs_deviation for i in diffs]

            print sorted(modified_z_score)

            outlier_bool = numpy.array(modified_z_score) > MAD_threshold
            outliers = numpy.array(ids_list)[outlier_bool]
        #write outliers to disk, if desired
        if outlier_fasta_output_dirpath:
            if outlier_fasta_output_dirpath[-1] != '/':
                outlier_fasta_output_dirpath += '/'
            if not os.path.exists(outlier_fasta_output_dirpath):
                os.makedirs(outlier_fasta_output_dirpath)
            fileout_dic = {}
            for i in outliers:
                seq_data = i.split('_')
                filename = seq_data[0]
                seq_id = seq_data[1]
                seq_sample_filepath = '%s%s.%s' % (self.dirpath, filename, self.file_suffix)
                sample = sequence_sample(filepath=seq_sample_filepath, count_attribute_name=self.count_attribute_name)
                seq_data = sample.get_seq_data(seq_id=seq_id)
                output_filepath = outlier_fasta_output_dirpath + filename + '.fasta'
                try:
                    fileout_dic[filename].write('>%s|%s=%s|%s\n%s\n' % (seq_id, self.count_attribute_name, seq_data['count'], '|'.join(['%s=%s' % (j, seq_data['other'][j]) for j in seq_data['other']]), seq_data['seq']))
                except KeyError:
                    fileout_dic[filename] = open(output_filepath, "w")
                    fileout_dic[filename].write('>%s|%s=%s|%s\n%s\n' % (seq_id, self.count_attribute_name, seq_data['count'], '|'.join(['%s=%s' % (j, seq_data['other'][j]) for j in seq_data['other']]), seq_data['seq']))
            for i in fileout_dic:
                fileout_dic[i].close()
        if remove_outliers_from_data:
            #group outliers by their input file
            outlier_dic = {}
            for i in outliers:
                seq_data = i.split('_')
                filename = seq_data[0]
                seq_id = seq_data[1]
                seq_sample_filepath = '%s%s.%s' % (self.dirpath, filename, self.file_suffix)
                try:
                    outlier_dic[seq_sample_filepath].append(seq_id)
                except KeyError:
                    outlier_dic[seq_sample_filepath] = [seq_id]
            #now overwrite data
            for i in outlier_dic:
                sample = sequence_sample(filepath=i, count_attribute_name=self.count_attribute_name)
                sample.remove_seq_entries(seq_indicators=outlier_dic[i], seq_ids=True)
                #sample.write_full_data_to_disk_fasta(output_filepath=i)

                output_filepath = '/Users/nstrauli/data/abr_hiv_coevo/temp_stuff/' + os.path.basename(self.dirpath[:-1]) + '/' + os.path.basename(i)
                if not os.path.exists(os.path.dirname(output_filepath)):
                    os.makedirs(os.path.dirname(output_filepath))
                sample.write_full_data_to_disk_fasta(output_filepath=output_filepath)

        return outliers

    def cluster_seqs_within_samples(self, output_dirpath, method='vsearch', min_ident_within_cluster=0.97, path_to_vsearch=None, temp_dirpath=None, max_edit_distance=1, use_comp_cluster=False, wait_till_cluster_done=False, path_to_needle=None, scale_freqs_by=None, dists_dstrb_output_dirpath=None, output_network_dirpath=None, output_node_attribute_dirpath=None, node_attributes_to_add=['count'], names_for_node_attributes_in_header=None, add_freq_prior_to_clustering=None, add_indiv_seq_attribute=None):
        """
        This method uses the 'cluster_seqs' method in the sequence_sample class to cluster the sequences in each of the samples in the data set.
        output_dirpath - This gives the path to the directory that the clustered sequences will be written to.
        method - This gives the method that will be used for clustering. Acceptable values are:
            'vsearch' - This means that the program Vsearch will be used to cluster based on a given identity threshold. If this method is chosen then centroids of the clusters will be returned.
            'by_edit_distance' - This means that the sequences will be clustered by a simple requirement that all seqs within a cluster must have an edit distance with at least one other seq in the cluster that is less than or equal to a provided max edit distance cutoff. This requires an all by all pairwise alignment, so be carefule with large datasets. This also used the program 'needle' to do the alignments, which needs to be in the PATH.
            'by_edit_distance_seqanpy' - This means that edit distance is still used as the metric to cluster seqs by, but the package 'seqanpy' is used to do the alignments, instead of 'needle'. This may be faster.
        min_ident_within_cluster - This is an important parameter that gives the distance allowed between the sequences within a cluster. It should be between 0 and 1. If a target sequence, when aligned to a cluster centroid, has an identity lower then this value, then it is not added to that cluster. This parameter is only considered if method='vsearch'.
        path_to_vsearch - This gives the path to the vsearch excecutable. If this is None (default) then the path is assumed to be 'vsearch' (i.e. in $PATH). This parameter is only considered if method='vsearch'.                                     
        'temp_dirpath' - This gives the path to the directory for which a temporary output directory will be created in. It will contain temp fasta files and output from vsearch. The temp directory that is created within 'temp_dirpath' will be removed at the end of the script, but the 'temp_dirpath' directory itself will not be removed. If this equals None, the the current working directory is used.
        max_edit_distance - Int. Default, 1. This gives the maximum edit distance between any two seqs, for them to be assigned to the same cluster. This parameter is only considered if method='by_edit_distance'
        use_comp_cluster - If True (default, False), this will instruct the method to submit an array job to a computational cluster using SGE. One job will be submitted per sample filepath in the sample set. Currently, only the method 'by_edit_distance' or 'by_edit_distance_seqanpy' options are supported when using a computational cluster.
        wait_till_cluster_done - If True (default, False), this will cause the method to wait until the jobs on the computational cluster had completed before exiting. This can be useful if one is string a bunch of methods together, and stuff down the pipe depends on the results of this step.
        path_to_needle - This gives the path to the needle executable. This is required if using the computational cluster for an array job because the PATH variable is all messed up when sending jobs to parallele nodes. One can also use this if needle is not in the PATH
        scale_freqs_by - If defined (default, None), this should be a dic that gives information as how to scale the frequencies calculated for each of the clusters in each of the samples. Specifically, this should be a dictionary, where the index gives the absolute filepath of each file in the sample set, and the definition will be a float that will be multiplied by the frequency that is calculated for all the clusters in that sample. For example if a certain cluster (in a certain sample) has a calculated freq of .5, and one is to scale this freq by .25, then the resulting freq will be .125. Note that if this is used, freqs will not sum to one!
        dists_dstrb_output_dirpath - This gives the directory path for which all the the edit distances will be written for each pairwise sequence comparison of seqs within a sample. One file per sample. This parameter is only considered if method=='by_edit_distance'.
        output_network_dirpath - If defined (default, None), this will, 1) instruct the method to creat a network file (as a simple interaction file ('.sif')), and 2) will be the path to this file. This is only considered if method='by_edit_distance'
        output_node_attribute_dirpath - This gives the path to the file that will contain the attributes for each of the unique nodes (i.e. unique seqs) in the data. This is only considered if 'output_network_dirpath' is defined.
        node_attributes_to_add - This is a list of sequence attributes to include in the node attribute file. This parameter is only considered if 'output_node_attribute_dirpath' is defined. The default is to only add the value for the 'count' attribute (i.e. the value for self.data[some_index]['count']). Acceptable values within this list of strings are:
            any string that is a key for one of the attributes for the sequences. So, could be 'count', 'total_freq', 'timepoint', etc. if those are keys to the dictionarys within self.data[some_index], or self.data[some_index]['other']
            'element_X_of_seq_id' - sometimes info about the seq is encoded in the seq ID, delimited by '_'. This means that the 'X'th element of the seq ID (when delimited by a '_') will be included in the node attributes file.
        names_for_node_attributes_in_header - This is an optional argument, which if defined (default, None), it will give the names for each of the node attributes included in the 'output_node_attribute_dirpath' header. This parameter is only considered if 'output_node_attribute_dirpath'. If 'output_node_attribute_dirpath' is defined and this is not, then the default is to simply use the values in 'node_attributes_to_add'. The names in this list need to be in the same order as the attributes in 'node_attributes_to_add'.
        add_freq_prior_to_clustering - If defined (default, None), this should be a string that gives the name of the frequency attribute that will be added to each of the seqs in the dataset. If freq information has not been calculated yet for the seq data, and one wants to include freq info in the network plots, then one can use this parameter to add freq info to each of the sample's data before clustering.
        add_indiv_seq_attribute - If defined (default, None), this should be a list of strings where each string is the same as an attribute for each of the individual seqs that should be included in the header for each of the clusters. Similar to the 'indiv_seq_freqs' attribute that is included for each of the seqs belonging to a cluster, this will give the provided attribute's values for each of the seqs in a cluster. This is only applied if 'overwrite_data' is set to True
        """
        if output_dirpath[-1] != '/':
            output_dirpath += '/'
        if not os.path.exists(output_dirpath):
            os.makedirs(output_dirpath)
        if temp_dirpath == None:
            temp_dirpath = os.getcwd() + '/'
        else:
            if temp_dirpath[-1] != '/':
                temp_dirpath += '/'
            if not os.path.exists(temp_dirpath):
                os.makedirs(temp_dirpath)
        if dists_dstrb_output_dirpath:
            if dists_dstrb_output_dirpath[-1] != '/':
                dists_dstrb_output_dirpath += '/'
            if not os.path.exists(dists_dstrb_output_dirpath):
                os.makedirs(dists_dstrb_output_dirpath)
        file_structure = master_directory(self.dirpath)
        file_structure.mirror_directory_tree_with_files(out_dirpath=output_dirpath, only_include_filetypes=['fasta', 'fa', 'fastq', 'fq'], include_file_suffix=False, avoid_files_with=self.avoid_files_with)
        output_filepaths = file_structure.mirrored_filepaths
        if dists_dstrb_output_dirpath:
            file_structure = master_directory(self.dirpath)
            file_structure.mirror_directory_tree_with_files(out_dirpath=dists_dstrb_output_dirpath, only_include_filetypes=['fasta', 'fa', 'fastq', 'fq'], include_file_suffix=False, avoid_files_with=self.avoid_files_with)
            output_filepaths_dists = sorted(file_structure.mirrored_filepaths)
        if output_network_dirpath:
            file_structure = master_directory(self.dirpath)
            file_structure.mirror_directory_tree_with_files(out_dirpath=output_network_dirpath, only_include_filetypes=['fasta', 'fa', 'fastq', 'fq'], include_file_suffix=False, avoid_files_with=self.avoid_files_with)
            output_filepaths_networks = sorted(file_structure.mirrored_filepaths)
        if output_node_attribute_dirpath:
            file_structure = master_directory(self.dirpath)
            file_structure.mirror_directory_tree_with_files(out_dirpath=output_node_attribute_dirpath, only_include_filetypes=['fasta', 'fa', 'fastq', 'fq'], include_file_suffix=False, avoid_files_with=self.avoid_files_with)
            output_filepaths_nodes = sorted(file_structure.mirrored_filepaths)

        if use_comp_cluster:
            random_suffix = str(random.random())
            temp_input_filepath_info_filepath = '%stemp_input_filepath_info_%s' % (temp_dirpath, random_suffix)
            fileout = open(temp_input_filepath_info_filepath, 'w')
            index = 0
            for file_out, file_in in itertools.izip(sorted(output_filepaths), sorted(self.sample_filepaths)):
                #check to make sure in and out files are right
                file_in_base = os.path.basename(file_in)
                file_out_base = os.path.basename(file_out)
                if file_in_base != '%s.%s' % (file_out_base, self.file_suffix):
                    print 'Input and output files are out of sync. Aborting.'
                    return
                fileout.write('%s\t%s' % (file_in, file_out))
                if dists_dstrb_output_dirpath:
                    fileout.write('\t%s.txt' % output_filepaths_dists[index])
                else:
                    fileout.write('\tNone')
                if output_network_dirpath:
                    fileout.write('\t%s.sif' % output_filepaths_networks[index])
                else:
                    fileout.write('\tNone')
                if output_node_attribute_dirpath:
                    fileout.write('\t%s_node_attributes.txt' % output_filepaths_nodes[index])
                else:
                    fileout.write('\tNone')
                fileout.write('\n')
                index += 1
            fileout.close()
            if scale_freqs_by:
                temp_scale_freqs_by_filepath = '%stemp_scale_freqs_by_%s' % (temp_dirpath, random_suffix)
                filein = open(temp_scale_freqs_by_filepath, "w")
                for i in scale_freqs_by:
                    filein.write('%s\t%s\n' % (i, scale_freqs_by[i]))
                filein.close()
            else:
                temp_scale_freqs_by_filepath = 'None'
            if names_for_node_attributes_in_header:
                names_for_node_attributes_in_header = ','.join(names_for_node_attributes_in_header)
            if add_indiv_seq_attribute:
                add_indiv_seq_attribute = ','.join(add_indiv_seq_attribute)
            job_name = 'clust_within_%s' % random_suffix
            p = Popen(['qsub', '-N', job_name, '-t', '1-%s' % len(self.sample_filepaths), 'sequence_sample_set_class.py', 'cluster_seqs_within_samples_compCluster', temp_input_filepath_info_filepath, self.count_attribute_name, '/scratch/', str(max_edit_distance), path_to_needle, temp_scale_freqs_by_filepath, str(dists_dstrb_output_dirpath), str(output_network_dirpath), str(output_node_attribute_dirpath), ','.join(node_attributes_to_add), str(names_for_node_attributes_in_header), str(add_freq_prior_to_clustering), str(add_indiv_seq_attribute), method], stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()
            print out
            print err
            if wait_till_cluster_done:
                submit_job_id = set([out.split()[2].split('.')[0]])
                #sleep till job is done
                good_to_go = False
                while good_to_go == False:
                    time.sleep(10)
                    p = Popen(['qstat'], stdout=PIPE, stderr=PIPE)
                    jobs = p.communicate()[0].split('\n')
                    job_ids = set()
                    if len(jobs) >= 1:
                        for i in jobs[2:-1]:
                            job_id = i.split()[0]
                            job_ids.update([job_id])
                        if len(job_ids.intersection(submit_job_id)) == 0:
                            good_to_go = True
                    else:
                        good_to_go = True
                for i in sorted(output_filepaths):
                    if not os.path.exists(i + '.fasta'):
                        print ''
                        print 'Output filepath does not exists and it should. Something went wrong babe:', i
                        print 'Offending job ID:', submit_job_id
                subprocess.call(['rm', temp_input_filepath_info_filepath])
                if scale_freqs_by:
                    subprocess.call(['rm', temp_scale_freqs_by_filepath])
            return job_name

        else:
            for file_out, file_in in itertools.izip(sorted(output_filepaths), sorted(self.sample_filepaths)):
                print file_in
                print file_out
                if scale_freqs_by != 'None':
                    scale_freqs_by = scale_freqs_by[file_in]
                sample = sequence_sample(filepath=file_in, count_attribute_name=self.count_attribute_name)
                sample.cluster_seqs(method=method, min_ident_within_cluster=min_ident_within_cluster, path_to_vsearch=path_to_vsearch, temp_dirpath=temp_dirpath, max_edit_distance=max_edit_distance, overwrite_data=True, path_to_needle=path_to_needle)
                sample.write_full_data_to_disk_fasta(output_filepath=file_out + '.fasta')
            return

    def cluster_seqs_across_samples(self, output_network_filepath=None, output_node_attribute_filepath=None, node_attributes_to_add=['count'], names_for_node_attributes_in_header=None, temp_dirpath=None, add_string_to_seq_ids=None, add_to_start_or_end='start', max_edit_distance=1, path_to_needle=None, id_outliers=False, outlier_def='clusters_with_dif_file', outlier_freq_cutoff=0.001, freq_attribute_name=None):
        """
        This method will pool all the seqs in the sample set, and cluster them. BEWARE, this requires all the data to be read into memory, and also requires that all the data be clustered. So, if the data is anywhere near large, then don't use this. It is recommended that one runs the method 'cluster_seqs_within_samples' and then run this method on the output of that.
        output_network_filepath - If defined (default, None), this will, 1) instruct the method to creat a network file (as a simple interaction file ('.sif')), and 2) will be the path to this file.
        output_node_attribute_filepath - This is the file that will contain all the node attributes for each of the seqs. This is for cytoscape. This parameter is only considered if 'output_network_filepath' is defined.
        node_attributes_to_add - This is a list of sequence attributes to include in the node attribute file. This parameter is only considered if 'output_node_attribute_filepath' is defined. The default is to only add the value for the 'count' attribute (i.e. the value for self.data[some_index]['count']). Acceptable values within this list of strings are:
            any string that is a key for one of the attributes for the sequences. So, could be 'count', 'total_freq', 'timepoint', etc. if those are keys to the dictionarys within self.data[some_index], or self.data[some_index]['other']
            'element_X_of_seq_id' - sometimes info about the seq is encoded in the seq ID, delimited by '_'. This means that the 'X'th element of the seq ID (when delimited by a '_') will be included in the node attributes file.
        names_for_node_attributes_in_header - This is an optional argument, which if defined (default, None), it will give the names for each of the node attributes included in the 'output_node_attribute_filepath' header. This parameter is only considered if 'output_node_attribute_filepath'. If 'output_node_attribute_filepath' is defined and this is not, then the default is to simply use the values in 'node_attributes_to_add'. The names in this list need to be in the same order as the attributes in 'node_attributes_to_add'.
        add_string_to_seq_ids - This tells whether or not to add strings to each seq ID for the purposes of making the seq IDs unique across the entire data set. Acceptable values are:
            None - Default. This means that nothing will be added to the seq IDs.
            'filename' - This means that the base filename (not including suffix (ex: '.fasta')) will be added to each seq ID, in that sample file.
            'dirname_and_filename' - This means that both the name of the containing directory and the filename (not including suffix) will be added to each seq ID for a given sample file. Ex: if the file is located in '~/stuff/seqs/blah/patient_1/1229.fasta' then 'patient_1_1229_' will be inserted at the beggining (or end) of each seq ID.
        add_to_start_or_end - This tells if the provided string will be added to the beginning or the end of each of the IDs. Only considered if 'add_string_to_seq_ids' is defined. Acceptable values are:
            'start' - Default. String will be added to the beggining.
            'end' - String will be added to the end.
        max_edit_distance - Int. Default, 1. This gives the maximum edit distance between any two seqs, for them to be assigned to the same cluster.
        path_to_needle - This gives the path to the needle executable. This is required if using the computational cluster for an array job because the PATH variable is all messed up when sending jobs to parallele nodes. One can also use this if needle is not in the PATH
        id_outliers - Boolean (default, False). Sometimes one does not want the seqs to cluster across samples, and instead wants to find those that do, and get rid of them. If this is defined it should be a string that gives the path to a fasta file that will contain any found outlier sequences. This will instruct the method to find each of the seqs in the data set that are outliers. An outlier is defined as a low frequency sequence that clusters with other seqs that it's not supposed to. Details on this below. If this is True, it will write the data to disk with this new boolean outlier attribute in the header of each seq.
        outlier_def - This gives the details about what defines an outlier seq. I.e. what determines if a seq is clustering with other seqs that it's not supposed to. Acceptable values are:
            'clusters_with_dif_file' - This means that if the seq is low frequency and clusters with others seqs from a different filepath, it is deamed an outlier.
            'clusters_with_dif_directory' - This means that if the seq is low frequency and clusters with others seqs from a different directory, it is deamed an outlier.
        outlier_freq_cutoff - One of the conditions of being an outlier is that it has to be of low frequency. This gives what the cutoff for being classified as low freq actually is. Default is 0.001.
        freq_attribute_name - If defined (default, None), this will give the name of the attribute that gives the frequency for each seq in the data set. If it is not defined then frequency is calculated by the seq 'count' / the total, and the freq_attribute_name is set to 'freq'.
        """
        if temp_dirpath == None:
            temp_dirpath = os.getcwd() + '/'
        else:
            if temp_dirpath[-1] != '/':
                temp_dirpath += '/'
            if not os.path.exists(temp_dirpath):
                os.makedirs(temp_dirpath)
        random_suffix = str(random.random())
        temp_fasta_filepath = '%stemp_seqfile_%s.fasta' % (temp_dirpath, random_suffix)
        if id_outliers:
            seq_filepath_dic = {}
        for i in self.sample_filepaths:
            sample = sequence_sample(filepath=i, count_attribute_name=self.count_attribute_name)
            if not freq_attribute_name:
                sample.add_freq_attribute(freq_attribute_name='freq')
            if add_string_to_seq_ids:
                if add_string_to_seq_ids == 'filename':
                    string_to_add = '.'.join(os.path.basename(i).split('.')[:-1]) + '_'
                elif add_string_to_seq_ids == 'dirname_and_filename':
                    directory_name = i.split('/')[-2]
                    file_name = '.'.join(os.path.basename(i).split('.')[:-1])
                    string_to_add = '%s_%s_' % (directory_name, file_name)
                sample.add_string_to_each_id(string_to_add, add_to_start_or_end=add_to_start_or_end)
            if id_outliers:
                for j in sample.data:
                    seq_filepath_dic[j['id']] = i
            sample.write_full_data_to_disk_fasta(output_filepath=temp_fasta_filepath, append_to_file=True)
        if not freq_attribute_name:
            freq_attribute_name = 'freq'
        sample = sequence_sample(filepath=temp_fasta_filepath, count_attribute_name=self.count_attribute_name)
        print 'Total unique sequences in pooled data:', len(sample.data)
        cluster_dic, connections_dic = sample.cluster_seqs(method='by_edit_distance', temp_dirpath=temp_dirpath, max_edit_distance=max_edit_distance, overwrite_data=False, output_network_filepath=output_network_filepath, output_node_attribute_filepath=output_node_attribute_filepath, node_attributes_to_add=node_attributes_to_add, names_for_node_attributes_in_header=names_for_node_attributes_in_header, path_to_needle=path_to_needle)
        if id_outliers:
            outlier_indices = []
            for index, i in enumerate(sample.data):
                is_outlier = False
                #seq must be below the outlier cutoff freq to be deemed an outlier
                if len(connections_dic[i['id']]) > 0 and float(i['other'][freq_attribute_name]) <= outlier_freq_cutoff:
                    if outlier_def == 'clusters_with_dif_file':
                        seq1_group = seq_filepath_dic[i['id']]
                    elif outlier_def == 'clusters_with_dif_directory':
                        seq1_group = os.path.dirname(seq_filepath_dic[i['id']])
                    for j in connections_dic[i['id']]:
                        if outlier_def == 'clusters_with_dif_file':
                            seq2_group = seq_filepath_dic[sample.data[j]['id']]
                        elif outlier_def == 'clusters_with_dif_directory':
                            seq2_group = os.path.dirname(seq_filepath_dic[sample.data[j]['id']])
                        #if the seq is connected to another seq of the group that is large (i.e. high freq), then it's not an outlier
                        if seq1_group == seq2_group and sample.data[j]['other'][freq_attribute_name] > outlier_freq_cutoff:
                            is_outlier = False
                            break
                        #else if the seq is conneted to another seq from a different group, then it is an outlier
                        elif seq1_group != seq2_group:
                            is_outlier = True
                if is_outlier:
                    outlier_indices.append(index)
            sample.write_subset_of_seqs_to_disk(seq_indices=outlier_indices, output_filepath=id_outliers, append_to_file=False)
        subprocess.call(['rm', temp_fasta_filepath])
        return


#############################################
# Below are scripts that are used for running
# jobs on the computational cluster. They
# should only be used by the methods above.
#############################################

def vsearch_pi_calculator_compCluster(min_sample_size, count_attribute_name, path_to_vsearch, temp_dirpath, num_parallel_cores, sge_override):
    """
    This script is to be used by 'down_sample_and_calc_pi_compCluster' when the 'method' parameter equals 'vsearch'. This script is what is called when the 'down_sample_and_calc_pi_compCluster' method submits an array job. So, each of the jobs run this script. Each job will downsample then calculate pi by way of using the vsearch software.
    min_sample_size - this is the value that each sequence sample will be downsampled to. If this equals 'Nope' then no downsampling occurs.
    count_attribute_name - this is the name of the count attribute in each of the headers in each of the sequences in the sequence file.
    path_to_vsearch - This is the path to the vsearch excecutable.
    temp_dirpath - This is the path to the temporary directory that will contain all the sequence sample filepaths. This directory will also contain the temporary output of vsearch, which can be quite large. So make sure there is enough space where ever this path leads.
    num_parallel_cores - This gives the number of parallel jobs that vsearch can use. vsearch's default is to use as many cores as it can, and this slows down the cluster for others. So, we need to set this to a fixed amount, and then let the cluster know when submitting array job. Default is 12.
    try_again - Boolean. If True, then the script will try to calculate pi once more if the output of vsearch is faulty. It seems that vsearch is a bit buggy and can write uninterpretable lines of output. This setting (when set to True) will simply try once more.
    """
    if sge_override != 'None':
        sge_task_id = int(sge_override)
    else:
        sge_task_id = int(os.environ['SGE_TASK_ID'])
    filein = open('%s%s_input_filepath' % (temp_dirpath, sge_task_id), "r")
    sample_filepath = filein.readline()
    filein.close()
    #read in sequence data
    sample = sequence_sample(sample_filepath, count_attribute_name)
    #downsample
    if min_sample_size != 'Nope':
        min_sample_size = float(min_sample_size)
        sample.down_sample(min_sample_size)
    #make temp dir for vsearch output on /scratch
    vsearch_temp_dir = tempfile.mkdtemp(dir='/scratch')
    #run vsearch
    pi = sample.calc_diversity_pi(method='vsearch', path_to_vsearch=path_to_vsearch, temp_dirpath=vsearch_temp_dir, num_parallel_cores=num_parallel_cores)
    subprocess.call(['rm', '-r', vsearch_temp_dir])
    fileout = open('%s%s_pi' % (temp_dirpath, sge_task_id), "w")
    fileout.write('%s\t%s' % (sample_filepath, pi))
    fileout.close()
    return

def get_immune_reads_with_changeo_compCluster(temp_input_filepath, count_attribute_name, path_to_changeo, path_to_igblast, temp_dirpath, imgt_ref_seqs, add_germline, add_selection, use_comp_cluster, add_divergence, overwrite_output, translate_VDJ, remove_seqs_with_stop):
    """
    This script is to be used by 'get_immune_reads_with_changeo_foreach_sample', when the 'use_comp_cluster' is set to True. See the 'get_immune_reads_with_changeo_foreach_sample' method above, or the 'get_immune_reads_with_changeo' method in the sequence_sample class for a descripttion of all the parameters here.
    """
    sge_task_id = int(os.environ['SGE_TASK_ID'])

    # if sge_task_id != 1:
    #     return

    #get input and output filepaths from the temp file
    filein = open(temp_input_filepath, "r")
    count = 0
    for i in filein:
        count += 1
        if count == sge_task_id:
            line = i[:-1].split('\t')
            input_filepath = line[0]
            output_filepath = line[1]
            break
    filein.close()
    if count_attribute_name == 'None':
        count_attribute_name = None
    if temp_dirpath == 'None':
        temp_dirpath = None
    if imgt_ref_seqs == 'None':
        imgt_ref_seqs = None
    if add_germline == 'None':
        add_germline = None
    else:
        add_germline = add_germline.split(',')
    if add_selection == 'False':
        add_selection = False
    elif add_selection == 'True':
        add_selection = True
    if use_comp_cluster == 'False':
        use_comp_cluster = False
    elif use_comp_cluster == 'True':
        use_comp_cluster = True
    if add_divergence == 'False':
        add_divergence = False
    if overwrite_output == 'True':
        overwrite_output = True
    elif overwrite_output == 'False':
        overwrite_output = False
    if translate_VDJ == 'True':
        translate_VDJ = True
    elif translate_VDJ == 'False':
        translate_VDJ = False
    if remove_seqs_with_stop == 'True':
        remove_seqs_with_stop = True
    elif remove_seqs_with_stop == 'False':
        remove_seqs_with_stop = False
    print input_filepath
    print output_filepath
    sample = sequence_sample(input_filepath, count_attribute_name=count_attribute_name)
    sample.get_immune_reads_with_changeo(output_filepath=output_filepath, path_to_changeo=path_to_changeo, path_to_igblast=path_to_igblast, temp_dirpath=temp_dirpath, imgt_ref_seqs=imgt_ref_seqs, add_germline=add_germline, add_selection=add_selection, use_comp_cluster=use_comp_cluster, add_divergence=add_divergence, overwrite_output=overwrite_output, translate_VDJ=translate_VDJ)
    return

def cluster_seqs_within_samples_compCluster(filepath_info, count_attribute_name, temp_dirpath, max_edit_distance, path_to_needle, temp_scale_freqs_by_filepath, dists_dstrb_output_dirpath, output_network_dirpath, output_node_attribute_dirpath, node_attributes_to_add, names_for_node_attributes_in_header, add_freq_prior_to_clustering, add_indiv_seq_attribute, method):
    """
    This script manages the array job that is submitted by the method 'cluster_seqs_within_samples' when the parameter 'use_comp_cluster' is set to True. It finds the correct input and output filepaths for a given job, and then uses the 'sequence_sample_class' to run 'cluster_seqs' to cluster those samples.
    """
    sge_task_id = int(os.environ['SGE_TASK_ID'])
    if node_attributes_to_add != 'None':
        node_attributes_to_add = node_attributes_to_add.split(',')
    else:
        node_attributes_to_add = None
    if names_for_node_attributes_in_header != 'None':
        names_for_node_attributes_in_header = names_for_node_attributes_in_header.split(',')
    else:
        names_for_node_attributes_in_header = None
    if add_freq_prior_to_clustering == 'None':
        add_freq_prior_to_clustering = None
    if add_indiv_seq_attribute == 'None':
        add_indiv_seq_attribute = None
    else:
        add_indiv_seq_attribute = add_indiv_seq_attribute.split(',')
    #try to open the file 3 times. The server can be quite slow, and can cause an IO error sometimes, so we try this.
    try:
        filein = open(filepath_info, "r")
    except IOError:
        time.sleep(10)
        try:
           filein = open(filepath_info, "r")
        except IOError:
            time.sleep(10)
            filein = open(filepath_info, "r")
    count = 1
    for i in filein:
        if count == sge_task_id:
            line = i[:-1].split('\t')
            file_in = line[0]
            file_out = line[1]
            if line[2] != 'None':
                file_out_dists = line[2]
            else:
                file_out_dists = None
            if line[3] != 'None':
                file_out_network = line[3]
            else:
                file_out_network = None
            if line[4] != 'None':
                file_out_node = line[4]
            else:
                file_out_node = None
            break
        count += 1
    filein.close()

    # if file_in != '/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/abr_seq_unique_VJ/changeo/2/IGHV1-2_IGHJ1/1002.fasta':
    #     return

    print file_in
    print file_out
    print file_out_dists
    print file_out_network
    print file_out_node
    print add_indiv_seq_attribute
    print method
    sys.stdout.flush()
    #convert scale_freqs_by to dic, if necessary
    if temp_scale_freqs_by_filepath != 'None':
        try:
            filein = open(temp_scale_freqs_by_filepath, "r")
        except IOError:
            time.sleep(10)
            try:
               filein = open(temp_scale_freqs_by_filepath, "r")
            except IOError:
                time.sleep(10)
                filein = open(temp_scale_freqs_by_filepath, "r")
        for i in filein:
            line = i[:-1].split('\t')
            input_filepath = line[0]
            if input_filepath == file_in:
                scale_freqs_by = float(line[1])
                break
    else:
        scale_freqs_by = None
    max_edit_distance = int(max_edit_distance)
    sample = sequence_sample(filepath=file_in, count_attribute_name=count_attribute_name)
    if add_freq_prior_to_clustering:
        sample.add_freq_attribute(freq_attribute_name=add_freq_prior_to_clustering)
    sample.cluster_seqs(method=method, temp_dirpath=temp_dirpath, max_edit_distance=max_edit_distance, overwrite_data=True, path_to_needle=path_to_needle, scale_freqs_by=scale_freqs_by, dists_dstrb_output_filepath=file_out_dists, output_network_filepath=file_out_network, output_node_attribute_filepath=file_out_node, node_attributes_to_add=node_attributes_to_add, names_for_node_attributes_in_header=names_for_node_attributes_in_header, add_indiv_seq_attribute=add_indiv_seq_attribute)
    sample.write_full_data_to_disk_fasta(output_filepath=file_out + '.fasta')
    return

if __name__ == '__main__':
    if sys.argv[1] == 'vsearch_pi_calculator_compCluster':
        vsearch_pi_calculator_compCluster(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
    elif sys.argv[1] == 'get_immune_reads_with_changeo_compCluster':
        get_immune_reads_with_changeo_compCluster(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13], sys.argv[14])
    elif sys.argv[1] == 'cluster_seqs_within_samples_compCluster':
        cluster_seqs_within_samples_compCluster(filepath_info=sys.argv[2], count_attribute_name=sys.argv[3], temp_dirpath=sys.argv[4], max_edit_distance=sys.argv[5], path_to_needle=sys.argv[6], temp_scale_freqs_by_filepath=sys.argv[7], dists_dstrb_output_dirpath=sys.argv[8], output_network_dirpath=sys.argv[9], output_node_attribute_dirpath=sys.argv[10], node_attributes_to_add=sys.argv[11], names_for_node_attributes_in_header=sys.argv[12], add_freq_prior_to_clustering=sys.argv[13], add_indiv_seq_attribute=sys.argv[14], method=sys.argv[15])
    elif sys.argv[1] == 'test':
        #here one can test out different methods in the class
        sample_set = sequence_sample_set(dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_fasta_files/10', count_attribute_name='DUPCOUNT')
        sample_set.cluster_seqs_within_samples(output_dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/sequence_clusters/hiv/clustered_by_edit_distance/max_edit_dist_1', method='by_edit_distance', min_ident_within_cluster=0.97, path_to_vsearch=None, temp_dirpath='/Users/nstrauli/data/abr_hiv_coevo/temp_stuff', max_edit_distance=1)
    else:
        print 'If running sequence_sample_set_class.py as a script, then must specify an appropriate sub-script to run.'
