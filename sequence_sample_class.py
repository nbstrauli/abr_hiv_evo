#!/usr/bin/python
#$ -S /usr/bin/python
#$ -e error
#$ -o out_2
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=16G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=336:00:00

import sys
sys.path.insert(0, './')
import random
from subprocess import Popen, PIPE
import subprocess
from scipy.misc import comb
import re
import os
import time
import itertools
import tempfile
import math
from Bio import Seq
geneticCode = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L",
               "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
               "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
               "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
                
               "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
               "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
               "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
               "GCT":"A","GCC":"A","GCA":"A","GCG":"A",

               "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
               "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
               "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
               "GAT":"D","GAC":"D","GAA":"E","GAG":"E",

               "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
               "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
               "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
               "GGT":"G","GGC":"G","GGA":"G","GGG":"G", '---':'-'}
import get_nonSyn_and_syn_sites_foreachCodon

class sequence_sample(object):
    """
    A class for sequence samples, where the 'sample' is a set of genetic sequences that are the result of targeted deep sequencing of some species. This species could be antibodies, viruses, etc. Input is a fasta file containing all the sequences from the sample (ideally, post QC).

    Attributes:
    filepath - path to the fasta or fastq formatted file. If fastq, however, the 3rd and 4th lines are ignored for each sequence entry.
    data - a list of dics, where each element has all the information for a given read (except the quality scores). Each read's info is the unique seq ID, the count of the sequence in the data, and the genetic sequence itself. It also contains any other information that is stored in the header. It stores this information as a dictionary called 'other' that is itself an element of the sequence's entry. This dictionary named 'other' has keys that are named after whatever other names are in the header line, and the keys values are whatever the value is listed as for that key. For example if there is "...|mut_freq=.012|..." in the header, than this information would be stored as 'other'={'mut_freq':'.012'}. Make sure that all headers in the sample have the same set of attributes listed in them.
    counts - a list of counts for each of the sequences.
    total - total number of seqs represented in the file (not just unique seqs).
    diversity - diversity calculated as the mean pairwise genetic distance in the seqs in the fastq file. Must run the calc_diversity method before this is has a value (other than None type)
    downsamp_indicator - Boolean. If True, then it means the data has been down sampled. If False, then no downsampling has occured.
    count_attribute_name - This is the name (as a string) of the attribute in the header that gives the count information for all the reads.
    filetype - This is either 'fasta' or 'fastq'. It is determined by the suffix of the input file name.
    """

    def __init__(self, filepath, count_attribute_name=None):
        self.filepath = filepath
        suffix = os.path.basename(self.filepath).split('.')[-1]
        if suffix == 'fasta' or suffix == 'fa':
            self.filetype = 'fasta'
        elif suffix == 'fastq' or suffix == 'fq':
            self.filetype = 'fastq'
        else:
            print 'Unrecognized file type. Must be either fasta of fastq. Check file suffix.'
            return
        self.count_attribute_name = count_attribute_name
        self.data = []
        self.counts = []
        self.pi = None
        self.downsamp_indicator = False
        filein = open(self.filepath, "r")
        count = 0
        missing_counts = False
        for i in filein:
            count += 1
            #get info in the header
            if count == 1:
                seq_entry_dic = {'id':'NA', 'count':1, 'seq':'NA', 'other':{}}
                line = i[1:-1].split('|')
                #the id is assumed to be the first entry in the header
                id = line[0]
                seq_entry_dic['id'] = id
                #get seq attribute info from header
                found_count = False
                for j in line[1:]:
                    attribute = j.split('=')
                    name = attribute[0]
                    value = attribute[1]
                    if name == self.count_attribute_name:
                        seq_entry_dic['count'] = int(value)
                        self.counts.append(int(value))
                        found_count = True
                    else:
                        seq_entry_dic['other'][name] = value
                if found_count == False:
                    self.counts.append(1)
                    if self.count_attribute_name:
                        missing_counts = True
            elif count == 2:
                seq = i[:-1]
                seq_entry_dic['seq'] = seq
                self.data.append(seq_entry_dic)
                if self.filetype == 'fasta':
                    count = 0
            elif count == 4:
                count = 0
        filein.close()
        if missing_counts:
            print "At least one sequence entry exists in the sample that does not have a count attribute, even though a count attribute was specified (i.e. not a None type). You might want to check this."
        self.total = float(sum(self.counts))
        return

    def num_entries(self):
        """Counts the number of reads in a fastq file and returns that number. Terribly simple"""
        return len(self.data)

    @staticmethod
    def convert_pmf2cmf(pmf):
        """This script takes the pdf and simply converts it to a cdf. Pretty straight forward."""
        total = 0
        cmf = []
        for i in pmf:
            cmf.append(total + i)
            total += i
        return cmf

    @staticmethod
    def inverse_transform_sampling(cmf):
        """This is a cool general way to sample from any distribution (especially empirical ones). Generates a random # b/t 0 and 1, then finds the largest element in the cdf (x-axis) that has a value (y-axis) that is greater or equal to the random number. Returns the index of the element that it finds ('chosen_cdr3_index')."""
        random_float = random.random()
        for index, i in enumerate(cmf):
            if i >= random_float:
                chosen_index = index
                break
        return chosen_index

    @staticmethod
    def update_dstrb(cmf, pmf, counts, total, chosen_index):
        """This updates the pmf and cmf according to the new counts that are contained in the 'counts' variable (which is just a list of counts for each unique sequence). The counts variable is only chaging by subtracting one from one of its elements, so there is probably a better way to do this, but this is what we have for now."""
        total -= 1
        total_prob = 0
        for i in xrange(len(cmf)):
            if i == chosen_index:
                counts[i] -= 1
            pmf[i] = counts[i] / total
            cmf[i] = pmf[i] + total_prob
            total_prob += pmf[i]
        return cmf, pmf, counts, total

    def down_sample(self, num_to_downsample_to, method='list_sampling'):
        """Down sample the sequences in the data. The idea is that each sequence represents a unique sequence in the data, and each unique sequence has counts that give the number of times it occurs. Often times one wants to downsample these sequences to a certain depth. For eaxample this needs to be done when calculating a diversity statistic. This method does the downsampling and permenantly saves the result to the self.data object.                                     
        'num_to_downsample_to' = The number of reads that the sequence data should be down sampled to.
        'method' = Can either be 'list_sampling' (default) or 'inverse_transform_sampling'. List sampling turns self.counts into a long list of indeces, and then randomly removes elements from this list, whereas 'inverse_transform_sampling' transforms the counts into a pmf and samples from this. 'list_sampling' is most likely fastest for most situations, but it might be advisable to experiment with both to see which is better for a specific case.
        """
        num_reads_to_remove = int(self.total - num_to_downsample_to)
        #perform inverse transform sampling if that is desired
        if method == 'inverse_transform_sampling':
            #now turn the counts into sampling probabilities
            pmf = [i / self.total for i in self.counts]
            cmf = self.convert_pmf2cmf(pmf)
            #now remove reads one-by-one
            for i in xrange(num_reads_to_remove):
                chosen_index = self.inverse_transform_sampling(cmf)
                cmf, pmf, self.counts, self.total = self.update_dstrb(cmf, pmf, self.counts, self.total, chosen_index)
        #else if the desired sampling method is list sampling then do that
        elif method == 'list_sampling':
            total_unique = len(self.counts)
            list_of_indeces = []
            for i, j in enumerate(self.counts):
                for k in xrange(j):
                    list_of_indeces.append(i)
            #now downsample the list of indeces
            for i in xrange(num_reads_to_remove):
                del list_of_indeces[random.randrange(0, self.total)]
                self.total -= 1
            #initialize the new list of counts
            self.counts = [0 for i in xrange(total_unique)]
            #fill in the new list of counts based upon the
            #down-sampled list_of_indices
            for i in list_of_indeces:
                self.counts[i] += 1

        #update the counts that are associated with each read in self.data      
        indeces_to_remove = []
        for i in xrange(len(self.counts)):
            self.data[i]['count'] = self.counts[i]
            if self.counts[i] == 0:
                indeces_to_remove.append(i)
        #remove the reads that have counts of 0
        self.data = [j for i, j in enumerate(self.data) if i not in indeces_to_remove]
        self.counts = [j for i, j in enumerate(self.counts) if i not in indeces_to_remove]
        self.downsamp_indicator = True
        return

    @staticmethod
    def parse_needle_output_and_get_sub_sum(alignment_filepath, counts_sub_set):
        """This script simultaneously parses the needle output in order to get the percent distance for each alignment, and scales that value by the count of cdr3 sequences that are of each corresponding type. Returns the sum of these scaled values."""
        filein = open(alignment_filepath, "r")
        sub_sum_of_distances = 0
        index = 0
        for i in filein:
            if i[:11] == '# Identity:':
                match = re.search('\(\s*([0-9]+\.[0-9])', i)
                percent_dist = 100 - float(match.group(1))
                sub_sum_of_distances += counts_sub_set[index] * percent_dist
                index += 1
        filein.close()
        return sub_sum_of_distances

    def calc_diversity_pi(self, method='vsearch', path_to_vsearch=None, path_to_needle=None, temp_dirpath=None, num_parallel_cores=12):
        """
        This method calculates the diversity statistic pi, which is the mean genetic distance between all possible pairs of sequnces in the sample.
        method - There are different aligning software to use in order to accomplish the all-by-all alignments. Acceptable values are:
            'vsearch' - This is a software package that can cluster genetic sequences, and also has an all pairwise alignment function.
            'needle' - This is the implementation of the Needleman-Wunsch global alignment algorithm in the EMBOSS package.
        path_to_vsearch - This gives the path to the vsearch excecutable. If this is None (default) then the path is assumed to be 'vsearch' (i.e. in $PATH).
        path_to_needle - This gives the path to the needle excecutable. If this is None (default) then the path is assumed to be 'needle' (i.e. in $PATH).
        temp_dirpath - This gives the path to the directory for which the temporary alignment files will be written. These files can get quite big, so make sure there is enough space where ever this path leads.
        num_parallel_cores - This gives the number of parallel jobs that vsearch can use. vsearch's default is to use as many cores as it can, and this slows down the cluster for others. So, if this method is being used by an array job on the cluster, we need to set this to a fixed amount, and then let the vsearch know what this is when running the program. Default is 12.
        """
        #make sure method parameter makes sense
        if method != 'vsearch' and method != 'needle':
            print 'method parameter does not have an acceptable value:', method
            return

        #if there is only one uniqe seq in the data then pi=0
        if len(self.data) == 1:
            self.pi = 0.
            return self.pi

        #if using needle
        if method == 'needle':
            #if the path to needle is not provided,
            #it is assumed to be already in the PATH
            if path_to_needle == None:
                path_to_needle = 'needle'
            total_unique_seqs = len(self.data)
            sum_of_distances = 0
            for i in xrange(len(self.data) - 1):
                #get 1st seq
                seq1 = self.data[i]['seq']
                random_suffix = str(random.random())
                #this file will contain the 2nd seqs
                temp_seq2_file = '2_%s.fasta' % random_suffix
                fileout = open(temp_seq2_file, "w")
                for j in xrange(i+1, total_unique_seqs):
                    #I think needle runs faster if the '-bsequence' is a fasta file
                    #with many sequences in it. This is why I save the sequences to
                    #disk in such a way. It makes for kind of wonky code, as I
                    #have to loop through each 'bsequence' to write them to disk,
                    #and then loop through each entry in the needle output to get
                    #the percent distance. But I think it ends up being faster.
                    fileout.write('>' + str(j) + '\n' + self.data[j]['seq'] + '\n')
                fileout.close()
                alignment_filepath = 'alignment_%s.txt' % random_suffix
                #do the alignment using the the Needleman-Wunsch algorithm
                #from EMBOSS
                p = Popen([path_to_needle, '-asequence', 'asis:'+seq1, '-bsequence', temp_seq2_file, '-gapopen', '10.0', '-gapextend', '0.5', '-outfile', alignment_filepath, '-brief', 'Y'], stdout=PIPE, stderr=PIPE)
                p.communicate()
                sub_sum_of_distances = self.parse_needle_output_and_get_sub_sum(alignment_filepath, self.counts[i+1:])
                subprocess.call(['rm', temp_seq2_file, alignment_filepath])
                sum_of_distances += self.counts[i] * sub_sum_of_distances
            self.pi = sum_of_distances / round(comb(self.total, 2))

        #else if using vsearch
        elif method == 'vsearch':
            if temp_dirpath == None:
                temp_dirpath = os.getcwd()
            #if downsampling has occured, then need to calc pi on downsampled data.
            #Which means we need to write the downsampled data to disk before
            #running vsearch
            random_suffix = str(random.random()) + 'temp'
            if self.downsamp_indicator == True:
                temp_input_fasta = '%s_%s.fasta' % (temp_dirpath, random_suffix)
                self.write_to_disk(temp_input_fasta)
            #make out filepath
            output_filepath = '%s/out_%s' % (temp_dirpath, random_suffix)
            #run vsearch all pairwise alignments
            if self.downsamp_indicator == True:
                p = Popen([path_to_vsearch, '--allpairs_global', temp_input_fasta, '--userout', output_filepath, '--acceptall', '--userfields', 'query+target+id', '--threads', str(num_parallel_cores)], stdout=PIPE, stderr=PIPE)
                out_err = p.communicate()
                subprocess.call(['rm', temp_input_fasta])
            elif self.downsamp_indicator == False:
                p = Popen([path_to_vsearch, '--allpairs_global', self.filepath, '--userout', output_filepath, '--acceptall', '--userfields', 'query+target+id', '--threads', str(num_parallel_cores)], stdout=PIPE, stderr=PIPE)
                out_err = p.communicate()
            #parse output
            filein = open(output_filepath, "r")
            sum_of_distances = 0.
            comparison_counts = 0
            for line_count, i in enumerate(filein):
                line = i[:-1].split('\t')

                try:
                    a_seq_header = line[0].split('|')
                    a_seq_count = [int(j.split('=')[1]) for j in a_seq_header if j.split('=')[0] == self.count_attribute_name][0]
                    b_seq_header = line[1].split('|')
                    b_seq_count = [int(j.split('=')[1]) for j in b_seq_header if j.split('=')[0] == self.count_attribute_name][0]
                except:
                    print self.filepath
                    print 'line count:', line_count
                    print i
                    print line
                    print '####################'
                    sys.stdout.flush()
                    return 'error'
                
                try:
                    percent_identity = float(line[2])
                except:
                    print self.filepath
                    print 'line count:', line_count
                    print i
                    print line
                    print '####################'
                    sys.stdout.flush()
                    return 'error'

                comparisons = a_seq_count * b_seq_count
                comparison_counts += comparisons
                sum_of_distances += comparisons * (100-percent_identity)
            filein.close()
            subprocess.call(['rm', output_filepath])
            comparison_count = comb(self.total, 2)
            self.pi = sum_of_distances / comparison_count

        return self.pi

    def calc_diversity_pi_compCluster(self, num_alignments_per_job=100000, method='needle', path_to_needle='needle', temp_dirpath=None):
        """
        This method calculates the diversity statistic pi for the set of sequences in the sample. It differs from the 'calc_diversity_pi' method in that it uses a computational cluster in order to handle sample with a whole lot of sequences. The results from these two methods should be identical however. This method submits an array job to the cluster, and calls the script 'pi_calculator_compCluster' (below). Each of the se array jobs calculate a portion of the genetic distances b/t seqs in the sample, then once they are finished, this method combines them all to get a mean. Returns this value.
        'num_alignments_per_job' - This gives the number of alignments that each job on the cluster will be responsible for completing. The higher the number, the less total jobs will be run. The lower, the more total jobs will be run. Adjust as needeed.
        'method' - This gives the method that will be used to get the genetic distance from the pairs of seqs. The exceptable values for this are:
            'needle' - This means the needleman-wunsch global alignment algorithm will be used in a program called 'needle' in the EMBOSS package
            'pairwise2' - This means the 'pairwise2' Biopython package will be used, which also implements the needleman-wunsch global alignment algorithm. Might be faster than needle though.
        'path_to_needle' - This is the path to the needle program by EMBOSS that impliments the needleman-wunsch global alignment algorithm. This is what we use to get the genetic distance between two seqs. This parameter is ignored unless method='needle'
        temp_dirpath - This gives the path to the temporary directory that will contain the fasta files to be submitted to the aligner (if the chosen aligner needs this), and the sub_sum_of_distances values outputed by 'pi_calculator_compCluster'. The temp fasta files can get quite large, so make sure that there is enough space where ever this path leads. If temp_dirpath=None (default) then the current working directory is used.
        """
        #check if there is only one seq. If so, pi=0
        if len(self.data) == 1:
            self.pi = 0.0
            return self.pi

        #check the method parameter
        if method != 'needle' and method != 'pairwise2':
            print '"method" parameter does not have an exceptable value:', method
            return

        #make temp output directory
        if temp_dirpath == None:
            temp_dir = tempfile.mkdtemp(dir=os.getcwd())
        else:
            temp_dir = tempfile.mkdtemp(dir=temp_dirpath)

        ######################################################
        # The below code make temp input fasta files for the #
        # desired aligner, if appropriate. It evenly divides #
        # all the pairwise seq comparisons based upon        #
        # 'num_alignments_per_job'                           #
        ######################################################
        if method == 'needle':
            print '\t\t\tCompiling temp alignment input files'
            start_time = time.time()
            comparison_count = 0
            job_num = 1
            #for all pairs of seqs in the sample
            for i in itertools.combinations(self.data, 2):
                comparison_count += 1
                a_seq = i[0]['seq']
                a_seq_count = i[0]['count']
                b_seq = i[1]['seq']
                b_seq_count = i[1]['count']
                #if at start of data for new job
                if comparison_count == 1:
                    os.makedirs('%s/%s/a_seqs' % (temp_dir, job_num))
                    os.makedirs('%s/%s/b_seqs' % (temp_dir, job_num))
                    seq_set = 0
                    a_seq_old = 'butt'
                #if the 1st seq has changed then make new set of
                #input files for aligner
                if a_seq_old != a_seq:
                    try:#use 'try' b/c if 1st seq pair then will throw NameError
                        fileout_a_seq.close()
                        fileout_b_seq.close()
                    except:
                        pass
                    seq_set += 1
                    fileout_a_seq = open('%s/%s/a_seqs/%s.fasta' % (temp_dir, job_num, seq_set), "w")
                    fileout_a_seq.write('>%s_acount=%s\n%s\n' % (comparison_count, a_seq_count, a_seq))
                    fileout_b_seq = open('%s/%s/b_seqs/%s.fasta' % (temp_dir, job_num, seq_set), "w")
                fileout_b_seq.write('>%s_bcount=%s\n%s\n' % (comparison_count, b_seq_count, b_seq))
                #if have reached the max number of alignments per job
                #start making data for new job
                if comparison_count == num_alignments_per_job:
                    comparison_count = 0
                    job_num += 1
                a_seq_old = a_seq

            job_num = int(math.ceil(float(comb(len(self.data), 2)) / num_alignments_per_job))
            fileout_a_seq.close()
            fileout_b_seq.close()
            end_time = time.time()
            print '\t\t\ttime elapsed:', end_time-start_time
        ######################################################

        elif method == 'pairwise2':
            #It seems that not all the nodes have BioPython installed
            #on them. So, we are only importing these modules if needed.
            #May need to revisit this problem.
            from Bio import pairwise2
            from Bio.SubsMat import MatrixInfo
            job_num = len(self.data) - 1
            #temporarily write seqs to disk so that they can be read
            #by each of the jobs on the cluster
            temp_seq_and_count_filepath = temp_dir + '/seqs_and_counts'
            fileout = open(temp_seq_and_count_filepath, "w")
            for i, j in itertools.izip(self.data, self.counts):
                fileout.write('%s\t%s\n' % (i['seq'].upper(), j))
            fileout.close()

        #submit array job
        print '\t\t\tsubmiting array job.'
        start_time = time.time()
        p = Popen(['qsub', '-t', '1-%s' % job_num, 'sequence_sample_class.py', 'pi_calculator_compCluster', temp_dir+'/', method, path_to_needle], stdout=PIPE, stderr=PIPE)
        out_err = p.communicate()
        job_id = out_err[0].split()[2].split('.')[0]
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
        end_time = time.time()
        print '\t\t\ttime elapsed:', end_time-start_time
        #when job is done sum up the sub_sum_of_distances
        #and get the mean.
        #Also, check to make sure there is the right number of
        #files (i.e. right number of comparisons) in the
        #directory. If not, probably some jobs failed
        sum_of_distances = 0.
        for i in xrange(job_num):
            input_filepath = '%s/%s_subsum' % (temp_dir, i+1)
            #check if necessary subsum file is there. If not
            #it could be because the server is quite slow,
            #so keep checking this for awhile.
            for j in xrange(1000):
                try:
                    filein = open(input_filepath, "r")
                    break
                except IOError:
                    time.sleep(1)
            sum_of_distances += float(filein.readline())
            filein.close()
        num_comparisons = comb(self.total, 2)
        self.pi = sum_of_distances / num_comparisons
        subprocess.call(['rm', '-r', temp_dir])
        return self.pi
        
    def cluster_seqs(self, method='vsearch', min_ident_within_cluster=0.97, path_to_vsearch=None, temp_dirpath=None, max_edit_distance=1, overwrite_data=False, output_network_filepath=None, output_node_attribute_filepath=None, node_attributes_to_add=['count'], names_for_node_attributes_in_header=None, path_to_needle=None, scale_freqs_by=None, dists_dstrb_output_filepath=None, add_indiv_seq_attribute=None):
        """
        This method will cluster the sequences in the data.
        method - This gives the method that will be used for clustering. Acceptable values are:
            'vsearch' - This means that the program Vsearch will be used to cluster based on a given identity threshold. If this method is chosen then centroids of the clusters will be returned.
            'by_edit_distance' - This means that the sequences will be clustered by a simple requirement that all seqs within a cluster must have an edit distance with at least one other seq in the cluster that is less than or equal to a provided max edit distance cutoff. This requires an all by all pairwise alignment, so be carefule with large datasets. This also used the program 'needle' to do the alignments, which needs to be in the PATH.
            'by_edit_distance_seqanpy' - This means that edit distance is still used as the metric to cluster seqs by, but the package 'seqanpy' is used to do the alignments, instead of 'needle', this may be faster.
        min_ident_within_cluster - This is an important parameter that gives the distance allowed between the sequences within a cluster. It should be between 0 and 1. If a target sequence, when aligned to a cluster centroid, has an identity lower then this value, then it is not added to that cluster. This parameter is only considered if method='vsearch'.
        path_to_vsearch - This gives the path to the vsearch excecutable. If this is None (default) then the path is assumed to be 'vsearch' (i.e. in $PATH). This parameter is only considered if method='vsearch'.                                     
        'temp_dirpath' = This gives the path to the directory for which a temporary output directory will be created in. It will contain temp fasta files and output from vsearch. The temp directory that is created within 'temp_dirpath' will be removed at the end of the script, but the 'temp_dirpath' directory itself will not be removed. If this equals None, the the current working directory is used.
        max_edit_distance - Int. Default, 1. This gives the maximum edit distance between any two seqs, for them to be assigned to the same cluster. This parameter is only considered if method='by_edit_distance'
        overwrite_data - Boolean. Default, False. If True, this will overwrite the data in self.data with the results of the clustering. So, each element in self.data will be a representative seq for the cluster, with its count (i.e. number of reads that map to that cluster), and frequence (i.e. proportion of reads that map to that cluster). If the method is 'vsearch' then the representative seq will be the centroid seq, and if method='by_edit_distance', then the representative seq will be the most numerous seq in the cluster.
        output_network_filepath - If defined (default, None), this will, 1) instruct the method to creat a network file (as a simple interaction file ('.sif')), and 2) will be the path to this file. This is only considered if method='by_edit_distance'
        output_node_attribute_filepath - This gives the path to the file that will contain the attributes for each of the unique nodes (i.e. unique seqs) in the data. This is only considered if 'output_network_filepath' is defined.
        node_attributes_to_add - This is a list of sequence attributes to include in the node attribute file. This parameter is only considered if 'output_node_attribute_filepath' is defined. The default is to only add the value for the 'count' attribute (i.e. the value for self.data[some_index]['count']). Acceptable values within this list of strings are:
            any string that is a key for one of the attributes for the sequences. So, could be 'count', 'total_freq', 'timepoint', etc. if those are keys to the dictionarys within self.data[some_index], or self.data[some_index]['other']
            'element_X_of_seq_id' - sometimes info about the seq is encoded in the seq ID, delimited by '_'. This means that the 'X'th element of the seq ID (when delimited by a '_') will be included in the node attributes file.
        names_for_node_attributes_in_header - This is an optional argument, which if defined (default, None), it will give the names for each of the node attributes included in the 'output_node_attribute_filepath' header. This parameter is only considered if 'output_node_attribute_filepath'. If 'output_node_attribute_filepath' is defined and this is not, then the default is to simply use the values in 'node_attributes_to_add'. The names in this list need to be in the same order as the attributes in 'node_attributes_to_add'.
        path_to_needle - This gives the path to the needle executable. This is required if using the computational cluster for an array job because the PATH variable is all messed up when sending jobs to parallele nodes. One can also use this if needle is not in the PATH
        scale_freqs_by - If defined (default, None), this should be a float that gives information as how to scale the frequencies calculated for each of the clusters in the sample. For example if a certain cluster has a calculated freq of .5, and one is to scale this freq by .25, then the resulting freq will be .125. Note that if this is used, freqs will not sum to one!
        dists_dstrb_output_filepath - This gives the filepath fr which all the the edit distances will be written for each pairwise sequence comparison. This parameter is only considered if method=='by_edit_distance'.
        add_indiv_seq_attribute - If defined (default, None), this should be a list of strings where are string is the same as an attribute for each of the individual seqs that should be included in the header for each of the clusters. Similar to the 'indiv_seq_freqs' attribute that is included for each of the seqs belonging to a cluster, this will give the provided attribute's values for each of the seqs in a cluster. This is only applied if 'overwrite_data' is set to True
        
        OUTPUT:
        If method='vsearch', returns a list called 'centroids'. This is a list where each element represents a cluster. In each element the cluster count (counts of seqs in the cluster), the cluster frequency (count / total number of seqs in data), and the cluster centroid (the sequence in the cluster that is closest to the center) are reported, as a list, in that order.
        If method='by_edit_distance', then adds a 'cluster_id' attribute to the self.data['other'] for each seq. The 'cluster_id' attribute is an int. where each unique value indicates a unique cluster.
        """

        time_start = time.time()

        if temp_dirpath == None:
            temp_dirpath = tempfile.mkdtemp(dir=os.getcwd()) + '/'
        temp_dirpath = tempfile.mkdtemp(dir=temp_dirpath) + '/'

        if method == 'vsearch':
            centroids = []
            rand_string = str(random.random())[2:]
            temp_fasta_filepath = '%s%s.fasta' % (temp_dirpath, rand_string)
            fileout = open(temp_fasta_filepath, "w")
            #for each seq
            count = 0
            for i in self.data:
                count += 1
                #write to temp fasta file
                fileout.write('>%s;size=%s;\n%s\n' % (count, i['count'], i['seq']))
            fileout.close()
            #now cluster seqs
            centroid_filepath = '%s%s_centroids.fasta' % (temp_dirpath, rand_string)
            p = Popen([path_to_vsearch, '--cluster_size', temp_fasta_filepath, '--id', str(min_ident_within_cluster), '--sizein', '--sizeout', '--centroids', centroid_filepath, '--fasta_width', str(0)], stdout=PIPE, stderr=PIPE)
            out_err = p.communicate()
            #retrieve centroid seqs for each cluster
            filein = open(centroid_filepath, "r")
            id = 0
            for i in filein:
                if i[0] == '>':
                    id += 1
                    line = i[:-1].split(';')
                    count = int(line[1].split('=')[1])
                    freq = count / self.total
                    if scale_freqs_by:
                        freq *= scale_freqs_by
                else:
                    centroids.append([count, freq, i[:-1]])
            filein.close()
            subprocess.call(['rm', centroid_filepath, temp_fasta_filepath])
            subprocess.call(['rm', '-r', temp_dirpath])
            if overwrite_data:
                self.data = []
                self.counts = []
                self.diversity = None
                for i in centroids:
                    self.data.append({'seq':i[2], 'count':i[0], 'other':{'freq':i[1]}})
                    self.counts.append(i[0])
            return centroids

        elif method == 'by_edit_distance' or method == 'by_edit_distance_seqanpy':

            if method == 'by_edit_distance':
                if not path_to_needle:
                    path_to_needle = 'needle'
                if output_network_filepath:
                    fileout_network = open(output_network_filepath, "w")
                if dists_dstrb_output_filepath:
                    fileout_dists = open(dists_dstrb_output_filepath, "w")
                    fileout_dists.write('genetic_distance\n')
                random_suffix = str(random.random())
                list_of_indices = range(len(self.data))
                #give all seqs a unique cluster ID.
                #cluster_dic is a dictinary where each Key is a unique cluster ID, and the def is the list of seq indices that belong to that cluster
                cluster_dic = {}
                #connections_dic is a dictionary where each Key is a seq ID and the def is a list of seq indices that have direct connections (i.e. edges) with that seq.
                connections_dic = {}
                for i in list_of_indices:
                    self.data[i]['other']['cluster_id'] = i+1
                    cluster_dic[i+1] = [i]
                    connections_dic[self.data[i]['id']] = []
                for seq1_index in list_of_indices[:-1]:
                    temp_fasta_filepath = '%stemp_fasta_%s_%s.fasta' % (temp_dirpath, random_suffix, seq1_index)
                    fileout_tempFasta = open(temp_fasta_filepath, "w")
                    alignment_filepath = '%stemp_alignment_%s_%s' % (temp_dirpath, random_suffix, seq1_index)
                    for seq2_index in list_of_indices[seq1_index+1:]:
                        fileout_tempFasta.write('>%s\n%s\n' % (self.data[seq2_index]['id'], self.data[seq2_index]['seq']))
                    fileout_tempFasta.close()
                    p = Popen([path_to_needle, '-asequence', 'asis:'+self.data[seq1_index]['seq'], '-bsequence', temp_fasta_filepath, '-gapopen', '10.0', '-gapextend', '0.5', '-outfile', alignment_filepath, '-brief', 'Y'], stdout=PIPE, stderr=PIPE)
                    out, err = p.communicate()
                    subprocess.call(['rm', temp_fasta_filepath])
                    filein = open(alignment_filepath, "r")
                    seq2_index = list_of_indices[seq1_index+1]
                    for i in filein:
                        if i[:11] == '# Identity:':
                            ratio = i.split()[2].split('/')
                            edit_dist = int(ratio[1]) - int(ratio[0])
                            if dists_dstrb_output_filepath:
                                fileout_dists.write('%s\n' % edit_dist)
                            #if seq1 and seq2 are close enough to be in same cluster
                            if edit_dist <= max_edit_distance:
                                connections_dic[self.data[seq1_index]['id']].append(seq2_index)
                                connections_dic[self.data[seq2_index]['id']].append(seq1_index)
                                seq1_cluster_id = self.data[seq1_index]['other']['cluster_id']
                                seq2_cluster_id = self.data[seq2_index]['other']['cluster_id']
                                if seq1_cluster_id != seq2_cluster_id:
                                    #change the cluster ID for all the seqs that map to the cluster that seq2 is in
                                    for j in cluster_dic[seq2_cluster_id]:
                                        self.data[j]['other']['cluster_id'] = seq1_cluster_id
                                    #add all the seq indeces in cluster_dic for the seq2 cluster to the seq1 cluster
                                    cluster_dic[seq1_cluster_id] += cluster_dic[seq2_cluster_id]
                                    #delete the seq2 cluster entry in cluster_dic
                                    del cluster_dic[seq2_cluster_id]
                                #add network connection
                                if output_network_filepath:
                                    fileout_network.write('%s\tpp\t%s\n' % (self.data[seq1_index]['id'], self.data[seq2_index]['id']))
                            seq2_index += 1
                    filein.close()
                    subprocess.call(['rm', alignment_filepath])

            elif method == 'by_edit_distance_seqanpy':
                import seqanpy
                if output_network_filepath:
                    fileout_network = open(output_network_filepath, "w")
                if dists_dstrb_output_filepath:
                    fileout_dists = open(dists_dstrb_output_filepath, "w")
                    fileout_dists.write('genetic_distance\n')
                list_of_indices = range(len(self.data))
                #give all seqs a unique cluster ID.
                #cluster_dic is a dictinary where each Key is a unique cluster ID, and the def is the list of seq indices that belong to that cluster
                cluster_dic = {}
                #connections_dic is a dictionary where each Key is a seq ID and the def is a list of seq indices that have direct connections (i.e. edges) with that seq.
                connections_dic = {}
                for i in list_of_indices:
                    self.data[i]['other']['cluster_id'] = i+1
                    cluster_dic[i+1] = [i]
                    connections_dic[self.data[i]['id']] = []
                for seq1_index in list_of_indices[:-1]:
                    for seq2_index in list_of_indices[seq1_index+1:]:
                        alignment = seqanpy.align_global(self.data[seq1_index]['seq'], self.data[seq2_index]['seq'])
                        edit_dist = 0
                        alignment_len = len(alignment[1])
                        for i in xrange(alignment_len):
                            if alignment[1][i] != alignment[2][i]:
                                edit_dist += 1
                        if dists_dstrb_output_filepath:
                            fileout_dists.write('%s\n' % edit_dist)
                        #if seq1 and seq2 are close enough to be in same cluster
                        if edit_dist <= max_edit_distance:
                            connections_dic[self.data[seq1_index]['id']].append(seq2_index)
                            connections_dic[self.data[seq2_index]['id']].append(seq1_index)
                            seq1_cluster_id = self.data[seq1_index]['other']['cluster_id']
                            seq2_cluster_id = self.data[seq2_index]['other']['cluster_id']
                            if seq1_cluster_id != seq2_cluster_id:
                                #change the cluster ID for all the seqs that map to the cluster that seq2 is in
                                for j in cluster_dic[seq2_cluster_id]:
                                    self.data[j]['other']['cluster_id'] = seq1_cluster_id
                                #add all the seq indeces in cluster_dic for the seq2 cluster to the seq1 cluster
                                cluster_dic[seq1_cluster_id] += cluster_dic[seq2_cluster_id]
                                #delete the seq2 cluster entry in cluster_dic
                                del cluster_dic[seq2_cluster_id]
                            #add network connection
                            if output_network_filepath:
                                fileout_network.write('%s\tpp\t%s\n' % (self.data[seq1_index]['id'], self.data[seq2_index]['id']))

            if dists_dstrb_output_filepath:
                fileout_dists.close()
            #make sure all seqs are represented in network file by adding them all to the end. Then close filehandle
            if output_network_filepath:
                for i in self.data:
                    fileout_network.write('%s\n' % i['id'])
                fileout_network.close()
            #make the ID of each cluster correspond to the seq ID that is most numerous in that cluster
            new_cluster_dic = {}
            for cluster_id in cluster_dic:
                #get counts for each seq in cluster
                seq_counts = [self.data[i]['count'] for i in cluster_dic[cluster_id]]
                #find max seq in terms of count
                max_seq_index = seq_counts.index(max(seq_counts))
                #get the seq ID that corresponds to the max seq
                max_seq_id = self.data[cluster_dic[cluster_id][max_seq_index]]['id']
                #change the cluster ID for each seq in the cluster
                for i in cluster_dic[cluster_id]:
                    self.data[i]['other']['cluster_id'] = str(max_seq_id)
                #give new cluster dic the correct cluster ID
                new_cluster_dic[max_seq_id] = cluster_dic[cluster_id][:]
            cluster_dic = new_cluster_dic

            if output_node_attribute_filepath:
                fileout = open(output_node_attribute_filepath, "w")
                if names_for_node_attributes_in_header:
                    header = 'node_ID\t%s\n' % '\t'.join(names_for_node_attributes_in_header)
                else:
                    header = 'node_ID\t%s\n' % '\t'.join(node_attributes_to_add)
                fileout.write(header)
                for i in self.data:
                    attribute_values_list = []
                    for attribute in node_attributes_to_add:
                        if attribute[:8] == 'element_':
                            element = int(attribute.split('_')[1])
                            attribute_values_list.append(i['id'].split('_')[element])
                        else:
                            try:
                                #look for attribute in self.data
                                attribute_values_list.append(i[attribute])
                            except KeyError:
                                #look for attributes in self.data[i]['other']
                                attribute_values_list.append(i['other'][attribute])
                    fileout.write('%s\t%s\n' % (i['id'], '\t'.join([str(j) for j in attribute_values_list])))
                fileout.close()

            if overwrite_data:
                #gather new data
                new_data = {}
                self.diversity = None#diversity statistic no longer applies (if calculated)
                if add_indiv_seq_attribute:
                    additional_seq_attrb_names = ['indiv_seq_%s' % seq_attrb for seq_attrb in add_indiv_seq_attribute]
                for i in self.data:
                    cluster_seq_id = i['other']['cluster_id']
                    freq = i['count']/self.total
                    if scale_freqs_by:
                        freq *= scale_freqs_by
                    try:
                        new_data[cluster_seq_id]['other']['indiv_seq_counts'].append(i['count'])
                        new_data[cluster_seq_id]['other']['indiv_seqs'].append(i['seq'])
                        new_data[cluster_seq_id]['other']['indiv_seq_ids'].append(i['id'])
                        new_data[cluster_seq_id]['count'] += i['count']
                        new_data[cluster_seq_id]['other']['total_freq'] += freq
                        new_data[cluster_seq_id]['other']['indiv_seq_freqs'].append(i['count']/self.total)
                        if add_indiv_seq_attribute:
                            for j in xrange(len(add_indiv_seq_attribute)):
                                new_data[cluster_seq_id]['other'][additional_seq_attrb_names[j]].append(i['other'][add_indiv_seq_attribute[j]])
                    except KeyError:
                        new_data[cluster_seq_id] = {'id':cluster_seq_id, 'count':i['count'], 'other':{'cluster_id':i['other']['cluster_id'], 'indiv_seq_ids':[i['id']], 'indiv_seqs':[i['seq']], 'indiv_seq_counts':[i['count']], 'total_freq':freq, 'indiv_seq_freqs':[i['count']/self.total]}}
                        if add_indiv_seq_attribute:
                            for j in xrange(len(add_indiv_seq_attribute)):
                                new_data[cluster_seq_id]['other'][additional_seq_attrb_names[j]] = [i['other'][add_indiv_seq_attribute[j]]]
                #now update data
                self.data = []
                self.counts = []
                for i in new_data:
                    self.data.append(new_data[i])
                    self.counts.append(new_data[i]['count'])
                    max_seq_index = new_data[i]['other']['indiv_seq_counts'].index(max(new_data[i]['other']['indiv_seq_counts']))
                    self.data[-1]['seq'] = new_data[i]['other']['indiv_seqs'][max_seq_index]
                    self.data[-1]['other']['indiv_seq_counts'] = ','.join([str(j) for j in self.data[-1]['other']['indiv_seq_counts']])
                    self.data[-1]['other']['indiv_seqs'] = ','.join(self.data[-1]['other']['indiv_seqs'])
                    self.data[-1]['other']['indiv_seq_ids'] = ','.join(self.data[-1]['other']['indiv_seq_ids'])
                    self.data[-1]['other']['indiv_seq_freqs'] = ','.join([str(j) for j in self.data[-1]['other']['indiv_seq_freqs']])
                    if add_indiv_seq_attribute:
                        for j in xrange(len(add_indiv_seq_attribute)):
                            self.data[-1]['other'][additional_seq_attrb_names[j]] = ','.join([str(k) for k in self.data[-1]['other'][additional_seq_attrb_names[j]]])
            subprocess.call(['rm', '-r', temp_dirpath])

            time_end = time.time()
            print time_end - time_start

            return cluster_dic, connections_dic

        else:
            print 'method dos not have a valid value: method =', method
            print 'Aborting.'
            return

    def write_to_disk(self, output_filepath, append_to_file=False):
        """
        This method simply writes the data in self.data to disk, in fasta format. Only writes seq ID and count information. Everything else in headers is left out.
        output_filepath - Path to the file that will be written.
        append_to_file - Boolean. If True (default, False), this will append to the provided filepath, rather than overwriting.
        """
        if append_to_file:
            fileout = open(output_filepath, "a")
        else:
            fileout = open(output_filepath, "w")
        for i in self.data:
            fileout.write('>%s|%s=%s\n%s\n' % (i['id'], self.count_attribute_name, i['count'], i['seq']))
        fileout.close()
        return

    @staticmethod
    def write_one_seq_entry_to_disk(self, data_entry, fileout, add_string_to_each_id_being_written=None, write_attribute_value_instead=None):
        """
        This is an internal method that is used by other methods to write one entry of a fasta file to an output file.
        add_string_to_each_id_being_written - If defined (default, None), this will be a string that is to be added to each of the seq IDs when they are written to disk. Does not permenantly change the data object in memory, just adds the string as they are written to disk.
        write_attribute_value_instead - If defined (default, None), this will write the value of the provided attribute. If defined this must be the name of an attribute in the seq data, and the value of this attribute should be a sequence of some sort. This is essentially a hack to allow the writing of CDR3 seq data.
        """
        if add_string_to_each_id_being_written:
            seq_id = add_string_to_each_id_being_written + data_entry['id']
        else:
            seq_id = data_entry['id']
        other_headers = []
        for i in data_entry['other']:
            if type(data_entry['other'][i]) is list:
                value = ','.join([str(j) for j in data_entry['other'][i]])
            else:
                value = data_entry['other'][i]
            other_headers.append('%s=%s' % (i, value))
        header = '>%s|%s=%s|%s\n' % (seq_id, self.count_attribute_name, data_entry['count'], '|'.join(other_headers))
        if write_attribute_value_instead:
            fileout.write('%s%s\n' % (header, data_entry['other'][write_attribute_value_instead]))
        else:
            fileout.write('%s%s\n' % (header, data_entry['seq']))
        return

    def write_full_data_to_disk_fasta(self, output_filepath=None, append_to_file=False, seq_id_fileout_dic=None):
        """
        This method writes all the data to disk. Including the seq data, and all the attributes that are defined for each sequence entry. Writes in fasta format.
        output_filepath - The path to the file for which all this data will be written.
        append_to_file - Boolean. If True (default, False), this will append to the provided filepath, rather than overwriting.
        seq_id_fileout_dic - If defined (default, None), then this will give a dictionary that gives different output filepaths for each of the seqs. Note that there is a limit on the number of open file handles in python for some reason. This parameter is only defined if output_filepath is None.
        """
        #check to make sure there are seqs in the data
        if len(self.data) == 0:
            print 'No data in file, so nothing to write to disk. Aborting.'
            return
        if not output_filepath and not seq_id_fileout_dic:
            print "Need to define either output filepath or an output filepath dictionary. Aborting."
            return
        if output_filepath:
            if append_to_file:
                fileout = open(output_filepath, "a")
            else:
                fileout = open(output_filepath, "w")
            for i in self.data:
                self.write_one_seq_entry_to_disk(self, i, fileout)
            fileout.close()
        elif seq_id_fileout_dic:
            fileout_dic = {}
            for i in seq_id_fileout_dic:
                try:
                    fileout_dic[seq_id_fileout_dic[i]]
                except KeyError:
                    fileout_dic[seq_id_fileout_dic[i]] = open(seq_id_fileout_dic[i], "w")
            for i in self.data:
                self.write_one_seq_entry_to_disk(self, i, fileout_dic[seq_id_fileout_dic[i['id']]])
            for i in fileout_dic:
                fileout_dic[i].close()
        return

    def write_subset_of_seqs_to_disk(self, seq_indices, output_filepath, append_to_file=False, seq_ids=None, add_string_to_each_id_being_written=None, write_attribute_value_instead=None):
        """
        This method is similar to 'write_full_data_to_disk_fasta', but it will only write a subset of the seqs in the data to disk.
        seq_indices - This is a list of ints, where each int is the index of a desired seq in self.data to write to file
        output_filepath - The path to the file for which all this data will be written.
        append_to_file - Boolean. If True (default, False), this will append to the provided filepath, rather than overwriting.
        seq_ids - If defined (default, None), this will be a list of strings, where each string is the ID of a seq that is to be written to disk. If this is defined, then the 'seq_indices' parameter will be ignored, and seqs will be written based on their ID value instead.
        add_string_to_each_id_being_written - If defined (default, None), this will be a string that is to be added to each of the seq IDs when they are written to disk. Does not permenantly change the data object in memory, just adds the string as they are written to disk.
        write_attribute_value_instead - If defined (default, None), this will write the value of the provided attribute. If defined this must be the name of an attribute in the seq data, and the value of this attribute should be a sequence of some sort. This is essentially a hack to allow the writing of CDR3 seq data.
        """
        if append_to_file:
            fileout = open(output_filepath, 'a')
        else:
            fileout = open(output_filepath, 'w')
        if not seq_ids:
            for i in seq_indices:
                self.write_one_seq_entry_to_disk(self, self.data[i], fileout, add_string_to_each_id_being_written, write_attribute_value_instead)
        else:
            no_id = True
            for i in self.data:
                if i['id'] in seq_ids:
                    no_id = False
                    self.write_one_seq_entry_to_disk(self, i, fileout, add_string_to_each_id_being_written, write_attribute_value_instead)
            if no_id:
                print 'Did not find any of the seq IDs in the data, so an empty file results. Output filepath:', output_filepath
        fileout.close()
        return

    def get_seqID_to_count_dic(self):
        """This method creates a dictionary where the index is the unique seq ID and the def for each index is the count for that unique seq."""
        id_count_dic = {}
        for i in self.data:
            id_count_dic[i['id']] = i['count']
        return id_count_dic

    def get_seqID_to_attribute_dic(self, query_attribute_name, only_for_seq_ids=None):
        """
        This method will map the unique seq ID to the value of whatever sequence attribute one is interested in, so long as it is not the 'count' attribute. This looks through the attributes that are stored in the 'other' dic that is within each sequences data dictionary.
        query_attribute_name - This gives the name of the attribute that is to be retrieved.
        only_for_seq_ids - If defined (default, None), this should be a list of seq IDs, and will instruct the script to only return the attribute values for these seq IDs. If the iequals None, then it returns the ID to attribute value dic for each seq in the data.
        """
        id_to_attribute_dic = {}
        for i in self.data:
            if only_for_seq_ids:
                if i['id'] in only_for_seq_ids:
                    id_to_attribute_dic[i['id']] = i['other'][query_attribute_name]
            else:
                id_to_attribute_dic[i['id']] = i['other'][query_attribute_name]
        return id_to_attribute_dic

    def get_mean_attribute_value(self, query_attribute_name, weight_by_counts=True, treat_NAs='ignore'):
        """This method will calculate the mean value for a given query sequence attribute. The sequence attribute must be one of the other attributes that are contained in the 'other' dic if the sequence data.
        weight_by_counts - If True (default), then this weights each of the sequence entries values by the 'count' of that sequence in the data (i.e. weighted average).
        treat_NAs - This informs the method how to treat 'NA' values. Acceptable values are:
            'ignore' - Default. This means that seq entries that have a value of 'NA' are simply ignored. Neither their value nor their count is included in the numerator nor the denominator of the mean calculation, respectively.
            'zero' - This means that a seq entry with a value of 'NA' will be treated as if this value were 0.
        """
        mean_attribute_value = 0.
        denominator = 0
        for i in self.data:
            value = i['other'][query_attribute_name]
            if value == 'NA':
                if treat_NAs == 'ignore':
                    continue
                elif treat_NAs == 'zero':
                    value = 0.
            else:
                value = float(value)
            if weight_by_counts == True:
                value *= i['count']
                denominator += i['count']
            else:
                denominator += 1
            mean_attribute_value += value
        if denominator == 0:
            return 'NA'
        else:
            mean_attribute_value = mean_attribute_value / denominator
            return mean_attribute_value

    def run_igblast(self, output_filepath, path_to_igblast='/Users/nstrauli/tools/ncbi-igblast-1.8.0/', preset_options='change-o', temp_dirpath=None):
        """
        This method will run IgBLAST for each of the sequences in the data.
        output_filepath - The path to the output file for which the output from IgBLAST will be written
        path_to_igblast - This gives the path to the directory that contains the excecutables for IgBLAST. That is, the IgBLAST excecutables should be in a directory called 'bin' in this directory. The other file such as 'database', 'internal_data', and 'optional_file' should be in this directory as well.
        preset_options - This gives the 'mode' for which IgBLAST will be run. This means it will instruct the script what set of parameters should be given to IgBLAST. For example, if using IgBLAST for change-o, then specific output parameters need to be used. Exceptable values are:
            'change-o' - This means that IgBLAST will be run with parameters that are compatible with using change-o on the output files.
        temp_dirpath - This gives the directory for which temporary files will be written. If None (default) then the temp dirpath is the parent directory to 'output_filepath'.
        """
        #if input filepath is in fastq format, then needs to be written to disk in fasta temporarily for IgBLAST input
        if path_to_igblast[-1] != '/':
            path_to_igblast += '/'
        if temp_dirpath == None:
            temp_dirpath = os.path.dirname(output_filepath) + '/'
        else:
            if temp_dirpath[-1] != '/':
                temp_dirpath += '/'
        if not os.path.exists(temp_dirpath):
            os.makedirs(temp_dirpath)
        if self.filetype == 'fastq':
            input_fasta_filepath = '%s%s_%s.fasta' % (temp_dirpath, random.random(), os.path.basename(self.filepath)[:-6])
            self.write_to_disk(output_filepath=input_fasta_filepath)
        elif self.filetype == 'fasta':
            input_fasta_filepath = self.filepath
        cur_dir = os.getcwd()
        os.chdir(path_to_igblast)
        if preset_options == 'change-o':
            p = Popen([path_to_igblast + 'bin/igblastn', '-germline_db_V', 'database/human_igh_v', '-germline_db_D', 'database/human_igh_d', '-germline_db_J', 'database/human_igh_j', '-auxiliary_data', 'optional_file/human_gl.aux', '-domain_system', 'imgt', '-ig_seqtype', 'Ig', '-organism', 'human', '-outfmt', '7 std qseq sseq btop', '-query', input_fasta_filepath, '-out', output_filepath], stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        print 'output for IgBLAST:'
        print out
        print err
        os.chdir(cur_dir)
        if self.filetype == 'fastq':
            subprocess.call(['rm', input_fasta_filepath])
        return

    def get_immune_reads_with_changeo(self, output_filepath, path_to_changeo='/Users/nstrauli/tools/changeo-0.3.9/', path_to_igblast='/Users/nstrauli/tools/ncbi-igblast-1.8.0/', temp_dirpath=None, imgt_ref_seqs=None, add_germline=None, add_selection=False, use_comp_cluster=False, add_divergence=False, overwrite_output=True, translate_VDJ=False, remove_seqs_with_stop=False):
        """
        This method will first use IgBLAST to map reads to immune loci (i.e. antibodies or TCR's), then it will use change-o to parse the results from IgBLAST to make a change-o database (i.e. a tab delimited file) of immune reads. Then it will write those immune reads in fasta format to an output file. Essentially, This method will take the sequence sample data and turn it into immune sample data.
        output_filepath - This is the path to the file that will contain the immune data. It will be fasta formatted.
        path_to_changeo - This is the path to the change-o excecutable. It should be a directory that contains 'bin', 'changeo/', etc.
        path_to_igblast - This gives the path to the directory that contains the excecutables for IgBLAST. That is, the IgBLAST excecutables should be in a directory called 'bin' in this directory. The other file such as 'database', 'internal_data', and 'optional_file' should be in this directory as well.
        temp_dirpath - This gives the directory for which temporary files will be written. If None (default) then the temp dirpath is the parent directory to 'output_filepath'.
        imgt_ref_seqs - This is the path to the files that contain the germline reference sequences for the genes that were used in the IgBLAST database. If this is None (default), then this will equal path_to_igblast + 'ref_seqs_from_imgt/'
        add_germline - This instructs if and what type of germline sequence to add for each sequence. Exceptable values are:
            'dmask' - This means that the germline seq created will have the D-region masked (i.e. all N's). This is needed in order to calculate selection values for each seq using Baseline.
            'full' - This gives the complete germline (i.e. includes V gene).
            'vonly' - This will return only the germline of the V gene.
            None - This will make it so the germline is not returned (default).
        add_selection - This indicates whether or not the method should invoke Baseline to add selection values to each of the seqs in the data. If this is True (default, False), then the method will run Baseline using the 'shazam' R package to get selection values.
        use_comp_cluster - This indicates if a computational cluster is being used via the SGE system. If it is, then special precautions need to be taken for using the correct version of python when invoking changeo.
        add_divergence - If defined (default, False), this will instruct the method to also include a divergence from the inferred unmutated sequence in the immune reads info. Acceptable values are:
            'shazam' - Causes divergence to be calculated using the 'calcObservedMutations' function to be run in shazam. This is an R script and is quite slow, but will give the number of replacement and silent mutations in addition to the total.
            'python' - Causes divergence to be calculated using a simple in house python script that gets the total number of muts between the germline and observed seq. These seqs must be the same length (i.e. does not do an alignment), but they should be because they are imgt numbered from changeo. It also ignores positions that are '.', '-', or 'N' in either the observed or germline seqs.
        overwrite_output - If True (default), then this will write the results to 'output_filepath' regardless of whether or not it exists. If False, the this method will first check to see if the output_filepath exists, and if it does, then it will abort with an error message.
        translate_VDJ - Boolean. If True (default, Flase), this will cause the method to add a field to the header that will be called 'AA_VDJ_sequence' that will give the translated amino acid sequence of the nucleotride sequence contained in the 'SEQUENCE_VDJ' field.
        remove_seqs_with_stop - Boolean. If True (default, False), this will cause the method to not include any seqs that were labled as containing stop codons (by looking at the 'STOP' field given by changeo.)
        """
        if not overwrite_output:
            if os.path.exists(output_filepath):
                print 'The "overwrite_output" parameter is set to false and the designated "output_filepath" already exists.'
                print "output_filepath:", output_filepath
                print 'aborting'
                return
        if path_to_changeo[-1] != '/':
            path_to_changeo += '/'
        if path_to_igblast[-1] != '/':
            path_to_igblast += '/'
        if temp_dirpath == None:
            temp_dirpath = os.path.dirname(output_filepath) + '/'
        else:
            if temp_dirpath[-1] != '/':
                temp_dirpath += '/'
        if not os.path.exists(temp_dirpath):
            os.makedirs(temp_dirpath)
        igblast_output_file = '%s%s_%s.fmt7' % (temp_dirpath, random.random(), os.path.basename(self.filepath)[:-6])
        print ''
        print '############# running IgBLAST #############'
        sys.stdout.flush()
        self.run_igblast(output_filepath=igblast_output_file, path_to_igblast=path_to_igblast, preset_options='change-o', temp_dirpath=temp_dirpath)
        cur_dir = os.getcwd()
        if imgt_ref_seqs == None:
            imgt_ref_seqs = path_to_igblast + 'ref_seqs_from_imgt/'
        print ''
        print '############# making changeo database #############'
        sys.stdout.flush()
        if use_comp_cluster == True:
            p = Popen(['bash', './call_changeo_makedb.bash', path_to_changeo, 'igblast', igblast_output_file, self.filepath, imgt_ref_seqs, os.path.dirname(output_filepath), os.path.basename(output_filepath)[:-6]], stdout=PIPE, stderr=PIPE)
        else:
            os.chdir(path_to_changeo)
            p = Popen(['./bin/MakeDb.py', 'igblast', '-i', igblast_output_file,  '-s', self.filepath, '-r', imgt_ref_seqs, '--regions', '--scores', '--outdir', os.path.dirname(output_filepath), '--outname', os.path.basename(output_filepath)[:-6], '--cdr3'], stdout=PIPE, stderr=PIPE)
            os.chdir(cur_dir)
        out, err = p.communicate()
        print 'output for MakeDb:'
        print out
        print err
        sys.stdout.flush()
        subprocess.call(['rm', igblast_output_file])
        changeo_db_filepath = output_filepath[:-6] + '_db-pass.tab'
        
        if add_germline:
            count = 0
            for germline_type in add_germline:
                count += 1
                print ''
                print '############# adding germline seqs: %s #############' % germline_type
                sys.stdout.flush()
                germline_output_filename = '%s_%s' % (os.path.basename(output_filepath)[:-6], germline_type)
                if use_comp_cluster == True:
                    p = Popen(['bash', './call_changeo_creategermlines.bash', path_to_changeo, changeo_db_filepath, imgt_ref_seqs, germline_type, germline_output_filename], stdout=PIPE, stderr=PIPE)
                else:
                    os.chdir(path_to_changeo)
                    p = Popen(['./bin/CreateGermlines.py', '-d', changeo_db_filepath, '-r', imgt_ref_seqs, '-g', germline_type, '--outname', germline_output_filename], stdout=PIPE, stderr=PIPE)
                    os.chdir(cur_dir)
                out, err = p.communicate()
                print 'output for CreateGermlines:'
                print out
                print err
                sys.stdout.flush()
                changeo_db_filepath_old = changeo_db_filepath[:]
                subprocess.call(['rm', changeo_db_filepath_old])
                changeo_db_filepath = '%s/%s_germ-pass.tab' % (os.path.dirname(output_filepath), germline_output_filename)
                # changeo_db_filepath = output_filepath[:-6] + '_db-pass_germ-pass.tab'

        if add_divergence:
            print ''
            print '############# calculating divergence #############'
            sys.stdout.flush()
            filein = open(changeo_db_filepath, "r")
            header = filein.readline()[:-2].split('\t')
            #get index for the full (i.e. non D-masked) germline seq
            non_D_masked_germline_index = header.index('GERMLINE_IMGT')
            div_stats = {}
            for i in filein:
                line = i[:-2].split('\t')
                seq_id = line[0]
                obs_seq = line[11]
                germline_seq = line[48]
                non_D_masked_germline_seq = line[non_D_masked_germline_index]
                #make sure the observed and germlines seqs are the same length
                if len(obs_seq) != len(germline_seq):
                    print 'the imgt formatted observed and germline seqs are not the same length.'
                    print 'observed seq:', repr(obs_seq)
                    print 'germline seq:', repr(germline_seq)
                    print 'aborting'
                    return
                if add_divergence == 'shazam':
                    divergence_filepath = output_filepath[:-6] + '_divergence_results.txt'
                    p = Popen(['Rscript', 'invoke_shazam_calcObservedMutations.R', obs_seq, germline_seq], stdout=PIPE, stderr=PIPE)
                    out_err = p.communicate()
                    values = out_err[0].split('\n')
                    if values[3] == 'NA':
                        div_stats[seq_id] = {'replacement_mut_count':'NA', 'synonymous_mut_count':'NA', 'total_mut_count':'NA', 'replacement_mut_freq':'NA', 'synonymous_mut_freq':'NA', 'total_mut_freq':'NA'}
                    else:
                        denominator = float(values[3])
                        div_stats[seq_id] = {'replacement_mut_count':values[0], 'synonymous_mut_count':values[1], 'total_mut_count':values[2], 'replacement_mut_freq':int(values[0])/denominator, 'synonymous_mut_freq':int(values[1])/denominator, 'total_mut_freq':(int(values[0])+int(values[1])) / denominator}
                elif add_divergence == 'python':
                    denominator = 0.
                    mut_count = 0
                    #get divergence relative to D-masked germline
                    for i, j in itertools.izip(obs_seq, germline_seq):
                        if i!='.' and j!='.' and i!='-' and j!='-' and i!='N' and j!='N':
                            denominator += 1
                            if i != j:
                                mut_count += 1
                    if denominator == 0:
                        div_stats[seq_id] = {'mut_count':'NA', 'mut_freq':'NA', 'mut_count_non_D_mask':'NA', 'mut_freq_non_D_mask':'NA'}
                    else:
                        div_stats[seq_id] = {'mut_count':mut_count, 'mut_freq':mut_count/denominator, 'mut_count_non_D_mask':'NA', 'mut_freq_non_D_mask':'NA'}
                    #get divergence relative to non-D-masked germline
                    denominator = 0.
                    mut_count = 0
                    for i, j in itertools.izip(obs_seq, non_D_masked_germline_seq):
                        if i!='.' and j!='.' and i!='-' and j!='-' and i!='N' and j!='N':
                            denominator += 1
                            if i != j:
                                mut_count += 1
                    if denominator > 0:
                        div_stats[seq_id]['mut_count_non_D_mask'] = mut_count
                        div_stats[seq_id]['mut_freq_non_D_mask'] = mut_count/denominator

        if add_selection:
            if 'dmask' not in add_germline:
                print '"add_germline" parameter must equal "dmask". Aborting.'
                return
            selection_filepath = output_filepath[:-6] + '_baseline_results.txt'
            print ''
            print '############# running baseline #############'
            sys.stdout.flush()
            p = Popen(['Rscript', 'invoke_baseline.R', changeo_db_filepath, selection_filepath])
            out, err = p.communicate()
            print 'output for baseline:'
            print out
            print err
            sys.stdout.flush()
            #get selection stats for each seq
            filein = open(selection_filepath, "r")
            filein.readline()
            sel_stats = {}
            for i in filein:
                line = i[:-1].split('\t')
                seq_id = line[0]
                region = line[1]
                sigma = line[2]
                ci_low = line[3]
                ci_high = line[4]
                p_val = line[5]
                try:
                    sel_stats[seq_id][region + '_sigma'] = sigma
                    sel_stats[seq_id][region + '_ci_low'] = ci_low
                    sel_stats[seq_id][region + '_ci_high'] = ci_high
                    sel_stats[seq_id][region + '_p_val'] = p_val
                except KeyError:
                    sel_stats[seq_id] = {'CDR_sigma':'NA', 'CDR_ci_low':'NA', 'CDR_ci_high':'NA', 'CDR_p_val':'NA', 'FWR_sigma':'NA', 'FWR_ci_low':'NA', 'FWR_ci_high':'NA', 'FWR_p_val':'NA'}
                    sel_stats[seq_id][region + '_sigma'] = sigma
                    sel_stats[seq_id][region + '_ci_low'] = ci_low
                    sel_stats[seq_id][region + '_ci_high'] = ci_high
                    sel_stats[seq_id][region + '_p_val'] = p_val
            filein.close()
            subprocess.call(['rm', selection_filepath])
        if translate_VDJ:
            print ''
            print '############# translating VDJ sequence #############'
            sys.stdout.flush()
            filein = open(changeo_db_filepath, "r")
            filein.readline()
            trans_stats = {}
            for i in filein:
                line = i[:-2].split('\t')
                seq_id = line[0]
                #find the coding frame for the seq
                imgt_nuc_seq = line[11]
                num_dots = 0
                for j in imgt_nuc_seq:
                    if j != '.':
                        break
                    num_dots += 1
                codon_remainder = num_dots % 3
                if codon_remainder == 0:
                    frame = 0
                elif codon_remainder == 1:
                    frame = 2
                else:
                    frame = 1
                nuc_seq = line[10].replace('-', '')#remove any dashes
                nuc_seq = nuc_seq[frame:]#start the seq in the right frame
                aa_seq = Seq.translate(nuc_seq)
                trans_stats[seq_id] = aa_seq
            filein.close()

        #write change-o database as fasta file
        fileout = open(output_filepath, "w")
        filein = open(changeo_db_filepath, "r")
        attribute_names = filein.readline()[:-2].split('\t')[2:]
        for i in filein:
            line = i[:-2].split('\t')
            seq_id = line[0]
            seq = line[1]
            if remove_seqs_with_stop:
                if line[4] == 'T':
                    continue
            header = ">" + seq_id
            for count, j in enumerate(line[2:]):
                header += '|%s=%s' % (attribute_names[count], j)
            if add_divergence:
                if add_divergence == 'shazam':
                    header += '|replacement_mut_count=%s|synonymous_mut_count=%s|total_mut_count=%s|replacement_mut_freq=%s|synonymous_mut_freq=%s|total_mut_freq=%s' % (div_stats[seq_id]['replacement_mut_count'], div_stats[seq_id]['synonymous_mut_count'], div_stats[seq_id]['total_mut_count'], div_stats[seq_id]['replacement_mut_freq'], div_stats[seq_id]['synonymous_mut_freq'], div_stats[seq_id]['total_mut_freq'])
                elif add_divergence == 'python':
                    header += '|total_mut_count=%s|total_mut_freq=%s|total_mut_count_non_D_masked=%s|total_mut_freq_non_D_masked=%s' % (div_stats[seq_id]['mut_count'], div_stats[seq_id]['mut_freq'], div_stats[seq_id]['mut_count_non_D_mask'], div_stats[seq_id]['mut_freq_non_D_mask'])
            if add_selection:
                header += '|CDR_baseline_selection=%s|CDR_baseline_selection_CI_upper=%s|CDR_baseline_selection_CI_lower=%s|CDR_baseline_selection_pval=%s|FWR_baseline_selection=%s|FWR_baseline_selection_CI_upper=%s|FWR_baseline_selection_CI_lower=%s|FWR_baseline_selection_pval=%s' % (sel_stats[seq_id]['CDR_sigma'], sel_stats[seq_id]['CDR_ci_high'], sel_stats[seq_id]['CDR_ci_low'], sel_stats[seq_id]['CDR_p_val'], sel_stats[seq_id]['FWR_sigma'], sel_stats[seq_id]['FWR_ci_high'], sel_stats[seq_id]['FWR_ci_low'], sel_stats[seq_id]['FWR_p_val'])
            if translate_VDJ:
                header += '|VDJ_AA_seq=%s' % trans_stats[seq_id]
            fileout.write('%s\n%s\n' % (header, seq))
        filein.close()
        fileout.close()
        subprocess.call(['rm', changeo_db_filepath])
        return

    def make_seqIDs_sequential_integers(self, start_int=1):
        """This method changes the sequence IDs so that they are sequential integers (i.e. 1, 2, 3, etc.). self.data['id'] will thus be changed by implementing this method.
        start_int - This should be an integer. It gives the starting point for the naming of the sequence IDs.
        """
        id = start_int
        for i in xrange(len(self.data)):
            self.data[i]['id'] = str(id)
            id += 1
        return

    def get_most_abundant_seq(self, output_filepath=None, translate_by_ref_seq=False, temp_dirpath=None, append_to_output=False, return_seq_index=False, path_to_needle='needle'):
        """
        This method finds the sequence that is the most abundant in the set. If there are ties, then it arbitrarily returns the sequence that is first (i.e. has lowest index) in the the 'data' attribute.
        output_filepath - If defined (default, None) then this is the path to the output file that will be written for the most abundant seq. Will be in fasta format.
        translate_by_ref_seq - If defined (default is False) this should be the path to a fasta formatted reference sequence for the seqs in the data. This reference seq needs to be annotated such that the coding frame is known. Meaning, it is known where to begin translating the ref seq. The annotation for the coding frame information must be in the header of the ref seq and should look like:'...|coding_frame_start=[X]|...', where 'X' gives the coding frame start position. The most abundant seq is then aligned to this ref seq so that the coding frame can be mapped from the ref seq to the most abundant seq. The alignment is done using 'needle' from the EMBOSS package. This needs to be in the PATH. The resulting translation is recorded in the 'output_filepath', so this should be defined too. If 'translate_by_ref_seq' equals False, then no translation occurs.
        temp_dirpath - Gives the path to a directory for which temp data can be written. If None (default), then uses the master directory to 'output_filepath', if there is not output_filepath given then uses the current working directory.
        append_to_output - If True (default, False) then this will cause the output file to be appended to, rather than over-written. This is only considered if 'output_filepath' is defined.
        return_seq_index - Boolean. If True (default, False), this will instruct the method to return the index of the most abundant seq, rather than the seq itself.
        path_to_needle - String. This give the path to the 'needle' executable. The default is simply 'needle', which means that it is already in your PATH.
        """
        if not self.count_attribute_name:
            print "No count information in the 'sequence_sample' object. Aborting."
            return
        if temp_dirpath:
            if temp_dirpath[-1] != '/':
                temp_dirpath += '/'
            if not os.path.exists(temp_dirpath):
                os.makedirs(temp_dirpath)
        else:
            if output_filepath:
                temp_dirpath = os.path.dirname(output_filepath) + '/'
            else:
                temp_dirpath = os.getcwd() + "/"
        random_suffix = str(random.random())

        #this finds the most abundant seq in the sample
        most_abun_seq_index = self.counts.index(max(self.counts))
        most_abun_seq = self.data[most_abun_seq_index]['seq']

        #if translating the seq by a reference
        if translate_by_ref_seq:
            filein = open(translate_by_ref_seq, "r")
            header = filein.readline()[1:-1].split('|')
            for i in header:
                attribute = i.split('=')
                if attribute[0] == 'coding_frame_start':
                    coding_frame_start = int(attribute[1])
            ref_seq = filein.readline()[:-1]
            filein.close()
            alignment_filepath = "%stemp_alignment_%s" % (temp_dirpath, random_suffix)
            p = Popen([path_to_needle, '-asequence', 'asis:'+ref_seq, '-bsequence', 'asis:'+most_abun_seq, '-gapopen', '10.0', '-gapextend', '0.5', '-outfile', alignment_filepath, '-brief', 'Y'], stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()
            print out
            print err
            filein = open(alignment_filepath, "r")
            ref_seq_alignment = ''
            most_abun_seq_alignment = ''
            first_seq = True
            for i in filein:
                if i[:4] == 'asis':
                    seq = i[:-1].split()[2]
                    if first_seq:
                        ref_seq_alignment += seq
                        first_seq = False
                    else:
                        most_abun_seq_alignment += seq
                        first_seq = True
            filein.close()
            subprocess.call(['rm', alignment_filepath])
            ref_seq_pos = 0
            most_abun_seq_pos = 0
            for i, j in itertools.izip(ref_seq_alignment, most_abun_seq_alignment):
                if i != '-':
                    ref_seq_pos += 1
                if j != '-':
                    most_abun_seq_pos += 1
                if ref_seq_pos == coding_frame_start:
                    most_abun_seq_coding_frame_start = most_abun_seq_pos
                    break
            most_abun_aa_seq = Seq.translate(most_abun_seq[most_abun_seq_coding_frame_start-1:])
        else:
            most_abun_aa_seq = None

        #if writing to disk
        if output_filepath:
            if append_to_output:
                fileout = open(output_filepath, "a")
            else:
                fileout = open(output_filepath, "w")
            fileout.write('>%s|%s=%s|amino_acid_seq=%s|coding_frame_start=%s\n%s\n' % (self.data[most_abun_seq_index]['id'], self.count_attribute_name, self.data[most_abun_seq_index]['count'], most_abun_aa_seq, most_abun_seq_coding_frame_start, most_abun_seq))
            fileout.close()

        if return_seq_index:
            return most_abun_seq_index
        else:
            return most_abun_seq

    def remove_seq_entries(self, seq_indicators, seq_ids=False):
        """
        This method will remove sequence entries from the data. It will adjust all of the data attributes to reflect this removal.
        seq_indicators - Gives a list of the indicator values that allows one to identify the entries to remove. Could be a list of indices (ints, indexed at zero), or a list of sequence IDs.
        seq_ids - Boolean (default False). If True then this indicates that the seq_indicators are a list of sequence IDs, otherwise it is assumed that they are a list of ideces.
        """
        if seq_ids:
            seq_indices = []
            for index, i in enumerate(self.data):
                if i['id'] in seq_indicators:
                    seq_indices.append(index)
        else:
            seq_indices = seq_indicators
        for i in sorted(seq_indices, reverse=True):
            count = self.data[i]['count']
            del self.counts[i]
            self.total -= count
            del self.data[i]
        return

    def translate_each_seq(self, coding_frame_start=1, ref_seq=None, temp_dirpath=None, remove_seqs_with_stop_codons=False, path_to_needle='needle'):
        """
        This method will translate each of the nucleotide seqs in the data. A challenge with this is to know where this coding frame is in the sequence. The coding frame can start at either the 1st, 2nd, or 3rd nucleotide in a given seq. The resulting translated amino acid seqs will be stored in memory as a string in self.data['other']['amino_acid_seq'].
        coding_frame_start - This gives the coding frame of each of the sequences. Acceptable values are:
            [1,2,3] - (read as 1, 2, or 3). This gives the location of the start of the coding frame in the sequence. For example if this is 1, then that means that the coding frame starts at the 1st nucleotide in the seq.
            'relative_to_reference' - This means that the start of the coding frame should be determined by mapping each sequence to a reference sequence, where the coding frame is known. This is done by using 'needle' from the emboss package to globally align each seq (pairwise) to the reference seq. Using this alignment, the coding frame start location from the reference is then mapped onto a given seq in the data. If this approach is used, then the 'ref_seq' parameter must be defined.
        ref_seq - If defined (default, None), then this will be the path to the fasta file that contains the reference sequence for each of the seqs in the data. This will be used as a reference for translation.
        temp_dirpath - Gives the path to a directory for which temp data can be written. If None (default), then uses the current working directory.
        remove_seqs_with_stop_codons - If True (default, False) this will remove seq entries in the data that, when translated, have stop codons.
        path_to_needle - String. This give the path to the 'needle' executable. The default is simply 'needle', which means that it is already in your PATH.
        """
        if coding_frame_start == 'relative_to_reference' and ref_seq == None:
            print "if coding frame is 'relative_to_reference' then need to provide the path to the reference sequence via the 'ref_seq' parameter. Aborting."
            return
        if temp_dirpath:
            if temp_dirpath[-1] != '/':
                temp_dirpath += '/'
            if not os.path.exists(temp_dirpath):
                os.makedirs(temp_dirpath)
        else:
            temp_dirpath = os.getcwd() + "/"
        random_suffix = str(random.random())
        if ref_seq:
            filein = open(ref_seq, "r")
            header = filein.readline()[1:-1].split('|')
            for i in header:
                attribute = i.split('=')
                if attribute[0] == 'coding_frame_start':
                    ref_frame_start = int(attribute[1])
                if attribute[0] == 'amino_acid_seq':
                    ref_aa_seq = attribute[1]
            ref_seq = filein.readline()[:-1]
            filein.close()

        if remove_seqs_with_stop_codons:
            seqs_to_remove = []
        for seq_index, i in enumerate(self.data):

            #if we are to find the coding frame relative to a reference seq
            if ref_seq:
                alignment_filepath = '%stemp_alignment_%s' % (temp_dirpath, random_suffix)
                p = Popen([path_to_needle, '-asequence', 'asis:'+ref_seq, '-bsequence', 'asis:'+i['seq'], '-gapopen', '10.0', '-gapextend', '0.5', '-outfile', alignment_filepath, '-brief', 'Y'], stdout=PIPE, stderr=PIPE)
                out, err = p.communicate()
                #print out
                #print err
                filein = open(alignment_filepath, "r")
                ref_seq_alignment = ''
                seq_alignment = ''
                first_seq = True
                for j in filein:
                    if j[:4] == 'asis':
                        seq = j[:-1].split()[2]
                        if first_seq:
                            ref_seq_alignment += seq
                            first_seq = False
                        else:
                            seq_alignment += seq
                            first_seq = True
                filein.close()
                subprocess.call(['rm', alignment_filepath])
                ref_seq_pos = 0
                seq_pos = 0
                for ref, seq in itertools.izip(ref_seq_alignment, seq_alignment):
                    if ref != '-':
                        ref_seq_pos += 1
                    if seq != '-':
                        seq_pos += 1
                    if (ref_seq_pos-ref_frame_start) % 3 == 0 and seq != '-':
                        coding_frame_start = seq_pos
                        break

            aa_seq = Seq.translate(i['seq'][coding_frame_start-1:])
            #if there is a stop codon in the seq
            if remove_seqs_with_stop_codons:
                if '*' in aa_seq:
                    seqs_to_remove.append(seq_index)
            self.data[seq_index]['other']['amino_acid_seq'] = aa_seq[:]
            self.data[seq_index]['other']['coding_frame_start'] = str(coding_frame_start)

        if remove_seqs_with_stop_codons:
            self.remove_seq_entries(seqs_to_remove)
        return

    @staticmethod
    def get_observed_syn_and_nonSyn_muts_equal_weights(self, ref_codon, query_codon, mut_positions):
        """
        This is an internal method that will return the observed number of non-syn and syn mutations when given a reference and query codon.
        ref_codon - 3 character string giving the reference (i.e. starting) codon.
        query_codon - 3 character string giving the query (i.e. ending, or observed) codon.
        mut_positions - List of ints (max length of 3) giving the positions of each mutation b/t the reference and query codons.
        """
        num_muts = len(mut_positions)
        Ns = [0 for i in xrange(num_muts)]
        Ss = [0 for i in xrange(num_muts)]
        for i in xrange(num_muts):
            intermediate_codon = ref_codon[:mut_positions[i]] + query_codon[mut_positions[i]] + ref_codon[mut_positions[i]+1:]
            if geneticCode[ref_codon] == geneticCode[intermediate_codon]:
                Ss[i] = 1
            else:
                Ns[i] = 1
            if num_muts > 1:
                N, S = self.get_observed_syn_and_nonSyn_muts_equal_weights(self, intermediate_codon, query_codon, mut_positions[:i] + mut_positions[i+1:])
                Ns[i] += N
                Ss[i] += S
        N = sum([i/float(len(Ns)) for i in Ns])
        S = sum([i/float(len(Ss)) for i in Ss])
        return N, S

    @staticmethod
    def get_observed_syn_and_nonSyn_muts_max_synMuts(self, ref_codon, query_codon, mut_positions):
        """
        This is an internaml method that will return the number of non-synonymous and synonymous mutations between a reference codon and a query codon. Specifically, it returns the number non-syn and syn mutations for the path that has the maximum number of synonymous muts. This way we get a more conservative estimate of the number of non-syn muts, as these mutations tend to be less likely in many contexts.
        ref_codon - 3 character string giving the reference (i.e. starting) codon.
        query_codon - 3 character string giving the query (i.e. ending, or observed) codon.
        mut_positions - List of ints (max length of 3) giving the positions of each mutation b/t the reference and query codons.
        """
        num_syn_muts = []
        for i in itertools.permutations(mut_positions):
            intermediate_codon = ref_codon[:]
            num_syn_muts.append(0)
            for j in i:
                new_codon = intermediate_codon[:j] + query_codon[j] + intermediate_codon[j+1:]
                if geneticCode[intermediate_codon] == geneticCode[new_codon]:
                    num_syn_muts[-1] += 1
                intermediate_codon = new_codon[:]
        S = max(num_syn_muts)
        N = len(mut_positions) - S
        return N, S

    def get_dN_dS_and_divergence(self, ref_seq, method='counting', temp_dirpath=None, method_counting='equal_weights', path_to_needle='needle'):
        """
        This method will get both the non-synonymous and synonymous divergence, as well as dN/dS, relative to some reference sequence. It does this for each of the seqs in the data. It will store the resulting values in memory as self.data['other']['divergence_non_syn'] and self.data['other']['divergence_syn']
        ref_seq - This gives the path to the fasta file that has the reference sequence that each of the sequences in the data will be compared to to get the divergence values. This fasta file should contain one sequence with no line breaks (two lines total in the whole file). It should have an attribute in the header of the sq that looks like this: '>blah_blah|amino_acid_seq=[X]|blah_blah_blah...', where [X] is a string of amino acid characters.
        method - This gives the method for which divergence is calculated. Acceptable values are:
            'counting' - This means that non-syn. and syn. divergence is calculated by first codon aligning the query seq to the ref seq and then counting all the non-syn and syn changes between the two seqs, and then normalizing by all the possible non-syn and syn changes that could have taken place between the two seqs. A tutorial that describes how this was implemented can be found here: http://bioinformatics.cvr.ac.uk/blog/calculating-dnds-for-ngs-datasets/ .
        temp_dirpath - Gives the path to a directory for which temp data can be written. If None (default), then uses the current working directory.
        method_counting - This gives the method for counting the number of non-synonymous and synonymous mutations between a ref and query seq. This paramter is only considered if method='counting'. Acceptable values are:
            'equal_weights' - Default. This means that the different possible mutational paths between two codons are weigted equally. I.e. non-syn and syn mutations are treated as equally likely.
            'max_synonymous' - This means that the mutational path with the maximum amount of syn muts will be chosen for each reference/query codon pair.
        path_to_needle - String. This give the path to the 'needle' executable. The default is simply 'needle', which means that it is already in your PATH.
        """
        try:
            self.data[0]['other']['amino_acid_seq']
        except KeyError:
            print "The sequences in the data must be translated before running this method. Run the method 'translate_each_seq' before running this method. Aborting."
            return
        if temp_dirpath:
            if temp_dirpath[-1] != '/':
                temp_dirpath += '/'
            if not os.path.exists(temp_dirpath):
                os.makedirs(temp_dirpath)
        else:
            temp_dirpath = os.getcwd() + "/"
        random_suffix = str(random.random())

        filein = open(ref_seq, "r")
        ref_aa_seq = []
        ref_seq = []
        ref_coding_frame = []
        ref_seq_codons_list = []
        for line in filein:
            if line[0] == '>':
                header = line[1:-1].split('|')
                for i in header:
                    attribute = i.split('=')
                    if attribute[0] == 'amino_acid_seq':
                        ref_aa_seq.append(attribute[1])
                    elif attribute[0] == 'coding_frame_start':
                        ref_coding_frame.append(int(attribute[1]))
            else:
                seq = line[:-1]
                if not seq.istitle():
                    seq = seq.upper()
                ref_seq.append(seq)
                ref_seq_codons_list.append([seq[i:i+3] for i in range(ref_coding_frame[-1]-1, len(seq), 3)])
        filein.close()
        if len(ref_seq) == 1:
            ref_aa_seq = ref_aa_seq[0]
            ref_seq = ref_seq[0]
            ref_coding_frame = ref_coding_frame[0]
            ref_seq_codons = ref_seq_codons_list[0]

        #this is a dic that has the number of non-syn and syn sites for each possible codon.
        nonSyn_syn_sites_dic = get_nonSyn_and_syn_sites_foreachCodon.get_num_nonSyn_and_syn_sites_foreach_codon()
        for seq_index, i in enumerate(self.data):
            query_seq_codons = [i['seq'][j:j+3] for j in range(int(i['other']['coding_frame_start'])-1, len(i['seq']), 3)]
            #make query seq upper case, if not
            if not query_seq_codons[0].istitle():
                query_seq_codons = [j.upper() for j in query_seq_codons]

            #if there are multiple ref seqs then find the one that is closest to query seq and use that one
            if isinstance(ref_seq, list):
                ref_aa_seq_alignments = []
                query_aa_seq_alignments = []
                idents = []
                for j in xrange(len(ref_seq)):
                    alignment_filepath = '%stemp_alignment_%s_%s' % (temp_dirpath, j, random_suffix)
                    p = Popen([path_to_needle, '-asequence', 'asis:'+ref_aa_seq[j], '-bsequence', 'asis:'+i['other']['amino_acid_seq'], '-gapopen', '10.0', '-gapextend', '0.5', '-outfile', alignment_filepath, '-brief', 'Y', '-sprotein1', 'Y', '-sprotein2', 'Y'], stdout=PIPE, stderr=PIPE)
                    out, err = p.communicate()
                    print out
                    print err
                    filein = open(alignment_filepath, "r")
                    ref_aa_seq_alignments.append('')
                    query_aa_seq_alignments.append('')
                    first_seq = True
                    for line in filein:
                        if line[:4] == 'asis':
                            seq = line[:-1].split()[2]
                            if first_seq:
                                ref_aa_seq_alignments[-1] += seq
                                first_seq = False
                            else:
                                query_aa_seq_alignments[-1] += seq
                                first_seq = True
                        elif line[:11] == '# Identity:':
                            ident = float(line.split('(')[-1].split(')')[0][:-1])
                            idents.append(ident)
                    filein.close()
                    subprocess.call(['rm', alignment_filepath])
                max_index = idents.index(max(idents))
                ref_aa_seq_alignment = ref_aa_seq_alignments[max_index]
                query_aa_seq_alignment = query_aa_seq_alignments[max_index]
                ref_seq_codons = ref_seq_codons_list[max_index][:]
            else:
                alignment_filepath = '%stemp_alignment_%s' % (temp_dirpath, random_suffix)
                p = Popen([path_to_needle, '-asequence', 'asis:'+ref_aa_seq, '-bsequence', 'asis:'+i['other']['amino_acid_seq'], '-gapopen', '10.0', '-gapextend', '0.5', '-outfile', alignment_filepath, '-brief', 'Y', '-sprotein1', 'Y', '-sprotein2', 'Y'], stdout=PIPE, stderr=PIPE)
                out, err = p.communicate()
                print out
                print err
                filein = open(alignment_filepath, "r")
                ref_aa_seq_alignment = ''
                query_aa_seq_alignment = ''
                first_seq = True
                for j in filein:
                    if j[:4] == 'asis':
                        seq = j[:-1].split()[2]
                        if first_seq:
                            ref_aa_seq_alignment += seq
                            first_seq = False
                        else:
                            query_aa_seq_alignment += seq
                            first_seq = True
                filein.close()
                subprocess.call(['rm', alignment_filepath])

            ref_codon_index = -1
            query_codon_index = -1
            N_num_sites = 0.
            S_num_sites = 0.
            N_observed = 0.
            S_observed = 0.
            for ref_aa, query_aa in itertools.izip(ref_aa_seq_alignment, query_aa_seq_alignment):
                if ref_aa != '-':
                    ref_codon_index += 1
                if query_aa != '-':
                    query_codon_index += 1
                if ref_aa != '-' and query_aa != '-':
                    ref_codon = ref_seq_codons[ref_codon_index]
                    query_codon = query_seq_codons[query_codon_index]
                    N_num_sites += nonSyn_syn_sites_dic[ref_codon]['non_syn']
                    S_num_sites += nonSyn_syn_sites_dic[ref_codon]['syn']
                    if ref_codon != query_codon:
                        mut_positions = []
                        index = 0
                        for ref_nuc, query_nuc in itertools.izip(ref_codon, query_codon):
                            if ref_nuc != query_nuc:
                                mut_positions.append(index)
                            index += 1
                        if method_counting == 'equal_weights':
                            N, S = self.get_observed_syn_and_nonSyn_muts_equal_weights(self, ref_codon, query_codon, mut_positions)
                        elif method_counting == 'max_synonymous':
                            N, S = self.get_observed_syn_and_nonSyn_muts_max_synMuts(self, ref_codon, query_codon, mut_positions)
                        N_observed += N
                        S_observed += S
            divergence_non_syn = N_observed / N_num_sites
            divergence_syn = S_observed / S_num_sites
            dN = (-3./4)*math.log(1 - (4*divergence_non_syn/3.))
            dS = (-3./4)*math.log(1 - (4*divergence_syn/3.))
            if dS == 0.:
                dN_dS = 'NA'
            else:
                dN_dS = dN / dS
            self.data[seq_index]['other']['divergence_non_syn'] = str(divergence_non_syn)
            self.data[seq_index]['other']['divergence_syn'] = str(divergence_syn)
            self.data[seq_index]['other']['dN_dS'] = str(dN_dS)
        return

    def add_string_to_each_id(self, string_to_add, add_to_start_or_end='start'):
        """
        This method adds the provided string to each of the sequence IDs. This info is then saved into memory.
        string_to_add - This should be a string, and will be added to each of the seq IDs.
        add_to_start_or_end - This tells if the provided string will be added to the beginning or the end of each of the IDs. Acceptable values are:
            'start' - Default. String will be added to the beggining.
            'end' - String will be added to the end.
        """
        for i in xrange(len(self.data)):
            if add_to_start_or_end == 'start':
                self.data[i]['id'] = string_to_add + self.data[i]['id']
            elif add_to_start_or_end == 'end':
                self.data[i]['id'] += string_to_add
            else:
                'parameter "add_to_start_or_end" not set correctly:', add_to_start_or_end
                return
        return

    def get_seq_data(self, seq_id):
        """
        When provided the 'ID' of a seq entry, this will return the sequence data for that entry.
        """
        for i in self.data:
            if i['id'] == seq_id:
                return i
        print "Didn't find the seq yo"
        return

    def add_freq_attribute(self, freq_attribute_name='freq'):
        """
        This method will add an attribute to self.data[some_index]['other'][freq_attribute_name] that gives the frequency for each of the seqs. Frequency is calculated as self.data[some_index]['count'] / self.total.
        """
        for i in xrange(len(self.data)):
            self.data[i]['other'][freq_attribute_name] = self.data[i]['count'] / self.total
        return

    def add_boolean_attribute_to_seqs(self, attribute_name, seq_ids_that_are_True, overwrite_existing=False):
        """
        This method will add a boolean attribute (i.e. the values of the attribute are either 'True' or 'False') to self.data[some_index]['other'].
        attribute_name - This give the name for the attribute that will be added (ex: 'is_outlier').
        seq_ids_that_are_True - List of strings. This gives the sequence ID for each of the sequences that will have a 'True' value for the new boolean attribute.
        overwrite_existing - If True (default, False), this will overwrite the data if there is already and existing attribute of 'attribute_name' in self.data[some_index]['other']. If False, then it will check for this and abort if 'attribute_name' already exists in the data.
        """
        if attribute_name in self.data[0]['other'] and not overwrite_existing:
            print attribute_name, 'already exists in sample data. Aborting.'
            return
        for index in xrange(len(self.data)):
            if self.data[index]['id'] in seq_ids_that_are_True:
                self.data[index]['other'][attribute_name] = 'True'
            else:
                self.data[index]['other'][attribute_name] = 'False'
        return



#############################################
# Below are scripts that are used for running
# jobs on the computational cluster. They
# should only be used by the methods above.
#############################################

def pi_calculator_compCluster(output_dirpath, method, path_to_needle):
    """
    This script is to only be used by the 'calc_diversity_pi_compCluster' method above. We needed that method be able to submit an array job on the computational cluster, and you can't run a method from within a class from the command line, so we made this script to accompany the method. This script calculates the sum of genetic differences for a portion of the sequences in a sample. It then writes this sum to an output file to be used later. The SGE_TASK_ID environmental variable gives which index of the sequences will be compared to all the seqs that follow it.
    'seqs' - This is a long string, where each unique sequence is separated by a ','
    'counts' - This is another long string of the counts that correspond to the sequences in 'seqs', and formatted the same way.
    'output_dirpath' - This is the path to directory that all the output files (sub sums) will be written.
    'method' - This gives the method that will be used to get the genetic distance from the pairs of seq\
s. The exceptable values for this are:
            'needle' - This means the needleman-wunsch global alignment algorithm will be used in a program called 'needle' in the EMBOSS package
            'pairwise2' - This means the 'pairwise2' Biopython package will be used, which also implements the needleman-wunsch global alignment algorithm. Might be faster than needle though.
    'path_to_needle' - This gives the path to the 'needle' program, which is the global aligner from EMBOSS that implements the needleman-wunsch algorithm. This can also just be 'needle', which would mean that the program's directory is in the PATH.
    """
    sge_task_id = int(os.environ['SGE_TASK_ID'])

    ##############################
    # if using the needle method #
    ##############################
    if method == 'needle':

        job_dirpath = '%s%s/' % (output_dirpath, sge_task_id)

        #1st make sure that the necessary input dirpath exists.
        #If it doesn't then it could be because the server is
        #slow, and hasn't finished writing yet
        good_to_go = False
        while good_to_go == False:
            try:
                a = os.listdir(job_dirpath + 'a_seqs')
                b = os.listdir(job_dirpath + 'b_seqs')
                good_to_go = True
            except OSError:
                print "Temp sequence set hasn't written yet. Will wait and try again:"
                print job_dirpath
                sys.stdout.flush()
                time.sleep(1)

        #cycle through the sequence sets that will be
        #submitted to needle
        sub_sum_of_distances = 0
        for a_seq_set, b_seq_set in itertools.izip(os.listdir(job_dirpath + 'a_seqs'), os.listdir(job_dirpath + 'b_seqs')):
            #setup the needle call. The first string gives the location
            #of needle. Change this as necessary
            p = Popen([path_to_needle, '-asequence', job_dirpath+'a_seqs/'+a_seq_set, '-bsequence', job_dirpath+'b_seqs/'+b_seq_set, '-gapopen', '10.0', '-gapextend', '0.5', '-brief', 'Y', '-stdout', '-auto'], stdout=PIPE, stderr=PIPE)
            out_err = p.communicate()
            #parse needle output and get seq counts too
            first_entry = True
            for i in out_err[0].split('\n'):
                if i[:5] == '# 1: ' and first_entry:
                    a_seq_count = int(i.split('_acount=')[1])
                    first_entry = False
                if i[:5] == '# 2: ':
                    b_seq_count = int(i.split('_bcount=')[1])
                if i[:11] == '# Identity:':
                    match = re.search('\(\s*([0-9]+\.[0-9])', i)
                    percent_dist = 100 - float(match.group(1))
                    sub_sum_of_distances += a_seq_count * b_seq_count * percent_dist
            subprocess.call(['rm', job_dirpath+'a_seqs/'+a_seq_set, job_dirpath+'b_seqs/'+b_seq_set])

    #################################
    # if using the pairwise2 method #
    #################################
    elif method == 'pairwise2':

        #read seqs and counts into memory
        filein = open(output_dirpath+'seqs_and_counts', "r")
        seqs = []
        counts = []
        for i in filein:
            line = i[:-1].split('\t')
            seqs.append(line[0])
            counts.append(int(line[1]))
        filein.close()

        query_seq_index = sge_task_id - 1
        query_seq = seqs[query_seq_index]
        sub_matrix = MatrixInfo.genetic
        sub_sum_of_distances = 0
        index = query_seq_index + 1
        for i in seqs[query_seq_index+1:]:
            alignment = pairwise2.align.globalds(query_seq, i, sub_matrix, -10, -0.5, one_alignment_only=True)
            length = 0
            percent_dist = 0.
            for seqa, seqb in itertools.izip(alignment[0][0], alignment[0][1]):
                length += 1
                if seqa == seqb:
                    percent_dist += 1
            percent_dist = (percent_dist/length) * 100
            sub_sum_of_distances += counts[index] * percent_dist
            index += 1
        
        sub_sum_of_distances = counts[query_seq_index] * sub_sum_of_distances

    #write the sub sum of distances to an output file
    output_filepath = '%s%s_subsum' % (output_dirpath, sge_task_id)
    fileout = open(output_filepath, "w")
    fileout.write(str(sub_sum_of_distances))
    fileout.close()
    return

if __name__ == '__main__':
    if sys.argv[1] == 'pi_calculator_compCluster':
        pi_calculator_compCluster(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == 'test':
        #here one can test out the different methods
        sample = sequence_sample('/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/QCed_data_clean_final/abr/1/1144.fastq', 'DUPCOUNT')
        sample.get_immune_reads_with_changeo(output_filepath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/seq_data/changeo_annotated_fasta_files/1/1144.fasta', path_to_changeo='/netapp/home/nstrauli/tools_c/changeo-0.3.9', path_to_igblast='/netapp/home/nstrauli/tools_c/ncbi-igblast-1.8.0', temp_dirpath='/hernandez/mandrill/users/nstrauli/data/abr_hiv_coevo/temp_stuff/', imgt_ref_seqs=None, add_germline='dmask', add_selection=True, use_comp_cluster=True, add_divergence=True, overwrite_output=False)
    else:
        print 'If running sequence_sample_class.py as a script, then must specify an appropriate sub-script to run.'
