#!/usr/bin/python
#$ -S /usr/bin/python
#$ -e error
#$ -o out_5
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=200G
#$ -l h_rt=336:00:00

import sys
sys.path.insert(0, './')
import os
from immune_sample_class import immune_sample
import itertools
import re
from subprocess import Popen, PIPE
import subprocess
import tempfile
import time
import random

class immune_time_series(object):
    """
    A class for a collection of immune sequence data time points.
    
    Input:
    The input for this class should be a directory that contains only fasta formatted files (files named README[.txt] will be ignored). Further the fasta files' headers should be formatted in a specific way. The header for each sequence should start with an ID (could be anything, ideally unique) and then should be followed by a series of attributes. Each attribute should be delimited by a '|'. Each attribute should have a name followed by and '=', followed by the value of the attribute. An example of a well formatted header: ">seq_id|vgene=IGHV1|dgene=NA|jgene=IGHJ3|count=65\n". Here the attributes are 'vgene', 'dgene', 'jgene', and 'count'. The first entry in the header is assumed to be the sequence unique ID.
    The fasta files should also have a specific naming scheme. Each of the filenames should end with a numeric value that indicates the time point, such that if these numeric time point indicators where sorted (numerically), the files would sort in chronological order. This numeric time point indicator in the file names should be delimited by an '_'. For example 'day_56.fasta' or 'month_34.fasta' are examples of properly formatted filenames. '110.fasta' would also work, however 'month_12_day_21.fasta' would not be advisable. This particular case would not crash the code, but would result in the time points not being sorted properly (because the time point is encoded by two numerics (month and day) as opposed to one).

    Attributes:
    dirpath - The path to the directory that contains the sequence data (in fasta format) for each time point in the time series.
    count_attribute_name - this gives the name that corresponds to the count of each clone/sequence entity.
    vgene_name - This gives the name of the attribute that gives the V gene name.                          
    dgene_name - This gives the name of the attribute that gives the D gene name.                          
    jgene_name - This gives the name of the attribute that gives the J gene name.                          
    cdr3_name - This gives the name of the attribute that gives the CDR3 sequence.
    sample_filepaths - A list of each of the time points filepaths, in chronological order (ascending).
    timepoints - a list of numerics for each of the timepoints (in ascending order)
    replace_slash_with - This gives the character that is used to replace a slash in an attribute value. Slashes can cause problems if attribute values are used in the creation of output files in some of the methods below. The default is None, which means that the attribute values are kept as is.
    avoid_files_with - If defined (defualt, None), this means that any file's basename that has this string in it (not including the file suffix, i.e. '.fasta') will be omitted from the sample set. This parameter is included because often there will be duplicates in the data set that aren't meant to be included for some methods.
    """
    def __init__(self, dirpath, count_attribute_name, vgene_name, dgene_name, jgene_name, cdr3_name, ignor_allele_info=True, replace_slash_with=None, avoid_files_with=None):
        if dirpath[-1] != '/':
            dirpath += '/'
        self.dirpath = dirpath
        self.count_attribute_name = count_attribute_name
        self.vgene_name = vgene_name
        self.dgene_name = dgene_name
        self.jgene_name = jgene_name
        self.cdr3_name = cdr3_name
        self.ignor_allele_info = ignor_allele_info
        self.replace_slash_with = replace_slash_with
        self.avoid_files_with = avoid_files_with
        #get filepaths for each timepoint
        timepoints = []
        for i in os.listdir(dirpath):
            if i[0] == '.' or i[:6] == 'README':
                continue
            timepoint = i.split('_')[-1][:-6]
            if avoid_files_with:
                if avoid_files_with in timepoint:
                    continue

            #below deals with the special case of the duplicate
            #samples. Not too elegent though...
            if '-' in timepoint:
                timepoint = re.sub('-', '.', timepoint)

            timepoint = float(timepoint)
            filepath = dirpath + i
            timepoints.append([timepoint, filepath])
        #now sort based on timepoints
        self.sample_filepaths = []
        self.timepoints = []
        for i in sorted(timepoints):
            self.sample_filepaths.append(i[1])
            self.timepoints.append(i[0])
        return

    def get_gene_expressions(self, gene_class, write_to_file=None, normailze_by_gene_length=False, drop_allele_info=True, sort_by='range', gene_segment_reference_dirpath='/Users/nstrauli/data/ref_seqs/Ig/from_imgt_appended/human/'):
        """
        This method gets the expression values for each possible gene segment of a provided gene segment class (IGHV, IGHJ, IGKJ, etc). It uses the 'get_gene_expression' method within 'immune_sample_class.py' to get the expression values for each timepoint.
        gene_class = This tells the script for what gene class expression levels should be calculated. If 'gene_class' is a list then the script will calculate expression levels for each gene class in the list. Exceptable values are 'IGHV', 'IGHD', IGHJ', 'IGLV', 'IGLJ', 'IGKV', 'IGKJ', 'TRAV', 'TRAJ', 'TRBV', 'TRBD', 'TRBJ', 'TRDV', 'TRDJ', 'TRDD', 'TRGV', or 'TRGJ'.
        normailze_by_gene_length = Boolean. Default=False. Instructs the script if each gene segments expression level should be normalized by gene length, in addition to the total number of reads.
        drop_allele_info = Boolean. Default=True. If true, the allele information in the gene names is disregarded. So IGHV4-7*01 will be identical to IGHV4-7*05.
        sort_by = This instructs how the genes will be sorted in the output. Acceptable values are:
            'range': max(expression value) - min(expression value)
            'mean': mean(expression value)
        write_to_file = Default=None. If defined (i.e. != None), then this should be the directory path for where the output should be written. Output is a tab delimited file, where each row is the expression trajectory for a gene segment. Each provided gene class will have a corresponding output file, with the same basename.
        gene_segment_reference_dirpath = the path to the directory that contains the fasta files for each of the gene classes. The fasta files will have all the germline genes (for the given class) and their sequences. However, BEWARE, the default directory has some pseudo genes listed that don't have accompanied sequences. They just have 'NA' where the sequence should be. This is fine if one simply needs a list of gene names, but not good if one needs to normalize by the length of the gene sequence (i.e. when using RNAseq data). So, if using RNAseq data then this filepath should be changed.
        """
        #check if gene_class is a string (i.e. singular)
        if isinstance(gene_class, basestring):
            #if it is, make it a list (one element long)
            gene_class = [gene_class]

        #initialize the gene expression dictionaries
        #(one for each queried gene class
        gene_expressions = [{} for i in gene_class]
        for i in self.sample_filepaths:
            sample = immune_sample(i, self.count_attribute_name, self.vgene_name, self.dgene_name, self.jgene_name, self.cdr3_name, self.ignor_allele_info)
            gene_expression = sample.get_gene_expression(gene_class, normailze_by_gene_length, drop_allele_info, gene_segment_reference_dirpath)
            #If only one query gene class
            if isinstance(gene_expression, dict):
                #then make it a list of len 1
                gene_expression = [gene_expression]
            #cycle through gene expression data for each class
            for j in xrange(len(gene_expression)):
                #cycle through each gene
                for k in gene_expression[j]:
                    gene_name = k
                    expression = gene_expression[j][k]
                    try:
                        gene_expressions[j][k].append(expression)
                    except KeyError:
                        gene_expressions[j][k] = [expression]
        #change gene_expressions to a list and sort
        gene_expressions_list = []
        for i in xrange(len(gene_expressions)):
            gene_expressions_list.append([])
            for j in gene_expressions[i]:
                if sort_by == 'range':
                    summary_stat = max(gene_expressions[i][j]) - min(gene_expressions[i][j])
                elif sort_by == 'mean':
                    summary_stat = sum(gene_expressions[i][j]) / float(len(gene_expressions[i][j]))
                gene_expressions_list[i].append([summary_stat, j] + gene_expressions[i][j])
            gene_expressions_list[i] = sorted(gene_expressions_list[i])
        #if instructed to, write sorted gene expressions
        #trajectories to output
        if write_to_file:
            if write_to_file[-1] != '/':
                output_dirpath = write_to_file + '/'
            else:
                output_dirpath = write_to_file
            for i in xrange(len(gene_class)):
                output_filepath = output_dirpath + gene_class[i]
                fileout = open(output_filepath, "w")
                #write the header, which lists the timepoint values
                fileout.write("\t%s\n" % '\t'.join([str(j) for j in self.timepoints]))
                for j in gene_expressions_list[i]:
                    gene_name = j[1]
                    expr_traj = [str(k) for k in j[2:]]
                    fileout.write("%s\t%s\n" % (gene_name, '\t'.join(expr_traj)))
                fileout.close()
        #if gene_expressions_list has only one element (i.e. one queried
        #gene class), then remove it's highest dimensionality
        if len(gene_expressions_list) == 1:
            gene_expressions_list = gene_expressions_list[0]
        return gene_expressions_list

    def get_clone_expressions(self, write_to_file=None, define='v_j_cdr3', sort_by='range'):
        """
        This method functions similarily to 'get_gene_expressions' in that it gets the expression level of antibody elements over time. Here, we are getting the expression levels of clones over time. In this method, a clone is defined as identical CDR3 sequence, and identical V and J inferred germline genes.
        write_to_file = Default=None. If defined (i.e. != None), then this should be the file path for where the output should be written. Output is a tab delimited file, where each row is the expression trajectory for a clone.
        define - This tells the script how a clone is defined. Exceptable values are:
            'v_j_cdr3' - Default. This means that the V and J gene segment identity and CDR3 sequence define a unique clone.
            'cdr3' - This means that only the CDR3 sequence defines a unique clone.
            'v_j' - This means that the V and J gene segment identity defines a clone.
        sort_by = This instructs how the genes will be sorted in the output. Acceptable values are:
            'range': max(expression value) - min(expression value)
            'mean': mean(expression value)
        """
        num_tpoints = len(self.timepoints)
        clone_expressions = {}
        for tpoint_count, i in enumerate(self.sample_filepaths):
            sample = immune_sample(i, self.count_attribute_name, self.vgene_name, self.dgene_name, self.jgene_name, self.cdr3_name, self.ignor_allele_info)
            clone_expression = sample.get_clone_expression(define)
            for j in clone_expression:
                try:
                    clone_expressions[j][tpoint_count] = clone_expression[j]
                except KeyError:
                    #add appropriate amount of 0's before entry if it is new
                    clone_expressions[j] = [0 for k in xrange(tpoint_count)]
                    #now add the expression level for the 'current' tpoint
                    clone_expressions[j].append(clone_expression[j])
                    #now add appropriate amount of 0's after the current tpoint
                    clone_expressions[j] = clone_expressions[j] + [0 for k in xrange(num_tpoints-(tpoint_count+1))]
        #turn clone_expressions into a list and sort
        clone_expressions_list = []
        for i in clone_expressions:
            if sort_by == 'range':
                sum_stat = max(clone_expressions[i]) - min(clone_expressions[i])
            elif sort_by == 'mean':
                sum_stat = sum(clone_expressions[i]) / len(clone_expressions[i])
            else:
                print 'ill defined value for "sort_by" parameter'
                return
            clone_expressions_list.append([sum_stat, i] + clone_expressions[i])
        clone_expressions_list = sorted(clone_expressions_list)
        if write_to_file:
            fileout = open(write_to_file, "w")
            fileout.write('\t%s\n' % '\t'.join([str(i) for i in self.timepoints]))
            for i in clone_expressions_list:
                fileout.write('%s\t%s\n' % (i[1], '\t'.join(str(j) for j in i[2:])))
            fileout.close()
        return clone_expressions_list

    def get_clonal_freq_spectra(self, define='v_j_cdr3', use_counts_or_freqs='counts', truncate_to='shortest', write_output=None, min_num_reads=None, freq_bin_width=0.00001):
        """
        This method uses 'get_clonal_freq_spectrum' from 'immune_sample_class.py' to get the frequency spectra of the clones in each of the time-points. 

        Input:
        define - This tells the script how a clone is defined. Exceptable values are:
            'v_j_cdr3' - Default. This means that the V and J gene segment identity and CDR3 sequence define a unique clone.
            'cdr3' - This means that only the CDR3 sequence defines a unique clone.
        use_counts_or_freqs - This tells the script if the bins of the frequency spectrum will be counts or frequencies. acceptable values are:
            'counts' - Default. This means that each element of the resulting freq spectrum can be interpereted as "the proportion of clones that have a count of X"
            'freqs' - This means that each element of the resulting freq spectrum can be interpereted as "the proportion of clones that have a frequency of X"
        truncate_to - This tell the script how each of the frequency spectra will be truncate (because they will have differing lengths, and we want them to have the same lengths). Exceptable values are:
            'shortest' - This means that the freq spectra will be shortened to be the length of the shortest one. The elements of each of the lists that are truncated will be summed, so that the final value in the truncated list reflects the cumulative value of the elements that were removed.
            None - This makes it so no truncation happens, and instead the lists are lengthened (by adding 0's) so that they are all the same length as the longest freq spectrum.
            integer - If truncate_to is an integer, then all the lists will be shortened to have a length of that integer. If a given freq spectrum is shorter than the given integer, then it will be lengthened so that it's length equals that if the given integer.
        write_output - Default is None, but if this is defined then it should be defined as a path to a directory for which the output shall be written. There is one output file per time-point/sample. Files are tab delimited, pretty straight forward.
        min_num_reads - If defined (default None), this will give a cutoff, where any freq spectrum that is from a data set that has less than this cutoff, it will be discarded. It will also discard the elements of 'self.sample_filepaths' and 'self.timepoints' that correspond to the samples that do not meet this threshold, so use with care.
        freq_bin_width - This gives the width of a frequency bin in the spectra. I.e. this determines how far apart bin edges are. This is ignored if use_counts_or_freqs='counts'.
        """
        #first get freq_spectra
        freq_spectra = []
        bin_edges_list = []
        num_reads = []
        for i in self.sample_filepaths:
            sample = immune_sample(i, self.count_attribute_name, self.vgene_name, self.dgene_name, self.jgene_name, self.cdr3_name, self.ignor_allele_info)
            freq_spectrum, bin_edges = sample.get_clonal_freq_spectrum(define=define, use_counts_or_freqs=use_counts_or_freqs, freq_bin_width=freq_bin_width)
            freq_spectra.append(freq_spectrum)
            bin_edges_list.append(bin_edges)
            num_reads.append(sample.total)
        #if applicable, remove the spectra that don't meet
        #'min_num_reads'
        if min_num_reads:
            new_freq_spectra = []
            new_self_sample_filepaths = []
            new_self_timepoints = []
            new_num_reads = []
            for i in xrange(len(num_reads)):
                if num_reads[i] >= min_num_reads:
                    new_freq_spectra.append(freq_spectra[i])
                    new_self_sample_filepaths.append(self.sample_filepaths[i])
                    new_self_timepoints.append(self.timepoints[i])
                    new_num_reads.append(num_reads[i])
            freq_spectra = new_freq_spectra
            self.sample_filepaths = new_self_sample_filepaths
            self.timepoints = new_self_timepoints
            num_reads = new_num_reads

        #now make all the freq_spectra the same length
        #in the appropriate way (given by 'truncate_to')
        if len(freq_spectra) == 0:
            print 'all samples were below "min_num_reads". Aborted'
            return
        freq_spectra_lens = [len(i) for i in freq_spectra]
        max_len = max(freq_spectra_lens)
        min_len = min(freq_spectra_lens)
        bin_edge_increment = bin_edges_list[0][1] - bin_edges_list[0][0]

        if truncate_to == None:
            truncate_to = max_len
        elif truncate_to == 'shortest':
            truncate_to = min_len
        for i in xrange(len(freq_spectra)):
            if freq_spectra_lens[i] < truncate_to:
                to_add = [0 for j in xrange(freq_spectra_lens[i], truncate_to)]
                freq_spectra[i] = freq_spectra[i] + to_add
                for j in xrange(freq_spectra_lens[i], truncate_to):
                    bin_edges_list[i].append(bin_edges_list[i][-1] + bin_edge_increment)
            elif freq_spectra_lens[i] == truncate_to:
                pass
            else:
                to_be_cut = freq_spectra[i][truncate_to-1:]
                del freq_spectra[i][truncate_to-1:]
                to_be_cut = sum(to_be_cut)
                freq_spectra[i].append(to_be_cut)
                del bin_edges_list[i][truncate_to:]
        if write_output:
            if write_output[-1] != '/':
                write_output += '/'
            for sample_count, i in enumerate(self.sample_filepaths):
                output_filepath = '%s%s' % (write_output, os.path.basename(i))
                if output_filepath[-6:] == '.fasta':
                    output_filepath = output_filepath[:-6]
                elif output_filepath[-3:] == '.fa':
                    output_filepath = output_filepath[:-3]
                fileout = open(output_filepath, "w")
                fileout.write('right_edge_of_freq_bin_exclusive\tproportion_of_clones\n')
                for freq_bin, freq_proportion in itertools.izip(bin_edges_list[sample_count], freq_spectra[sample_count]):
                    fileout.write('%s\t%s\n' % (freq_bin, freq_proportion))
                fileout.close()
        return freq_spectra, bin_edges_list

    def group_clones_foreach_sample(self, define='v_j_cdr3', write=False, output_dirpath=None, sort_by_count=False):
        """
        This method uses the 'group_clones' method in the immune_sample_class to get the groups of clones for each of the time-points in the set.
        define - This tells the script how a clone is defined. Acceptable values are:
            'v_j_cdr3' - Default. This means that the V and J gene segment identity and CDR3 sequence define a unique clone.
            'cdr3' - This means that only the CDR3 sequence defines a unique clone.
        write - If true (default, False) then the information for each of the clones will be written to disk.
        output_dirpath - Gives the path to the directory that the clones will be written to. The file names for the individual time-points will be the same as the basename of the corresponding sequence filepaths. Disregarded if write=False.
        sort_by_count - If True (default False), then will sort the clones by their counts in decending order before writing. Disregarded if write=False.
        """
        if output_dirpath:
            if output_dirpath[-1] != '/':
                output_dirpath += '/'
            if not os.path.exists(output_dirpath):
                os.makedirs(output_dirpath)
        #this will be a list, where each element in the list
        #contains the clonal sequence information for a given
        #time-point. Elements of this list are in chron.
        #order (increasing).
        clone_list = []
        for i in xrange(len(self.sample_filepaths)):
            sample = immune_sample(self.sample_filepaths[i], self.count_attribute_name, self.vgene_name, self.dgene_name, self.jgene_name, self.cdr3_name, self.ignor_allele_info)
            sample.group_clones(define)
            clone_list.append(sample.clones)
            if write:
                output_filepath = '%s%s.fasta' % (output_dirpath, os.path.basename(self.sample_filepaths[i])[:-6])
                sample.write_clone_seqs(output_filepath, sort_by_count)
        return clone_list

    def get_lineage_expression_trajectories(self, min_ident_within_cluster=0.97, min_ident_across_tpoints=.67, output_filepath=None, path_to_vsearch=None, temp_dirpath=None, freqs_or_counts='freqs', sort_lineages_by='range'):
        """
        This method will cluster the immune seqs within each single timepoint, and then cluster the centroid sequences that result from this, across time-points. In this way we can estimate lineages, and track their relative expression over time.
        min_ident_within_cluster - This gives the clustering parameter for clustering seqs within a tpoint. See the immune_sample_class.cluster_clones for explanation
        min_ident_across_tpoints - This gives the minimun identity for centroid sequences across time-points to be considered in the same cluster (i.e. lineage). This value will most likely be lower than 'min_ident_within_cluster'.
        output_filepath - This give the path to the directory for which the lineage expression trajectories will be written. The format is tab delimited, where the first column gives the identity of the lineage (i.e. Vgene identity, Jgene identity, and centroid sequence, separted by '_'), and the following columns give the expression level of each lineage at each of the respective time-points. If equals None, then no output will be written.
        path_to_vsearch - This gives the path to the vsearch excecutable. If this is None (default) then the path is assumed to be 'vsearch' (i.e. in $PATH).
        temp_dirpath - This gives the path to the directory for which temporary directories will be made within. The temporary directories will be deleted at the end of the script, but 'temp_dirpath' will not.
        freqs_or_counts - This indicates whether the relative frequency (default) or the absolute counts for each lineage should be reported at each of the time-points. If 'freq' then gives relative frequency, if 'count' then gives absolute counts.
        sort_lineages_by - This tells how the lineages should be sorted when written to file. Exceptable vlaues are:
            'range' - This means the absolute range is used to sort them (in descending order). Range is defined as the tpoint with max expression minus the tpoint with min expression.
            'sum' - This means the lineages will be sorted by the summation of their expression values over the time-course.
        """
        #make sure freqs_or_counts has exceptable value
        if freqs_or_counts != 'freqs' and freqs_or_counts != 'counts':
            print 'freqs_or_counts parameter does not have an exceptable value:', freqs_or_counts
            return
        #make sure sort_lineages_by has exceptable value
        if sort_lineages_by != 'range' and sort_lineages_by != 'sum':
            print 'sort_lineages_by parameter does not have an exceptable value:', sort_lineages_by
            return
        if temp_dirpath == None:
            temp_dirpath = tempfile.mkdtemp(dir=os.getcwd()) + '/'
        else:
            temp_dirpath = tempfile.mkdtemp(dir=temp_dirpath) + '/'
        #collect all the centroids over the time-course
        master_centroid_dic = {}
        for count, i in enumerate(self.sample_filepaths):
            print '\tclustering sequences for sample:', os.path.basename(i)
            start_time = time.time()
            sample = immune_sample(i, self.count_attribute_name, self.vgene_name, self.dgene_name, self.jgene_name, self.cdr3_name, self.ignor_allele_info)
            centroid_dic = sample.cluster_clones(min_ident_within_cluster, path_to_vsearch, temp_dirpath)
            for j in centroid_dic:
                try:
                    master_centroid_dic[j]
                except KeyError:
                    master_centroid_dic[j] = []
                for k in centroid_dic[j]:
                    #'count' contains the tpoint info
                    master_centroid_dic[j].append(k + [count])
            end_time = time.time()
            print '\ttime elapsed:', end_time-start_time
        #This will be a dictionary where each index represents a
        #unique lineage. This index is defined by the V and J gene
        #pair, and the centroid seq for the lineage. Each index is
        #defined by the counts of sequences belonging to that lineage
        #in each timepoint
        lineage_expr_trajs = {}
        #Cluster the centroids across all the tpoints, that are within
        #a V and J gene pair group
        print '\tclustering centroids across samples'
        start_time = time.time()
        for i in master_centroid_dic:
            #write all centroids within this group as a temp
            #fasta file
            temp_fasta_filepath = '%s%s.fasta' % (temp_dirpath, i)
            fileout = open(temp_fasta_filepath, "w")
            for count, j in enumerate(master_centroid_dic[i]):
                fileout.write('>%s;size=%s;tpoint=%s;freq=%s;\n%s\n' % (count, j[0], j[3], j[1], j[2]))
            fileout.close()
            #cluster the centroids within this group
            msa_output_filepath = '%s%s_msa.fasta' % (temp_dirpath, i)
            p = Popen([path_to_vsearch, '--cluster_size', temp_fasta_filepath, '--id', str(min_ident_across_tpoints), '--sizein', '--sizeout', '--msaout', msa_output_filepath, '--fasta_width', str(0)], stdout=PIPE, stderr=PIPE)
            out_err = p.communicate()
            #parse MSA output to get centroids and counts for each timepoint
            filein = open(msa_output_filepath, "r")
            while True:
                line = filein.readline()
                #if EOF, break loop
                if not line:
                    break
                elif line[0] == '>':
                    #if this is the centroid sequence
                    if line[1] == '*':
                        line = line[:-1].split(';')
                        tpoint = int(line[2].split('=')[1])
                        #get centroid seq for index
                        centroid_seq = filein.readline()[:-1]
                        #remove gaps
                        centroid_seq = re.sub('-', '', centroid_seq)
                        index = i + '_' + centroid_seq
                        #initialize the dic entry
                        lineage_expr_trajs[index] = [0 for j in self.timepoints]
                        if freqs_or_counts == 'freqs':
                            freq = float(line[3].split('=')[1])
                            lineage_expr_trajs[index][tpoint] += freq
                        else:
                            count = int(line[1].split('=')[1])
                            lineage_expr_trajs[index][tpoint] += count
                    elif line[:-1] != '>consensus':
                        line = line[:-1].split(';')
                        tpoint = int(line[2].split('=')[1])
                        if freqs_or_counts == 'freqs':
                            freq = float(line[3].split('=')[1])
                            lineage_expr_trajs[index][tpoint] += freq
                        else:
                            count = int(line[1].split('=')[1])
                            lineage_expr_trajs[index][tpoint] += count
            filein.close()
            #delete temp files
            subprocess.call(['rm', temp_fasta_filepath, msa_output_filepath])
        subprocess.call(['rm', '-r', temp_dirpath])
        end_time = time.time()
        print '\ttime elapsed:', end_time-start_time
        if output_filepath:
            #turn lineage_expr_trajs into a list for sorting
            if sort_lineages_by == 'sum':
                lineage_expr_trajs_list = [[sum(lineage_expr_trajs[i]), i]+lineage_expr_trajs[i] for i in lineage_expr_trajs]
            elif sort_lineages_by == 'range':
                lineage_expr_trajs_list = [[max(lineage_expr_trajs[i])-min(lineage_expr_trajs[i]), i]+lineage_expr_trajs[i] for i in lineage_expr_trajs]
            fileout = open(output_filepath, "w")
            #write header
            fileout.write('\t' + '\t'.join([str(i) for i in self.timepoints]) + '\n')
            for i in sorted(lineage_expr_trajs_list, reverse=True):
                fileout.write(i[1]+'\t'+'\t'.join([str(j) for j in i[2:]]) + '\n')
            fileout.close()
        return lineage_expr_trajs

    def simulate_time_points(self, element_to_simulate='clone', trials=1, serial_samp_pop_counts=None, update_freqs_with_serial_samp=True, output_dirpath=None, clone_def='v_j_cdr3', use_comp_cluster=False, temp_dirpath=None, time_per_generation=None, store_full_data_in_ram=True):
        """
        This method will simulate the time-course, based upon the frequencies of the 1st time-point. It largely relies upon the 'simulate_new_samples' method from 'immune_sample_class.py'. See the method for further details.
        element_to_simulate - Signifies which component of the sequence data will be simulated. Essentially any component of the immune population that has a frequency could be simulated (i.e. V gene, clones, etc), so this tells the script what that component is. Exceptable values are:
            'clone' - (default) This means that the frequencies of clones will be calculated and used to simulate new samples of clones.
            'v_gene' - This means that the frequencies of V genes will be calculated and used to simulate new samples of V genes.
            'j_gene' - Same as 'v_gene', but for J genes.
            'entries' - This means that each entry in the in input fasta file will be simulated. For each entry, its frequency is determined by the count attribute for that entry, and then simulated to make new samples of sequence entries.
        trials - This gives the number (int) of independent simulations to run
        serial_samp_pop_counts - If the immune population size changes with your time-points (i.e. differing number of B-cells in different samples), then this can be incorporated into the simulation with this variable. This should be a numeric list of population sizes, ordered by time. The length of this list should be equal to the number of time-points in the sample set. If this parameter equals None (default) then the popsize for each time-point is assumed to be equal to the number of reads (i.e. self.total) for each of the samples in the set.
        update_freqs_with_serial_samp - If True (default) then the frequencies of elements of the previous time-point are used to determine the probabiliteis of selecting the elements for the next time-point (similar to a Wright-Fisher sim). If False, then the frequencies of the 1st (the true, observed time-point) are always used for selecting immune elements in all time-points.
        output_dirpath - If defined (default None), then this should be defined as a path (str) to the directory for which all of the simulations will be written. Each output file corresponds to one 'trial'. The format of the output is a tab delimited list, where the 1st column is the names (IDs) of each simulated element, the 2nd column is the observed counts for each element in the sample (i.e. 1st time-point), the following columns are the simulated counts for each element, and columns are ordered with serial samples.
        clone_def - If element_to_simulate='clone' then this provides the information for how a clone is defined. This parameter is then passed to the 'group_clones' method in the immune_sample_class. If element_to_simulate does not equal 'clone' then this parameter is ignored. Exceptable values are:
            'v_j_cdr3' - Default. This means that the V and J gene segment identity and CDR3 sequence define a unique clone.
            'cdr3' - This means that only the CDR3 sequence defines a unique clone.
        use_comp_cluster - If True (default, False) then this will instruct the method to use an SGE interface to a computational cluster. Each job will be the computation for one of the trials. So, the number of jobs will equal 'trials'. It will still loop through the sample/time-points as before, but within a sample the trials will computed in parallele
        temp_dirpath - If defined (default, None), then this should be the path to a directory where temp files will be written. This parameter is only considered if 'use_comp_cluster' is True.
        time_per_generation - This gives how many units of time per 1 generation (i.e. per 1 sampling event). A unit of time is assumed to be the same as that used in 'serial_samp_tpoints'. This tells the method how many serial sampling events should take place between simulated time-points. If this equals None (default) then there is only one sampling event between each serial sample. if 'serial_samp_tpoints' is not defined then it is assumed the there is 1 unit of time between each consequetive sample. So, this would need to be set to <= 0.5 in order to sample multiple times between time-points, in that case. This can also be a list, in which case multiple values will be cycled through. This means that the entire program will be run for each of the parameter values in this list. So, all the trials, time-points, and elements will be simulated for each of the 'time_per_generation' values given in the list. Use with caution. If this is set to be a list than the dimensionality of the output is increased by one.
        store_full_data_in_ram - If True (default), then the method will store all the results of the sims in the output variable 'sim_element_counts_record', and this will be returned as output. If this is False, then the method will still return sim_elements_count_record, but instead of this variable containing the counts for each of the elements in each of the sims, it will contain the (absolute) paths to the files that contain this information. These paths will lead to files located in the 'temp_dirpath', so this variable must be defined as well.

        OUTPUT:
        Returns 'element_names' and 'sim_element_counts_record'. 'element_names' is a list that gives some unique name, or ID, for each of the simulated elements of the immune population. For example, if simulating V genes, this would be an ordered list of V gene names. 'sim_element_counts_record' is a record of the counts for each simulated element for each of the trials/serial samples. This is a 3 dimensional list. The 1st level of the list corresponds to trial number, the 2nd level corresponds to serial sample number, and the 3rd corresponds unique simulated element (i.e. unique V gene, if that is what one is simulating). The elements of the 3rd level are given in the same order as 'element_names'.
        """

        #Check to make sure there is a pop count (if defined) for each time-point
        if serial_samp_pop_counts:
            if len(serial_samp_pop_counts) != len(self.timepoints):
                print "'serial_samp_pop_counts' variable does not have same length as the number of time-points in the sample set. len(serial_samp_pop_counts):", len(serial_samp_pop_counts)
                return
        #if pop count is not defined then get those values from the read count for each sample
        else:
            serial_samp_pop_counts = []
            for i in self.sample_filepaths:
                sample = immune_sample(filepath=i, count_attribute_name=self.count_attribute_name, vgene_name=self.vgene_name, dgene_name=self.dgene_name, jgene_name=self.jgene_name, cdr3_name=self.cdr3_name, ignor_allele_info=self.ignor_allele_info)
                serial_samp_pop_counts.append(sample.total)
        #make immune sample class of 1st time-point
        sample = immune_sample(filepath=self.sample_filepaths[0], count_attribute_name=self.count_attribute_name, vgene_name=self.vgene_name, dgene_name=self.dgene_name, jgene_name=self.jgene_name, cdr3_name=self.cdr3_name, ignor_allele_info=self.ignor_allele_info)
        #if simulating clones, then must define and gather clonal seqs first
        if element_to_simulate=='clone':
            sample.group_clones(define=clone_def)
        #run simulations
        element_names, sim_element_counts_record = sample.simulate_new_samples(element_to_simulate=element_to_simulate, trials=trials, serial_samples_to_sim=len(self.timepoints)-1, serial_samp_pop_counts=serial_samp_pop_counts, serial_samp_tpoints=self.timepoints, update_freqs_with_serial_samp=update_freqs_with_serial_samp, output_dirpath=output_dirpath, include_observed_sample=True, use_comp_cluster=use_comp_cluster, temp_dirpath=temp_dirpath, time_per_generation=time_per_generation, store_full_data_in_ram=store_full_data_in_ram)
        return element_names, sim_element_counts_record

    def remove_time_points(self, tpoints_to_remove):
        """
        This method will remove the time-points given from the time-series.
        tpoints_to_remove - This is a numeric list that gives the indeces of the time-points that should be removed. So, if the time-points are [0,5,7,10,15], and you want to remove the '5' and '10' time-points, then tpoints_to_remove would equal [1,3].
        """
        new_sample_filepaths = [self.sample_filepaths[i] for i in xrange(len(self.sample_filepaths)) if not i in tpoints_to_remove]
        self.sample_filepaths = new_sample_filepaths[:]
        new_timepoints = [self.timepoints[i] for i in xrange(len(self.timepoints)) if not i in tpoints_to_remove]
        self.timepoints = new_timepoints[:]
        return

    def add_annotations_using_mixcr_foreach_sample(self, output_dirpath, fasta_output_dirpath, annotations_to_get, delete_output=False, debug=False):
        """
        This method uses the method in 'add_annotations_using_mixcr' to first run MIXCR align, and then use the 'mixcr_align_output' class to parse desired information out of the MIXCR output. The MIXCR executables MUST be in your PATH for this to work.
        output_dirpath - The path to the directory that will contain the output from running MIXCR align on all the samples.
        fasta_output_dirpath - The path to the directory for which the fasta sequence files will be written after adding the annotations from MIXCR.
        annotations_to_get - A list of strings. Each element represents an annotation that should be retrieved from the data. Exceptable values (as strings) within this list are:
            'fwr1' - will retrieve the nucleotide sequence of FWR1
            'cdr1' - will retrieve the nucleotide sequence of CDR1
            'fwr2' - will retrieve the nucleotide sequence of FWR2
            'cdr2' - will retrieve the nucleotide sequence of CDR2
            'fwr3' - will retrieve the nucleotide sequence of FWR3
            'cdr3' - will retrieve the nucleotide sequence of CDR3
            'fwr4' - will retrieve the nucleotide sequence of FWR4
            'region_bounds' - will retieve the bounds (coordinants from beggining to end, indexed at 0) for each of the regions listed above
            'fwr1aa' - will retrieve the amino acid sequence of FWR1
            'cdr1aa' - will retrieve the amino acid sequence of CDR1
            'fwr2aa' - will retrieve the amino acid sequence of FWR2
            'cdr2aa' - will retrieve the amino acid sequence of CDR2
            'fwr3aa' - will retrieve the amino acid sequence of FWR3
            'cdr3aa' - will retrieve the amino acid sequence of CDR3
            'fwr4aa' - will retrieve the amino acid sequence of FWR4
            'seqVDJaa' - will retrieve the amino acid sequence of the VDJ region of the input seq
        delete_output - Boolean. If True (default=False), then will delete the output files generated by MIXCR after extracting the desirable info from them.
        debug - Boolean. If True (default, False) will run through some checks at the end of the script to make sure things add up.
        """
        if output_dirpath[-1] != '/':
            output_dirpath += '/'
        if not os.path.exists(output_dirpath):
            os.makedirs(output_dirpath)
        if fasta_output_dirpath[-1] != '/':
            fasta_output_dirpath += '/'
        if not os.path.exists(fasta_output_dirpath):
            os.makedirs(fasta_output_dirpath)
        for i in xrange(len(self.sample_filepaths)):
            sample = immune_sample(filepath=self.sample_filepaths[i], count_attribute_name=self.count_attribute_name, vgene_name=self.vgene_name, dgene_name=self.dgene_name, jgene_name=self.jgene_name, cdr3_name=self.cdr3_name, ignor_allele_info=self.ignor_allele_info)
            output_filepath = output_dirpath + str(self.timepoints[i])
            fasta_output_filepath = fasta_output_dirpath + os.path.basename(self.sample_filepaths[i]).split('.')[0] + '.fasta'
            sample.add_annotations_using_mixcr(output_filepath=output_filepath, annotations_to_get=annotations_to_get, delete_output=delete_output, debug=debug)
            sample.write_full_data_to_disk(output_filepath=fasta_output_filepath)
        if delete_output:
            subprocess.call(['rm', '-r', output_dirpath])
        return

    def get_selection_values_using_baseline_forech_sample(self, output_dirpath, temp_data_dirpath=None, group_sel_stats=False, use_comp_cluster=False):
        """
        This method uses the method in 'get_selection_values_using_baseline' in the immune_sample_class to get the selection values for each of the samples in the time series.
        output_dirpath - path to directory for which the output from Baseline for each of the samples will be written. Writes the selection values.
        temp_data_dirpath - Path to a directory for which all the temporary input and output data for baseline will be stored. This directory is deleted upon completion.
        group_sel_stats - If not False (default False), then this will instruct the the method to group the selection values for all the seqs in a sample ('convolve') to get one summary stat for the entire sample regarding selection. If defined this should equal the directory path to the directory that will contain this sample wide selection summary stat info for each sample.
        use_comp_cluster - This indicates if the computational cluster should be used. If True (default, False), then the SGE system will be used to submit jobs to a computational cluster. One job will be submitted for each sample (regardless of the number of seqs in each file).
        """
        if output_dirpath[-1] != '/':
            output_dirpath += '/'
        if not os.path.exists(output_dirpath):
            os.makedirs(output_dirpath)
        if not temp_data_dirpath:
            temp_data_dirpath = output_dirpath
        else:
            if temp_data_dirpath[-1] != '/':
                temp_data_dirpath += '/'
            if not os.path.exists(temp_data_dirpath):
                os.makedirs(temp_data_dirpath)
        if group_sel_stats:
            if group_sel_stats[-1] != '/':
                group_sel_stats += '/'
            if not os.path.exists(group_sel_stats):
                os.makedirs(group_sel_stats)
        for i in xrange(len(self.sample_filepaths)):
            output_filepath = output_dirpath + '.'.join(os.path.basename(self.sample_filepaths[i]).split('.')[:-1])
            temp_data_sub_dirpath = temp_data_dirpath + '.'.join(os.path.basename(self.sample_filepaths[i]).split('.')[:-1]) + '_' + str(random.random()) + '/'
            if group_sel_stats:
                group_sel_stats_filepath = group_sel_stats + '.'.join(os.path.basename(self.sample_filepaths[i]).split('.')[:-1])
            if use_comp_cluster:
                p = Popen(['qsub', 'immune_time_series_class.py', 'get_baseline_selection_values_compcluster', output_filepath, str(temp_data_sub_dirpath), str(group_sel_stats_filepath), self.sample_filepaths[i], str(self.count_attribute_name), str(self.vgene_name), str(self.dgene_name), str(self.jgene_name), str(self.cdr3_name), str(self.ignor_allele_info)], stdout=PIPE, stderr=PIPE)
                out_err = p.communicate()
            else:
                sample = immune_sample(filepath=self.sample_filepaths[i], count_attribute_name=self.count_attribute_name, vgene_name=self.vgene_name, dgene_name=self.dgene_name, jgene_name=self.jgene_name, cdr3_name=self.cdr3_name, ignor_allele_info=self.ignor_allele_info)
                sample.get_selection_values_using_baseline(output_filepath=output_filepath, temp_data_dirpath=temp_data_sub_dirpath, group_sel_stats=group_sel_stats_filepath)
                subprocess.call(['rm', '-r', temp_data_sub_dirpath])
        return

    def divide_data_by_VJ_combo(self, output_dirpath, use_cluster=False):
        """
        Uses 'divide_data_by_VJ_combo' from the 'immune_sample_class' to divide the time-series by which V and J gene segment combo a given seq has. Does this for all seqs in the data.
        output_dirpath - This is the path to the directory that the method writes the divided data to.
        use_cluster - If True (defualt False), this will cause the method to implement using a computation cluster using SGE interface.
        """
        if output_dirpath[-1] != '/':
            output_dirpath += '/'
        if not use_cluster:
            for i in self.sample_filepaths:
                sample = immune_sample(filepath=i, count_attribute_name=self.count_attribute_name, vgene_name=self.vgene_name, dgene_name=self.dgene_name, jgene_name=self.jgene_name, cdr3_name=self.cdr3_name, ignor_allele_info=self.ignor_allele_info, replace_slash_with=self.replace_slash_with)
                VJ_data_dic = sample.divide_data_by_VJ_combo()
                filename_base = os.path.basename(i)
                for j in VJ_data_dic:
                    out_dir = '%s%s/' % (output_dirpath, j)
                    if not os.path.exists(out_dir):
                        os.makedirs(out_dir)
                    output_filepath = '%s%s' % (out_dir, filename_base)
                    sample.write_subset_of_seqs_to_disk_immune(seq_indices=VJ_data_dic[j], output_filepath=output_filepath)
        else:
            p = Popen(['qsub', 'immune_time_series_class.py', 'divide_data_by_VJ_combo', self.dirpath,  output_dirpath, self.count_attribute_name, self.vgene_name, self.dgene_name, self.jgene_name, self.cdr3_name, str(self.ignor_allele_info), str(self.replace_slash_with)], stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()
        return

    def convert_to_cdr3_seqs_and_write_to_disk(self, output_dirpath):
        """
        This method will first combine all the seqs that have the same CDR3 sequence, and the same V and J genes, and combine them into the same seq entry. It will then write this data to disk. It will also chnage all necessary attributes, such as 'counts' in the data. For each CDR3 sequence entry: the seq is the CDR3, the V and J genes are determined by the original seqs that make up the entry, and the CDR3 attribute will be the same the sequence entry. The rest of the attibutes will be retained, but may not make a whole lot of sense. For example, if divergence is an attribute, then it probably shouldn't be considered since the CDR3 is quite difficult to align to germline seqs. However, for compatability reasons with down stream code, we keep all the original attributes in the CDR3 seq entries. When multiple seqs are combined into one CDR3 seq, we take the attributes from the original seq that has the highest count to be represented in the new CDR3 seq entry.
        output_dirpath - the path to the directory that will contain the CDR3 seq data, post processing
        """
        if output_dirpath[-1] != '/':
            output_dirpath += '/'
        if not os.path.exists(output_dirpath):
            os.makedirs(output_dirpath)
        output_filepaths = []
        for i in self.sample_filepaths:
            filename = '.'.join(os.path.basename(i).split('.')[:-1])
            output_filepath = '%s%s.fasta' % (output_dirpath, filename)
            output_filepaths.append(output_filepath)
        num_jobs = len(output_filepaths)
        p = Popen(['qsub', '-t', '1-%s' % num_jobs, 'immune_time_series_class.py', 'convert_to_cdr3_seqs_and_write_to_disk_compcluster', ','.join(self.sample_filepaths), ','.join(output_filepaths), self.count_attribute_name, self.vgene_name, self.dgene_name, self.jgene_name, self.cdr3_name, str(self.ignor_allele_info), str(self.replace_slash_with)], stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        print out
        print err
        return


#############################################
# Below are scripts that are used for running
# jobs on the computational cluster. They
# should only be used by the methods above.
#############################################

def get_baseline_selection_values_compcluster(output_filepath, temp_data_dirpath, group_sel_stats, sample_filepath, count_attribute_name, vgene_name, dgene_name, jgene_name, cdr3_name, ignor_allele_info):
    """This uses qsub to call the method 'get_selection_values_using_baseline' from the 'immune_sample_class' in a computational cluster environment using SGE. See the 'get_selection_values_using_baseline' method for descriptions of each of these parameters."""
    if group_sel_stats == 'False':
        group_sel_stats = False
    if count_attribute_name == 'None':
        count_attribute_name = None
    if ignor_allele_info == 'True':
        ignor_allele_info = True
    else:
        ignor_allele_info = False
    sample = immune_sample(filepath=sample_filepath, count_attribute_name=count_attribute_name, vgene_name=vgene_name, dgene_name=dgene_name, jgene_name=jgene_name, cdr3_name=cdr3_name, ignor_allele_info=ignor_allele_info)
    selection_values, selection_sum_stats = sample.get_selection_values_using_baseline(output_filepath=output_filepath, temp_data_dirpath=temp_data_dirpath, group_sel_stats=group_sel_stats)
    subprocess.call(['rm', '-r', temp_data_dirpath])
    return

def divide_data_by_VJ_combo_compcluster(input_dirpath, output_dirpath, count_attribute_name, vgene_name, dgene_name, jgene_name, cdr3_name, ignor_allele_info, replace_slash_with):

    print input_dirpath
    sys.stdout.flush()

    if ignor_allele_info == 'True':
        ignor_allele_info = True
    else:
        ignor_allele_info = False
    if replace_slash_with == 'None':
        replace_slash_with = None
    for i in os.listdir(input_dirpath):
        if i[0] == '.' or i[:6] == 'README':
            continue

        print '#######'
        print i
        print '#######'
        sys.stdout.flush()

        sample_filepath = input_dirpath + i
        sample = immune_sample(filepath=sample_filepath, count_attribute_name=count_attribute_name, vgene_name=vgene_name, dgene_name=dgene_name, jgene_name=jgene_name, cdr3_name=cdr3_name, ignor_allele_info=ignor_allele_info, replace_slash_with=replace_slash_with)
        VJ_data_dic = sample.divide_data_by_VJ_combo()
        filename_base = os.path.basename(i)
        for j in VJ_data_dic:

            print '\t', j
            sys.stdout.flush()

            out_dir = '%s%s/' % (output_dirpath, j)
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            output_filepath = '%s%s' % (out_dir, filename_base)
            sample.write_subset_of_seqs_to_disk_immune(seq_indices=VJ_data_dic[j], output_filepath=output_filepath)
    return

def convert_to_cdr3_seqs_and_write_to_disk_compcluster(input_filepaths, output_filepaths, count_attribute_name, vgene_name, dgene_name, jgene_name, cdr3_name, ignor_allele_info, replace_slash_with):
    if ignor_allele_info == 'True':
        ignor_allele_info = True
    else:
        ignor_allele_info = False
    if replace_slash_with == 'None':
        replace_slash_with = None
    sge_task_id = int(os.environ['SGE_TASK_ID'])
    input_filepath = input_filepaths.split(',')[sge_task_id-1]
    output_filepath = output_filepaths.split(',')[sge_task_id-1]

    print input_filepath
    print output_filepath
    sys.stdout.flush()

    imm_samp = immune_sample(filepath=input_filepath, count_attribute_name=count_attribute_name, vgene_name=vgene_name, dgene_name=dgene_name, jgene_name=jgene_name, cdr3_name=cdr3_name, ignor_allele_info=ignor_allele_info, replace_slash_with=replace_slash_with)
    imm_samp.convert_to_cdr3_seqs_and_write_to_disk(output_filepath=output_filepath)
    return

if __name__ == '__main__':
    if sys.argv[1] == 'get_baseline_selection_values_compcluster':
        get_baseline_selection_values_compcluster(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11])
    elif sys.argv[1] == 'divide_data_by_VJ_combo':
        divide_data_by_VJ_combo_compcluster(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10])
    elif sys.argv[1] == 'convert_to_cdr3_seqs_and_write_to_disk_compcluster':
        convert_to_cdr3_seqs_and_write_to_disk_compcluster(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10])
    elif sys.argv[1] == 'test':
        #here one can test out different methods in the class
        #samples = immune_time_series('/Users/nstrauli/data/ms/clone_seqs_fasta_format/bcr/opera/IgG/10311', 'count', 'V_gene', None, 'J_gene', 'CDR3')
        #element_names, sim_element_counts_record = samples.simulate_time_points(element_to_simulate='clone', trials=3, serial_samp_pop_counts=[948000,3000,4000,2000,6000,12000,5000,9000], update_freqs_with_serial_samp=True, output_dirpath=None, clone_def='v_j_cdr3', use_comp_cluster=False, temp_dirpath='/Users/nstrauli/data/ms/temp_stuff', time_per_generation=[4,5], store_full_data_in_ram=False)
        samples = immune_time_series('/hernandez/mandrill/users/nstrauli/data/ms/clone_seqs_fasta_format/bcr/opera/IgG/10311', 'count', 'V_gene', None, 'J_gene', 'CDR3')
        element_names, sim_element_counts_record = samples.simulate_time_points(element_to_simulate='clone', trials=3, serial_samp_pop_counts=[948000,3000,4000,2000,6000,12000,5000,9000], update_freqs_with_serial_samp=True, output_dirpath=None, clone_def='v_j_cdr3', use_comp_cluster=True, temp_dirpath='/hernandez/mandrill/users/nstrauli/data/ms/temp_stuff', time_per_generation=[4,5], store_full_data_in_ram=False)
        print sim_element_counts_record
    else:
        print 'If running immune_time_series_class.py as a script, then must specify an appropriate sub-script to run.'
        print sys.argv[1]
