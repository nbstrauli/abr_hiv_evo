#!/usr/bin/python
#$ -S /usr/bin/python
#$ -e error
#$ -o out_3
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
from sequence_time_series_class import sequence_time_series
from sequence_clusters_class import sequence_clusters
from sequence_sample_class import sequence_sample
import random
from subprocess import Popen, PIPE
import subprocess
import re
from colour import Color
import seqanpy

class sequence_clusters_time_series(sequence_time_series):
    """
    This class inherits from the 'sequence_time_series' class. It takes in a dataset that consists of a time series of clustered sequences.

	Input:
    The input for this class should be a directory that contains only fasta formatted files (files named README[.txt] will be ignored). Further the fasta files' headers should be formatted in a specific way. The header for each sequence should start with an ID (could be anything, ideally unique) and then should be followed by a series of attributes. Each attribute should be delimited by a '|'. Each attribute should have a name followed by and '=', followed by the value of the attribute. An example of a well formatted header: '>seq_id|count=65\n'. Here, there is only one attribute, which is 'count'. The first entry in the header is assumed to be the sequence unique ID.                                                                         
    The fasta files should also have a specific naming scheme. Each of the filenames should end with a numeric value that indicates the time point, such that if these numeric time point indicators where sorted (numerically), the files would sort in chronological order. This numeric time point indicator in the file names should be delimited by an '_'. For example 'day_56.fasta' or 'month_34.fasta' are examples of properly formatted filenames. '110.fasta' would also work, however 'month_12_day_21.fasta' would not be advisable. This particular case would not crash the code, but would result in the time points not being sorted properly (because the time point is encoded by two numerics (month and day) as opposed to one).

    Attributes:
    dirpath - The path to the directory that contains all the sequence samples. There could be other sub-directories within this master-directory that contain other sequence samples. All ending in a fasta or fastq suffix will be included in the set. If this parameter equals None then 'list_of_filepaths' is considered
    sample_filepaths - A list of the absolute filepaths for each of the samples in the set. In no particular order.
    count_attribute_name - this gives the name that corresponds to the count of each sequence entity.
    file_suffix - This gives the suffix that is used in the input seq file data. This is determined by the first filepath that is listed in self.sample_filepaths
    freq_attribute_name - This gives the name of the seq attribute that gives the relative frequency of each of the clusters in the dataset
    indiv_seqs_attribute_name - This gives the name of the seq attribute that lists individual sequences of the seqs for each cluster in the dataset.
    indiv_seq_ids_attribute_name - This gives the name of the seq attribute that lists individual sequence IDs of the seqs for each cluster in the dataset
    indiv_seq_counts_attribute_name - This gives the name of the seq attribute that lists individual sequence counts of the seqs for each cluster in the dataset
    indiv_seq_freqs_attribute_name - This gives the name of the seq attribute that lists individual sequence frequencies of the seqs for each cluster in the dataset
    """
    
    def __init__(self, dirpath, count_attribute_name=None, freq_attribute_name=None, indiv_seqs_attribute_name='indiv_seqs', indiv_seq_ids_attribute_name='indiv_seq_ids', indiv_seq_counts_attribute_name='indiv_seq_counts', indiv_seq_freqs_attribute_name='indiv_seq_freqs'):
        if dirpath[-1] != '/':
            dirpath += '/'
        self.dirpath = dirpath
        self.count_attribute_name = count_attribute_name
        self.freq_attribute_name = freq_attribute_name
        self.indiv_seqs_attribute_name = indiv_seqs_attribute_name
        self.indiv_seq_ids_attribute_name = indiv_seq_ids_attribute_name
        self.indiv_seq_counts_attribute_name = indiv_seq_counts_attribute_name
        self.indiv_seq_freqs_attribute_name = indiv_seq_freqs_attribute_name
        #get filepaths for each timepoint
        timepoints = []
        for i in os.listdir(dirpath):
            if i[0] == '.' or i[:-4] == 'README':
                continue
            timepoint = i.split('_')[-1][:-6]
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

    @staticmethod
    def get_genetic_dists_from_needle_output(needle_output_filepath, distance_units='edit_distance'):
        """
        This is an internal method that is to be used by other methods within this class. When provided the path to the output needle (the pairwise gloabla seq alignment tool), it will parse the alignments in the output file and get the genetic distances for each of the alignments in the file. Returns a list of floats, that gives the genetic distances in the order that they occur in the file.
        distance_units - This gives how distance between a pair of seqs is actually calculated from an alignment. Acceptable values are:
            'edit_distance' - This means that the edit distance (i.e. number of changes between the sequences) will be used.
            'percent_distance' - This means that the number of changes in the alignment, normalized by the length of the alignment is used, then turned into a percent. Also known as 'percent genetic distance'.
        """
        filein = open(needle_output_filepath, "r")
        genetic_dists = []
        index = 0
        for i in filein:
            if i[:11] == '# Identity:':
                if distance_units == 'percent_distance':
                    match = re.search('\(\s*([0-9]+\.[0-9])', i)
                    dist = 100 - float(match.group(1))
                    genetic_dists.append(percent_dist)
                elif distance_units == 'edit_distance':
                    ratio = i.split()[2].split('/')
                    dist = int(ratio[1]) - int(ratio[0])
                genetic_dists.append(dist)
        filein.close()
        return genetic_dists

    @staticmethod
    def get_closest_seq_pair_dist(self, seq_list_1, seq_list_2, temp_dirpath, path_to_needle, distance_units='edit_distance'):
    	"""
		This is an internal method that is to be used by other methods within this class. It does an all-by-all comparison between the seqs in 'seq_list_1' and 'seq_list_2' to find the seq pair that has the lowest genetic distance. Returns that distance.
        distance_units - This gives how distance between a pair of seqs is actually calculated from an alignment. Acceptable values are:
            'edit_distance' - This means that the edit distance (i.e. number of changes between the sequences) will be used.
            'percent_distance' - This means that the number of changes in the alignment, normalized by the length of the alignment is used, then turned into a percent. Also known as 'percent genetic distance'.
    	"""
        random_suffix = str(random.random())
        #make temp list if seq2s in a fasta file
        temp_seq2_fasta_filepath = '%stemp_fasta_%s.fasta' % (temp_dirpath, random_suffix)
        fileout = open(temp_seq2_fasta_filepath, "w")
        for index2, seq2 in enumerate(seq_list_2):
            fileout.write(">%s\n%s\n" % (index2, seq2))
        fileout.close()
        min_genetic_dists = []
        for index1, seq1 in enumerate(seq_list_1):
            alignment_filepath = '%stemp_alignment_%s_%s' % (temp_dirpath, random_suffix, index1)
            p = Popen([path_to_needle, '-asequence', 'asis:'+seq1, '-bsequence', temp_seq2_fasta_filepath, '-gapopen', '10.0', '-gapextend', '0.5', '-outfile', alignment_filepath, '-brief', 'Y'], stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()
            genetic_dists = self.get_genetic_dists_from_needle_output(needle_output_filepath=alignment_filepath, distance_units=distance_units)
            subprocess.call(['rm', alignment_filepath])
            min_genetic_dists.append(min(genetic_dists))
        subprocess.call(['rm', temp_seq2_fasta_filepath])
        min_genetic_dist = min(min_genetic_dists)
        return min_genetic_dist

    @staticmethod
    def get_dist_between_rep_seqs(self, rep_seq_1, rep_seq_2, temp_dirpath, distance_units, path_to_needle):
        """
        This is an internal method that is to be used by other methods within this class. It compares the representative sequences of two clusters to one-another in order to get a distance between the two clusters.
        """
        random_suffix = str(random.random())
        alignment_filepath = '%stemp_alignment_%s' % (temp_dirpath, random_suffix)
        p = Popen([path_to_needle, '-asequence', 'asis:'+rep_seq_1, '-bsequence', 'asis:'+rep_seq_2, '-gapopen', '10.0', '-gapextend', '0.5', '-outfile', alignment_filepath, '-brief', 'Y'], stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        genetic_dists = self.get_genetic_dists_from_needle_output(needle_output_filepath=alignment_filepath, distance_units=distance_units)
        subprocess.call(['rm', alignment_filepath])
        return genetic_dists[0]

    @staticmethod
    def get_dist_between_rep_seqs_seqan(rep_seq_1, rep_seq_2, distance_units):
        """
        This is and internal method that uses the python package 'seqanpy' to find the distance between two sequences.
        """
        alignment = seqanpy.align_global(rep_seq_1, rep_seq_2)
        num_diffs = 0
        alignment_len = len(alignment[1])
        for i in xrange(alignment_len):
            if alignment[1][i] != alignment[2][i]:
                num_diffs += 1
        if distance_units == 'edit_distance':
            return num_diffs
        elif distance_units == 'percent_distance':
            return (num_diffs/float(alignment_len)) * 100

    @staticmethod
    def get_closest_seq_pair_dist_seqan(self, seq_list_1, seq_list_2, distance_units='edit_distance'):
        """
        This is an internal method that is to be used by other methods within this class. It does the same thing as 'get_closest_seq_pair_dist' but uses the python package 'seqanpy'.
        """
        distances = []
        for seq1 in seq_list_1:
            for seq2 in seq_list_2:
                distance = self.get_dist_between_rep_seqs_seqan(seq1, seq2, distance_units=distance_units)
                distances.append(distance)
        return min(distances)

    def create_lineages(self, distance_metric='closest_seq_pair', temp_dirpath=None, path_to_needle=None, distance_units='edit_distance', compare_tpoint_to_all_previous=False, genetic_dist_filepath=None, output_fasta_dirpath=None, write_full_original_seqs=None, muller_plot_output_dirpath=None, start_end_colors=None, max_distance=None, min_cluster_freq=None, use_comp_cluster=False, lineage_min_freq_cutoff=None, make_master_network_plot=None, indiv_sample_network_dirpath=None, mut_count_attribute_name=None, alignment_method='needle', sim_null_lineages_fasta_dirpath=None, sim_null_lineages_muller_plot_dirpath=None, multiply_num_sim_lineages_by=1, only_simulate=False, prob_of_new_lineage_filepath=None, hold_for_job_id=None):
    	"""
		This method compare each of the clusters in a given time-point to the clusters before it, to see which of them are connected as a lineage. A cluster is connected to one previous to it if the previous cluster has the closest distance to it, out of all the other previous clusters.
		distance_metric - This tells the method how 'distance' between two clusters is calculated. Acceptable values are:
			'closest_seq_pair' - This means that there is an all pairwise comparison of the seqs between the two clusters, and the gentic distance between the most similar seq pair is used.
            'dist_between_rep_seqs' - This means that the distance between the representative seqs for each of the clusters is used to calculate the distance between the cluster pair.
		temp_dirpath - Gives the path to a directory for which temp data can be written. If None (default), then uses the current working directory. If the parameter 'use_comp_cluster' is set to True, then this is automatically set to '/scratch/'
        path_to_needle - This gives the path to the needle executable. This is required if using the computational cluster for an array job because the PATH variable is all messed up when sending jobs to parallele nodes. One can also use this if needle is not in the PATH. If this is none, then it is assumed that 'needle' is in the PATH.
        distance_units - This gives how distance between a pair of seqs is actually calculated from an alignment. Acceptable values are:
            'edit_distance' - Default. This means that the edit distance (i.e. number of changes between the sequences) will be used.
            'percent_distance' - This means that the number of changes in the alignment, normalized by the length of the alignment is used, then turned into a percent. Also known as 'percent genetic distance'.
        compare_tpoint_to_all_previous - Boolean, default False. If True, this will instruct the method to compare the clusters in a given time-point to the clusters within all previous timepoints. Thus a cluster could be linked to another cluster that is 2, or more time-points back in time, rather than just the previous time-point (as is default, when True).
        genetic_dist_filepath - If defined (default, None), this gives the path to an output file that will contain a list of the genetic distances for each of the cluster comparisons performed. This is done for the purposes of plotting the distribution of these values.
        output_fasta_dirpath - If defined (default, None), this will be the path to the master directory that will contain sets of fasta files that represent each lineage. The 1st level of the directory tree will be one folder for each lineage, followed by one fasta file for each time-point.
        write_full_original_seqs - If defined (default, None), this should be the path to the directory that has the fully annotated original seqs that make up the clusters. The method will write these fully annotated seqs to disk, rather than the seq clusters (which are made up of these seqs). Importantly, these seq need to be the same seq data that was used to make the sequence cluster files. This parameter is only considered if 'output_fasta_dirpath' is defined.
        muller_plot_output_dirpath - If defined (default, None), then this will give the path to the directory that will contain the information needed to make Muller plots by the muller.plot function in R. If this is None, then this information will not be written.
        start_end_colors - If defined (default, None), this is a two element list of strings, that give the common names of the starting and ending colors for the range of lineages in the data. So, if there are 10 lineages, then there will be 10 colors and they will be evenly spaced on the color spectrum between the color named in the 1st element of this list to the color given by the last element of the list. If this is None, then colors will be randomly assigned by the muller.plot R script when plotting.
        max_distance - If defined (default, None), this will inform the method on thresholds for when two clusters will be joined across time-points. Acceptable values are:
            [float] - this will be a float that gives the maxumum distance allowed for two clusters to be joined across time-points, and thus considered a lineage. If a cluster in a time-point has no other clusters in a previous time-point that have a distance below (or equal to) this value, then it will be considered a unique lineage that starts at that time-point.
            'parsimony' - This is a max-parsimony based approach, where a cluster in the current time-point will be joined with the cluster in the previous time-point that has the min distance between the two, if and only if, the distance b/t the cluster pair is less than the distance between the cluster in the current time-point and its inferred starting point. The distance (or number of mutations) between the current cluster and a starting point is given by an attribute in the header of each cluster entry. The name of this attribute is given by the parameter 'mut_count_attribute_name'. This parameter MUST be defined if 'parsimony' is used.
            None - Default. All min distance cluster pairs joined.
        min_cluster_freq - None (default), or Float. This gives the minimum frequency allowed for a sequence cluster to be included in the data. If None, then all clusters are included. This parameter is only considered if 'self.freq_attribute_name' is defined.
        use_comp_cluster - If True (default, False), this will cause the method to be submitted as a single job to a computational cluster (using an SGE interface)
        lineage_min_freq_cutoff - If defined (default, None), this will be a float that gives the cutoff, where each lineage must have a cluster that exceeds this cutoff freq in order to be included in the data. Importantly, this only considers the clusters of a given lineage once it is separated from its parent lineage (if applicable).
        make_master_network_plot - If defined (default, None), this should be a string that gives the path to the directory that will contain a large network file (.sif) and node attribute file that will show all the connections between all the seqs in the time series. This will inform the method to make a large network file (and accompaning node attribute file) that gives all the seqs as nodes, with there connections, and then connections b/t clusters, as edges between each of the clusters' representative seqs. This is largely accomplished by appending each of the .sif and node attribute files for the sampleset together (see below).
        indiv_sample_network_dirpath - If defined (default, None), this should be a string that gives the path to the directory that contains all the network files for each of the time-points in the data. These network files should be the data for the connections b/t each of the uniqe seqs in each of the samples in the dataset (i.e. not connections b/t clusters, but connections b/t individual seqs). The network files should be in simple interaction file format (.sif). They should be named identical to their corresponding time-point sequence_cluster files. The node attribute files for each of the time-points should be named as follows: '[timepoint]_node_attributes.txt', where [timepoint] is identicle to the name in the corresponding sequence_cluster file. All node attribute files need to have the same attributes listed. Node attribute files and .sif files for all time-points should be in this directory. If make_master_network_plot is defined, this MUST be defined.
        mut_count_attribute_name - This gives the name of the attribute in the header of each cluster entry that lists the number of mutations in each seq that belongs to that cluster.
        alignment_method - This gives the method that will be used to globally align pairs of sequences. Acceptable values are:
            'needle' - The 'needle' program from the EMBL package will be used.
            'seqanpy' - The python package 'seqanpy' will be used. Should be faster.
        sim_null_lineages_fasta_dirpath - If defined (default, None), this means that null lineages will be simulated and then written in fasta format (similar to those written in 'output_fasta_dirpath'). The lineages are simulated by randomly choosing an ancestor for a given cluster, including the possibility of no ancestor, and building this out moving down the time-series.
        sim_null_lineages_muller_plot_dirpath - If defined (default, False), this will be the path to the directory that will contain the information to make Muller plots of the simulated null lineages. This is only considered if 'sim_null_lineages_fasta_dirpath' is defined.
        multiply_num_sim_lineages_by - Int. This gives the fold increase for the number of clusters per time-point in the null lineages. This effectively increases the size of the null lineages data by X fold
        only_simulate - If True (default, False), this will cause the method to only create and record simulated lineages (see 'sim_null_lineages_fasta_dirpath'). If this is True, then no 'true' lineages will be created and none of that information will be recorded, so the parameters like 'muller_plot_output_dirpath' and 'output_fasta_dirpath' will be ignored. If this is true, then the parameter 'prob_of_new_lineage_filepath' MUST be defined. This is because these probabilities are learned from the observed data.
        prob_of_new_lineage_filepath - If defined (default, None), then this is the path to the file that contains the information for what the probability for a new lineage occurring in each of the time-points (the 1st time-point, it is of course, 1). If 'only_simulate' is True, then it is assumed that this file already exists, and if it is False, then this file is made from the data.
        hold_for_job_id - If defined (default, None), this should be a string, and will hold starting the array job until a job with the provided ID is finished.
    	"""
        if not sim_null_lineages_fasta_dirpath and only_simulate:
            print 'If one wants to only create simulated lineages then the "sim_null_lineages_fasta_dirpath" must be defined. Aborting.'
            return

        if use_comp_cluster:
            if start_end_colors:
                start_end_colors = ','.join(start_end_colors)
            random_suffix = str(random.random())
            job_name = 'clust_across_%s' % random_suffix
            if hold_for_job_id:
                p = Popen(['qsub', '-hold_jid', hold_for_job_id, '-N', job_name, 'sequence_clusters_time_series_class.py', 'create_lineages_comp_cluster', distance_metric, '/scratch/', str(path_to_needle), distance_units, str(compare_tpoint_to_all_previous), str(genetic_dist_filepath), str(output_fasta_dirpath), str(write_full_original_seqs), str(muller_plot_output_dirpath), str(start_end_colors), str(max_distance), str(min_cluster_freq), self.dirpath, str(self.count_attribute_name), str(self.freq_attribute_name), self.indiv_seqs_attribute_name, self.indiv_seq_ids_attribute_name, self.indiv_seq_counts_attribute_name, self.indiv_seq_freqs_attribute_name, str(make_master_network_plot), str(indiv_sample_network_dirpath), str(lineage_min_freq_cutoff), str(mut_count_attribute_name), alignment_method, str(sim_null_lineages_fasta_dirpath), str(multiply_num_sim_lineages_by), str(only_simulate), str(prob_of_new_lineage_filepath)], stdout=PIPE, stderr=PIPE)
            else:
                p = Popen(['qsub', '-N', job_name, 'sequence_clusters_time_series_class.py', 'create_lineages_comp_cluster', distance_metric, '/scratch/', str(path_to_needle), distance_units, str(compare_tpoint_to_all_previous), str(genetic_dist_filepath), str(output_fasta_dirpath), str(write_full_original_seqs), str(muller_plot_output_dirpath), str(start_end_colors), str(max_distance), str(min_cluster_freq), self.dirpath, str(self.count_attribute_name), str(self.freq_attribute_name), self.indiv_seqs_attribute_name, self.indiv_seq_ids_attribute_name, self.indiv_seq_counts_attribute_name, self.indiv_seq_freqs_attribute_name, str(make_master_network_plot), str(indiv_sample_network_dirpath), str(lineage_min_freq_cutoff), str(mut_count_attribute_name), alignment_method, str(sim_null_lineages_fasta_dirpath), str(multiply_num_sim_lineages_by), str(only_simulate), str(prob_of_new_lineage_filepath)], stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()
            return job_name

        if len(self.sample_filepaths) <= 1:
            print 'Number of time-points is:', len(self.sample_filepaths)
            print 'This needs to be 2 or more. Aborting.'
            return

        if compare_tpoint_to_all_previous:
            "'compare_tpoint_to_all_previous' is currently not supported. Aborting."
            return

    	if temp_dirpath:
            if temp_dirpath[-1] != '/':
                temp_dirpath += '/'
            if not os.path.exists(temp_dirpath):
                os.makedirs(temp_dirpath)
        else:
            temp_dirpath = os.getcwd() + "/"
        if not path_to_needle:
            path_to_needle = 'needle'
        #load all data into memory as sequence_cluster objects
        #this will be a list of chron. ordered sequence_cluster objects
    	seq_clusts_time_series = []
        sample_clust_id_dic = {}
        filenames = []
        indeces_to_remove = []
        num_clusters_time_series = []#num clusts for each time-point
        for index, i in enumerate(self.sample_filepaths):
            cluster_sample = sequence_clusters(filepath=i, count_attribute_name=self.count_attribute_name, freq_attribute_name=self.freq_attribute_name, indiv_seqs_attribute_name=self.indiv_seqs_attribute_name, indiv_seq_ids_attribute_name=self.indiv_seq_ids_attribute_name, indiv_seq_counts_attribute_name=self.indiv_seq_counts_attribute_name, indiv_seq_freqs_attribute_name=self.indiv_seq_freqs_attribute_name, min_freq=min_cluster_freq)
            if len(cluster_sample.data) == 0:
                indeces_to_remove.append(index)
                continue
            seq_clusts_time_series.append(cluster_sample)
            filenames.append('.'.join(os.path.basename(i).split('.')[:-1]))
            if make_master_network_plot:
                sample_clust_id_dic[i] = set()
                for j in cluster_sample.data:
                    sample_clust_id_dic[i].update([j['id']])
            num_clusters_time_series.append(len(cluster_sample.data))
        self.remove_timepoints(timepoint_indeces=indeces_to_remove)

        list_of_timepoint_indices = range(len(seq_clusts_time_series))
        if genetic_dist_filepath and not only_simulate:
            fileout_dist_dstrb = open(genetic_dist_filepath, "w")
            fileout_dist_dstrb.write('genetic_distance\n')
            min_dist_filepath = '%s/min_%s' % (os.path.dirname(genetic_dist_filepath), os.path.basename(genetic_dist_filepath))
            fileout_min_dist_dstrb = open(min_dist_filepath, "w")
            fileout_min_dist_dstrb.write('min_genetic_distance\n')

        if make_master_network_plot and not only_simulate:
            if make_master_network_plot[-1] != '/':
                make_master_network_plot += '/'
            if not os.path.exists(make_master_network_plot):
                os.makedirs(make_master_network_plot)
            if indiv_sample_network_dirpath[-1] != '/':
                indiv_sample_network_dirpath += '/'
            fileout_sif = open('%snetwork_file.sif' % make_master_network_plot, "w")
            fileout_node_attrbs = open('%snode_attributes.txt' % make_master_network_plot, "w")
            first_file = True
            for sample in self.sample_filepaths:
                tpoint = '.'.join(os.path.basename(sample).split('.')[:-1])
                filein_node = open('%s%s_node_attributes.txt' % (indiv_sample_network_dirpath, tpoint), "r")
                if first_file:
                    fileout_node_attrbs.write('%s\ttime_point\tcluster_representative\n' % filein_node.readline()[:-1])
                    first_file = False
                else:
                    filein_node.readline()
                for i in filein_node:
                    line = i[:-1].split('\t')
                    seq_id = line[0]
                    if seq_id in sample_clust_id_dic[sample]:
                        cluster_main_seq = 'True'
                        sample_clust_id_dic[sample].remove(seq_id)
                    else:
                        cluster_main_seq = 'False'
                    tpoint_seq_id = '%s_%s' % (tpoint, seq_id)
                    fileout_node_attrbs.write('%s\t%s\t%s\t%s\n' % (tpoint_seq_id, '\t'.join(line[1:]), tpoint, cluster_main_seq))
                filein_node.close()
                filein_sif = open('%s%s.sif' % (indiv_sample_network_dirpath, tpoint), 'r')
                for i in filein_sif:
                    line = i[:-1].split('\t')
                    tpoint_seq_id_1 = '%s_%s' % (tpoint, line[0])
                    if len(line) == 3:
                        tpoint_seq_id_2 = '%s_%s' % (tpoint, line[2])
                        fileout_sif.write('%s\tpp\t%s\n' % (tpoint_seq_id_1, tpoint_seq_id_2))
                    else:
                        fileout_sif.write('%s\n' % (tpoint_seq_id_1))
                filein_sif.close()
            fileout_node_attrbs.close()

        #init lineages dic with clusters from 1st tpoint
        #This will be a dic with a numeric lineage ID as key and a list of of clusters that belong to that lineage as def. A cluster is represented as a two element list of [time-point_index, cluster_index].
        lineages = {}
        #this will be a dic where the index is the lineage ID and the def is that lineage, except the first element gives the lineage ID of the parent of the given lineage. If there is no parent then this is given as 'NA'. This dic of lineages is formatted for Muller plots from the Muller.plot R package. Importantly, the lineages that are the defs of this dic will be a dic that is structured that same as 'lineages' (i.e. a dic, with lineage ID as key and list of clusters as def), however it will not include the ancestral clusters in a lineage, once a new lineage is born. In other words, if a new lineage is born from another ancestral lineage, then this new lineage will not include the clusters from the ancestral lineage prior to the current time-point.
        muller_lineages = {}
        #this is a dic that maps each cluster to the lineage IDs that it belongs to.
        cluster_to_lineage_dic = {}
        if not only_simulate:
            for index, i in enumerate(seq_clusts_time_series[0].data):
                lineage_id = index + 1
                lineages[lineage_id] = [[0, index]]
                cluster_to_lineage_dic['0_%s' % index] = lineage_id
                muller_lineages[lineage_id] = ['NA', [0, index]]
            new_lineage = lineage_id + 1
        ###### if making null lineages ######
        if sim_null_lineages_fasta_dirpath:
            lineages_null = {}
            muller_lineages_null = {}
            cluster_to_lineage_dic_null = {}
            lineage_id_null = 0
            for index, i in enumerate(seq_clusts_time_series[0].data):
                for i in xrange(multiply_num_sim_lineages_by):
                    lineage_id_null += 1
                    lineages_null[lineage_id_null] = [[0, index]]
                    cluster_to_lineage_dic_null['0_%s' % index] = lineage_id_null
                    muller_lineages_null[lineage_id_null] = ['NA', [0, index]]
            new_lineage_null = lineage_id_null + 1
        #####################################

        if prob_of_new_lineage_filepath and not only_simulate:
            fileout_new_lin_prob = open(prob_of_new_lineage_filepath, "w")
            fileout_new_lin_prob.write('timepoint_index\tprob_of_new_lineage\n0\t1.0\n')
        elif only_simulate:
            filein_new_lin_prob = open(prob_of_new_lineage_filepath, "r")
            filein_new_lin_prob.readline()
            tpoint_index_to_new_lin_prob_dic = {}
            for i in filein_new_lin_prob:
                line = i[:-1].split('\t')
                tpoint_index = int(line[0])
                new_lin_prob = float(line[1])
                tpoint_index_to_new_lin_prob_dic[tpoint_index] = new_lin_prob
            filein_new_lin_prob.close()
    	for index_1 in list_of_timepoint_indices[1:]:
            if not compare_tpoint_to_all_previous:
                start_index = index_1 - 1
            else:
                start_index = 0
            for index_2 in list_of_timepoint_indices[start_index:index_1]:
                print '###########################'
                print 'time-point 1:', index_1
                print 'time-point 2:', index_2
                print '###########################'
                sys.stdout.flush()

                if not only_simulate:
                    #this will be a dic that has the previous time-point cluster_index's as the key and the def is the cluster_index's in the current time-point that they are connected to
                    connections_with_previous_tpoint = {}
                    #this is used for simulating null lineages
                    prob_of_new_lineage = 0.
                    for data_obj_1_index, data_obj_1 in enumerate(seq_clusts_time_series[index_1].data):
                        print '\tcluster 1 ID:', data_obj_1['id']
                        genetic_dists = []
                        for data_obj_2_index, data_obj_2 in enumerate(seq_clusts_time_series[index_2].data):
                            if distance_metric == 'closest_seq_pair':
                                if alignment_method == 'needle':
                                    genetic_dist = self.get_closest_seq_pair_dist(self, seq_list_1=data_obj_1['other'][self.indiv_seqs_attribute_name], seq_list_2=data_obj_2['other'][self.indiv_seqs_attribute_name], temp_dirpath=temp_dirpath, path_to_needle=path_to_needle, distance_units=distance_units)
                                elif alignment_method == 'seqanpy':
                                    genetic_dist = self.get_closest_seq_pair_dist_seqan(self, seq_list_1=data_obj_1['other'][self.indiv_seqs_attribute_name], seq_list_2=data_obj_2['other'][self.indiv_seqs_attribute_name], distance_units=distance_units)
                            elif distance_metric == 'dist_between_rep_seqs':
                                if alignment_method == 'needle':
                                    genetic_dist = self.get_dist_between_rep_seqs(self, rep_seq_1=data_obj_1['seq'], rep_seq_2=data_obj_2['seq'], temp_dirpath=temp_dirpath, distance_units=distance_units, path_to_needle=path_to_needle)
                                elif alignment_method == 'seqanpy':
                                    genetic_dist = self.get_dist_between_rep_seqs_seqan(rep_seq_1=data_obj_1['seq'], rep_seq_2=data_obj_2['seq'], distance_units=distance_units)
                            if genetic_dist_filepath:
                                fileout_dist_dstrb.write('%s\n' % genetic_dist)
                            genetic_dists.append([genetic_dist, data_obj_2_index])
                        connected_cluster = min(genetic_dists)
                        connection_cluster_index = connected_cluster[1]
                        connection_cluster_dist = connected_cluster[0]
                        if genetic_dist_filepath:
                            fileout_min_dist_dstrb.write('%s\n' % connection_cluster_dist)
                        if max_distance:
                            if max_distance == 'parsimony':
                                #get num muts for representative cluster
                                rep_seq_id = data_obj_1['id']
                                rep_seq_index = data_obj_1['other']['indiv_seq_ids'].index(rep_seq_id)
                                num_muts = float(data_obj_1['other'][mut_count_attribute_name].split(',')[rep_seq_index])
                                #if num muts less then distance to nearest cluster in previous tpoint
                                if connection_cluster_dist > num_muts:
                                    connection_cluster_index = 'NA'
                                    prob_of_new_lineage += 1
                            elif connection_cluster_dist > max_distance:
                                connection_cluster_index = 'NA'
                                prob_of_new_lineage += 1
                        try:
                            connections_with_previous_tpoint[connection_cluster_index].append([connection_cluster_dist, data_obj_1_index])
                        except KeyError:
                            connections_with_previous_tpoint[connection_cluster_index] = [[connection_cluster_dist, data_obj_1_index]]
                    prob_of_new_lineage /= data_obj_1_index + 1
                    if prob_of_new_lineage_filepath:
                        fileout_new_lin_prob.write('%s\t%s\n' % (index_1, prob_of_new_lineage))
                
                ###### if making null lineages ######
                if sim_null_lineages_fasta_dirpath:
                    if only_simulate:
                        prob_of_new_lineage = tpoint_index_to_new_lin_prob_dic[index_1]
                    print 'prob of new lineage:', prob_of_new_lineage
                    for i in xrange(multiply_num_sim_lineages_by):
                        connections_with_previous_tpoint_null = {}
                        num_previous_clusters = num_clusters_time_series[index_2]
                        for data_obj_1_index in xrange(num_clusters_time_series[index_1]):
                            if random.random() < prob_of_new_lineage:#prob of new lineage
                                connection_cluster_index = 'NA'
                            else:#or choose a random parent lineage
                                connection_cluster_index = random.randint(0, num_previous_clusters-1)
                            try:
                                connections_with_previous_tpoint_null[connection_cluster_index].append([0, data_obj_1_index])
                            except KeyError:
                                connections_with_previous_tpoint_null[connection_cluster_index] = [[0, data_obj_1_index]]
                        for connected_cluster_index in connections_with_previous_tpoint_null:
                            if connected_cluster_index == 'NA':
                                for i in connections_with_previous_tpoint_null[connected_cluster_index]:
                                    muller_lineages_null[new_lineage_null] = ['NA', [index_1, i[1]]]
                                    lineages_null[new_lineage_null] = [[index_1, i[1]]]
                                    cluster_to_lineage_dic_null['%s_%s' % (index_1, i[1])] = new_lineage_null
                                    new_lineage_null += 1
                                continue
                            lineage = cluster_to_lineage_dic_null['%s_%s' % (index_2, connected_cluster_index)]
                            muller_lineages_null[lineage].append([index_1, connections_with_previous_tpoint_null[connected_cluster_index][0][1]])
                            lineages_null[lineage].append([index_1, connections_with_previous_tpoint_null[connected_cluster_index][0][1]])
                            cluster_to_lineage_dic_null['%s_%s' % (index_1, connections_with_previous_tpoint_null[connected_cluster_index][0][1])] = lineage
                            for i in connections_with_previous_tpoint_null[connected_cluster_index][1:]:
                                muller_lineages_null[new_lineage_null] = [lineage, [index_1, i[1]]]
                                lineages_null[new_lineage_null] = lineages_null[lineage][:-1]
                                lineages_null[new_lineage_null].append([index_1, i[1]])
                                cluster_to_lineage_dic_null['%s_%s' % (index_1, i[1])] = new_lineage_null
                                new_lineage_null += 1
                #####################################

                if not only_simulate:
                    for connected_cluster_index in connections_with_previous_tpoint:
                        
                        #if there are clusters in current time-point that did not meet the max distance criterion
                        if max_distance and connected_cluster_index == 'NA':
                            for i in connections_with_previous_tpoint[connected_cluster_index]:
                                muller_lineages[new_lineage] = ['NA', [index_1, i[1]]]
                                lineages[new_lineage] = [[index_1, i[1]]]
                                cluster_to_lineage_dic['%s_%s' % (index_1, i[1])] = new_lineage
                                new_lineage += 1
                            continue
                        
                        lineage = cluster_to_lineage_dic['%s_%s' % (index_2, connected_cluster_index)]
                        #if the cluster from the previous time-point had multiple children
                        if len(connections_with_previous_tpoint[connected_cluster_index]) > 1:
                            #sort connections by genetic distance
                            sorted_connections = sorted(connections_with_previous_tpoint[connected_cluster_index])
                            #the child cluster that is closest in genetic distance to the parent cluster gets to keep the same lineage ID
                            ancestral_cluster = sorted_connections[0]
                            muller_lineages[lineage].append([index_1, sorted_connections[0][1]])
                            lineages[lineage].append([index_1, sorted_connections[0][1]])
                            cluster_to_lineage_dic['%s_%s' % (index_1, sorted_connections[0][1])] = lineage
                            
                            #add an edge to the sif file, if desired.
                            if make_master_network_plot:
                                previous_tpoint_seq_id = '%s_%s' % (filenames[index_2], seq_clusts_time_series[index_2].data[connected_cluster_index]['id'])
                                current_tpoint_seq_id = '%s_%s' % (filenames[index_1], seq_clusts_time_series[index_1].data[sorted_connections[0][1]]['id'])
                                fileout_sif.write('%s\tpp\t%s\n' % (previous_tpoint_seq_id, current_tpoint_seq_id))

                            for i in sorted_connections[1:]:
                                muller_lineages[new_lineage] = [lineage, [index_1, i[1]]]
                                lineages[new_lineage] = lineages[lineage][:-1]
                                lineages[new_lineage].append([index_1, i[1]])
                                cluster_to_lineage_dic['%s_%s' % (index_1, i[1])] = new_lineage
                                new_lineage += 1
                                #add an edge to the sif file, if desired.
                                if make_master_network_plot:
                                    current_tpoint_seq_id = '%s_%s' % (filenames[index_1], seq_clusts_time_series[index_1].data[i[1]]['id'])
                                    fileout_sif.write('%s\tpp\t%s\n' % (previous_tpoint_seq_id, current_tpoint_seq_id))
                        #else the cluster from the previous time-point has one child
                        else:
                            muller_lineages[lineage].append([index_1, connections_with_previous_tpoint[connected_cluster_index][0][1]])
                            lineages[lineage].append([index_1, connections_with_previous_tpoint[connected_cluster_index][0][1]])
                            cluster_to_lineage_dic['%s_%s' % (index_1, connections_with_previous_tpoint[connected_cluster_index][0][1])] = lineage
                            #add an edge to the sif file, if desired.
                            if make_master_network_plot:
                                previous_tpoint_seq_id = '%s_%s' % (filenames[index_2], seq_clusts_time_series[index_2].data[connected_cluster_index]['id'])
                                current_tpoint_seq_id = '%s_%s' % (filenames[index_1], seq_clusts_time_series[index_1].data[ connections_with_previous_tpoint[connected_cluster_index][0][1]]['id'])
                                fileout_sif.write('%s\tpp\t%s\n' % (previous_tpoint_seq_id, current_tpoint_seq_id))

        if make_master_network_plot:
            fileout_sif.close()

        if genetic_dist_filepath:
            fileout_dist_dstrb.close()
            fileout_min_dist_dstrb.close()

        if prob_of_new_lineage_filepath and not only_simulate:
            fileout_new_lin_prob.close()

        ##################################################################
        #Below is code for writing, then reading the lineage data to disk.
        #Use this when trouble shooting and you don't want to wait for
        #lineage structures to be sorted out by the code above.
        ##################################################################
        # fileout = open('/netapp/home/nstrauli/data/abr_hiv_coevo/temp_stuff/lineages', 'w')
        # for i in lineages:
        #     lineage = ','.join(['%s:%s' % (j[0],j[1]) for j in lineages[i]])
        #     fileout.write('%s,%s\n' % (i, lineage))
        # fileout.close()
        # fileout = open('/netapp/home/nstrauli/data/abr_hiv_coevo/temp_stuff/muller_lineages', 'w')
        # for i in muller_lineages:
        #     lineage = ','.join(['%s:%s' % (j[0],j[1]) for j in muller_lineages[i][1:]])
        #     fileout.write('%s,%s,%s\n' % (i, muller_lineages[i][0], lineage))
        # fileout.close()

        # lineages = {}
        # filein = open('/netapp/home/nstrauli/data/abr_hiv_coevo/temp_stuff/lineages', 'r')
        # for i in filein:
        #     line = i[:-1].split(',')
        #     lin_id = int(line[0])
        #     lineage = [j.split(':') for j in line[1:]]
        #     for j in xrange(len(lineage)):
        #         for k in xrange(len(lineage[j])):
        #             lineage[j][k] = int(lineage[j][k])
        #     lineages[lin_id] = lineage
        # filein.close()
        # muller_lineages = {}
        # filein = open('/netapp/home/nstrauli/data/abr_hiv_coevo/temp_stuff/muller_lineages', 'r')
        # for i in filein:
        #     line = i[:-1].split(',')
        #     lin_id = int(line[0])
        #     parent_id = line[1]
        #     if parent_id != 'NA':
        #         parent_id = int(parent_id)
        #     lineage = [j.split(':') for j in line[2:]]
        #     for j in xrange(len(lineage)):
        #         for k in xrange(len(lineage[j])):
        #             lineage[j][k] = int(lineage[j][k])
        #     muller_lineages[lin_id] = [parent_id] + lineage
        # filein.close()
        ##################################################################
        ##################################################################

        if lineage_min_freq_cutoff:
            lineages_to_remove = []
            if not only_simulate:
                for lineage_id in sorted(muller_lineages, reverse=True):
                    freq_too_low = True
                    for cluster in muller_lineages[lineage_id][1:]:
                        tpoint_index = cluster[0]
                        cluster_index = cluster[1]
                        abundance = seq_clusts_time_series[tpoint_index].data[cluster_index]['other'][self.freq_attribute_name]
                        if abundance >= lineage_min_freq_cutoff:
                            freq_too_low = False
                            break
                    if freq_too_low:
                        lineages_to_remove.append(lineage_id)
                for lineage_id in lineages_to_remove:
                    #need to check if any lineages descend from this lineage, and if so, need to fix the decsendancy of the lineages so that they now descend from an existing lineage, or none at all
                    child_lineage_lens = []
                    for other_lineage_id in muller_lineages:
                        if muller_lineages[other_lineage_id][0] == lineage_id:
                            child_lineage_lens.append([len(muller_lineages[other_lineage_id][1:]), other_lineage_id])
                    if len(child_lineage_lens) > 0:
                        child_lineage_lens = sorted(child_lineage_lens)
                        new_parent_lineage_id = child_lineage_lens[-1][1]
                        #shift the parent of the lineage to be deleted to the new chosen parent lineage
                        muller_lineages[new_parent_lineage_id][0] = muller_lineages[lineage_id][0]
                        #shift the earlier (low freq.) time-points of the lineage that will be deleted to the new parent lineage
                        first_tpoint = muller_lineages[new_parent_lineage_id][1][0]
                        clusters_to_add = []
                        for cluster in muller_lineages[lineage_id][1:]:
                            if cluster[0] < first_tpoint:
                                clusters_to_add.append(cluster)
                        muller_lineages[new_parent_lineage_id] = muller_lineages[new_parent_lineage_id][:1] + clusters_to_add + muller_lineages[new_parent_lineage_id][1:]
                        #change the rest of the child lineages to decsend from the new parent lineage
                        for child_lineage_len in child_lineage_lens[:-1]:
                            child_lineage_id = child_lineage_len[1]
                            muller_lineages[child_lineage_id][0] = new_parent_lineage_id
                    del muller_lineages[lineage_id]
                    del lineages[lineage_id]
            ###### if making null lineages ######
            if sim_null_lineages_fasta_dirpath:
                lineages_to_remove = []
                for lineage_id in sorted(muller_lineages_null, reverse=True):
                    freq_too_low = True
                    for cluster in muller_lineages_null[lineage_id][1:]:
                        tpoint_index = cluster[0]
                        cluster_index = cluster[1]
                        abundance = seq_clusts_time_series[tpoint_index].data[cluster_index]['other'][self.freq_attribute_name]
                        if abundance >= lineage_min_freq_cutoff:
                            freq_too_low = False
                            break
                    if freq_too_low:
                        lineages_to_remove.append(lineage_id)
                for lineage_id in lineages_to_remove:
                    #need to check if any lineages descend from this lineage, and if so, need to fix the decsendancy of the lineages so that they now descend from an existing lineage, or none at all
                    child_lineage_lens = []
                    for other_lineage_id in muller_lineages_null:
                        if muller_lineages_null[other_lineage_id][0] == lineage_id:
                            child_lineage_lens.append([len(muller_lineages_null[other_lineage_id][1:]), other_lineage_id])
                    if len(child_lineage_lens) > 0:
                        child_lineage_lens = sorted(child_lineage_lens)
                        new_parent_lineage_id = child_lineage_lens[-1][1]
                        #shift the parent of the lineage to be deleted to the new chosen parent lineage
                        muller_lineages_null[new_parent_lineage_id][0] = muller_lineages_null[lineage_id][0]
                        #shift the earlier (low freq.) time-points of the lineage that will be deleted to the new parent lineage
                        first_tpoint = muller_lineages_null[new_parent_lineage_id][1][0]
                        clusters_to_add = []
                        for cluster in muller_lineages_null[lineage_id][1:]:
                            if cluster[0] < first_tpoint:
                                clusters_to_add.append(cluster)
                        muller_lineages_null[new_parent_lineage_id] = muller_lineages_null[new_parent_lineage_id][:1] + clusters_to_add + muller_lineages_null[new_parent_lineage_id][1:]
                        #change the rest of the child lineages to decsend from the new parent lineage
                        for child_lineage_len in child_lineage_lens[:-1]:
                            child_lineage_id = child_lineage_len[1]
                            muller_lineages_null[child_lineage_id][0] = new_parent_lineage_id
                    del muller_lineages_null[lineage_id]
                    del lineages_null[lineage_id]
            #####################################

        if output_fasta_dirpath and not only_simulate:
            print 'Writing lineage fasta files.'
            sys.stdout.flush()
            if output_fasta_dirpath[-1] != '/':
                output_fasta_dirpath += '/'
            if not os.path.exists(output_fasta_dirpath):
                os.makedirs(output_fasta_dirpath)
            if write_full_original_seqs:
                if write_full_original_seqs[-1] != '/':
                    write_full_original_seqs += '/'
            for lineage in lineages:
                lineage_dirpath = '%s%s/' % (output_fasta_dirpath, lineage)
                if not os.path.exists(lineage_dirpath):
                    os.makedirs(lineage_dirpath)
                for cluster in lineages[lineage]:
                    tpoint_index = cluster[0]
                    output_file_basename = os.path.basename(self.sample_filepaths[tpoint_index])
                    cluster_index = cluster[1]
                    output_filepath = '%s%s' % (lineage_dirpath, output_file_basename)
                    if write_full_original_seqs:
                        orig_seqs_filepath = '%s%s' % (write_full_original_seqs, output_file_basename)
                        sample = sequence_sample(filepath=orig_seqs_filepath, count_attribute_name=self.count_attribute_name)
                        seq_ids = seq_clusts_time_series[tpoint_index].data[cluster_index]['other'][self.indiv_seq_ids_attribute_name]
                        sample.write_subset_of_seqs_to_disk(seq_indices=None, output_filepath=output_filepath, append_to_file=False, seq_ids=seq_ids, add_string_to_each_id_being_written=None)
                    else:
                        seq_clusts_time_series[tpoint_index].write_subset_of_seqs_to_disk(seq_indices=[cluster_index], output_filepath=output_filepath, append_to_file=False, seq_ids=None, add_string_to_each_id_being_written=None)

        ###### if making null lineages ######
        if sim_null_lineages_fasta_dirpath:
            print 'Writing simulated lineage fasta files.'
            sys.stdout.flush()
            if sim_null_lineages_fasta_dirpath[-1] != '/':
                sim_null_lineages_fasta_dirpath += '/'
            if not os.path.exists(sim_null_lineages_fasta_dirpath):
                os.makedirs(sim_null_lineages_fasta_dirpath)
            if write_full_original_seqs:
                if write_full_original_seqs[-1] != '/':
                    write_full_original_seqs += '/'
            for lineage in lineages_null:
                lineage_dirpath = '%s%s/' % (sim_null_lineages_fasta_dirpath, lineage)
                if not os.path.exists(lineage_dirpath):
                    os.makedirs(lineage_dirpath)
                for cluster in lineages_null[lineage]:
                    tpoint_index = cluster[0]
                    output_file_basename = os.path.basename(self.sample_filepaths[tpoint_index])
                    cluster_index = cluster[1]
                    output_filepath = '%s%s' % (lineage_dirpath, output_file_basename)
                    if write_full_original_seqs:
                        orig_seqs_filepath = '%s%s' % (write_full_original_seqs, output_file_basename)
                        sample = sequence_sample(filepath=orig_seqs_filepath, count_attribute_name=self.count_attribute_name)
                        seq_ids = seq_clusts_time_series[tpoint_index].data[cluster_index]['other'][self.indiv_seq_ids_attribute_name]
                        sample.write_subset_of_seqs_to_disk(seq_indices=None, output_filepath=output_filepath, append_to_file=False, seq_ids=seq_ids, add_string_to_each_id_being_written=None)

                        filein = open(output_filepath, "r")
                        l = filein.readline()
                        filein.close()
                        if l == '':
                            print seq_clusts_time_series[tpoint_index].filepath
                            print orig_seqs_filepath
                            print lineage
                            print lineages_null[lineage]
                            print cluster
                            print seq_ids
                            print ''

                    else:
                        seq_clusts_time_series[tpoint_index].write_subset_of_seqs_to_disk(seq_indices=[cluster_index], output_filepath=output_filepath, append_to_file=False, seq_ids=None, add_string_to_each_id_being_written=None)
        #####################################

        if muller_plot_output_dirpath and not only_simulate:
            print 'Writing Muller plot data.'
            if muller_plot_output_dirpath[-1] != '/':
                muller_plot_output_dirpath += '/'
            if not os.path.exists(muller_plot_output_dirpath):
                os.makedirs(muller_plot_output_dirpath)
            attributes_filepath = muller_plot_output_dirpath + 'attributes.txt'
            fileout_attributes = open(attributes_filepath, "w")
            if start_end_colors:
                fileout_attributes.write('names\tparents\tcolors\n')
                start_color = Color(start_end_colors[0])
                end_color = Color(start_end_colors[1])
                colors = list(start_color.range_to(end_color, len(muller_lineages)))
                hex_colors = [i.hex_l for i in colors]
            else:
                fileout_attributes.write('names\tparents\n')
            population_data_table_filepath = muller_plot_output_dirpath + 'population_data_table.txt'
            fileout_pop_data = open(population_data_table_filepath, "w")
            fileout_pop_data.write('\t%s\n' % '\t'.join([str(i) for i in self.timepoints]))
            for color_index, lineage in enumerate(muller_lineages):
                fileout_attributes.write('%s\t%s' % (lineage, muller_lineages[lineage][0]))
                if start_end_colors:
                    fileout_attributes.write('\t%s' % hex_colors[color_index])
                fileout_attributes.write('\n')
                #add 0's to beginning of lineage if need be
                full_lineage = [0. for i in xrange(muller_lineages[lineage][1][0])]
                for i in muller_lineages[lineage][1:]:
                    tpoint_index = i[0]
                    cluster_index = i[1]
                    abundance = seq_clusts_time_series[tpoint_index].data[cluster_index]['other'][self.freq_attribute_name]
                    full_lineage.append(abundance)
                #add 0's to end of lineage if need be
                full_lineage += [0 for i in xrange(list_of_timepoint_indices[-1] - muller_lineages[lineage][-1][0])]
                fileout_pop_data.write('%s\t%s\n' % (lineage, '\t'.join([str(i) for i in full_lineage])))
            fileout_attributes.close()
            fileout_pop_data.close()

        ###### if making null lineages ######
        if sim_null_lineages_fasta_dirpath:
            if sim_null_lineages_muller_plot_dirpath:
                if sim_null_lineages_muller_plot_dirpath[-1] != '/':
                    sim_null_lineages_muller_plot_dirpath += '/'
                if not os.path.exists(sim_null_lineages_muller_plot_dirpath):
                    os.makedirs(sim_null_lineages_muller_plot_dirpath)
                attributes_filepath = sim_null_lineages_muller_plot_dirpath + 'attributes.txt'
                fileout_attributes = open(attributes_filepath, "w")
                if start_end_colors:
                    fileout_attributes.write('names\tparents\tcolors\n')
                    start_color = Color(start_end_colors[0])
                    end_color = Color(start_end_colors[1])
                    colors = list(start_color.range_to(end_color, len(muller_lineages_null)))
                    hex_colors = [i.hex_l for i in colors]
                else:
                    fileout_attributes.write('names\tparents\n')
                population_data_table_filepath = sim_null_lineages_muller_plot_dirpath + 'population_data_table.txt'
                fileout_pop_data = open(population_data_table_filepath, "w")
                fileout_pop_data.write('\t%s\n' % '\t'.join([str(i) for i in self.timepoints]))
                for color_index, lineage in enumerate(muller_lineages_null):
                    fileout_attributes.write('%s\t%s' % (lineage, muller_lineages_null[lineage][0]))
                    if start_end_colors:
                        fileout_attributes.write('\t%s' % hex_colors[color_index])
                    fileout_attributes.write('\n')
                    #add 0's to beginning of lineage if need be
                    full_lineage = [0. for i in xrange(muller_lineages_null[lineage][1][0])]
                    for i in muller_lineages_null[lineage][1:]:
                        tpoint_index = i[0]
                        cluster_index = i[1]
                        abundance = seq_clusts_time_series[tpoint_index].data[cluster_index]['other'][self.freq_attribute_name]
                        full_lineage.append(abundance)
                    #add 0's to end of lineage if need be
                    full_lineage += [0 for i in xrange(list_of_timepoint_indices[-1] - muller_lineages_null[lineage][-1][0])]
                    fileout_pop_data.write('%s\t%s\n' % (lineage, '\t'.join([str(i) for i in full_lineage])))
                fileout_attributes.close()
                fileout_pop_data.close()
        #####################################

        return


#############################################
# Below are scripts that are used for running
# jobs on the computational cluster. They
# should only be used by the methods above.
#############################################

def create_lineages_comp_cluster(distance_metric, temp_dirpath, path_to_needle, distance_units, compare_tpoint_to_all_previous, genetic_dist_filepath, output_fasta_dirpath, write_full_original_seqs, muller_plot_output_dirpath, start_end_colors, max_distance, min_cluster_freq, dirpath, count_attribute_name, freq_attribute_name, indiv_seqs_attribute_name, indiv_seq_ids_attribute_name, indiv_seq_counts_attribute_name, indiv_seq_freqs_attribute_name, make_master_network_plot, indiv_sample_network_dirpath, lineage_min_freq_cutoff, mut_count_attribute_name, alignment_method, sim_null_lineages_fasta_dirpath, multiply_num_sim_lineages_by, only_simulate, prob_of_new_lineage_filepath):

    if temp_dirpath == 'None':
        temp_dirpath = None
    if path_to_needle == 'None':
        path_to_needle = None
    if compare_tpoint_to_all_previous == 'False':
        compare_tpoint_to_all_previous = False
    else:
        compare_tpoint_to_all_previous = True
    if genetic_dist_filepath == 'None':
        genetic_dist_filepath = None
    if output_fasta_dirpath == 'None':
        output_fasta_dirpath = None
    if write_full_original_seqs == 'None':
        write_full_original_seqs = None
    if muller_plot_output_dirpath == 'None':
        muller_plot_output_dirpath = None
    if start_end_colors == 'None':
        start_end_colors = None
    else:
        start_end_colors = start_end_colors.split(",")
    if max_distance == 'None':
        max_distance = None
    elif max_distance == 'parsimony':
        pass
    else:
        max_distance = float(max_distance)
    if min_cluster_freq == 'None':
        min_cluster_freq = None
    else:
        min_cluster_freq = float(min_cluster_freq)
    if count_attribute_name == 'None':
        count_attribute_name = None
    if freq_attribute_name == 'None':
        freq_attribute_name = None
    if make_master_network_plot == 'None':
        make_master_network_plot = None
    if indiv_sample_network_dirpath == 'None':
        indiv_sample_network_dirpath = None
    if lineage_min_freq_cutoff == 'None':
        lineage_min_freq_cutoff = None
    else:
        lineage_min_freq_cutoff = float(lineage_min_freq_cutoff)
    if mut_count_attribute_name == 'None':
        mut_count_attribute_name = None
    if sim_null_lineages_fasta_dirpath == 'None':
        sim_null_lineages_fasta_dirpath = None
    multiply_num_sim_lineages_by = int(multiply_num_sim_lineages_by)
    if only_simulate == 'True':
        only_simulate = True
    elif only_simulate == 'False':
        only_simulate = False
    if prob_of_new_lineage_filepath == 'None':
        prob_of_new_lineage_filepath = None

    seq_clusters_set = sequence_clusters_time_series(dirpath=dirpath, count_attribute_name=count_attribute_name, freq_attribute_name=freq_attribute_name, indiv_seqs_attribute_name=indiv_seqs_attribute_name, indiv_seq_ids_attribute_name=indiv_seq_ids_attribute_name, indiv_seq_counts_attribute_name=indiv_seq_counts_attribute_name, indiv_seq_freqs_attribute_name=indiv_seq_freqs_attribute_name)
    seq_clusters_set.create_lineages(distance_metric=distance_metric, temp_dirpath=temp_dirpath, path_to_needle=path_to_needle, distance_units=distance_units, compare_tpoint_to_all_previous=compare_tpoint_to_all_previous, genetic_dist_filepath=genetic_dist_filepath, output_fasta_dirpath=output_fasta_dirpath, write_full_original_seqs=write_full_original_seqs, muller_plot_output_dirpath=muller_plot_output_dirpath, start_end_colors=start_end_colors, max_distance=max_distance, min_cluster_freq=min_cluster_freq, use_comp_cluster=False, lineage_min_freq_cutoff=lineage_min_freq_cutoff, make_master_network_plot=make_master_network_plot, indiv_sample_network_dirpath=indiv_sample_network_dirpath, mut_count_attribute_name=mut_count_attribute_name, alignment_method=alignment_method, sim_null_lineages_fasta_dirpath=sim_null_lineages_fasta_dirpath, multiply_num_sim_lineages_by=multiply_num_sim_lineages_by, only_simulate=only_simulate, prob_of_new_lineage_filepath=prob_of_new_lineage_filepath)

    return


if __name__ == '__main__':
    if sys.argv[1] == 'create_lineages_comp_cluster':
        create_lineages_comp_cluster(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13], sys.argv[14], sys.argv[15], sys.argv[16], sys.argv[17], sys.argv[18], sys.argv[19], sys.argv[20], sys.argv[21], sys.argv[22], sys.argv[23], sys.argv[24], sys.argv[25], sys.argv[26], sys.argv[27], sys.argv[28], sys.argv[29])


#clust_time_series = sequence_clusters_time_series(dirpath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/sequence_clusters/hiv/clustered_by_edit_distance/most_numerous_seq_foreach_cluster/for_lineage_assignment/max_edit_dist_1/7', count_attribute_name='DUPCOUNT', freq_attribute_name='total_freq', indiv_seqs_attribute_name='indiv_seqs', indiv_seq_ids_attribute_name='indiv_seq_ids', indiv_seq_counts_attribute_name='indiv_seq_counts', indiv_seq_freqs_attribute_name='indiv_seq_freqs')
#clust_time_series.create_lineages(distance_metric='dist_between_rep_seqs', temp_dirpath='/Users/nstrauli/data/abr_hiv_coevo/temp_stuff', path_to_needle='/Users/nstrauli/tools/EMBOSS-6.6.0/emboss/needle', distance_units='edit_distance', compare_tpoint_to_all_previous=False, genetic_dist_filepath='/Users/nstrauli/Desktop/gen_dists.txt', output_fasta_dirpath='/Users/nstrauli/Desktop/output_fasta_tmp', write_full_original_seqs='/Users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_fasta_files/7', muller_plot_output_dirpath='/Users/nstrauli/Desktop/muller_plot', start_end_colors=None)#['red', 'purple'])
