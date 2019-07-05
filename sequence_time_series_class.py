import os
from sequence_sample_set_class import sequence_sample_set
from subprocess import Popen, PIPE
import subprocess
import random
from sequence_sample_class import sequence_sample
import re
import tempfile
import numpy
import itertools

class sequence_time_series(sequence_sample_set):
    """
    This class inherits from sequence_sample_set_class.py. It's defining difference is that each of the sequence samples is a time-point, so there is a natural ordering by time, similar to 'immune_time_series_class.py'.

    Input:                                                                                          
    The input for this class should be a directory that contains only fasta formatted files (files named README[.txt] will be ignored). Further the fasta files' headers should be formatted in a specific way. The header for each sequence should start with an ID (could be anything, ideally unique) and then should be followed by a series of attributes. Each attribute should be delimited by a '|'. Each attribute should have a name followed by and '=', followed by the value of the attribute. An example of a well formatted header: ">seq_id|count=65\n". Here, there is only one attribute, which is 'count'. The first entry in the header is assumed to be the sequence unique ID.                                                                         
    The fasta files should also have a specific naming scheme. Each of the filenames should end with a numeric value that indicates the time point, such that if these numeric time point indicators where sorted (numerically), the files would sort in chronological order. This numeric time point indicator in the file names should be delimited by an '_'. For example 'day_56.fasta' or 'month_34.fasta' are examples of properly formatted filenames. '110.fasta' would also work, however 'month_12_day_21.fasta' would not be advisable. This particular case would not crash the code, but would result in the time points not being sorted properly (because the time point is encoded by two numerics (month and day) as opposed to one).
    """
    def __init__(self, dirpath, count_attribute_name):
        if dirpath[-1] != '/':
            dirpath += '/'
        self.dirpath = dirpath
        self.count_attribute_name = count_attribute_name
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

    def cluster_sequences(self, method='vsearch', min_ident_within_cluster=0.90, output_filepath=None, path_to_vsearch=None, temp_dirpath=None, freqs_or_counts='freqs', sort_lineages_by='range', output_network_filepath=None, output_node_attribute_filepath=None, node_attributes_to_add=['count'], names_for_node_attributes_in_header=None, max_edit_distance=1, path_to_needle=None, freq_attribute_name=None):
        """
        This method will cluster all the sequences in all the time-points. BEWARE, in order to do this, all the data must be read into memory, so make sure they aren't super huge files. It clusters the sequences using Vsearch, and also gets the relative size of each cluster as well. Can write this information as expression trajectories for each of the clusters found in the data.
        method - This gives the method that will be used for clustering. Acceptable values are:
            'vsearch' - Default. This means that vsearch will be used to cluster the seqs
            'by_edit_distance' - This means that an in house algorithm of joining seqs to be in the same cluster if they have an edit distance at or below some given threshold. The 'needle' program is used for determining the edit distance between pairs of seqs.
        output_filepath - This give the path to the directory for which the lineage expression trajectories will be written. The format is tab delimited, where the first column gives the identity of the lineage (i.e. Vgene identity, Jgene identity, and centroid sequence, separted by '_'), and the following columns give the expression level of each lineage at each of the respective time-points. If equals None, then no output will be written.
        sort_lineages_by - This tells how the lineages should be sorted when written to file. Exceptable vlaues are:
            'range' - This means the absolute range is used to sort them (in descending order). Range is defined as the tpoint with max expression minus the tpoint with min expression.
            'sum' - This means the lineages will be sorted by the summation of their expression values over the time-course.

        #### Below are parameters only considered if method=='vsearch' ####

        min_ident_within_cluster - This gives the clustering parameter for clustering seqs within a tpoint. See the immune_sample_class.cluster_clones for explanation.
        path_to_vsearch - This gives the path to the vsearch excecutable. If this is None (default) then the path is assumed to be 'vsearch' (i.e. in $PATH).                                          
        temp_dirpath - This gives the path to the directory for which temporary directories will be made within. The temporary directories will be deleted at the end of the script, but 'temp_dirpath' will not.
        freqs_or_counts - This indicates whether the relative frequency (default) or the absolute counts for each lineage should be reported at each of the time-points. If 'freq' then gives relative frequency, if 'count' then gives absolute counts.

        #### Below are parameters only considered if method=='by_edit_distance' ####

        output_network_filepath - If defined (default, None), this will, 1) instruct the method to creat a network file (as a simple interaction file ('.sif')), and 2) will be the path to this file.
        output_node_attribute_filepath - This is the file that will contain all the node attributes for each of the seqs. This is for cytoscape. This parameter is only considered if 'output_network_filepath' is defined.
        node_attributes_to_add - This is a list of sequence attributes to include in the node attribute file. This parameter is only considered if 'output_node_attribute_filepath' is defined. The default is to only add the value for the 'count' attribute (i.e. the value for self.data[some_index]['count']). Acceptable values within this list of strings are:
            any string that is a key for one of the attributes for the sequences. So, could be 'count', 'total_freq', 'timepoint', etc. if those are keys to the dictionarys within self.data[some_index], or self.data[some_index]['other']
            'element_X_of_seq_id' - sometimes info about the seq is encoded in the seq ID, delimited by '_'. This means that the 'X'th element of the seq ID (when delimited by a '_') will be included in the node attributes file.
        names_for_node_attributes_in_header - This is an optional argument, which if defined (default, None), it will give the names for each of the node attributes included in the 'output_node_attribute_filepath' header. This parameter is only considered if 'output_node_attribute_filepath'. If 'output_node_attribute_filepath' is defined and this is not, then the default is to simply use the values in 'node_attributes_to_add'. The names in this list need to be in the same order as the attributes in 'node_attributes_to_add'.
        max_edit_distance - Int. Default, 1. This gives the maximum edit distance between any two seqs, for them to be assigned to the same cluster.
        path_to_needle - This gives the path to the needle executable. This is required if using the computational cluster for an array job because the PATH variable is all messed up when sending jobs to parallele nodes. One can also use this if needle is not in the PATH
        freq_attribute_name - If defined (default, None), this will give the name of the attribute that gives the frequency for each seq in the data set. If it is not defined then frequency is calculated by the seq 'count' / the total, and the freq_attribute_name is set to 'freq'.
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
        
        if method == 'vsearch':
            #write all the data to a temp fasta file
            random_string = str(random.random())[2:]
            temp_fasta_filepath = '%s%s.fasta' % (temp_dirpath, random_string)
            fileout = open(temp_fasta_filepath, "w")
            count = 0
            for i in xrange(len(self.sample_filepaths)):
                sample = sequence_sample(self.sample_filepaths[i], self.count_attribute_name)
                for j in sample.data:
                    count += 1
                    fileout.write('>%s;size=%s;tpoint=%s;freq=%s;\n%s\n' % (count, j['count'], i, j['count']/sample.total, j['seq']))
            fileout.close()
            #cluster all this data
            msa_output_filepath = '%s%s_msa.fasta' % (temp_dirpath, random_string)
            p = Popen([path_to_vsearch, '--cluster_size', temp_fasta_filepath, '--id', str(min_ident_within_cluster), '--sizein', '--sizeout', '--msaout', msa_output_filepath, '--fasta_width', str(0)], stdout=PIPE, stderr=PIPE)
            out_err = p.communicate()
            #parse MSA output to get centroids and counts for each timepoint
            lineage_expr_trajs = {}
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
                        index = centroid_seq
                        #initialize the dic entry
                        lineage_expr_trajs[index] = [0.0 for j in self.timepoints]
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
        
        elif method == 'by_edit_distance':
            random_suffix = str(random.random())
            temp_fasta_filepath = '%stemp_seqfile_%s.fasta' % (temp_dirpath, random_suffix)
            for i in self.sample_filepaths:
                sample = sequence_sample(filepath=i, count_attribute_name=self.count_attribute_name)
                if not freq_attribute_name:
                    sample.add_freq_attribute(freq_attribute_name='freq')
                tpoint = float('.'.join(os.path.basename(i).split('.')[:-1]))
                sample.add_string_to_each_id(str(tpoint)+'_', add_to_start_or_end='start')
                sample.write_full_data_to_disk_fasta(output_filepath=temp_fasta_filepath, append_to_file=True)
            if not freq_attribute_name:
                freq_attribute_name = 'freq'
            sample = sequence_sample(filepath=temp_fasta_filepath, count_attribute_name=self.count_attribute_name)
            print 'Total unique sequences in pooled data:', len(sample.data)
            cluster_dic, connections_dic = sample.cluster_seqs(method='by_edit_distance', temp_dirpath=temp_dirpath, max_edit_distance=max_edit_distance, overwrite_data=False, output_network_filepath=output_network_filepath, output_node_attribute_filepath=output_node_attribute_filepath, node_attributes_to_add=node_attributes_to_add, names_for_node_attributes_in_header=names_for_node_attributes_in_header, path_to_needle=path_to_needle)
            tpoint_to_index_dic = {}
            for index, i in enumerate(self.timepoints):
                tpoint_to_index_dic[i] = index
            lineage_expr_trajs = {}
            for cluster_id in cluster_dic:
                lineage_expr_trajs[cluster_id] = [0.0 for j in self.timepoints]
                for seq_index in cluster_dic[cluster_id]:
                    tpoint = float(sample.data[seq_index]['id'].split('_')[0])
                    traj_index = tpoint_to_index_dic[tpoint]
                    lineage_expr_trajs[cluster_id][traj_index] += float(sample.data[seq_index]['other'][freq_attribute_name])
            subprocess.call(['rm', temp_fasta_filepath])

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

        subprocess.call(['rm', '-r', temp_dirpath])

        return lineage_expr_trajs

    def cluster_over_grid_of_dists(self, path_to_vsearch=None, temp_dirpath=None):
        """
        This method will do the same clustering on the seqs in the time-series data as 'cluster_sequences', but it will do this over a grid of 'min_ident_within_cluster' values. The purpose of this is to figure out what a good distance is for a given sequence set. It will return the number of clusters for each of the distance parameters. The idea is that if there is some natural structure to the population then there will be a plateu where the number of cluster stays relatively consistent over a range of distance values.
        path_to_vsearch - This gives the path to the vsearch excecutable. If this is None (default) then the path is assumed to be 'vsearch' (i.e. in $PATH).
        temp_dirpath - This gives the path to the directory for which temporary directories will be made within. The temporary directories will be deleted at the end of the script, but 'temp_dirpath' will not.

        OUTPUT:
        Returns 'num_clusters' which is a list that represents the number of clusters found for each of the distance parameters. Also returns 'ident_values' which is a list of each of the identity values that were used to cluster, in ascending order, and order corresponds to 'num_clusters'. Also returns 'mean_cluster_sizes' which is the mean number of unique sequences across the clusters, for each of the ident_values.
        """

        ######## important parameters ########
        #this gives the grid for ident values
        ident_values = numpy.linspace(.8, 1, num=101)
        ######################################

        if temp_dirpath == None:
            temp_dirpath = tempfile.mkdtemp(dir=os.getcwd()) + '/'
        else:
            temp_dirpath = tempfile.mkdtemp(dir=temp_dirpath) + '/'
        num_clusters = []
        mean_cluster_sizes = []
        for i in ident_values:
            print 'clustering with min identity within cluster parameter:', i
            #write all the data to a temp fasta file
            random_string = str(random.random())[2:]
            temp_fasta_filepath = '%s%s.fasta' % (temp_dirpath, random_string)
            fileout = open(temp_fasta_filepath, "w")
            count = 0
            for j in xrange(len(self.sample_filepaths)):
                sample = sequence_sample(self.sample_filepaths[j], self.count_attribute_name)
                for k in sample.data:
                    count += 1
                    fileout.write('>%s;size=%s;tpoint=%s;freq=%s;\n%s\n' % (count, k['count'], j, k['count']/sample.total, k['seq']))
            fileout.close()
            #cluster all this data
            msa_output_filepath = '%s%s_msa.fasta' % (temp_dirpath, random_string)
            p = Popen([path_to_vsearch, '--cluster_size', temp_fasta_filepath, '--id', str(i), '--sizein', '--sizeout', '--msaout', msa_output_filepath, '--fasta_width', str(0)], stdout=PIPE, stderr=PIPE)
            out_err = p.communicate()
            #parse MSA output to get centroids and counts for each timepoint
            cluster_sizes = []
            cluster_count = 0
            filein = open(msa_output_filepath, "r")
            while True:
                line = filein.readline()
                #if EOF, break loop
                if not line:
                    break
                elif line[0] == '>':
                    #if this is the centroid sequence
                    if line[1] == '*':
                        cluster_count += 1
                        cluster_size = 1
                    #else if this is a member seq of the cluster
                    elif line[:-1] != '>consensus':
                        cluster_size += 1
                    #else if this is the end of the cluster seqs
                    elif line[:-1] == '>consensus':
                        cluster_sizes.append(cluster_size)
            filein.close()
            #delete temp files
            subprocess.call(['rm', temp_fasta_filepath, msa_output_filepath])
            num_clusters.append(cluster_count)
            mean_cluster_sizes.append(float(sum(cluster_sizes)) / len(cluster_sizes))
        subprocess.call(['rm', '-r', temp_dirpath])
        return num_clusters, mean_cluster_sizes, ident_values

    def make_MSA(self, output_filepath, temp_dirpath=None, method='mafft', path_to_mafft=None, add_ref_seq=None, ref_seq_name='ref_seq', write_seq_id_first=None):
        """
        This method takes all the sequences in the time series and creates one large multiple sequence alignment (MSA) of all the seqs in the data
        output_filepath - This gives the path for which the output MSA will be saved to.
        temp_dirpath - This gives the path to the directory that temporary data will be written. If None, then the temp_dirpath will be the directory above output_filepath.
        method - This gives the alignment tool that will be used to make the MSA's. Acceptable values are:
            'mafft' - This means that MAFFT will be used.
        path_to_mafft - If defined (default = None), then this gives the path to the MAFFT excecutable. If None, then this method will assume that mafft is in the PATH.
        add_ref_seq - If defined (default, None), this should be the sequence of a reference sequence that will be added to the top of the alignment. for example, if the data is HIV seq data, then one might want to add the HXB2 sequence to the alignment. This ref seq will be included with the input to the aligner. If this is defined and has a '/' in it (i.e. it is a filepath), then this will be treated as a filepath. If this is the case then the reference sequence and ref seq name will be retrieved from this file. Needs to be in fasta format.
        ref_seq_name - This give the name that will be used in the fasta header for the reference seq. This parameter is only considered if 'add_ref_seq' is defined. Default is 'ref_seq'
        write_seq_id_first - Sometimes we want a certain seq to be 1st in an MSA. If defined (default, None), this will be a string that gives the sequence ID for the seq that we want at the top of the MSA. The seq ID needs to equal the seq ID that will be in the MSA, so if one wants the 3rd seq of timepoint 47 to be first in the MSA, then write_seq_id_first='47_3' (assuming the seq IDs in the input fasta files are sequential ints).

        NOTE:
        This method has to load all the sequence data (for the entire time series) into memory, and then run MAFFT on that. So, caution must be excercised if the dataset is quite large.
        """
        if temp_dirpath:
            if temp_dirpath[-1] != '/':
                temp_dirpath += '/'
            if not os.path.exists(temp_dirpath):
                os.makedirs(temp_dirpath)
        else:
            temp_dirpath = os.path.dirname(output_filepath) + '/'
        random_suffix = str(random.random())
        if not path_to_mafft:
            path_to_mafft = 'mafft'
        temp_allseqs_fasta = "%sall_seqs_%s.fasta" % (temp_dirpath, random_suffix)
        fileout = open(temp_allseqs_fasta, "w")
        if add_ref_seq:
            if '/' in add_ref_seq:
                filein = open(add_ref_seq, "r")
                ref_seq_name = filein.readline()[1:-1]
                add_ref_seq = filein.readline()[:-1]
                filein.close()
            fileout.write('>%s\n%s\n' % (ref_seq_name, add_ref_seq))
        total_num_seqs = 0
        for tpoint, file in itertools.izip(self.timepoints, self.sample_filepaths):
            filein = open(file, "r")
            for j in filein:
                if j[0] == '>':
                    total_num_seqs += 1
                    fileout.write('>%s_%s' % (tpoint, j[1:]))
                    if write_seq_id_first:
                        seq_id = '%s_%s' % (tpoint, j[1:].split('|')[0])
                        if seq_id == write_seq_id_first:
                            first_seq_header = '>%s_%s' % (tpoint, j[1:])
                else:
                    fileout.write(j)
            filein.close()
        fileout.close()
        #make sure there is more than one seq in the input fasta file
        if total_num_seqs < 2:
            print 'The input data only has one sequence. Aborting.'
            return
        p = Popen(['bash', 'call_mafft.bash', temp_allseqs_fasta, output_filepath+'_unformatted', path_to_mafft], stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        print out
        print err
        #if a certain seq is desired to be first, get that seq
        if write_seq_id_first:
            first_seq_found = False
            first_seq = ''
            filein = open(output_filepath+'_unformatted', "r")
            for i in filein:
                if i[0] == '>':
                    if not first_seq_found:
                        seq_id = i[1:-1].split('|')[0]
                        if seq_id == write_seq_id_first:
                            first_seq_found = True
                    else:
                        break
                elif first_seq_found:
                    first_seq += i[:-1]
            filein.close()
        #format the fasta file correctly
        filein = open(output_filepath+'_unformatted', "r")
        fileout = open(output_filepath, "w")
        filein_temp_allseqs = open(temp_allseqs_fasta, "r")

        if write_seq_id_first:
            fileout.write('%s%s' % (first_seq_header, first_seq))
            on_desired_first_seq = False
            for i in filein:
                if i[0] == '>':
                    header = filein_temp_allseqs.readline()
                    filein_temp_allseqs.readline()#burn non aligned seq line
                    seq_id = header[1:].split('|')[0]
                    if seq_id == write_seq_id_first:
                        on_desired_first_seq = True
                    else:
                        if on_desired_first_seq:
                            on_desired_first_seq = False
                        fileout.write('\n' + header)
                else:
                    if on_desired_first_seq:
                        pass
                    else:
                        fileout.write(i[:-1])
        else:
            fileout.write(filein_temp_allseqs.readline())
            filein.readline()
            for i in filein:
                if i[0] == '>':
                    filein_temp_allseqs.readline()
                    #mafft cuts off long seq headers, so need the original seq headers
                    fileout.write('\n' + filein_temp_allseqs.readline())
                else:
                    fileout.write(i[:-1])
        fileout.write('\n')
        filein.close()
        filein_temp_allseqs.close()
        fileout.close()
        subprocess.call(['rm', temp_allseqs_fasta])
        subprocess.call(['rm', output_filepath+'_unformatted'])
        return

    def get_first_tpoint_most_abundant_seq(self, output_filepath=None, translate_by_ref_seq=False, temp_dirpath=None, append_to_output=False):
        """
        This method returns the most abundant seq from the first time-point.
        output_filepath - If defined (default, None) then this is the path to the output file that will be written for the most abundant seq. Will be in fasta format.
        translate_by_ref_seq - If defined (default is False) this should be the path to a fasta formatted reference sequence for the seqs in the alignment. This reference seq needs to be annotated such that the coding frame is known. Meaning, it is known where to begin translating the ref seq. The annotation for the coding frame information must be in the header of the ref seq and should look like:'...|coding_frame_start=[X]|...', where 'X' gives the coding frame start position. The most abundant seq is then aligned to this ref seq so that the coding frame can be mapped from the ref seq to the most abundant seq. The alignment is done using 'needle' from the EMBOSS package. This needs to be in the PATH. The resulting translation is recorded in the 'output_filepath', so this should be defined too. If 'translate_by_ref_seq' equals False, then no translation occurs.
        temp_dirpath - Gives the path to a directory for which temp data can be written. If None (default), then uses the master directory to 'output_filepath', if there is not output_filepath given then uses the current working directory.
        append_to_output - If True (default, False) then this will cause the output file to be appended to, rather than over-written. This is only considered if 'output_filepath' is defined.
        """
        sample = sequence_sample(filepath=self.sample_filepaths[0], count_attribute_name=self.count_attribute_name)
        most_abun_seq = sample.get_most_abundant_seq(output_filepath=output_filepath, translate_by_ref_seq=translate_by_ref_seq, temp_dirpath=temp_dirpath, append_to_output=append_to_output)
        return most_abun_seq

    def get_most_abundant_seq_from_tpoint(self, tpoint_index, output_filepath=None, translate_by_ref_seq=False, temp_dirpath=None, append_to_output=False, path_to_needle='needle'):
        """
        This method does the same thing as that of 'get_first_tpoint_most_abundant_seq', but it will get the most abundant seq for any desired time-point.
        tpoint_index - Int. This gives the index (indexed at 0) for the query time-point for which the mst abundant seq is desired. For example, if one wants the most abun. seq from the 2nd time-point this would be 1.
        append_to_output - If True (default, False) then this will cause the output file to be appended to, rather than over-written. This is only considered if 'output_filepath' is defined.
        path_to_needle - String. This give the path to the 'needle' executable. The default is simply 'needle', which means that it is already in your PATH.
        """
        sample = sequence_sample(filepath=self.sample_filepaths[tpoint_index], count_attribute_name=self.count_attribute_name)
        most_abun_seq = sample.get_most_abundant_seq(output_filepath=output_filepath, translate_by_ref_seq=translate_by_ref_seq, temp_dirpath=temp_dirpath, append_to_output=append_to_output, path_to_needle=path_to_needle)
        return most_abun_seq

    def get_time_series_sample_counts(self, return_dic=False):
        """
        This gets the sample counts for each of the sequence samples. Returns a list of floats where each float gives the total seqs in a given sample (in chronological order).
        return_dic - If True (default, False), this will return a dictionary, where the index a the time-point (as a float) and the def is the sample count (as a float).
        """
        sample_counts = []
        sample_counts_dic = {}
        for index, i in enumerate(self.sample_filepaths):
            sample = sequence_sample(filepath=i, count_attribute_name=self.count_attribute_name)
            sample_counts.append(sample.total)
            sample_counts_dic[self.timepoints[index]] = sample.total
        if return_dic:
            return sample_counts_dic
        else:
            return sample_counts

    def write_most_abundant_seq_foreach_tpoint(self, output_dirpath, write_attribute_value_instead=None):
        """
        This method finds the most abundant sequence in each of the time-points, and writes them (as separate files) to disk.
        output_dirpath - This give the path to the directory that will contain the most abundant seqs.
        write_attribute_value_instead - If defined (default, None), this will write the value of the provided attribute. If defined this must be the name of an attribute in the seq data, and the value of this attribute should be a sequence of some sort. This is essentially a hack to allow the writing of CDR3 seq data.
        """
        for tpoint_filepath in self.sample_filepaths:
            filename = '.'.join(os.path.basename(tpoint_filepath).split('.')[:-1])
            output_filepath = '%s%s.fasta' % (output_dirpath, filename)
            sample = sequence_sample(filepath=tpoint_filepath, count_attribute_name=self.count_attribute_name)
            most_abun_index = sample.get_most_abundant_seq(output_filepath=None, translate_by_ref_seq=False, temp_dirpath=None, append_to_output=False, return_seq_index=True)
            sample.write_subset_of_seqs_to_disk(seq_indices=[most_abun_index], output_filepath=output_filepath, append_to_file=False, seq_ids=None, add_string_to_each_id_being_written=None, write_attribute_value_instead=write_attribute_value_instead)
        return

    def remove_timepoints(self, timepoint_indeces):
        """
        This script will remove a given list of time-points from the data set.
        timepoint_indeces - List of ints. This gives a list of the indeces of the time-points that one wants to remove.
        """
        for i in sorted(timepoint_indeces, reverse=True):
            del self.sample_filepaths[i]
            del self.timepoints[i]
        return

        
#sample_set = sequence_time_series('/Users/nstrauli/data/abr_hiv_coevo/seq_data/hiv_fasta_files/1', 'DUPCOUNT')
#sample_set.make_MSA('/Users/nstrauli/data/abr_hiv_coevo/seq_data/alignments/hiv/1.fasta')
