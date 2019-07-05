#!/usr/bin/python
#$ -S /usr/bin/python
#$ -e error
#$ -o out3
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=200G
#$ -l h_rt=336:00:00

import sys
sys.path.insert(0, './')
import numpy
import os
from sequence_sample_class import sequence_sample
import tempfile
from subprocess import Popen, PIPE
import subprocess
from mixcr_align_output_class import mixcr_align_output
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import random
import time
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

class immune_sample(sequence_sample):
    """
    A class for the sequence data in a single sample of immune repertoire data.

    Attributes:
    filepath - path to the fasta formatted file that contains the sequence information for each clone/sequence entity. Can also be fastq formatted.
    count_attribute_name - this gives the name that corresponds to the count of each clone/sequence entity. This count should be given in the header of each sequence entry in the fasta file. The format of the header should be such that each attribute given in the header is delimited by a '|', and should consist of the name of the attribute followed by and equals sign ('=') and then the value of the attribute. An example of a well formatted header: ">seq_id|vgene=IGHV1|dgene=NA|jgene=IGHJ3|count=65\n". Here the attributes are 'vgene', 'dgene', 'jgene', and 'count'. The first entry in the header is assumed to be the sequence unique ID.
    vgene_name - This gives the name of the attribute that gives the V gene name.
    dgene_name - This gives the name of the attribute that gives the D gene name.
    jgene_name - This gives the name of the attribute that gives the J gene name.
    cdr3_name - This gives the name of the attribute that gives the CDR3 sequence.
    data - A list that contains all the data in the fasta file. Each entry in the list corresponds to one sequence entry in the fasta file. This is a list of dics, where each dic contains a number of essential bits of annotation for the sequence. The essential annotations are: 'vgene', 'dgene', 'jgene', 'cdr3', 'count', 'seq', and 'other'. The 'other' key in the dictionary contains all the remaining information (if any) that are in the headers of the sequences. The definition of 'other' is itself a dictionary that contains the remaining annotation information.
    total = gives the total number of sequences represented by the data. That is, taking into account 'count' information for each of the seq entries.
    clones - This is a list of clonal seqs. Each entry contains the information for one unique clone. Clones can be defined in a variety of different ways. This attribute is not defined until the 'group_clones' method is run. See this method for more details.
    clone_define - This tells how a clone is defined. Initially this attributes equals None. This attribute can be defined by setting the 'define' parameter in the 'group_clones' method. Acceptable values are:
        'v_j_cdr3' - Default. This means that the V and J gene segment identity and CDR3 sequence define a unique clone.
        'cdr3' - This means that only the CDR3 sequence defines a unique clone.
    ignor_allele_info - Boolean. Tells whether or not the allele information in the gene names is ignored. Default=True
    counts - a list of counts for each of the sequences.
    vgene_counts - This is a dictionary, where each index is an observed V gene and the data, and it is defined by the count (number of reads that mapped to that V gene).
    dgene_counts - This is a dictionary, where each index is an observed D gene and the data, and it is defined by the count (number of reads that mapped to that D gene).
    jgene_counts - This is a dictionary, where each index is an observed J gene and the data, and it is defined by the count (number of reads that mapped to that J gene).
    filetype - This is either 'fasta' or 'fastq'. It is determined by the suffix of the input file name.
    replace_slash_with - This gives the character that is used to replace a slash in an attribute value. Slashes can cause problems if attribute values are used in the creation of output files in some of the methods below. The default is None, which means that the attribute values are kept as is. This only applies to the values of V, D, and J gene segment names.
    """
    def __init__(self, filepath, count_attribute_name, vgene_name, dgene_name, jgene_name, cdr3_name, ignor_allele_info=True, replace_slash_with=None):
        self.filepath = filepath
        suffix = os.path.basename(self.filepath).split('.')[-1]
        if suffix == 'fasta' or suffix == 'fa':
            self.filetype = 'fasta'
        elif suffix == 'fastq' or suffix == 'fq':
            self.filetype = 'fastq'
        else:
            print 'Unrecognized file type. Must be either fasta of fastq. Check file suffix.'
            return
        self.filepath = filepath
        self.count_attribute_name = count_attribute_name
        self.vgene_name = vgene_name
        self.dgene_name = dgene_name
        self.jgene_name = jgene_name
        self.cdr3_name = cdr3_name
        filein = open(filepath, "r")
        self.data = []
        self.counts = []
        self.total = 0.
        self.clones = None
        self.clone_define = None
        self.ignor_allele_info = ignor_allele_info
        self.vgene_counts = {}
        self.dgene_counts = {}
        self.jgene_counts = {}
        self.replace_slash_with = replace_slash_with
        count = 0
        for i in filein:
            count += 1
            #if the line is the seq header, then fill in the header info
            if count == 1:
                seq_entry_dic = {'vgene':'NA', 'dgene':'NA', 'jgene':'NA', 'cdr3':'NA', 'count':1, 'seq':'NA', 'other':{}}
                line = i[:-1].split('|')
                seq_entry_dic['id'] = line[0][1:]
                found_count = False
                for j in line[1:]:
                    attribute = j.split('=')
                    name = attribute[0]
                    value = attribute[1]
                    if name == vgene_name:
                        if self.ignor_allele_info == True:
                            value = value.split('*')[0]
                        if self.replace_slash_with:
                            if '/' in value:
                                value = value.replace('/', self.replace_slash_with)
                        seq_entry_dic['vgene'] = value
                    elif name == dgene_name:
                        if self.ignor_allele_info == True:
                            value = value.split('*')[0]
                        if self.replace_slash_with:
                            if '/' in value:
                                value = value.replace('/', self.replace_slash_with)
                        seq_entry_dic['dgene'] = value
                    elif name == jgene_name:
                        if self.ignor_allele_info == True:
                            value = value.split('*')[0]
                        if self.replace_slash_with:
                            if '/' in value:
                                value = value.replace('/', self.replace_slash_with)
                        seq_entry_dic['jgene'] = value
                    elif name == cdr3_name:
                        seq_entry_dic['cdr3'] = value
                    elif name == count_attribute_name:
                        seq_entry_dic['count'] = int(value)
                        self.total += int(value)
                        self.counts.append(int(value))
                        found_count = True
                    else:
                        seq_entry_dic['other'][name] = value
                if found_count == False:
                    self.counts.append(1)
                    self.total += 1
                #fill in the gene counts
                try:
                    self.vgene_counts[seq_entry_dic['vgene']] += seq_entry_dic['count']
                except KeyError:
                    self.vgene_counts[seq_entry_dic['vgene']] = seq_entry_dic['count']
                try:
                    self.dgene_counts[seq_entry_dic['dgene']] += seq_entry_dic['count']
                except KeyError:
                    self.dgene_counts[seq_entry_dic['dgene']] = seq_entry_dic['count']
                try:
                    self.jgene_counts[seq_entry_dic['jgene']] += seq_entry_dic['count']
                except KeyError:
                    self.jgene_counts[seq_entry_dic['jgene']] = seq_entry_dic['count']
            #now fill in the seq info
            elif count == 2:
                seq_entry_dic['seq'] = i[:-1]
                self.data.append(seq_entry_dic)
                if self.filetype == 'fasta':
                    count = 0
            elif count == 4:
                count = 0
        filein.close()
        return
    
    def get_gene_expression(self, gene_class, normailze_by_gene_length=False, drop_allele_info=True, gene_segment_reference_dirpath='/Users/nstrauli/data/ref_seqs/Ig/from_imgt_appended/human/'):
        """
        This method uses the data in 'self.data' to get the expression level of each possible gene segment. That includes gene segments that aren't in the data, so a path to a file that contains all the reference gene segments is hard coded in this script. Epression is defined as the counts for each gene segment, normalized by the total counts in the data. User also the option to normalize by gene length as well (this should be done if the data comes from RNAseq). If this option is selected, then the value is multiplied by 100 as well (so that number aren't so small, similar to FPKM)

        Input:
        gene_class = This tells the script for what gene class expression levels should be calculated. If 'gene_class' is a list then the script will calculate expression levels for each gene class in the list. Exceptable values are 'IGHV', 'IGHD', IGHJ', 'IGLV', 'IGLJ', 'IGKV', 'IGKJ', 'TRAV', 'TRAJ', 'TRBV', 'TRBD', 'TRBJ', 'TRDV', 'TRDJ', 'TRDD', 'TRGV', or 'TRGJ'.
        normailze_by_gene_length = Boolean. Default=False. Instructs the script if each gene segments expression level should be normalized by gene length, in addition to the total number of reads.
        drop_allele_info = Boolean. Default=True. If true, the allele information in the gene names is disregarded. So IGHV4-7*01 will be identical to IGHV4-7*05.
        gene_segment_reference_dirpath = the path to the directory that contains the fasta files for each of the gene classes. The fasta files will have all the germline genes (for the given class) and their sequences. However, BEWARE, the default directory has some pseudo genes listed that don't have accompanied sequences. They just have 'NA' where the sequence should be. This is fine if one simply needs a list of gene names, but not good if one needs to normalize by the length of the gene sequence (i.e. when using RNAseq data). So, if using RNAseq data then this filepath should be changed.
        Output:
        gene_expression = If 'gene_class' is a list then this is a list of dictionarys. If 'gene_class' is a singular string, then this is a dictionary. Each element of the dictionaries are the names of gene segments, and their definitions are that respective gene's expression level. If 'gene_class' is a list, each dictionary in 'gene_expression' corresponds to the each gene class listed in 'gene_class'. Further, the dictionaries are in the same respective order as 'gene_class'.
        """

        #check if gene_class is a string (i.e. singular)
        if isinstance(gene_class, basestring):
            #if it is, make it a list (one element long)
            gene_class = [gene_class]

        #get reference sequence info
        reference_info = []
        for i in gene_class:
            filein = open("%s%s.fasta" % (gene_segment_reference_dirpath, i), "r")
            gene_lens_dic = {}
            for j in filein:
                if j[0] == '>':
                    gene_name = j[1:-1]
                    if drop_allele_info:
                        gene_name = gene_name.split('*')[0]
                else:
                    length = float(len(j[:-1]))
                    try:
                        gene_lens_dic[gene_name].append(length)
                    except KeyError:
                        gene_lens_dic[gene_name] = [length]
            filein.close()
            for j in gene_lens_dic:
                if drop_allele_info:
                    gene_lens_dic[j] = round(sum(gene_lens_dic[j]) / len(gene_lens_dic[j]), 0)
                else:
                    gene_lens_dic[j] = gene_lens_dic[j][0]
            reference_info.append(gene_lens_dic)

        #get counts for each gene segment in the data
        #first initialize the dic
        gene_expression = []
        for i in reference_info:
            gene_expression.append({})
            for j in i:
                gene_expression[-1][j] = 0.
        #now update the gene expression dic with counts
        #while simultaneously normalizing by self.total
        gene_class_dic = {'IGHV':'vgene', 'IGHD':'dgene', 'IGHJ':'jgene', 'IGLV':'vgene', 'IGLJ':'jgene', 'IGKV':'vgene', 'IGKJ':'jgene', 'TRAV':'vgene', 'TRAJ':'jgene', 'TRBV':'vgene', 'TRBD':'dgene', 'TRBJ':'jgene', 'TRDV':'vgene', 'TRDJ':'jgene', 'TRDD':'dgene', 'TRGV':'vgene', 'TRGJ':'jgene'}
        for i in self.data:
            for j,k in enumerate(gene_class):
                query_gene_class = gene_class_dic[k]
                query_gene_name = i[query_gene_class]
                if drop_allele_info:
                    query_gene_name = query_gene_name.split('*')[0]
                query_read_counts = i['count']
                query_expression = query_read_counts / self.total
                if normailze_by_gene_length:
                    query_expression = (query_expression / reference_info[j][query_gene_name]) * 1000
                gene_expression[j][query_gene_name] += query_expression
        if len(gene_expression) == 1:
            gene_expression = gene_expression[0]
        return gene_expression

    def get_gene_counts(self, drop_allele_info=True):
        """
        This method is similar to 'get_gene_expression', excpet it ignores the reference info for gene segments. That is it does not reort the expression for all genes possible, it only gets the counts (number of reads) that map to a given germline gene. If a gene is not observed in the data, then it is not reported. Returns this information for both V and J genes. Running this method will define the attributes: 'v_gene_counts' and 'j_gene_counts'.
        'drop_allele_info'
        """

    def get_clone_expression(self, define='v_j_cdr3'):
        """
        This method functions similarily to 'get_gene_expression' in that it gets the expression level of antibody elements over time. Here, we are getting the expression levels of clones over time. The expression level of a given clone is defined as the counts for that clone divided by the total number of counts (i.e. read counts) in the sample.
        define - This tells the script how a clone is defined. Exceptable values are:
            'v_j_cdr3' - Default. This means that the V and J gene segment identity and CDR3 sequence define a unique clone.
            'cdr3' - This means that only the CDR3 sequence defines a unique clone.
            'v_j' - This means that the V and J gene segment identity defines a clone.

        Output:
        'clone_expression' - A dictionary. Each index represents a unique cloneand is coded as the V gene name, the J gene name, and the CDR3 seq (delimited by a '_'). The definition of that index is the expression level of the clone.
        """
        clone_expression = {}
        for i in self.data:
            vgene = i['vgene']
            jgene = i['jgene']
            cdr3 = i['cdr3']
            count = i['count']
            expression = count / self.total
            if define == 'v_j_cdr3':
                index = '_'.join([vgene, jgene, cdr3])
            elif define == 'cdr3':
                index = cdr3
            elif define == 'v_j':
                index = '%s_%s' % (vgene, jgene)
            try:
                clone_expression[index] += expression
            except KeyError:
                clone_expression[index] = expression
        return clone_expression

    def get_clonal_freq_spectrum(self, define='v_j_cdr3', use_counts_or_freqs='counts', freq_bin_width=0.00001):
        """
        This method gets the frequency spectrum of clones in the immune sample. A frequency spectrum gives the proportion of clones that are singletons (i.e. only occur once), the proportion of clones that are doubletons (i.e. only occur twice), tripletons, etc...
        
        Input:
        define - This tells the script how a clone is defined. Acceptable values are:
            'v_j_cdr3' - Default. This means that the V and J gene segment identity and CDR3 sequence define a unique clone.
            'cdr3' - This means that only the CDR3 sequence defines a unique clone.
        
        use_counts_or_freqs - This tells the script if the bins of the frequency spectrum will be counts or frequencies. acceptable values are:
            'counts' - Default. This means that each element of the resulting freq spectrum can be interpereted as "the proportion of clones that have a count of X"
            'freqs' - This means that each element of the resulting freq spectrum can be interpereted as "the proportion of clones that have a frequency of X"
        freq_bin_width - This gives the width of a frequency bin in the spectrum. I.e. this determines how far apart bin edges are. This is ignored if use_counts_or_freqs='counts'.
        
        Output:
        freq_spectrum - A list. Each element corresponds to a frequency bin. The first element in the list is the proportion of singleton clones, the second is the proportion of doubleton clones, ..., the last element is the proportion of N-ton clones (where N is the total count of immune reads in the sample). So, freq_spectrum has a length of N.
        """

        #make sure 'use_counts_or_freqs' variable is defined correctly
        if use_counts_or_freqs != 'counts' and use_counts_or_freqs != 'freqs':
            print 'use_counts_or_freqs variable mast be defined as either "counts" or "freqs":', use_counts_or_freqs
            return
        #now get the counts for each unique clone in the data
        clone_counts = {}
        for i in self.data:
            vgene = i['vgene']
            jgene = i['jgene']
            cdr3 = i['cdr3']
            count = i['count']
            if define == 'v_j_cdr3':
                index = '_'.join([vgene, jgene, cdr3])
            elif define == 'cdr3':
                index = cdr3
            try:
                clone_counts[index] += count
            except KeyError:
                clone_counts[index] = count
        #transform clone_counts into a numeric list of counts.
        #don't need clone identity info anymore
        if use_counts_or_freqs == 'freqs':
            clone_counts = [clone_counts[i]/self.total for i in clone_counts]
            #this gives the bin edges for the distribution
            bin_edges = numpy.arange(0, max(clone_counts)+freq_bin_width, freq_bin_width)
            freq_spectrum = numpy.histogram(clone_counts, bin_edges)
        else:
            clone_counts = [clone_counts[i] for i in clone_counts]
            bin_edges = range(1, max(clone_counts)+2)
            freq_spectrum = numpy.histogram(clone_counts, bin_edges)
        #convert array types to list types
        #remove the first bin edge, as this is the 
        #beginning of the first bin. We only report the
        #right side of each edge, exclusive.
        bin_edges = freq_spectrum[1][1:].tolist()
        freq_spectrum = freq_spectrum[0].tolist()
        #was having trouble with the 'density' parameter in numpy.histogram, so
        #decided to normalize manually
        total_unique_clones = float(sum(freq_spectrum))
        freq_spectrum = [i/total_unique_clones for i in freq_spectrum]
        return freq_spectrum, bin_edges

    def group_clones(self, define='v_j_cdr3'):
        """
        This method groups all the sequences in the sample based on some definition of 'clonality'.
        define - This tells the script how a clone is defined. Acceptable values are:
            'v_j_cdr3' - Default. This means that the V and J gene segment identity and CDR3 sequence define a unique clone.
            'cdr3' - This means that only the CDR3 sequence defines a unique clone.
        """
        #make sure 'define' is defined correctly
        if define != 'v_j_cdr3' and define != 'cdr3':
            print 'The "define" parameter does not have an acceptable value:', define
            return

        self.clone_define = define

        #first gather unique clones, based on 'define'
        clone_counts = {}
        for i in self.data:
            #set indices of dic based on the value of 'define'
            if define == 'v_j_cdr3':
                if i['vgene'] == 'NA' or i['jgene'] == 'NA' or i['cdr3'] == 'NA':
                    continue
                index = '_'.join([i['vgene'], i['jgene'], i['cdr3']])
            elif define == 'cdr3':
                if i['cdr3'] == 'NA':
                    continue
                index = i['cdr3']
            #now populate clone_counts dic
            try:
                clone_counts[index]['count'] += i['count']
            except KeyError:
                if define == 'v_j_cdr3':
                    clone_counts[index] = {'count':i['count'], 'vgene':i['vgene'], 'jgene':i['jgene'], 'cdr3':i['cdr3'], 'id':index}
                elif define == 'cdr3':
                    clone_counts[index] = {'count':i['count'], 'cdr3':i['cdr3'], 'id':index}
        #now turn 'clone_counts' into a list of dics
        self.clones = [clone_counts[i] for i in clone_counts]
        return

    def cluster_clones(self, min_ident_within_cluster=0.97, path_to_vsearch=None, temp_dirpath=None):
        """
        This method will first group all the sequences by V and J gene identity, and then cluster the sequences within each of those groups. It uses the vsearch program to do the clustering.
        min_ident_within_cluster - This is an important parameter that gives the distance allowed between the sequences within a cluster. It should be between 0 and 1. If a target sequence, when aligned to a cluster centroid, has an identity lower then this value, then it is not added to that cluster.
        path_to_vsearch - This gives the path to the vsearch excecutable. If this is None (default) then the path is assumed to be 'vsearch' (i.e. in $PATH).
        'temp_dirpath' = This gives the path to the directory for which a temporary output directory will be created in. It will contain temp fasta files and output from vsearch. The temp directory that is created within 'temp_dirpath' will be removed at the end of the script, but the 'temp_dirpath' directory itself will not be removed. If this equals None, the the current working directory is used.

        OUTPUT:
        Returns 'centroid_dic'. This is a dictionary where its index is a unique V and J gene pair found in the data. These V and J gene pairs form the criteria for grouping together the sequences within the data. Withiin these groups, clustering takes place. So, in each of the indeces of 'centroid_dic' are list of centroids, where each centroid is the representative sequence (the centroid) of a cluster, along with the counts of sequences in that cluster. Each index in these lists of centroids is a two element list, where the first element corresponds to the counts of seq in that cluster, and the second element gives the actual sequence of the centroid.
        """
        if temp_dirpath == None:
            temp_dirpath = tempfile.mkdtemp(dir=os.getcwd()) + '/'
        else:
            temp_dirpath = tempfile.mkdtemp(dir=temp_dirpath) + '/'

        #fisrt group sequences by V and J gene identity
        seq_groups = {}
        for i in self.data:
            index = '%s_%s' % (i['vgene'], i['jgene'])
            try:
                seq_groups[index].append([i['count'], i['seq']])
            except KeyError:
                seq_groups[index] = [[i['count'], i['seq']]]
        #this will be a dictionary, where each index is a sequence group.
        #The sequence groups are defined by unique V and J gene identity.
        #The definition of each index are lists of centroid sequences
        #with there respective counts, and relative freqs
        #(i.e. three element lists).
        centroid_dic = {}
        #for each sequence group
        for i in seq_groups:
            centroid_dic[i] = []
            #make a temp fasta file
            temp_fasta_filepath = '%s%s.fasta' % (temp_dirpath, i)
            fileout = open(temp_fasta_filepath, "w")
            #for each sequence in the group
            count = 0
            for j in seq_groups[i]:
                count += 1
                #write to temp fasta file
                fileout.write('>%s;size=%s;\n%s\n' % (count, j[0], j[1]))
            fileout.close()
            #now cluster the seqs
            centroid_filepath = '%s%s_centroids.fasta' % (temp_dirpath, i)
            p = Popen([path_to_vsearch, '--cluster_size', temp_fasta_filepath, '--id', str(min_ident_within_cluster), '--sizein', '--sizeout', '--centroids', centroid_filepath, '--fasta_width', str(0)], stdout=PIPE, stderr=PIPE)
            out_err = p.communicate()
            #retrieve centroid seqs for each cluster
            filein = open(centroid_filepath, "r")
            id = 0
            for j in filein:
                if j[0] == '>':
                    id += 1
                    line = j[:-1].split(';')
                    count = int(line[1].split('=')[1])
                    freq = count / self.total
                else:
                    centroid_dic[i].append([count, freq, j[:-1]])
            filein.close()
            subprocess.call(['rm', centroid_filepath, temp_fasta_filepath])
        subprocess.call(['rm', '-r', temp_dirpath])
        return centroid_dic

    def write_clone_seqs(self, output_filepath, sort_by_count=False):
        """
        This method writes the clonal sequences (as well as all info necessary to define them) to the provided filepath ('output_filepath'), in fasta format.
        'output_filepath' - The path to the file that will be written
        'sort_by_count' - If True (default False), then will sort the clones by their counts in decending order before writing.
        """
        #make sure clones have been compiled 1st
        if self.clones == None:
            print 'The clonal sequences have not been gathered yet from the sample. Run the "group_clones" method first'
            return

        fileout = open(output_filepath, "w")

        if sort_by_count:
            #add count to beginning of each element in the list
            clones_with_count = []
            for i in self.clones:
                clones_with_count.append([i['count'], i])
            #now sort based on count and write to fileout
            seq_id = 0
            for i in sorted(clones_with_count, reverse=True):
                seq_id += 1
                header = '|'.join(['%s=%s' % (j, i[1][j]) for j in i[1] if j!='cdr3'])
                fileout.write('>%s|%s\n%s\n' % (seq_id, header, i[1]['cdr3']))

        else:
            seq_id = 0
            for i in self.clones:
                seq_id += 1
                header = '|'.join(['%s=%s' % (j, i[j]) for j in i if j!='cdr3'])
                fileout.write('>%s|%s\n%s\n' % (seq_id, header, i['cdr3']))

        fileout.close()

        return

    def simulate_new_samples(self, element_to_simulate='clone', trials=1, serial_samples_to_sim=1, serial_samp_pop_counts=None, serial_samp_tpoints=None, update_freqs_with_serial_samp=True, output_dirpath=None, include_observed_sample=False, use_comp_cluster=False, temp_dirpath=None, time_per_generation=None, store_full_data_in_ram=True):
        """
        This method will use multinomial sampling to use the frequencies of elements in the immune population (an 'element' could be clones, V genes, J genes, unique sequences, etc.) to simulate new samples with those elements randomly picked to be in them. The frequencies of the elements are used to get the probability that they will be sampled (with replacement) to fill the new simulated sample. Thus, in these simulated samples, there is no chance of getting novel elements (i.e. not observed in the real immune population), an important thing to keep in mind. One can simulate a time-course of serially sampled immune samples, where the probability of selecting an immune element for the next time-point is determined by the element frequencies in the previous time-point. One can also run multiple 'trials' of simulations, where the whole process is repeated a given number of times.
        element_to_simulate - Signifies which component of the sequence data will be simulated. Essentially any component of the immune population that has a frequency could be simulated (i.e. V gene, clones, etc), so this tells the script what that component is. Exceptable values are:
            'clone' - (default) This means that the frequencies of clones will be calculated and used to simulate new samples of clones.
            'v_gene' - This means that the frequencies of V genes will be calculated and used to simulate new samples of V genes.
            'j_gene' - Same as 'v_gene', but for J genes.
            'entries' - This means that each entry in the in input fasta file will be simulated. For each entry, its frequency is determined by the count attribute for that entry, and then simulated to make new samples of sequence entries.
        trials - This gives the number (int) of independent simulations to run
        serial_samples_to_sim - This gives the number (int) of serial samples to simulate within each 'trial'. For example, one wishes to simulate a timecourse of 5 time-points, and do this 10 times then trials=10, and serial_samples_to_sim=5. If one does not wish to simulate a time course then simply set serial_samples_to_sim=1 (default).
        serial_samp_pop_counts - If the immune population size changes with your time-points (i.e. differing number of B-cells in different samples), then this can be incorporated into the simulation with this variable. This should be a numeric list of population sizes, ordered by time. The length of this list should be equal to serial_samples_to_sim. If this parameter equals None (default) then the popsize is assumed to be consistent and equal to 'self.total' (i.e. number of reads) of the immune sample.
        update_freqs_with_serial_samp - If True (default) then the frequencies of elements of the previous time-point are used to determine the probabiliteis of selecting the elements for the next time-point (similar to a Wright-Fisher sim). If False, then the frequencies of the 1st (the true, observed time-point) are always used for selecting immune elements in all time-points. 
        output_dirpath - If defined (default None), then this should be defined as a path (str) to the directory for which all of the simulations will be written. Each output file corresponds to one 'trial'. The format of the output is a tab delimited list, where the 1st column is the names (IDs) of each simulated element, the 2nd column is the observed counts for each element in the sample (i.e. 1st time-point), the following columns are the simulated counts for each element, and columns are ordered with serial samples.
        serial_samp_tpoints - If defined (default None), then this is a list of floats that gives the time-point of each of the simulated serial samples, as well as the observed sample. This is only used for labeling the the columns of the output files. If output_dirpath=None then this parameter is ignored. If output_dirpath is defined and this varialbe is undefined then the time-point labels are consequetive integers starting at 0.0. The length of the list should be equal to serial_samples_to_sim+1 (+1 to account for the time-point for the observed sample)
        include_observed_sample - If True (default False) then the count information for each element in the observed data (i.e. the 1st tpoint) will be include in the output (i.e. 'sim_element_counts_record')
        use_comp_cluster - If True (default, False) then this will instruct the method to use an SGE interface to a computational cluster. Each job will be the computation for one of the trials. So, the number of jobs will equal 'trials'.
        temp_dirpath - If defined (default, None), then this should be the path to a directory where temp files will be written. This parameter is only considered if 'use_comp_cluster' is True.
        time_per_generation - This gives how many units of time per 1 generation (i.e. per 1 sampling event). A unit of time is assumed to be the same as that used in 'serial_samp_tpoints'. This tells the method how many serial sampling events should take place between simulated time-points. If this equals None (default) then there is only one sampling event between each serial sample. if 'serial_samp_tpoints' is not defined then it is assumed the there is 1 unit of time between each consequetive sample. So, this would need to be set to <= 0.5 in order to sample multiple times between time-points, in that case. This can also be a list, in which case multiple values will be cycled through. This means that the entire program will be run for each of the parameter values in this list. So, all the trials, time-points, and elements will be simulated for each of the 'time_per_generation' values given in the list. Use with caution. If this is set to be a list than the dimensionality of the output is increased by one.
        store_full_data_in_ram - If True (default), then the method will store all the results of the sims in the output variable 'sim_element_counts_record', and this will be returned as output. If this is False, then the method will still return sim_elements_count_record, but instead of this variable containing the counts for each of the elements in each of the sims, it will contain the (absolute) paths to the files that contain this information. These paths will lead to files located in the 'temp_dirpath', so this variable must be defined as well.

        OUTPUT:
        Returns 'element_names' and 'sim_element_counts_record'. 'element_names' is a list that gives some unique name, or ID, for each of the simulated elements of the immune population. For example, if simulating V genes, this would be an ordered list of V gene names. 'sim_element_counts_record' is a record of the counts for each simulated element for each of the trials/serial samples. This is a 3 dimensional list. The 1st level of the list corresponds to trial number, the 2nd level corresponds to serial sample number, and the 3rd corresponds unique simulated element (i.e. unique V gene, if that is what one is simulating). The elements of the 3rd level are given in the same order as 'element_names'.
        """
        #make sure we're not adding sampling events if were not updating the frequencies of the elements with each event (i.e. generation).
        if update_freqs_with_serial_samp == False and time_per_generation:
            print "There is no point in adding sampling events (i.e. multiple generations between time-points) if the 'update_freqs_with_serial_samp' variable is False. Aborting."
            return None, None
        #make sure clones have been compiled 1st
        if self.clones == None and element_to_simulate=='clone':
            print 'The clonal sequences have not been gathered yet from the sample. Run the "group_clones" method first'
            return None, None
        #if using comp cluster make sure temp dirpath is defined
        if use_comp_cluster and temp_dirpath == None:
            print 'If using a computational cluster, the "temp_dirpath" parameter must be defined. Aborting'
            print 'temp_dirpath:', temp_dirpath
            return None, None
        if temp_dirpath:
            if temp_dirpath[-1] != '/':
                temp_dirpath += '/'
            if not os.path.exists(temp_dirpath):
                os.makedirs(temp_dirpath)
        #if NOT storing data in RAM, make sure temp_dirpath is defined
        if not store_full_data_in_ram and not temp_dirpath:
            print 'If not "store_full_data_in_ram", then "temp_dirpath" must be defined. Aborting'
            return None, None
        #get names for each element
        if element_to_simulate == 'clone':
            element_names = [i['id'] for i in self.clones]
        elif element_to_simulate == 'v_gene':
            element_names = [i for i in self.vgene_counts]
        elif element_to_simulate == 'j_gene':
            element_names = [i for i in self.jgene_counts]
        elif element_to_simulate == 'entries':
            element_names = [i['id'] for i in self.data]
        #set initial probabilities for each element based on observed data
        if element_to_simulate == 'clone':
            init_element_probs = [i['count']/self.total for i in self.clones]
        elif element_to_simulate == 'v_gene':
            init_element_probs = [self.vgene_counts[i]/self.total for i in self.vgene_counts]
        elif element_to_simulate == 'j_gene':
            init_element_probs = [self.jgene_counts[i]/self.total for i in self.jgene_counts]
        elif element_to_simulate == 'entries':
            init_element_probs = [i['count']/self.total for i in self.data]
        if not serial_samp_tpoints:
            serial_samp_tpoints = range(serial_samples_to_sim+1)
        if not serial_samp_pop_counts:
            serial_samp_pop_counts = [self.total for i in xrange(serial_samples_to_sim+1)]
        #make sure pop sizes are floats, if not force them to be.
        if not isinstance(serial_samp_pop_counts[0], float):
            serial_samp_pop_counts = [float(i) for i in serial_samp_pop_counts]
        if use_comp_cluster or not store_full_data_in_ram:
            random_suffix = str(random.random())
            base_filepath = temp_dirpath + 'input_' + random_suffix + '_'
        if not isinstance(time_per_generation, list):
            time_per_generation = [time_per_generation]

        #if using computational cluster
        if use_comp_cluster:
            #create input files foreach job
            for i in xrange(trials):
                #this will be a temp file that contains all the data needed for a single trial on the cluster
                input_parameter_filepath = base_filepath + str(i+1)
                fileout = open(input_parameter_filepath, "w")
                #write element_probs to temp input file
                fileout.write(','.join([str(j) for j in init_element_probs]) + '\n')
                #write pop sizes
                fileout.write(','.join([str(j) for j in serial_samp_pop_counts]) + '\n')
                #write time-points
                fileout.write(','.join([str(j) for j in serial_samp_tpoints]) + '\n')
                #write whether or not to include observed sample in output
                if include_observed_sample:
                    fileout.write('True\n')
                else:
                    fileout.write('False\n')
                #write whether or not to update serial sample freqs
                if update_freqs_with_serial_samp:
                    fileout.write('True\n')
                else:
                    fileout.write('False\n')
                #write how many gens per unit time
                fileout.write(','.join([str(j) for j in time_per_generation]) + '\n')
                fileout.close()
            p = Popen(['qsub', '-t', '1-'+str(trials), 'immune_sample_class.py', 'simulate_samples_using_comp_cluster', base_filepath], stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()
            job_id = set([out.split()[2].split('.')[0]])
            #wait for jobs to complete
            good_to_go = False
            while good_to_go == False:
                p = Popen(['qstat'], stdout=PIPE, stderr=PIPE)
                jobs = p.communicate()[0].split('\n')
                jobs = set([i.split()[0] for i in jobs[2:-1]])
                if job_id & jobs: #if the the intersection of the job sets is non-empty
                    time.sleep(10)
                else:
                    good_to_go = True
            #when they are done, read in the sim results
            sim_element_counts_record = []
            for gen_time_index, gen_time in enumerate(time_per_generation):
                sim_element_counts_record.append([])
                for i in xrange(trials):
                    #delete the file that had the input parameters for a job
                    trial_input_filepath = '%s%s' % (base_filepath, i+1)
                    if os.path.exists(trial_input_filepath):
                        subprocess.call(['rm', trial_input_filepath])
                    trial_output_filepath = '%soutput_%s_%s_%s' % (temp_dirpath, random_suffix, gen_time, i+1)
                    if store_full_data_in_ram:
                        filein = open(trial_output_filepath, "r")
                        one_sim = [[int(k) for k in j[:-1].split(',')] for j in filein]
                        sim_element_counts_record[gen_time_index].append(one_sim)
                        filein.close()
                        subprocess.call(['rm', trial_output_filepath])
                    else:
                        sim_element_counts_record[-1].append(trial_output_filepath)

        #if not using computational cluster
        else:
            sim_element_counts_record = []
            for gen_time in time_per_generation:
                sim_element_counts_record.append([])
                #if there is a time_per_generation set, figure out time-points and pop_sizes over time
                if gen_time:
                    serial_samp_tpoints_full = []
                    serial_samp_pop_counts_full = []
                    for i in xrange(len(serial_samp_tpoints)-1): #for each space between time points
                        hidden_tpoints = numpy.arange(start=serial_samp_tpoints[i], stop=serial_samp_tpoints[i+1], step=gen_time).tolist()
                        serial_samp_tpoints_full += hidden_tpoints #add extra 'hidden' time points between observed time points
                        serial_samp_pop_counts_full += [serial_samp_pop_counts[i] for j in hidden_tpoints]#add pop sizes for each of these hidden time points as well
                    #add final values to new time-point and pop size lists
                    serial_samp_tpoints_full.append(serial_samp_tpoints[-1])
                    serial_samp_pop_counts_full.append(serial_samp_pop_counts[-1])
                else:
                    serial_samp_pop_counts_full = serial_samp_pop_counts[:]
                for i in xrange(trials):
                    element_probs = init_element_probs[:]
                    if not store_full_data_in_ram:
                        output_filepath = '%soutput_%s_%s_%s' % (temp_dirpath, random_suffix, gen_time, i+1)
                        sim_element_counts_record[-1].append(output_filepath)
                        fileout = open(output_filepath, "w")
                    else:
                        sim_element_counts_record[-1].append([])
                    if include_observed_sample:#if the observed count data is desired as well
                        if store_full_data_in_ram:
                            sim_element_counts_record[-1][-1].append([int(round(j * serial_samp_pop_counts[0])) for j in element_probs])
                        else:
                            element_counts = [int(round(j * serial_samp_pop_counts[0])) for j in element_probs]
                            fileout.write(','.join([str(j) for j in element_counts]) + '\n')
                    for index, pop_size in enumerate(serial_samp_pop_counts_full[1:]):
                        tpoint_index = index + 1
                        sim_element_counts = numpy.random.multinomial(pop_size, element_probs)
                        if gen_time:
                            if serial_samp_tpoints_full[tpoint_index] in serial_samp_tpoints:
                                if store_full_data_in_ram:
                                    sim_element_counts_record[-1][-1].append(sim_element_counts.tolist())
                                else:
                                    fileout.write(','.join([str(j) for j in sim_element_counts]) + '\n')
                        else:
                            if store_full_data_in_ram:
                                sim_element_counts_record[-1][-1].append(sim_element_counts.tolist())
                            else:
                                fileout.write(','.join([str(j) for j in sim_element_counts]) + '\n')
                        if update_freqs_with_serial_samp:
                            element_probs = [k/pop_size for k in sim_element_counts]
                    if not store_full_data_in_ram:
                        fileout.close()

        #if desired, write sim results to disk
        if output_dirpath:
            if output_dirpath[-1] != '/':
                output_dirpath += '/'
            if not os.path.exists(output_dirpath):
                os.makedirs(output_dirpath)
            for gen_time_index, gen_time in enumerate(time_per_generation):
                gen_time_output_dirpath = output_dirpath + str(gen_time) + '/'
                if not os.path.exists(gen_time_output_dirpath):
                    os.makedirs(gen_time_output_dirpath)
                for i in xrange(trials):
                    fileout = open(gen_time_output_dirpath+str(i), "w")
                    fileout.write('\t' + '\t'.join([str(k) for k in serial_samp_tpoints]) + '\n')
                    for j in xrange(len(element_names)):
                        fileout.write(element_names[j] + '\t')
                        #if we don't include the observed value in the output, then we need to add it to the written output here
                        if include_observed_sample == False:
                            if element_to_simulate == 'clone':
                                observed_value = int(round((self.clones[j]['count']/self.total) * serial_samp_pop_counts[0]))
                            elif element_to_simulate == 'v_gene':
                                observed_value = int(round((self.vgene_counts[element_names[j]]/self.total) * serial_samp_pop_counts[0]))
                            elif element_to_simulate == 'j_gene':
                                observed_value = int(round((self.jgene_counts[element_names[j]]/self.total) * serial_samp_pop_counts[0]))
                            elif element_to_simulate == 'entries':
                                observed_value = int(round((self.data[j]['count']/self.total) * serial_samp_pop_counts[0]))
                            fileout.write(str(observed_value) + '\t')
                        sim_values = [str(sim_element_counts_record[gen_time_index][i][k][j]) for k in xrange(serial_samples_to_sim)]
                        fileout.write('\t'.join(sim_values) + '\n')
                    fileout.close()

        if len(time_per_generation) == 1:
            sim_element_counts_record = sim_element_counts_record[0]
        return element_names, sim_element_counts_record

    def add_annotations_using_mixcr(self, output_filepath, annotations_to_get, delete_output=False, debug=False):
        """
        This method first runs MIXCR align, and then uses the 'mixcr_align_output' class to parse desired information out of the MIXCR output. The MIXCR executables MUST be in your PATH for this to work.
        output_filepath - The path to the file that will contain the output from running MIXCR align. This should not contain a suffix such as '.txt' or '.csv'
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

        NOTE:
        This method will change the data. If a sequence has no output from MIXCR align then it will be removed. Also, if an annotation to get from the MIXCR output already exists in the data (ex: CDR3 seq), then it will be replaced by the values found by MIXCR.
        """

        #run MIXCR
        print 'running MIXCR...'
        p = Popen(['mixcr', 'align', '-f', '--write-all', self.filepath, output_filepath + '.vdjca', ], stdout=PIPE, stderr=PIPE)
        out_err = p.communicate()
        p = Popen(['mixcr', 'exportAlignments', '-f', '-readId', '-sequence', '-defaultAnchorPoints', '-nFeature', 'FR1', '-nFeature', 'CDR1', '-nFeature', 'FR2', '-nFeature', 'CDR2', '-nFeature', 'FR3', '-nFeature', 'CDR3', '-nFeature', 'FR4', '-aaFeature', 'FR1', '-aaFeature', 'CDR1', '-aaFeature', 'FR2', '-aaFeature', 'CDR2', '-aaFeature', 'FR3', '-aaFeature', 'CDR3', '-aaFeature', 'FR4', '-aaFeature', 'VDJRegion', output_filepath + '.vdjca', output_filepath + '.txt'])
        out_err = p.communicate()

        #get annotations from MIXCR align output
        alignment = mixcr_align_output(output_filepath + '.txt')
        annotations = alignment.extract_annotations(annotations_to_get=annotations_to_get)
        if delete_output:
            subprocess.call(['rm', output_filepath + '.vdjca', output_filepath + '.txt'])

        #cycle through seq data and add new annotations from MIXCR
        indices_to_remove = []
        for i in xrange(len(self.data)):
            if not i in annotations:
                indices_to_remove.append(i)
                continue
            if 'cdr3' in annotations_to_get:
                self.data[i]['cdr3'] = annotations[i]['cdr3'][:]
                del annotations[i]['cdr3']
            self.data[i]['other'].update(annotations[i])

        #remove elements in all parts of the data that did not have output fom MIXCR alignment
        for i in sorted(indices_to_remove, reverse=True):
            seq_count = self.counts[i]
            vgene_name = self.data[i]['vgene']
            dgene_name = self.data[i]['dgene']
            jgene_name = self.data[i]['jgene']
            del self.data[i]
            self.total -= seq_count
            self.vgene_counts[vgene_name] -= seq_count
            self.dgene_counts[dgene_name] -= seq_count
            self.jgene_counts[jgene_name] -= seq_count
        #need to redo the grouping of clones after changing the underlying data
        if not self.clones == None:
            self.group_clones(define=self.clone_define)

        #check to make sure things have been properly added, and that everything adds up
        if debug:
            for i in self.data:
                try:
                    for j in annotations_to_get:
                        if j != 'cdr3':
                            i['other'][j]
                except:
                    print 'found entry that does not have updated attributes:'
                    print j
                    print i
            total_data = 0
            for i in self.data:
                total_data += i['count']
            if total_data != self.total:
                print 'The "count" attribute in self.data does not add up to self.total'
                print 'Sum of counts:', total_data
                print 'self.total:', self.total
            total_vgenes = 0
            for i in self.vgene_counts:
                total_vgenes += self.vgene_counts[i]
            if total_vgenes != self.total:
                print 'All of the counts of V genes in self.vgene_counts do not add up to self.total'
                print 'Sum of counts:', total_vgenes
                print 'self.total:', self.total

        return

    def convert_to_tab_delimited_data(self, output_filepath, seq_column_name='SEQUENCE_INPUT', seq_id_column_name='SEQUENCE_ID', add_cols=None):
        """
        This method writes the data to disk such that it is tab delimited. This format is compatible with change-o database formatting (conditional on he names of the attributes being apporopriate).
        output_filepath - The path to the file that will contain the new tab delimited data. 
        seq_column_name - This gives the name of the column that will house the 'input sequence' (i.e. the sequence that occupies of seq space of each fasta entry). The default is set to be compatible with change-o database formatting.
        seq_id_column_name - The name of the column that gives the ID for each of the seqs. The default is set ot be compatable with change-o database formatting.
        add_cols - If defined (default None), then this will add column(s) that will be defined identically for each row (seq) in the data. This should be a list of two-element lists. Each two element list represents an added column, and the first element in a two element list gives the name of an added column, and the second gives the value that will be added for each row in the column.
        """
        fileout = open(output_filepath, "w")
        header = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (seq_id_column_name, seq_column_name, self.vgene_name, self.dgene_name, self.jgene_name, self.cdr3_name, self.count_attribute_name)
        header += '\t'.join([i for i in sorted(self.data[0]['other'])])
        if add_cols:
            header += '\t' + '\t'.join([str(i[0]) for i in add_cols])
        fileout.write(header + '\n')
        for i in self.data:
            line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (i['id'], i['seq'], i['vgene'], i['dgene'], i['jgene'], i['cdr3'], i['count'])
            line += '\t'.join([str(i['other'][j]) for j in sorted(i['other'])])
            if add_cols:
                line += '\t' + '\t'.join([str(j[1]) for j in add_cols])
            fileout.write(line + '\n')
        fileout.close()
        return

    def get_selection_values_using_baseline(self, output_filepath=None, temp_data_dirpath=None, group_sel_stats=False):
        """
        This method uses the program baseline to assess the selection level for each of the antibody seqs in the data. It uses baseline via the shazam R package, so this package must be installed on whatever default R environment is used. The data in the sample must have a few attributes that are interpretable by baseline. Namely, 'SEQUENCE_IMGT', AND 'GERMLINE_IMGT_D_MASK'. These attributes make it possible for baseline to make the correct amount of counts for each of the different categories (i.e. CDR, FWR, silent, and replacement). This method also gives summary statistics of selection for all the seqs in the sample.
        output_filepath - path to file for which the output from Baseline will be written. Writes the selection values for each seq in the data, as well some stats like upper and lower confidence interval, and p values. If this equals None (default) then no output is written to disk, and selection values are only returned in memory.
        temp_data_dirpath - Path to a directory for which all the temporary input and output data for baseline will be stored. This directory is deleted upon completion. If equals None (default) then this temp dir is made in the dir that contain the output filepath.
        group_sel_stats - If not False (default False), then this will instruct the the method to group the selection values for all the seqs in the sample ('convolve') to get one summary stat for the entire sample regarding selection. If defined this should equal the filepath to the file that will contain this sample wide selection summary stat info.

        OUPUT:
        Returns 'selection_values', which is a dic, where each key is the unique ID of the sequence, and its definition is the selection values reported by baseline (as a dictionary).
        Also returns 'selection_sum_stats' which is a dic that has summary statistics the selection values of all the seqs in the sample (i.e. mean, upper/lower CI, p value). If, however, 'group_sel_stats' equals False then this returned variable will also equal None
        """
        if temp_data_dirpath == None and output_filepath == None:
            print 'Both "temp_data_dirpath" and "output_filepath" cannont equal None. Aborting'
            return
        if not temp_data_dirpath:
            temp_data_dirpath = os.path.dirname(output_filepath) + '/temp_' + str(random.random()) + '/'
        else:
            if temp_data_dirpath[-1] != '/':
                temp_data_dirpath += '/'
            if not os.path.exists(temp_data_dirpath):
                os.makedirs(temp_data_dirpath)
        if not output_filepath:
            output_filepath = temp_data_dirpath + os.path.basename(self.filepath)[:-6] + str(random.random())

        #check to make sure data has the needed attributes
        if not 'SEQUENCE_IMGT' in self.data[0]['other'] or not 'GERMLINE_IMGT_D_MASK' in self.data[0]['other']:
            print "The data must have the attributes 'SEQUENCE_IMGT' and 'GERMLINE_IMGT_D_MASK' in it. Aborting."
            return

        #make temp change-o tab delimited database out of the data
        changeo_db_filepath = temp_data_dirpath + os.path.basename(self.filepath)[:-6] + '.tab'
        self.convert_to_tab_delimited_data(changeo_db_filepath, add_cols=[['GROUP', 1]])

        #run baseline
        if group_sel_stats:
            p = Popen(['Rscript', './invoke_baseline.R', changeo_db_filepath, output_filepath, 'GROUP', group_sel_stats], stdout=PIPE, stderr=PIPE)
        else:
            p = Popen(['Rscript', './invoke_baseline.R', changeo_db_filepath, output_filepath], stdout=PIPE, stderr=PIPE)
        out_err = p.communicate()
        for i in out_err[0].split('\n'):
            print i
        print out_err[1]
        subprocess.call(['rm', changeo_db_filepath])

        #read selection values into memory
        filein = open(output_filepath, "r")
        filein.readline()
        selection_values = {}
        for i in filein:
            line = i[:-1].split('\t')
            seq_id = line[0]
            region = line[1]
            sigma = line[2]
            ci_low = line[3]
            ci_high = line[4]
            p_val = line[5]
            try:
                selection_values[seq_id][region + '_sigma'] = sigma
                selection_values[seq_id][region + '_ci_low'] = ci_low
                selection_values[seq_id][region + '_ci_high'] = ci_high
                selection_values[seq_id][region + '_p_val'] = p_val
            except KeyError:
                selection_values[seq_id] = {'CDR_sigma':'NA', 'CDR_ci_low':'NA', 'CDR_ci_high':'NA', 'CDR_p_val':'NA', 'FWR_sigma':'NA', 'FWR_ci_low':'NA', 'FWR_ci_high':'NA', 'FWR_p_val':'NA'}
                selection_values[seq_id][region + '_sigma'] = sigma
                selection_values[seq_id][region + '_ci_low'] = ci_low
                selection_values[seq_id][region + '_ci_high'] = ci_high
                selection_values[seq_id][region + '_p_val'] = p_val
        filein.close()
        #If it is not desired to keep the baseline output
        if os.path.dirname(output_filepath) == temp_data_dirpath[:-1]:
            subprocess.call(['rm', output_filepath])

        #if summary stats are desired, then get those too
        if group_sel_stats:
            selection_sum_stats = {}
            filein = open(group_sel_stats, "r")
            filein.readline()
            for i in filein:
                line = i[:-1].split('\t')
                region = line[1]
                sigma = line[2]
                ci_low = line[3]
                ci_high = line[4]
                p_val = line[5]
                selection_sum_stats[region + '_sigma'] = sigma
                selection_sum_stats[region + '_ci_low'] = ci_low
                selection_sum_stats[region + '_ci_high'] = ci_high
                selection_sum_stats[region + '_p_val'] = p_val
        else:
            selection_sum_stats = None

        return selection_values, selection_sum_stats

    @staticmethod
    def write_one_seq_entry_to_disk_immune(self, data_entry, fileout):
        """
        This is an internal method that is used by other methods to write one entry of a fasta file to an output file.
        """
        other_headers = []
        for i in data_entry['other']:
            if type(data_entry['other'][i]) is list:
                value = ','.join([str(j) for j in data_entry['other'][i]])
            else:
                value = data_entry['other'][i]
            other_headers.append('%s=%s' % (i, value))
        header = '>%s|%s=%s|%s=%s|%s=%s|%s=%s|%s=%s|%s\n' % (data_entry['id'], self.vgene_name, data_entry['vgene'], self.dgene_name, data_entry['dgene'], self.jgene_name, data_entry['jgene'], self.cdr3_name, data_entry['cdr3'], self.count_attribute_name, data_entry['count'], '|'.join(other_headers))
        fileout.write(header + data_entry['seq'] + '\n')
        return

    def write_full_data_to_disk(self, output_filepath):
        """
        This method writes all the data to disk. Including the seq data, and all the attributes that are defined for each sequence entry. Writes in fasta format.
        output_filepath - The path to the file for which all this data will be written.
        """
        fileout = open(output_filepath, "w")
        for i in self.data:
            self.write_one_seq_entry_to_disk_immune(self, i, fileout)
        fileout.close()
        return

    def write_subset_of_seqs_to_disk_immune(self, seq_indices, output_filepath):
        """
        This method only writes a subset of the seqs to disk, given by the 'seq_indices', which is a list of indices for the seqs that are desired to be written.
        """
        fileout = open(output_filepath, 'w')
        for i in seq_indices:
            self.write_one_seq_entry_to_disk_immune(self, self.data[i], fileout)
        fileout.close()
        return

    def divide_data_by_VJ_combo(self):
        """
        This method divides the self.data object based upon the combination of V and J gene segments that a given sequence has. It returns a dictionary, indexed by V gene name and J gene name, and def is a list of the indeces of the seqs that have that given V and J combo.
        """
        VJ_data_dic = {}
        for i in xrange(len(self.data)):
            vgene = self.data[i]['vgene']
            jgene = self.data[i]['jgene']
            try:
                VJ_data_dic['%s_%s' % (vgene, jgene)].append(i)
            except KeyError:
                VJ_data_dic['%s_%s' % (vgene, jgene)] = [i]
        return VJ_data_dic

    def convert_to_cdr3_seqs_and_write_to_disk(self, output_filepath):
        """
        This method will first combine all the seqs that have the same CDR3 sequence, and the same V and J genes, and combine them into the same seq entry. It will then write this data to disk. It will also chnage all necessary attributes, such as 'counts' in the data. For each CDR3 sequence entry: the seq is the CDR3, the V and J genes are determined by the original seqs that make up the entry, and the CDR3 attribute will be the same the sequence entry. The rest of the attibutes will be retained, but may not make a whole lot of sense. For example, if divergence is an attribute, then it probably shouldn't be considered since the CDR3 is quite difficult to align to germline seqs. However, for compatability reasons with down stream code, we keep all the original attributes in the CDR3 seq entries. When multiple seqs are combined into one CDR3 seq, we take the attributes from the original seq that has the highest count to be represented in the new CDR3 seq entry.
        output_filepath - the path to the file that will contain the CDR3 seq data, post processing
        """
        new_data = {}
        for seq_entry in self.data:
            #CDR3, V gene, and J gene must all be defined. If not, skip
            if seq_entry['cdr3'] == 'NA' or seq_entry['vgene'] == 'NA' or seq_entry['jgene'] == 'NA':
                continue
            cdr3_v_j = '%s_%s_%s' % (seq_entry['cdr3'], seq_entry['vgene'], seq_entry['jgene'])
            try:
                new_data[cdr3_v_j].append([seq_entry['count'], seq_entry])
            except KeyError:
                new_data[cdr3_v_j] = [[seq_entry['count'], seq_entry]]
        self.data = []
        self.counts = []
        self.total = 0.
        for cdr3_entry in new_data:
            new_data[cdr3_entry] = sorted(new_data[cdr3_entry])
            self.data.append(new_data[cdr3_entry][-1][1])
            self.data[-1]['seq'] = self.data[-1]['cdr3']
            count = 0
            for seq_entry in new_data[cdr3_entry]:
                count += seq_entry[0]
            self.data[-1]['count'] = count
            self.counts.append(count)
            self.total += count
        self.write_full_data_to_disk(output_filepath)
        return


#############################################
# Below are scripts that are used for running
# jobs on the computational cluster. They
# should only be used by the methods above.
#############################################

def simulate_samples_using_comp_cluster(input_parameter_filepath_base):
    """
    This script runs each of the trials in the 'simulate_new_samples' method as job on a computational cluster. 
    """
    sge_task_id = int(os.environ['SGE_TASK_ID'])
    #read in input parameters
    filein = open(input_parameter_filepath_base+str(sge_task_id), "r")
    element_probs_init = [float(i) for i in filein.readline()[:-1].split(',')]
    serial_samp_pop_counts = [float(i) for i in filein.readline()[:-1].split(',')]
    serial_samp_tpoints = [float(i) for i in filein.readline()[:-1].split(',')]
    include_observed_sample = filein.readline()[:-1]
    if include_observed_sample == "True":
        include_observed_sample = True
    else:
        include_observed_sample = False
    update_freqs_with_serial_samp = filein.readline()[:-1]
    if update_freqs_with_serial_samp == "True":
        update_freqs_with_serial_samp = True
    else:
        update_freqs_with_serial_samp = False
    time_per_generation = filein.readline()[:-1].split(',')
    time_per_generation_new = []
    for i in time_per_generation:
        if i != 'None':
            time_per_generation_new.append(float(i))
        else:
            time_per_generation_new.append(None)
    time_per_generation = time_per_generation_new[:]
    filein.close()

    #set up output filepaths
    temp_dirpath = os.path.dirname(input_parameter_filepath_base)
    random_suffix = os.path.basename(input_parameter_filepath_base).split('_')[1]

    #run sims and write to output file throughout
    for gen_time in time_per_generation:
        element_probs = element_probs_init[:]
        output_filepath = '%s/output_%s_%s_%s' % (temp_dirpath, random_suffix, gen_time, sge_task_id)
        fileout = open(output_filepath, "w")
        #if there is a time_per_generation set, figure out time-points and pop_sizes over time
        if gen_time:
            serial_samp_tpoints_full = []
            serial_samp_pop_counts_full = []
            for i in xrange(len(serial_samp_tpoints)-1): #for each space between time points
                hidden_tpoints = numpy.arange(start=serial_samp_tpoints[i], stop=serial_samp_tpoints[i+1], step=gen_time).tolist()
                serial_samp_tpoints_full += hidden_tpoints #add extra 'hidden' time points between observed time points
                serial_samp_pop_counts_full += [serial_samp_pop_counts[i] for j in hidden_tpoints]#add pop sizes for each of these hidden time points as well
            #add final values to new time-point and pop size lists
            serial_samp_tpoints_full.append(serial_samp_tpoints[-1])
            serial_samp_pop_counts_full.append(serial_samp_pop_counts[-1])
        else:
            serial_samp_pop_counts_full = serial_samp_pop_counts[:]
        if include_observed_sample:#if the observed count data is desired as well
            element_counts = [int(round(i * serial_samp_pop_counts_full[0])) for i in element_probs]
            fileout.write(','.join([str(i) for i in element_counts]) + '\n')
        for i, pop_size in enumerate(serial_samp_pop_counts_full[1:]):
            tpoint_index = i+1
            sim_element_counts = numpy.random.multinomial(pop_size, element_probs)
            if gen_time:
                if serial_samp_tpoints_full[tpoint_index] in serial_samp_tpoints:
                    fileout.write(','.join([str(j) for j in sim_element_counts.tolist()]) + '\n')
            else:
                fileout.write(','.join([str(j) for j in sim_element_counts.tolist()]) + '\n')
            if update_freqs_with_serial_samp:
                element_probs = [j/pop_size for j in sim_element_counts]
        fileout.close()

    return

if __name__ == '__main__':
    if sys.argv[1] == 'simulate_samples_using_comp_cluster':
        simulate_samples_using_comp_cluster(input_parameter_filepath_base=sys.argv[2])
    elif sys.argv[1] == 'test':
        f = immune_sample('/Users/nstrauli/data/abr_hiv_coevo/seq_data/changeo_annotated_fasta_files/1/1144.fasta', 'DUPCOUNT', 'V_CALL', 'D_CALL', 'J_CALL', 'CDR3_IMGT')
        f.convert_to_cdr3_seqs_and_write_to_disk(output_filepath='/Users/nstrauli/data/abr_hiv_coevo/seq_data/cdr3_seqs/changeo/1/1144.fasta')
    else:
        print 'huh?'
