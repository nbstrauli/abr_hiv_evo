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
               "GGT":"G","GGC":"G","GGA":"G","GGG":"G"}

def get_num_nonSyn_and_syn_sites_foreach_codon():
     nonSyn_syn_sites_dic = {}
     for codon in geneticCode:
          nonSyn_syn_sites_dic[codon] = {'non_syn':0, 'syn':0}
          for position in xrange(3):
               for mutation in ['A', 'C', 'T', 'G']:
                    if mutation != codon[position]:
                         mut_codon = codon[:position] + mutation + codon[position+1:]
                         if geneticCode[codon] == geneticCode[mut_codon]:
                              nonSyn_syn_sites_dic[codon]['syn'] += 1/3.
                         else:
                              nonSyn_syn_sites_dic[codon]['non_syn'] += 1/3.
     return nonSyn_syn_sites_dic
