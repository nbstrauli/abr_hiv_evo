class mixcr_align_output(object):
	"""
	A class to manipulate and extract information out of the output from MIXCR align.

	Attributes:
	filepath - path to the file that contains the MIXCR align output.
	annotations - A list of the annotations that are offered in the MIXCR align output file (i.e. those that are in the header of the tab delimited file)
	seq_annotations - A list of dics. Each element represent the info for a single seq in the data (i.e. a single line in the MIXCR align output file). The keys in each of the dics are equal to the strings that populate self.annotations. The definition for each of the keys are the values that are in the file for that line (or seq).
	"""
	def __init__(self, filepath):
		self.filepath = filepath
		self.seq_annotations = []
		self.annotations = []
		filein = open(filepath, 'r')
		header = filein.readline()
		self.annotations = [i for i in header[:-1].split('\t')]
		for i in filein:
			line = i[:-1].split('\t')
			entry = {}
			for count, j in enumerate(line):
				entry[self.annotations[count]] = j
			self.seq_annotations.append(entry)
		filein.close()
		return

	def extract_annotations(self, annotations_to_get):
		"""
		This method will retrieve the desired information (or annotations) from each of the sequences (or lines) in the data.
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

		Output:
		annotations - A dic of dics, where each dic (within the main dic) contains the annotations, and values for each of the annotations, for each of the seqs in the data. Each element in the dic contains the desired data for a given seq. The keys in each of the dics are the same as the strings that populate 'annotations_to_get'. The main dic is indexed by the seq ID given to each sequence by MIXCR. This ID is simply the order of the input sequence (starting at 0), so the first seq in the input file would have an ID of 0, the 2nd seq would have an ID of 1, etc.
		"""
		annotations = {}
		for i in self.seq_annotations:
			seq_id = int(i['readId'])
			annotations[seq_id] = {}
			annotations[seq_id]['seq'] = i['readSequence']
			if 'fwr1' in annotations_to_get:
				if i['nSeqFR1'] == '':
					value = 'NA'
				else:
					value = i['nSeqFR1']
				annotations[seq_id]['fwr1'] = value
			if 'cdr1' in annotations_to_get:
				if i['nSeqCDR1'] == '':
					value = 'NA'
				else:
					value = i['nSeqCDR1']
				annotations[seq_id]['cdr1'] = value
			if 'fwr2' in annotations_to_get:
				if i['nSeqFR2'] == '':
					value = 'NA'
				else:
					value = i['nSeqFR2']
				annotations[seq_id]['fwr2'] = value
			if 'cdr2' in annotations_to_get:
				if i['nSeqCDR2'] == '':
					value = 'NA'
				else:
					value = i['nSeqCDR2']
				annotations[seq_id]['cdr2'] = value
			if 'fwr3' in annotations_to_get:
				if i['nSeqFR3'] == '':
					value = 'NA'
				else:
					value = i['nSeqFR3']
				annotations[seq_id]['fwr3'] = value
			if 'cdr3' in annotations_to_get:
				if i['nSeqCDR3'] == '':
					value = 'NA'
				else:
					value = i['nSeqCDR3']
				annotations[seq_id]['cdr3'] = value
			if 'fwr4' in annotations_to_get:
				if i['nSeqFR4'] == '':
					value = 'NA'
				else:
					value = i['nSeqFR4']
				annotations[seq_id]['fwr4'] = value
			if 'region_bounds' in annotations_to_get:
				annotations[seq_id]['region_bounds'] = {}
				value = i['refPoints'].split(':')
				if value[4]=='' or value[5]=='':
					annotations[seq_id]['region_bounds']['fwr1'] = ['NA', 'NA']
				else:
					annotations[seq_id]['region_bounds']['fwr1'] = [int(value[4]), int(value[5])]
				if value[5]=='' or value[6]=='':
					annotations[seq_id]['region_bounds']['cdr1'] = ['NA', 'NA']
				else:
					annotations[seq_id]['region_bounds']['cdr1'] = [int(value[5]), int(value[6])]
				if value[6]=='' or value[7]=='':
					annotations[seq_id]['region_bounds']['fwr2'] = ['NA', 'NA']
				else:
					annotations[seq_id]['region_bounds']['fwr2'] = [int(value[6]), int(value[7])]
				if value[7]=='' or value[8]=='':
					annotations[seq_id]['region_bounds']['cdr2'] = ['NA', 'NA']
				else:
					annotations[seq_id]['region_bounds']['cdr2'] = [int(value[7]), int(value[8])]
				if value[8]=='' or value[9]=='':
					annotations[seq_id]['region_bounds']['fwr3'] = ['NA', 'NA']
				else:
					annotations[seq_id]['region_bounds']['fwr3'] = [int(value[8]), int(value[9])]
				if value[9]=='' or value[18]=='':
					annotations[seq_id]['region_bounds']['cdr3'] = ['NA', 'NA']
				else:
					annotations[seq_id]['region_bounds']['cdr3'] = [int(value[9]), int(value[18])]
				if value[18]=='' or value[19]=='':
					annotations[seq_id]['region_bounds']['fwr4'] = ['NA', 'NA']
				else:
					annotations[seq_id]['region_bounds']['fwr4'] = [int(value[18]), int(value[19])]
			if 'fwr1aa' in annotations_to_get:
				if i['aaSeqFR1'] == '':
					value = 'NA'
				else:
					value = i['aaSeqFR1']
				annotations[seq_id]['fwr1aa'] = value
			if 'cdr1aa' in annotations_to_get:
				if i['aaSeqCDR1'] == '':
					value = 'NA'
				else:
					value = i['aaSeqCDR1']
				annotations[seq_id]['cdr1aa'] = value
			if 'fwr2aa' in annotations_to_get:
				if i['aaSeqFR2'] == '':
					value = 'NA'
				else:
					value = i['aaSeqFR2']
				annotations[seq_id]['fwr2aa'] = value
			if 'cdr2aa' in annotations_to_get:
				if i['aaSeqCDR2'] == '':
					value = 'NA'
				else:
					value = i['aaSeqCDR2']
				annotations[seq_id]['cdr2aa'] = value
			if 'fwr3aa' in annotations_to_get:
				if i['aaSeqFR3'] == '':
					value = 'NA'
				else:
					value = i['aaSeqFR3']
				annotations[seq_id]['fwr3aa'] = value
			if 'cdr3aa' in annotations_to_get:
				if i['aaSeqCDR3'] == '':
					value = 'NA'
				else:
					value = i['aaSeqCDR3']
				annotations[seq_id]['cdr3aa'] = value
			if 'fwr4aa' in annotations_to_get:
				if i['aaSeqFR4'] == '':
					value = 'NA'
				else:
					value = i['aaSeqFR4']
				annotations[seq_id]['fwr4aa'] = value
			if 'seqVDJaa' in annotations_to_get:
				if i['aaSeqVDJRegion'] == '':
					value = 'NA'
				else:
					value = i['aaSeqVDJRegion']
				annotations[seq_id]['seqVDJaa'] = value
		return annotations
