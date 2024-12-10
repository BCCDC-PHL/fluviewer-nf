#%%
import os, sys, re
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
from Bio.SeqRecord import SeqRecord

BASE_DICT = dict(zip(list("ACGT"), range(4)))

def init_aligner(mode='global', match=2, mismatch=-1, open_gap=-10, x_gap=-0.5):
	aligner = Align.PairwiseAligner()
	aligner.mode = mode
	aligner.match_score = match
	aligner.mismatch_score = mismatch
	aligner.open_gap_score = open_gap
	aligner.extend_gap_score = x_gap
	aligner.target_end_gap_score = 0.0
	aligner.query_end_gap_score = 0.0
	return aligner

def get_boundaries(str):
    gap_prefix = re.compile('^[-]+')
    gap_suffix = re.compile('[-]+$')
    # return a tuple giving indices of subsequence without gap prefix and suffix
    start, stop = 0, len(str)
    left = gap_prefix.findall(str)
    right = gap_suffix.findall(str)
    if left:
        start = len(left[0])
    if right:
        stop = len(str) - len(right[0])
    return start, stop

def flex_translate(nt_seq, debug=False):
	"""
	Function to find the most appropriate amino acid translation by testing all three reading frames.
	Designed to receive Biopython SeqRecord objects.
	Returns three outputs as a tuple:
		- A SeqRecord object containing the best amino acid translation 
		- An integer indicating the best reading frame that was used for translation (0, 1, or 2)
		- An integer indicating the number of stop codons in this best translation candidate
	"""
	min_count = np.inf
	best_seq = None
	best_frame = -1

	def pad(seq):
		"""
		Supporting function to prevent Biopython translate() function from complaining about non-multiples-of-3
		"""
		mod = len(seq) % 3
		if mod != 0:
			toadd = 3 - mod
			seq += "N" * toadd
		return seq

	for i in range(3):
		nt_tmp = pad(nt_seq[i:])
		aa_seq = nt_tmp.translate()
		aa_count = aa_seq.count("*")
		
		if debug: 
			print(aa_seq)
	
		if aa_count < min_count:
			min_count = aa_count
			best_seq = aa_seq
			best_frame = i

	best_seq.id = nt_seq.id
	best_seq.name = nt_seq.name
	best_seq.description = nt_seq.description
	return best_seq, best_frame, min_count


# performs a pairwise alignment between two Biopython SeqRecord objects
def pairwise_alignment(aligner, ref, qry):
	"""
	Arguments
	- `aligner`: an object of the Biopython alignment tool to be used for alignment
	- `ref`: a Bio.SeqRecord object representing the reference sequence
	- `qry`: a Bio.SeqRecord object representing the query sequence

	Returns
	- `ref_aln`: a string representing the aligned reference sequence
	- `qry_aln`: a string representing the aligned query sequence

	The function performs a pairwise alignment between the reference and query sequences using the provided `aligner` object from the Biopython alignment tool. 
	The resulting alignment is then cut down to only the region of interest that overlaps with the reference sequence, as determined by the `get_boundaries()` function. 
	The aligned reference and query sequences are returned as strings.
	"""
	ref_aln, qry_aln = next(aligner.align(str(ref.seq), str(qry.seq)))

	start, stop = get_boundaries(ref_aln)
	ref_aln = ref_aln[start:stop]
	qry_aln = qry_aln[start:stop]
	return ref_aln, qry_aln


def get_mutations(ref, qry):
	mutations = []
	insertion = ''
	deletion = ''

	ref_pos = 0
	qry_pos = 0
	for ref_char, qry_char in zip(ref, qry):
		# increment the ref_pos at any valid reference position
		if ref_char != '-':  
			ref_pos += 1
		if qry_pos != '-':
			qry_pos += 1

		if ref_char == '-':     # insertion
			if deletion: 
				mutations += [(deletion, ref_pos-1, "-")]
				deletion = ""	
			insertion += qry_char
		elif qry_char == '-':   # deletion
			if insertion:
				mutations += [("-", ref_pos-1, insertion)]
				insertion = ""
			deletion += ref_char

		else:					# neither insertion nor deletion
			if insertion:
				mutations += [("-", ref_pos-1, insertion)]
				insertion = ""
			if deletion: 
				mutations += [(deletion, ref_pos-1, "-")]
				deletion = ""

			if ref_char != qry_char and qry_char not in ['X', 'N'] and ref_char not in ['X', 'N']: # standard snp mismatch 
				mutations += [(ref_char, ref_pos, qry_char)]

	return mutations

def get_nt_mutations(ref, qry):
	mutations = []
	insertion = ''
	deletion = ''

	ref_pos = 0
	qry_pos = 0
	for ref_char, qry_char in zip(ref, qry):
		# increment the ref_pos at any valid reference position
		if ref_char != '-':  
			ref_pos += 1
		if qry_pos != '-':
			qry_pos += 1

		if ref_char == '-':     # insertion
			if deletion: 
				mutations += [(deletion, ref_pos-1, "-", qry_pos)]
				deletion = ""	
			insertion += qry_char
		elif qry_char == '-':   # deletion
			if insertion:
				mutations += [("-", ref_pos-1, insertion, qry_pos)]
				insertion = ""
			deletion += ref_char

		else:					# neither insertion nor deletion
			if insertion:
				mutations += [("-", ref_pos-1, insertion, qry_pos)]
				insertion = ""
			if deletion: 
				mutations += [(deletion, ref_pos-1, "-", qry_pos)]
				deletion = ""

			if ref_char != qry_char and qry_char not in ['X', 'N'] and ref_char not in ['X', 'N']: # standard snp mismatch 
				mutations += [(ref_char, ref_pos, qry_char, qry_pos)]

	return mutations

# def get_nt_mutations(ref, qry, pileup_dict):
# 	mutations = []
# 	insertion = ''
# 	deletion = ''

# 	ref_pos = 0
# 	qry_pos = 0
# 	for ref_char, qry_char in zip(ref, qry):
# 		# increment the ref_pos at any valid reference position
# 		if ref_char != '-':  
# 			ref_pos += 1
# 		if qry_pos != '-':
# 			qry_pos += 1

# 		if ref_char == '-':     # insertion
# 			if deletion: 
# 				mutations += [(deletion, ref_pos-1, "-", '.')]
# 				deletion = ""	
# 			insertion += qry_char
# 		elif qry_char == '-':   # deletion
# 			if insertion:
# 				mutations += [("-", ref_pos-1, insertion, '.')]
# 				insertion = ""
# 			deletion += ref_char

# 		else:					# neither insertion nor deletion
# 			if insertion:
# 				mutations += [("-", ref_pos-1, insertion, '.')]
# 				insertion = ""
# 			if deletion: 
# 				mutations += [(deletion, ref_pos-1, "-", '.')]
# 				deletion = ""

# 			if ref_char != qry_char and qry_char not in ['X', 'N'] and ref_char not in ['X', 'N']: # standard snp mismatch 
# 				pileup = pileup_dict[qry_pos]
# 				vaf = pileup[BASE_DICT[qry_char]] / pileup[4] if pileup[4] != 0 else 0
# 				mutations += [(ref_char, ref_pos, qry_char, round(vaf, 2), pileup[4], pileup[5])]

# 	return mutations