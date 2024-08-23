import os, sys
from Bio import Align
import numpy as np
import re

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

def init_aligner(mode='global', open_gap=-1.0, x_gap=-0.1):
	aligner = Align.PairwiseAligner()
	aligner.mode = mode
	aligner.open_gap_score = open_gap
	aligner.extend_gap_score = x_gap
	aligner.target_end_gap_score = 0.0
	aligner.query_end_gap_score = 0.0
	return aligner

def get_boundaries(str):
    gap_prefix = re.compile('^[-]+')
    gap_suffix = re.compile('[-]+$')
    # return a tuple giving indices of subsequence without gap prefix and suffix
    res = [0,len(str)]
    left = gap_prefix.findall(str)
    right = gap_suffix.findall(str)
    if left:
        res[0] = len(left[0])
    if right:
        res[1] = len(str) - len(right[0])
    return res

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
	# cut down the alignment to only the region of interest (that which overlaps with reference)
	start, stop = get_boundaries(ref_aln)
	ref_aln = ref_aln[start:stop]
	qry_aln = qry_aln[start:stop]
	return ref_aln, qry_aln
