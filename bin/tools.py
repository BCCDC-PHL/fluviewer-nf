#%%
import os, sys, re
import numpy as np
import pandas as pd 
from typing import Tuple, List
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
from Bio.SeqRecord import SeqRecord
from functools import partial
import gzip 

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


SEGMENTS = 'PB2 PB1 PA HA NP NA M NS'.split()

def read_vcf(file_path: str , max_header_lines=300) -> Tuple[pd.DataFrame, List[str]]:
    """
    Parse a Variant Call Format (VCF) file into a pandas DataFrame and header list.
    
    Args:
        file_path: Path to the VCF file
        
    Returns:
        Tuple containing:
            - pandas DataFrame with VCF data
            - List of header lines (including metadata)
    """
    # Store header lines
    header_lines = []
    
    # Read header lines
    read_function = partial(open, mode='r') if not file_path.endswith('.gz') else partial(gzip.open, mode='rt')

    with read_function(file_path) as infile:
        line_count = 0 
        line = infile.readline()
        while not line.lstrip('#').startswith("CHROM"):
            line_count += 1
            header_lines.append(line.strip())
            line = infile.readline()

            if line_count > max_header_lines:
                print(f"ERROR: Exceeded {max_header_lines} lines to search for the CHROM start point.")
                sys.exit(1)

        header_lines.append(line.strip())
        
        column_names = header_lines[-1].lstrip("#").split("\t")

        df = pd.read_csv(infile, sep='\t', names=column_names)
    
    return df, header_lines


def write_vcf(df, header_lines, output_path):
    """
    Write out a Variant Call Format (VCF) file using a Pandas df, list of headers, and output path
    
    Args:
        file_path: Path to the VCF file
        
    Returns:
        Tuple containing:
            - pandas DataFrame with VCF data
            - List of header lines (including metadata)
    """

    # Read header lines
    with open(output_path, 'w') as outfile:
        for line in header_lines:
            outfile.write(line + '\n')

        df.to_csv(outfile, sep='\t', index=False, header=False)