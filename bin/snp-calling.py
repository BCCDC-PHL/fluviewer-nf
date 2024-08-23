#!/usr/bin/env python3
#%%
import argparse
import os, sys, re
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from tools import flex_translate, init_aligner, pairwise_alignment
import yaml
from subprocess import run, PIPE
from glob import glob

script_dir = os.path.dirname(__file__)

EXPR_SEGMENT = re.compile("\{[^\}]+\}")

def path_check(path):
	if EXPR_SEGMENT.search(path):
		return path
	else:
		raise argparse.ArgumentTypeError('ERROR: Output TSV path must include {segment} to indicate the segment variable.')

def init_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-q', '--query', required=True, help='FASTA consensus sequence containing all query sequences. Segment MUST be included in the header.')
	parser.add_argument('-r', '--reference_db', required=True, help='FASTA reference database containing references for all segments x subtypes')
	# parser.add_argument('-p', '--positions', required=True, help='YAML file containing gene positions for all segments and subtypes')
	parser.add_argument('-b', '--blast', required=True, help='Results from running BLASTx of consensus seqs (queries) against reference typing database')
	parser.add_argument('-s', '--header_pos_seg',  default=1, help='Zero-indexed position of the segment in the consensus headers')
	parser.add_argument('-t', '--header_pos_st', default=2, help='Zero-indexed position of the subtype in the consensus headers')
	parser.add_argument('-d', '--header_delim', default='|', help='FASTA reference database containing references for all segments x subtypes')
	parser.add_argument('-o', '--output', type=path_check, help='Mutation outputs in TSV format. Requires {segment} marker to indicate the gene name in the output.') 
	parser.add_argument('-O', '--outaln', type=path_check, help='Pairwise alignments used to call mutations. Requires {segment} marker to indicate the gene name in the output.') 
	return parser

class MissingSubtypeException(Exception):
    def __init__(self, message):            
        super().__init__(message)

def get_splice_product(ntseq, product_name, pos_list):
	# convert this to a 0-indexed position list, (top index stays the same)
	pos_list = [[x[0]-1, x[1]] for x in pos_list]   

	# initialize a final SeqRecord object 
	final_seq = SeqRecord(
		Seq(''),
		id=ntseq.id + '_' + product_name,
		name=ntseq.id + '_' + product_name,
		description=ntseq.id + '_' + product_name
	)

	# append each translated segment to the final output sequence 
	for start, end in pos_list:
		segment = ntseq[start: end]
		# assert len(segment) % 3 == 0, print("ERROR: not a multiple of three")  # for testing purposes 

		aaseq, frame, stop_count = flex_translate(segment)
		final_seq += aaseq

	return final_seq

def get_position_dict(gff_path):
	"""
	A function to read all GFF3 files in a given path and return a dictionary in the following format:
	pos_dict[SUBTYPE][GENE_PRODUCT] = [(START, END),(START, END)]
	"""
	cols = 'seqid source type start end score strand frame attributes'.split()
	
	pos_dict = defaultdict(dict)

	for filename in glob(os.path.join(gff_path, "*.gff3")):
		
		subtype = os.path.basename(filename).split("_")[0]

		df = pd.read_csv(filename,sep='\t',comment="#", header=None, names=cols)
		df = df.loc[df['type']=='CDS'].reset_index(drop=True)

		attrib_dict = df['attributes'].apply(lambda x : dict([i.split("=") for i in x.split(";")]))
		attrib_df = pd.DataFrame.from_records(attrib_dict)
		final = pd.concat([df.drop('attributes',axis=1), attrib_df],axis=1)

		positions = final.groupby("gene")[['start','end']].agg(list)

		for name, start, end in positions.itertuples():
			pos_dict[subtype][name] = list(zip(start, end))

	return pos_dict

def parse_blast(blast_path):
	cols = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

	blastdf = pd.read_csv(blast_path, sep='\t', names=cols.split())

	idxmax = blastdf.groupby(['qseqid'])['bitscore'].idxmax()
	filtered = blastdf.loc[idxmax]

	return filtered.set_index("qseqid")['sseqid'].to_dict()

def get_mutations(ref, qry):
	mutations = []
	insertion = []
	deletion = []

	ref_pos = 0
	for ref_char, qry_char in zip(ref, qry):
		# increment the ref_pos at any valid reference position
		if ref_char != '-':  
			ref_pos += 1

		if ref_char == '-':     # insertion
			if deletion: 
				mutations += [("".join(deletion), ref_pos-1, "-")]
				deletion = []	
			insertion += qry_char
		elif qry_char == '-':   # deletion
			if insertion:
				mutations += [("-", ref_pos-1, "".join(insertion))]
				insertion = []
			deletion += ref_char

		else:					# neither insertion nor deletion
			if insertion:
				mutations += [("-", ref_pos-1, "".join(insertion))]
				insertion = []
			if deletion: 
				mutations += [("".join(deletion), ref_pos-1, "-")]
				deletion = []

			if ref_char != qry_char and qry_char != "X": 	# standard snp mismatch 
				mutations += [(ref_char, ref_pos, qry_char)]

	return mutations


def exit(outpath):
	mutations_df = pd.DataFrame([], columns=['CHROM','POS','REF','ALT'])
	mutations_df.to_csv(outpath, sep='\t', index=False)

#%%
def main():
	parser = init_parser()
	args = parser.parse_args()	

	# load the query sequences into a dict structure
	# sequence = qry_seqs[SEGMENT]
	qry_seqs = list(SeqIO.parse(args.query, 'fasta'))

	ref_database = SeqIO.to_dict(SeqIO.parse(args.reference_db, 'fasta'))

	aligner = init_aligner()

	ref_dict = parse_blast(args.blast)

	for seq in qry_seqs:

		# something went wrong in this case; skip for now 
		if seq.id not in ref_dict:
			print("No BLAST result found for ", seq.id)
			continue
		
		# retrieve the name of the most appropriate reference (based on a BLAST alignment)
		ref_name = ref_dict[seq.id]
		
		segment = ref_name.split("|")[0]

		print("Segment : ", segment)

		aa_qry, _ , _ = flex_translate(seq)

		# grab the appropriate reference sequence for the subtype
		aa_ref = ref_database[ref_name]

		# compute all mutations using the pairwise alignment
		ref_aln, qry_aln = pairwise_alignment(aligner, aa_ref, aa_qry)

		if args.outaln:
			outaln = EXPR_SEGMENT.sub(segment, args.outaln)
			with open(outaln, 'w') as outfile:
				outfile.write(f">{ref_name}\n{ref_aln}\n>{seq.id}\n{qry_aln}\n")

		mutations = get_mutations(ref_aln, qry_aln)

		# convert list of mutations to data frame
		mutations_df = pd.DataFrame(mutations, columns=['POS','REF','ALT'])
		mutations_df.insert(0, "CHROM", aa_ref.name)

		# format output path 
		# outname = seq.id + "_" + segment + '_mutations.tsv'

		outpath = EXPR_SEGMENT.sub(segment, args.output)

		# output mutations to CSV file 
		mutations_df.to_csv(outpath, sep='\t', index=False)

		print(f'Completed segment {segment}')



if __name__ == '__main__':
	main()
# %%
