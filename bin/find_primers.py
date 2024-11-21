#!/usr/bin/env python3
# 4-space indented, v0.0.1
# File name: find_primers.py
# Description: Use biopython to filter out reads that don't have the forward and reverse primer sequences.
# Author: Robert Linder
# Date: 2024-10-09

import argparse
from subprocess import Popen, PIPE, check_output
import os
import glob
from Bio import SeqIO
from Bio import Align
from Bio.Seq import Seq
from Bio.Seq import reverse_complement
from concurrent.futures import ProcessPoolExecutor


def parse_args():
	"""this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
	parser = argparse.ArgumentParser(description="Filter out reads without both the forward and reverse primer sequences")
	parser.add_argument("fasta", type=str, help="path to the aligned fasta file")
	parser.add_argument('-p1', '--primer1', type=str, help="The primer1 sequence")
	parser.add_argument('-p2', '--primer2', type=str, help="The primer2 sequence")
	parser.add_argument('-t', '--threads', type=int, default=8, help="The number of threads to use to parallelize processing")
	parser.add_argument('-m', '--mismatches', type=int, default=2, help="The number of mismatches and/or gaps allowed between each primer and a read")
	args = parser.parse_args()
	return args

def primer_filter(fasta_chunk, primer1, primer2, mismatches, file_name):
    """Filter out reads without both primers (allowing at most -m mismatches) using pairwise local alignments"""
    with open(file_name, 'w') as outfile:
        print("got here")
        primer1_seq = Seq(primer1)
        primer2_seq = Seq(primer2)
        primer1_score_threshold = len(primer1_seq) - 2*mismatches
        primer2_score_threshold = len(primer2_seq) - 2*mismatches
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        # the scoring system ensures that only two gaps and/or mismatches are allowed between each primer sequence and the read
        aligner.match_score = 1.0
        aligner.mismatch_score = -1
        aligner.gap_score = -1
        for record in fasta_chunk:
            find_primer1_score = aligner.score(record, primer1_seq)
            find_primer1_score_rc = aligner.score(record, reverse_complement(primer1_seq))
            if find_primer1_score >= primer1_score_threshold or find_primer1_score_rc >= primer1_score_threshold:
                find_primer2_score = aligner.score(record, primer2_seq)
                find_primer2_score_rc = aligner.score(record, reverse_complement(primer2_seq))
                if find_primer2_score >= primer2_score_threshold or find_primer2_score_rc >= primer2_score_threshold:
                    outfile.write(f">{str(record.id)}\n{str(record.seq)}\n")

def chunk_iterator(fasta_file, chunk_size):
    """Yield chunks of 100,000 reads from the FASTA file."""
    chunk = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        chunk.append(record)
        if len(chunk) == chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk

def main():
    inputs = parse_args()
    fasta =  inputs.fasta
    primer1 = inputs.primer1
    primer2 = inputs.primer2
    mismatches = inputs.mismatches
    num_processes = inputs.threads
    merged_file_name = fasta.split('/')[-1].split('.')[0]
    file_counter = 0
    # create a generator object to iterate through below
    cmd = ["grep", "-c",  ">", fasta]
    res = check_output(cmd, universal_newlines=True)
    calculate_chunks = int(res) // int(num_processes)
    chunk_maker = chunk_iterator(fasta, calculate_chunks)
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        for chunk in chunk_maker:
            file_counter += 1
            temp_file = f"temp_output_{file_counter}.fasta"
            print(temp_file)
            executor.submit(primer_filter, chunk, primer1, primer2, mismatches, temp_file)

    print("INFO: Merging temp files")
    with open(f"{merged_file_name}_filtered.fasta", 'w') as outfile:
        for fname in glob.glob('temp*.fasta'):
            with open(fname, 'r') as infile:
                for line in infile:
                    outfile.write(line)

    for fname in glob.glob('temp*.fasta'):
        os.remove(fname)       

if __name__ =="__main__":
    main()