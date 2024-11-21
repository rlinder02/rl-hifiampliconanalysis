#!/usr/bin/env python3
# 4-space indented, v0.0.1
# File name: primer_boundaries.py
# Description: Find the primer boundaries on the reference fasta sequence for subsequent plotting.
# Author: Robert Linder
# Date: 2024-10-25

import argparse
from subprocess import Popen, PIPE, check_output
import os
import re
import glob
from Bio import SeqIO
from Bio import Align
from Bio.Seq import Seq
from Bio.Seq import reverse_complement


def parse_args():
	"""this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
	parser = argparse.ArgumentParser(description="Filter out reads without both the forward and reverse primer sequences")
	parser.add_argument("ref_fasta", type=str, help="path to the reference fasta file")
	parser.add_argument('-p1', '--primer1', type=str, help="The primer1 sequence")
	parser.add_argument('-p2', '--primer2', type=str, help="The primer2 sequence")
	args = parser.parse_args()
	return args

def primer_bounds(fasta, primer1, primer2, file_name):
    """Find the primer boundaries on the reference fasta sequence"""
    with open(f"{file_name}.txt", 'w') as outfile:
        primer1_seq = Seq(primer1)
        primer2_seq = Seq(primer2)
        primer1_rc = primer1_seq.reverse_complement()
        primer2_rc = primer2_seq.reverse_complement()
        ref_fasta = SeqIO.parse(fasta, "fasta")
        boundaries = []
        for record in ref_fasta:
            if primer1_seq in record.seq:
                for match in re.finditer(str(primer1_seq), str(record.seq)):
                    start_pos = match.start() + 1
            elif primer1_rc in record.seq:
                for match in re.finditer(str(primer1_rc), str(record.seq)):
                    start_pos = match.end()
            if primer2_seq in record.seq:
                for match in re.finditer(str(primer2_seq), str(record.seq)):
                    end_pos = match.start() + 1
            elif primer2_rc in record.seq:
                for match in re.finditer(str(primer2_rc), str(record.seq)):
                    end_pos = match.end()
            boundaries.append([start_pos, end_pos])
        outfile.write(f"{boundaries[0][0]}\n{boundaries[0][1]}")

def main():
    inputs = parse_args()
    fasta =  inputs.ref_fasta
    primer1 = inputs.primer1
    primer2 = inputs.primer2
    file_name = fasta.split('/')[-1].split('.')[0]
    primer_bounds(fasta, primer1, primer2, file_name)

if __name__ =="__main__":
    main()