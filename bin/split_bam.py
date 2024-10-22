#!/usr/bin/env python3 
# 4-space indented, v0.0.1
# File name: split_bam.py
# Description: Split bam file into a separate bam file for each cluster.
# Author: Robert Linder
# Date: 2024-10-21

import re
import os
import argparse
import pysam

def parse_args():
    """this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
    parser = argparse.ArgumentParser(description="Load bam file to separate by cluster")
    parser.add_argument("bam_file", type=str, help="Sorted bam file to separate aligned cluster from")
    args = parser.parse_args()
    return args
    
def find_reads(bam):
    read_dict = {}
    with pysam.AlignmentFile(bam, "rb") as bamfile:
        for read in bamfile.fetch():
            cluster_id = read.query_name.split('_')[3]
            if cluster_id in read_dict:
                read_dict[cluster_id].append(read.query_name)
            else:
                read_dict[cluster_id] = [read.query_name]
        return read_dict

def extract_reads(cluster_dict, bam, file_name):
    with pysam.AlignmentFile(bam, "rb") as bamfile:
        for read in bamfile.fetch():
            cluster_id = read.query_name.split('_')[3]
            if read.query_name in cluster_dict.get(cluster_id):
                with pysam.AlignmentFile(f"{file_name}_{cluster_id}.bam", "wb", template=bamfile) as bam_file:
                    bam_file.write(read)
def main():
    inputs = parse_args()
    bam =  inputs.bam_file
    pysam.index(bam)
    file_name = bam.split('/')[-1].split('.')[0]
    cluster_dict = find_reads(bam)
    extract_reads(cluster_dict, bam, file_name)

if __name__ =="__main__":
    main()
