#!/usr/bin/env python3 
# 4-space indented, v0.0.1
# File name: tag_clusters.py
# Description: Tag reads for each cluster.
# Author: Robert Linder
# Date: 2024-10-22

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
    
def add_tag(bam, file_name):
    with pysam.AlignmentFile(bam, "rb") as bamfile:
        with pysam.AlignmentFile(f"{file_name}_modified.bam", "wb", template=bamfile) as out_file:
            for read in bamfile.fetch():
                cluster_id = read.query_name.split('_')[3]
                read.set_tag("CL", cluster_id, "Z")
                out_file.write(read)
            
def main():
    inputs = parse_args()
    bam =  inputs.bam_file
    pysam.index(bam)
    file_name = bam.split('/')[-1].split('.')[0]
    add_tag(bam, file_name)

if __name__ =="__main__":
    main()