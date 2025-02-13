#!/usr/bin/env Rscript

##############################################################################

# Project:     Generate consensus clusters of amplicons
# Author:      Robert Linder
# Date:        2024-10-14
# Title:       generate_clusters_consensus
# Description: Generates clusters of amplicons based on sequence similarity, aligns sequences within clusters, then uses the alignments to find a consensus sequence for each cluster using the Decipher package. 
# Version:     1.0.0

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Usage: generate_clusters_consensus.R <fasta> <threads> <bounds>", call.=FALSE)
}

fasta <-  args[1]
threads <- as.numeric(args[2])
bounds <- args[3]

# ============================================================================
# Load packages and sourced files
library(DECIPHER)
library(data.table)

# ============================================================================
# Set global options

options(digits = 10)
projectDir <- getwd()

# ============================================================================
# Custom functions

clusterize_recurse <- function(dna, cutoff, threads, minCov) {
  c1 <- Clusterize(dna, cutoff=cutoff, processors=threads, penalizeGapLetterMatches = TRUE, minCoverage = -minCov)
  cluster_list <- lapply(sort(unique(c1$cluster)), function(clust) {
    rownames(c1)[c1$cluster==clust]
  })
  # cap at 370 clusters per amplicon sequenced by recursively calling the Clusterize function
  if(length(cluster_list) > 70 & cutoff < 0.5) {
    cutoff <- cutoff + 0.05
    print("Next iteration")
    print(length(cluster_list))
    print(cutoff)
    flush.console()
    clusterize_recurse(dna, cutoff, threads, minCov)
  } else {
    print(paste0("new length is ", length(cluster_list)))
    flush.console()
    return(cluster_list)
  }
}
# ============================================================================
# Load data
dna <- readDNAStringSet(fasta)
bounds_dt <- fread(bounds)
amplicon_length <- as.numeric(bounds_dt[2]) - as.numeric(bounds_dt[1])
# filter out sequences longer than the length of the PCR product (calculated from primer boundaries), as these are likely concatamers formed during PCR
dna <- dna[width(dna) <= amplicon_length]
output_name <- strsplit(fasta, "/")[[1]]
output_name <- output_name[[length(output_name)]]
output_name <- strsplit(output_name, "\\.")[[1]][1]

# ============================================================================
# use clusterize to cluster similar sequences, with a similarity cutoff based off the median sequence length
set.seed(123)
# cluster sequences that are at least 90% or more similar to one another 
cutoff <- 0.1
# query sequences must overlap the cluster representative with at most 100 bases that don't align to cluster together
min_cov <- round((amplicon_length - 100)/amplicon_length, 2)

cluster_list <- clusterize_recurse(dna, cutoff, threads, min_cov)
print(paste0("Final length of clusters is ", length(cluster_list)))
flush.console()
# ============================================================================
# Iterate through clusters, write out each cluster to a fasta file if there are at least 5 reads in that cluster

counter <- 0
denote_clusters <- unlist(DNAStringSetList(lapply(cluster_list, function(clust) {
  counter <<- counter + 1
  print(counter)
  flush.console()
  cluster_seqs <- dna[clust]
  if(length(cluster_seqs) >= 5) {
    names(cluster_seqs) <- paste0(names(cluster_seqs), "_cluster", counter)
    print(counter)
    flush.console()
    cluster_seqs
  }
})))
writeXStringSet(denote_clusters, file = paste0(output_name, "_consensus.fasta"))