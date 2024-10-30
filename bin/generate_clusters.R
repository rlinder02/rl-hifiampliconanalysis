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
if (length(args) < 2) {
  stop("Usage: generate_clusters_consensus.R <fasta> <threads>", call.=FALSE)
}

fasta <-  args[1]
threads <- as.numeric(args[2])

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

clusterize_recurse <- function(dna, cutoff, threads) {
  c1 <- Clusterize(dna, cutoff=cutoff, processors=threads, penalizeGapLetterMatches = TRUE)
  cluster_list <- lapply(sort(unique(c1$cluster)), function(clust) {
    rownames(c1)[c1$cluster==clust]
  })
  # cap at 10 clusters per amplicon sequenced by recursively calling the Clusterize function
  if(length(cluster_list) > 10 & cutoff <= 0.9) {
    cutoff <- cutoff + 0.05
    print("Next iteration")
    print(cutoff)
    flush.console()
    clusterize_recurse(dna, cutoff, threads)
  } else {
    return(cluster_list)
  }
  print("Finished clustering")
  flush.console()
  print(length(cluster_list))
  flush.console()
  cluster_list
}
# ============================================================================
# Load data

#fasta <- "/Users/rlinder/Library/CloudStorage/OneDrive-SanfordBurnhamPrebysMedicalDiscoveryInstitute/Chun_lab/Pipelines/HiFi_PCR_analysis/tests/HU_PCR_MAPT_AMPLICON_FTLD_CBD_POSITIVE_filtered.fasta"
dna <- readDNAStringSet(fasta)
output_name <- strsplit(fasta, "/")[[1]]
output_name <- output_name[[length(output_name)]]
output_name <- strsplit(output_name, "\\.")[[1]][1]

# ============================================================================
# use clusterize to cluster similar sequences, with a similarity cutoff based off the median sequence length
set.seed(123)

median_width <- median(width(dna)) # round to the nearest hundredth
#cutoff_dt <- data.table(bp_start = c(0, 501, seq(1001, 15001, 1000)), bp_end = c(500, 1000, seq(2000, 16000, 1000)), cutoff = c(0.1, seq(0.2, 0.95, 0.05)))
cutoff_dt <- data.table(bp_start = c(0, 501, seq(1001, 5001, 1000)), bp_end = c(500, 1000, seq(2000, 6000, 1000)), cutoff = c(0.1, seq(0.2, 0.95, 0.15)))
if(median_width > cutoff_dt$bp_end[nrow(cutoff_dt)]) {
  cutoff <- 0.95
} else {
  find_cutoff <- cutoff_dt[median_width %between% list(bp_start, bp_end)]
  cutoff <- find_cutoff$cutoff
}
# cluster sequences that are at least x% or higher similar to one another based on the median length of the input sequences
cluster_list <- clusterize_recurse(dna, cutoff, threads)

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
    cluster_seqs
  }
})))
writeXStringSet(denote_clusters, file = paste0(output_name, "_consensus.fasta"))