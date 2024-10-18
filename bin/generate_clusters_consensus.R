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
library(msa)


# ============================================================================
# Set global options

options(digits = 10)
projectDir <- getwd()

# ============================================================================
# Custom functions

# ============================================================================
# Load data

#fasta <- "/Users/rlinder/Library/CloudStorage/OneDrive-SanfordBurnhamPrebysMedicalDiscoveryInstitute/Chun_lab/Pipelines/HiFi_PCR_analysis/tests/HU_PCR_MAPT_AMPLICON_FTLD_CBD_POSITIVE_filtered.fasta"
dna <- readDNAStringSet(fasta)
output_name <- strsplit(fasta, "/")[[1]]
output_name <- output_name[[length(output_name)]]
output_name <- strsplit(output_name, "\\.")[[1]][1]

# ============================================================================
# use clusterize to cluster similar sequences
set.seed(123)
# cluster sequences that are at least 90% or higher (cutoff of 0.1) similar to one another and make the centers negative to find them later
c1 <- Clusterize(dna, cutoff=0.1, processors=threads, penalizeGapLetterMatches = TRUE, invertCenters = TRUE) # use invertCenters=TRUE to find a representative for each cluster as below
w <- which(c1 < 0 & !duplicated(c1))
c1$reads <- rownames(c1)
rownames(c1) <- 1:nrow(c1)
print(c1)
print(w)
selected <- c1[w,]
print(c1[w,])
b <- match(c1$cluster, abs(selected$cluster))
print(c1$cluster)
print(abs(selected$cluster))
print(b)
flush.console()

# ============================================================================
# Iterate through clusters, align each cluster, then find the consensus sequence within each aligned cluster 
cluster_list <- lapply(sort(unique(c1$cluster)), function(clust) {
  rownames(c1)[c1$cluster==clust]
})

counter <- 0
align_seqs <- unlist(DNAStringSetList(lapply(cluster_list, function(clust) {
  counter <<- counter + 1
  print(counter)
  flush.console()
  cluster_seqs <- dna[clust]
  unique_seqs <- unique(cluster_seqs)
  index <- match(cluster_seqs, unique_seqs)
  if(length(unique_seqs) == 1) {
    unique_seqs[2] <- unique_seqs[1]
    names(unique_seqs[2]) <- names(unique_seqs[1])
  }
  if(length(unique_seqs) > 1000) {
    # recluster a cluster that is too large for quick multiple sequence alignment, using a 1% cutoff to split the cluster and then pick each sub-cluster center for alignment 
    recluster <- Clusterize(unique_seqs, cutoff=0.05, processors=threads, penalizeGapLetterMatches = TRUE, invertCenters = TRUE)
    # reclustered_lengths <- lapply(unique(recluster$cluster), function(reclust) {
    #   clust_members <- length(which(recluster$cluster == reclust))
    #   data.frame(cluster = reclust, members = clust_members)
    # })
    # reclustered_lengths_df <- do.call('rbind', reclustered_lengths)
    # print(reclustered_lengths_df)
    # flush.console()
    centers <- which(recluster < 0 & !duplicated(recluster))
    reduced_seqs <- unique_seqs[centers]
    # all members of a cluster will be collapsed to having a single representative sequence here
    index2 <- match(reduced_seqs, unique_seqs)
    print(unique_seqs)
    print(reduced_seqs)
    print(index2)
    print(length(reduced_seqs))
    flush.console()
    aligned_seqs <- AlignSeqs(reduced_seqs, verbose=FALSE, processors=threads)
    all_unique_seqs <- aligned_seqs[index2]
    names(all_unique_seqs) <- names(unique_seqs)
    all_seqs <- all_unique_seqs[index]
    names(all_seqs) <- names(cluster_seqs)
  } else {
    aligned_seqs <- AlignSeqs(unique_seqs, verbose=FALSE, processors=threads)
    all_seqs <- aligned_seqs[index]
    print(aligned_seqs)
    names(all_seqs) <- names(cluster_seqs)
    print(all_seqs)
    flush.console()
  }
  print("Completed alignment")
  flush.console()
  find_consensus <- ConsensusSequence(all_seqs)
  #print(find_consensus)
  find_consensus
})))
#print(str(align_seqs))
print(align_seqs)
writeXStringSet(align_seqs, file = paste0(output_name, "_consensus.fasta"))