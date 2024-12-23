#!/usr/bin/env Rscript

##############################################################################

# Project:     Make linear plot of consensus amplicon sequences
# Author:      Robert Linder
# Date:        2024-12-16
# Title:       generate_linear_plots
# Description: Generates linear plots of each consensus amplicon species using the ggtranscript package. 
# Version:     1.0.1

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
  stop("Usage: generate_circos_plots.R <vcfs> <bed> <bounds> <total_reads> <file_name> <orfs>", call.=FALSE)
}

vcfs <- args[1]
bed <- args[2]
bounds <- args[3]
total_reads <- args[4]
gene_name <- args[5]
orfs <- args[6]

# ============================================================================
# For trouble-shooting locally

# vcfs <- "vcf_fofn.txt"
# bed <- "hTARDBP_cDNA_full.bed"
# bounds <- "hTARDBP_cDNA.txt"
# total_reads <- "total_reads_fofn.txt"
# gene_name <- "TARDBP"
# orfs <- "orf_fofn.txt"


# ============================================================================
# Load packages and sourced files

library(data.table)
library(ggtranscript)
library(tidyverse)
library(Polychrome)

# ============================================================================
# Set global options

options(digits = 10)
projectDir <- getwd()

#setwd("/Users/rlinder/Library/CloudStorage/OneDrive-SanfordBurnhamPrebysMedicalDiscoveryInstitute/Chun_lab/Projects/gencDNA/PCR_Southerns/Human/TARDBP/2024-12-13_run")


# ============================================================================
# Custom functions

pre.process.bed <- function(bed_file, bounds_file) {
  ref_bounds_dt <- fread(bounds_file)
  ref_bed_dt <- fread(bed_file)
  coord_shift <- ref_bed_dt$V2[1]
  ref_bed_dt$V1 <- ref_bed_dt$V4
  ref_bed_dt <- ref_bed_dt[, V4 := NULL]
  ref_bed_dt[, V2 := V2 - coord_shift]
  ref_bed_dt[, V3 := V3 - coord_shift]
  ref_bed_dt[, V4 := V3 - V2]
  ref_bed_dt[, V5 := cumsum(V4)]
  ref_bed_dt[, start := c(0, V5[-length(V5)])]
  ref_bed_dt[, end := V5]
  ref_bed_dt[, feature := gsub("CDS_", "", V1)]
  columns <- c("feature", "start", "end")
  ref_bed_dt <- ref_bed_dt[, ..columns]
  # limit plotting to primer bounds; change so limit plotting to CDS only if ORFs not included in primer bounds
  if(ref_bounds_dt$V1[1] < ref_bed_dt$end[1]) {
    start_row <- 1
    ref_bed_dt$start[1] <- ref_bounds_dt$V1[1]
  } else {
    start_row <- which(ref_bed_dt$feature == 1)
  }
  #ref_bed_dt$start[start_row] <- ref_bounds_dt$V1[1]
  if(ref_bounds_dt$V1[2] > ref_bed_dt$start[nrow(ref_bed_dt)]) {
    end_row <- nrow(ref_bed_dt)
    ref_bed_dt$end[nrow(ref_bed_dt)] <- ref_bounds_dt$V1[2]
  } else {
    end_row <- which(ref_bed_dt$feature == max(as.numeric(ref_bed_dt$feature), na.rm = TRUE))
  }
  #ref_bed_dt$end[end_row] <- ref_bounds_dt$V1[2]
  ref_bed_dt <- ref_bed_dt[start_row:end_row]
}

pre.process.bed_wt <- function(bed_file) {
  ref_bed_dt <- fread(bed_file)
  coord_shift <- ref_bed_dt$V2[1]
  ref_bed_dt$V1 <- ref_bed_dt$V4
  ref_bed_dt <- ref_bed_dt[, V4 := NULL]
  ref_bed_dt[, V2 := V2 - coord_shift]
  ref_bed_dt[, V3 := V3 - coord_shift]
  ref_bed_dt[, V4 := V3 - V2]
  ref_bed_dt[, V5 := cumsum(V4)]
  ref_bed_dt[, start := c(0, V5[-length(V5)]+1)]
  ref_bed_dt[, end := V5]
  ref_bed_dt[, feature := gsub("CDS_", "", V1)]
  columns <- c("feature", "start", "end")
  ref_bed_dt <- ref_bed_dt[, ..columns]
  ref_bed_dt[, featureType := ifelse(grepl("UTR", feature), "UTR", "CDS")]
  ref_bed_dt
}

pre.process.vcf.structure <- function(vcf_file, ref_bed_dt) {
  vcf_dt <- fread(vcf_file)
  vcf_dt[, depth := as.numeric(str_match(vcf_dt$INFO, 'DP=(\\d+)')[,2])]
  #vcf_dt <- vcf_dt[depth > 4]
  # subset for regions with coverage - remove sites with 0 coverage for preprocessed pseudogenes, for all other samples, remove sites with less than 5x coverage
  if(grepl("cluster", vcf_file)) {
    vcf_dt <- vcf_dt[depth > 4]
  } else {
    vcf_dt <- vcf_dt[depth > 0]
  }
  vcf_max_depth <- max(as.numeric(str_match(vcf_dt$INFO, 'DP=(\\d+)')[,2]), na.rm = TRUE)
  names(vcf_dt)[10] <- "SAMPLE"
  # subtract one from the vcf file coordinates so is in bed coordinate space (0-based)
  vcf_dt[, POS := POS - 1 ]
  vcf_dt[ref_bed_dt, on=.(POS >= start, POS <= end), feature := i.feature]
  vcf_dt <- vcf_dt[!is.na(feature)]
  # eliminate duplicate rows as, for some reason, indels cause the position to be duplicated (likely has to do with the ref allele being called differently in the indel TCC vs just T, for instance)
  vcf_dt <- vcf_dt[!duplicated(POS, fromLast=TRUE)]
  # add a new column that delineates covered sites
  vcf_dt[, diffs := diff(c(min(POS)-1, POS))]
  rles <- rle(vcf_dt$diffs)
  rles$values <- 1:length(rles$values)
  rles$values[rles$lengths == 1] <- rles$values[rles$lengths == 1] + 1 
  vcf_dt[, runs := rep(rles$values, rles$lengths)]
  vcf_dt[, c("start", "end") := .(min(POS), max(POS)), by = c("runs", "feature")]
  # ends block of code delineating covered sites
  #vcf_dt[, c("start", "end") := .(min(POS), max(POS)), by = feature]
  struct_columns <- c("feature", "start", "end")
  vcf_dt_structure <- unique(vcf_dt[, ..struct_columns])
  vcf_dt_structure[, maxDepth := vcf_max_depth]
  cluster_id <- strsplit(gsub("_modified.*", "", vcf_file), "_")[[1]]
  cluster_id <- cluster_id[length(cluster_id)]
  vcf_dt_structure[, c("sample", "cluster") := .(gsub("_cluster.*|_pp.*", "", vcf_file), cluster_id)]
  vcf_dt_structure
}

pre.process.vcf.mutations <- function(vcf_file, ref_bed_dt) {
  vcf_dt <- fread(vcf_file)
  names(vcf_dt)[10] <- "SAMPLE"
  # subtract one from the vcf file coordinates so is in bed coordinate space (0-based)
  vcf_dt[, POS := POS - 1 ]
  vcf_dt[ref_bed_dt, on=.(POS >= start, POS <= end), feature := i.feature]
  vcf_dt <- vcf_dt[!is.na(feature)]
  vcf_dt <- vcf_dt[!duplicated(POS, fromLast=TRUE)]
  vcf_dt[, c("start", "end") := .(min(POS), max(POS)), by = feature]
  vcf_dt[, total_reads := as.numeric(str_match(INFO, 'DP=(\\d+)')[,2])]
  # keep only positions with called mutations 
  vcf_dt_muts <- vcf_dt[grepl("AC=", INFO) & grepl("PASS", FILTER)]
  vcf_dt_muts[, TYPE := ifelse(grepl("INDEL", INFO), "INDEL", "SNV")]
  vcf_dt_muts[, symbol := ifelse(grepl("SNV", TYPE), 16, 17)]
  vcf_dt_muts[, alt_count := as.numeric(gsub(".*,", "", SAMPLE))]
  vcf_dt_muts[, value := alt_count/total_reads]
  vcf_dt_muts[, POS_END := POS +1]
  muts_columns <- c("feature", "POS", "POS_END", "value", "symbol", "REF", "ALT")
  vcf_dt_muts_dt <- vcf_dt_muts[, ..muts_columns]
  setnames(vcf_dt_muts_dt, old = c("POS", "POS_END"), new = c("start", "end"))
  cluster_id <- strsplit(gsub("_modified.*", "", vcf_file), "_")[[1]]
  cluster_id <- cluster_id[length(cluster_id)]
  vcf_dt_muts_dt[, c("sample", "cluster") := .(gsub("_cluster.*|_pp.*", "", vcf_file), cluster_id)]
  vcf_dt_muts_dt
}

pre.process.orf <- function(orf_file, vcf_file, ref_bed_dt) {
  vcf_dt <- fread(vcf_file)
  vcf_dt[, depth := as.numeric(str_match(vcf_dt$INFO, 'DP=(\\d+)')[,2])]
  # subset for regions with coverage - remove sites with 0 coverage 
  #vcf_dt <- vcf_dt[depth > 0]
  if(grepl("cluster", vcf_file)) {
    vcf_dt <- vcf_dt[depth > 4]
  } else {
    vcf_dt <- vcf_dt[depth > 0]
  }
  # subtract one from the vcf file coordinates so is in bed coordinate space (0-based)
  vcf_dt[, POS := POS - 1 ]
  vcf_dt[ref_bed_dt, on=.(POS >= start, POS <= end), feature := i.feature]
  vcf_dt <- vcf_dt[!is.na(feature)]
  # eliminate duplicate rows as, for some reason, indels cause the position to be duplicated (likely has to do with the ref allele being called differently in the indel TCC vs just T, for instance)
  vcf_dt <- vcf_dt[!duplicated(POS, fromLast=TRUE)]
  # add a new column that delineates covered sites
  vcf_dt[, diffs := diff(c(min(POS)-1, POS))]
  rles <- rle(vcf_dt$diffs)
  rles$values <- 1:length(rles$values)
  rles$values[rles$lengths == 1] <- rles$values[rles$lengths == 1] + 1 
  vcf_dt[, runs := rep(rles$values, rles$lengths)]
  vcf_dt[, c("start", "end") := .(min(POS), max(POS)), by = c("runs", "feature")]
  struct_columns <- c("feature", "start", "end", "runs")
  vcf_dt_structure <- unique(vcf_dt[, ..struct_columns])
  cds_dt <- vcf_dt_structure[!grepl("UTR", feature)]
  orf_dt <- fread(orf_file, header = F, sep = "\t")
  orf_dt_cds <- orf_dt[V2 >= cds_dt$start[1] & V3 <= cds_dt$end[nrow(cds_dt)]]
  orf_dt_cds$length <- as.numeric(str_match(orf_dt_cds$V4, 'ORF_len=(\\d+)')[,2])
  orf_dt_cds_longest <- orf_dt_cds[order(-length)][1]
  # find out if longest CDS orf is in-frame
  in_frame <- FALSE
  if((orf_dt_cds_longest$V2 - (ref_bed_dt$start[ref_bed_dt$feature == "1"][1])) %% 3 == 0) {
    in_frame <- TRUE
  } 
  expanded_dt <- data.table(POS = orf_dt_cds_longest$V2:orf_dt_cds_longest$V3, strand = orf_dt_cds_longest$V6, in_frame = in_frame)
  expanded_dt[vcf_dt_structure, on=.(POS >= start, POS <= end), c("feature", "runs", "in_frame") := .(i.feature, i.runs, in_frame)]
  expanded_dt <- expanded_dt[!is.na(feature)]
  expanded_dt[, c("start", "end") := .(min(POS), max(POS)), by = c("runs", "feature")]
  struct_columns <- c("feature", "start", "end", "strand", "in_frame")
  expanded_dt_struct <- unique(expanded_dt[, ..struct_columns])
  cluster_id <- strsplit(gsub(".bed", "", orf_file), "_")[[1]]
  cluster_id <- cluster_id[length(cluster_id)-2]
  expanded_dt_struct[, c("sample", "cluster") := .(gsub("_cluster.*|_pp.*", "", orf_file), cluster_id)]
  expanded_dt_struct
}

sample.coverage.calc <- function(vcf_file, ref_bed_dt) {
  vcf_dt <- fread(vcf_file)
  vcf_dt[, start_pos := POS - 1 ]
  vcf_dt[, end_pos := start_pos + 1]
  vcf_dt[ref_bed_dt, on=.(POS >= start, POS <= end), feature := i.feature]
  vcf_dt <- vcf_dt[!is.na(feature)]
  vcf_dt[, total_reads := as.numeric(str_match(INFO, 'DP=(\\d+)')[,2])]
  columns <- c("feature", "start_pos", "end_pos", "total_reads")
  vcf_dt <- vcf_dt[, ..columns]
  cluster_id <- strsplit(gsub("_modified.*", "", vcf_file), "_")[[1]]
  cluster_id <- cluster_id[length(cluster_id)]
  vcf_dt[, c("sample", "cluster") := .(gsub("_cluster.*|_pp.*", "", vcf_file), cluster_id)]
  vcf_dt
}

identify.identical.dfs <- function(df_list, columns) {
  df_list <- lapply(df_list, function(df) {df[, ..columns]})
  # Create a list of indices of identical data frames
  identical_groups <- list()
  checked <- rep(FALSE, length(df_list))
  for (i in seq_along(df_list)) {
    if (!checked[i]) {
      identical_indices <- which(map_lgl(df_list, ~ identical(df_list[[i]], .)))
      identical_groups <- append(identical_groups, list(identical_indices))
      checked[identical_indices] <- TRUE
    }
  }
  identical_groups
}

find.common.values <- function(lst1, lst2) {
  common_values <- list()
  # Iterate over each element in the first list
  for (element1 in lst1) {
    # Iterate over each element in the second list
    for (element2 in lst2) {
      find_intersection <- intersect(element1, element2)
      if (length(find_intersection) > 1) {
        common_values <- append(common_values, list(find_intersection))
      }
    }
  }
  return(unique(common_values))
}

vcf.read.depth <- function(vcf_file) {
  vcf_dt <- fread(vcf_file)
  vcf_dt_depth <- max(as.numeric(str_match(vcf_dt$INFO, 'DP=(\\d+)')[,2]), na.rm = TRUE)
  vcf_dt_depth
}

rescale <- function(x, old_min, new_min, old_max, new_max) {
  rescaled_value <- new_min + ((x - old_min) / (old_max - old_min)) * (new_max - new_min)
  rescaled_value
}

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha)) 
}   

keep.unique.vals <- function(column, delimiter) {
  sapply(strsplit(column, delimiter, fixed = TRUE), function(x) paste(unique(x), collapse = ";"))
}
# ============================================================================
# Load data

vcf_list <- fread(vcfs, header = F)
orf_list <- fread(orfs, header = F)
total_reads_list <- fread(total_reads, header = F)
#total_reads_num <- as.numeric(total_reads_dt$V1[1])
# ============================================================================
# Preprocess bed file
ref_bed_dt <- pre.process.bed(bed, bounds)

## assign names in loop below as go through each sample 
# base_name <- paste(strsplit(file_name, "_")[[1]][c(3,4)], collapse = "_")
# gene_name <- strsplit(file_name, "_")[[1]]
# gene_name <- gene_name[length(gene_name)]
# ============================================================================
# Combine consensus sequences if they end up having the same structure and mutation profile after filtering in the callconsensus module
total_reads_dfs <- lapply(total_reads_list$V1, function(depth) {
  reads <- fread(depth, header = F)
  reads[, sample := gsub("_total.*", "", depth)]
  reads
})
total_reads_dfs <- do.call('rbind', total_reads_dfs)

orf_dfs <- Map(pre.process.orf, orf_list$V1, vcf_list$V1, rep(list(ref_bed_dt), length(vcf_list$V1)))

vcf_structs <- lapply(vcf_list$V1, function(vcf) {
  vcf_struct_df <- pre.process.vcf.structure(vcf, ref_bed_dt)
  vcf_struct_df
} ) 

vcf_muts <- lapply(vcf_list$V1, function(vcf) {
  vcf_muts_df <- pre.process.vcf.mutations(vcf, ref_bed_dt)
  vcf_muts_df 
} )

# Calculate coverage across all aligned sites 
sample_covs <- lapply(vcf_list$V1, function(vcf) {
  cov_df <- sample.coverage.calc(vcf, ref_bed_dt)
  cov_df
})
combined_coverage_dt <- do.call('rbind', sample_covs)
sample_coverage_dt <- combined_coverage_dt[, .(read_depth = sum(total_reads)), by = c("feature", "start_pos", "end_pos", "sample")]

# Identify groups of identical data frames
identical_mut_groups <- identify.identical.dfs(vcf_muts, columns = c("feature", "start", "end", "symbol", "REF", "ALT"))
identical_struct_groups <- identify.identical.dfs(vcf_structs, columns = c("feature", "start", "end"))
# Find the intersection of data frames that are identical between both the structural and mutation data frames, outputting a list
common_values <- find.common.values(identical_struct_groups, identical_mut_groups)

# make an initial dataframe of sample/cluster names and corresponding genc ID and row index
base_names <- gsub(".vcf.gz", "", vcf_list$V1)
sample_names <- gsub("_cluster.*|_pp.*", "", vcf_list$V1)
cluster_id <- unlist(lapply(gsub("_modified.*", "", vcf_list$V1), function(clust) {
  splitting <- strsplit(clust, "_")[[1]]
  splitting[length(splitting)]
} ) )

id_dt <- data.table("sample" = c(sample_names, "wildtype"), "sample_cluster" = c(base_names, "wildtype"), "cluster_id" = c(cluster_id, "wildtype"), "genc_id" = c(paste0(gene_name, "_", seq(1, length(base_names))), "wildtype_0"), "struct_id" = seq(1, length(base_names)+1))
id_dt[, IDX := .I ]

modify_by_struct <- lapply(identical_struct_groups, function(struct) {
  id_dt[struct, struct_id := struct_id[1]]
})

if(length(common_values) > 0) {
  # modify the dataframe in place such that identical sample/clusters have identical ids
  updated_dt <- lapply(common_values, function(common) {
    id_dt[common, genc_id := genc_id[1]]
    id_dt
  })
  # look for within-sample duplicates that can be merged for plotting purposes 
  remove_idces <- c()
  same_within_sample <- lapply(unique(id_dt$genc_id), function(genc) {
    genc_dt <- id_dt[genc_id == genc]
    if(nrow(genc_dt) > 1) {
      same_sample_dt <- lapply(unique(genc_dt$sample), function(samp) {
        same_samp_dt <- genc_dt[sample == samp]
        if(nrow(same_samp_dt) > 1) {
          # subset the list of dataframes for just the unique ones within samples
          remove_idces <- append(remove_idces, same_samp_dt$IDX)
          # merge mutation frequency dataframes
          same_vcf_muts <- do.call('rbind', vcf_muts[c(same_samp_dt$IDX)])
          mean_vcf_muts <- same_vcf_muts[, .(value = mean(value)), by = .(feature, start, end, symbol, REF, ALT, sample)]
          setcolorder(mean_vcf_muts, c("feature", "start", "end", "value", "symbol", "REF", "ALT", "sample"))
          cluster_id <- same_vcf_muts$cluster[1]
          mean_vcf_muts[, cluster := cluster_id]
          vcf_muts <- append(vcf_muts, list(mean_vcf_muts))
          # merge vcf structure dataframes
          same_vcf_structs <- do.call('rbind', vcf_structs[c(same_samp_dt$IDX)])
          combined_depth <- same_vcf_structs[, .(maxDepth = sum(maxDepth)), by = .(feature, start, end, sample)]
          setcolorder(combined_depth, c("feature", "start", "end", "maxDepth", "sample"))
          combined_depth[, cluster := cluster_id]
          vcf_structs <- append(vcf_structs, list(combined_depth))
          # merge orf structure dataframes
          same_orf_structs <- do.call('rbind', orf_dfs[c(same_samp_dt$IDX)])
          combined_strand <- same_orf_structs[, .(strand = unique(strand)), by = .(feature, start, end, sample)]
          setcolorder(combined_strand, c("feature", "start", "end", "strand"))
          combined_strand[, cluster := cluster_id]
          orf_dfs <- append(orf_dfs, list(combined_strand))
        }
    } ) 
  }
 } )

  if(length(remove_idces) > 0) {
    vcf_muts <- vcf_muts[-remove_idces]
    vcf_structs <- vcf_structs[-remove_idces]
    orf_dfs <- orf_dfs[-remove_idces]
  }

}  

vcf_muts <- do.call('rbind', vcf_muts)
vcf_structs <- do.call('rbind', vcf_structs)
ref_struct <- data.table(feature = ref_bed_dt$feature, start = ref_bed_dt$start, end = ref_bed_dt$end, maxDepth = 1, sample = "wildtype", cluster = "wildtype")
vcf_structs <- rbind(vcf_structs, ref_struct)
orf_dfs <- do.call('rbind', orf_dfs)
ref_orf_dt <- ref_bed_dt[!(grepl("UTR", feature))]
ref_orf <- data.table(feature = ref_orf_dt$feature, start = ref_orf_dt$start, end = ref_orf_dt$end, strand = "+", in_frame = TRUE, sample = "wildtype", cluster = "wildtype")
orf_dfs <- rbind(orf_dfs, ref_orf)
setnames(id_dt, old = "cluster_id", new = "cluster")
# ============================================================================
# Generate ggtranscript plot
wt_dt <- id_dt[.N]
id_dt2 <- rbind(wt_dt, id_dt)
id_dt2 <- id_dt2[1:(nrow(id_dt)-1)]

id_vcf_structs <- vcf_structs[id_dt2, on=.(sample, cluster)]
id_vcf_structs[, c("seqnames", "strand", "type", "gene_name", "transcript_name") := .(1, "+", ifelse(grepl("UTR", feature), "UTR", "CDS"), gene_name, genc_id)]
setcolorder(id_vcf_structs, neworder = c("seqnames", "start", "end", "strand", "type", "gene_name", "transcript_name"))
id_vcf_structs[, genc_id_order := as.numeric(gsub(".*_", "", transcript_name))]
id_vcf_structs$transcript_name <- factor(id_vcf_structs$transcript_name)
id_vcf_structs[, features := gsub("3'|5'", "", feature)]
dup_samples <- id_vcf_structs[, unique(sample_cluster), by = genc_id]
dup_rows <- dup_samples[, which(duplicated(genc_id))]
dup_sample_clusters <- dup_samples$V1[dup_rows]
id_vcf_structs <- id_vcf_structs[!sample_cluster %in% c(dup_sample_clusters)]
cds <- id_vcf_structs[type == "CDS"]


id_orf_dfs <- orf_dfs[id_dt2, on = .(sample, cluster)]
id_orf_dfs <- id_orf_dfs[, genc_id_order := as.numeric(gsub(".*_", "", genc_id))][order(genc_id_order)]
correct_order_orfs <- data.table(genc_id_order = unique(id_orf_dfs$genc_id_order), correct_order = match(unique(id_orf_dfs$genc_id_order),unique(id_vcf_structs$genc_id_order)))
id_orf_dfs <- id_orf_dfs[correct_order_orfs, on = "genc_id_order"]
id_orf_dfs[, transcript_name := factor(genc_id)]
id_orf_dfs[, features := gsub("3'|5'", "", feature)]
id_orf_dfs[, orf_frame := ifelse(in_frame == TRUE, "In frame ORF", "Out of frame\nORF")]

# from the id_orf_dfs, get the info needed for labelling samples as brackets 
unique_positions <- id_orf_dfs[,unique(correct_order), by = sample]
min_positions <- unique_positions[, min(V1), by = sample]
min_positions$features <- 1
#min_positions$sample[3] <- "E13-27\nCTL"
min_positions[, sample := sub("_", "\n", sample)]
min_positions[, sample := sub("_.*$", "", sample)]
max_positions <- unique_positions[, max(V1), by = sample]


id_vcf_muts <- vcf_muts[id_dt2, on = .(sample, cluster)]
id_vcf_muts <- na.omit(id_vcf_muts)
id_vcf_muts <- id_vcf_muts[, genc_id_order := as.numeric(gsub(".*_", "", genc_id))][order(genc_id_order)]
id_vcf_muts[, `Mutation type` := ifelse(symbol == 16, "SNV", "INDEL")]
correct_order_muts <- data.table(genc_id_order = unique(id_vcf_muts$genc_id_order), correct_order = match(unique(id_vcf_muts$genc_id_order),unique(id_vcf_structs$genc_id_order)))
id_vcf_muts <- id_vcf_muts[correct_order_muts, on = "genc_id_order"]
id_vcf_muts[, features := gsub("3'|5'", "", feature)]

color_pal <- createPalette(length(unique(cds$feature)),  c("#ff0000", "#00ff00", "#0000ff"))
feature_colors <- data.table(features = c("In frame ORF", "Out of frame\nORF", as.character(sort(unique(as.numeric(cds$feature)))), "UTR"), color = c("#008000","#808080", color_pal, "#FFFFFF")) 

struct_plot <- ggplot() + 
  geom_range(data = id_vcf_structs, aes(xstart = start, xend = end, y = reorder(transcript_name, genc_id_order)), fill = "white", height = 0.25) +
  geom_range(data = cds, aes(xstart = start, xend = end, y = reorder(transcript_name, genc_id_order), fill = features), alpha = 0.5) +
  geom_point(data = id_vcf_muts, aes(start, correct_order, colour = `Mutation type`), shape = 16, size = 0.5) +
  geom_rect(data = id_orf_dfs, aes(xmin = start, xmax = end, ymin = correct_order + 0.25, ymax = correct_order + 0.5, fill = orf_frame)) +
  geom_segment(data = min_positions, aes(x = I(-0.2), xend = I(-0.2), y = V1-0.25, yend = max_positions$V1 +0.25), linetype=1, linewidth=0.5) +
  geom_text(data=min_positions, aes(x=I(-0.24), y=(V1-0.25+max_positions$V1+0.25)/2,  label = sample), angle = 90, fontface = "plain", size = 3) +
  scale_fill_manual(values = feature_colors$color, limits = c(feature_colors$features)) + 
  scale_x_continuous(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  coord_cartesian(xlim = c(0, ceiling(ref_bed_dt$end[nrow(ref_bed_dt)])), clip="off") +
  ylab("") +
  xlab("position (bp)") +
  theme_bw() +
  theme(plot.margin = unit(c(0.75,0.75,0.25,1), "cm")) +
  theme(legend.key=element_rect(colour="black"),legend.background=element_blank()) + 
  theme(aspect.ratio = (0.017544 + 0.081579*length(unique(id_orf_dfs$genc_id)))) +
  guides(fill = guide_legend(override.aes = list(shape = NA, border = NA)), colour = guide_legend(override.aes = list(size = 2)))


ggsave(file = paste0(gene_name, "_transcript_plot.png"), struct_plot, width = 8, height = 9, units = "in", dpi = 350)

# , ylim = c(0,20), expand = FALSE # add to coord_cartesian
# ============================================================================

# Generate a bed file of wild-type and all amplicons' exon structures; may want to just keep individual gencDNAs by id, then can have separate plot showing in which samples they were detected in
# vcf_struct_gsds <- data.table(sample = sample, genc_id = genc, struct_id = genc_dt$struct_id[1], start = vcf_struct_df$start, end = vcf_struct_df$end, featureType = vcf_struct_df$feature)
# struct_df <- do.call('rbind', sample_loop)
# repeated_features <- struct_df[, rle(featureType), by = genc_id]
# repeat_id <- unlist(sapply(repeated_features$lengths, function(x) seq(1, x)))
# struct_df[, featureType := paste0(featureType, ".", repeat_id)]
# struct_same_dt <- struct_df[, lapply(.SD, function(x) paste(unique(x), collapse = ";")), by = struct_id]
# struct_same_df <- struct_same_dt %>% separate_longer_delim(c(start, end, featureType), delim = ";")
# 
# struct_bed <- struct_same_df[, c("genc_id", "start", "end", "featureType")]
# struct_bed <- setDT(struct_bed)[, featureType := gsub("\\..*", "", featureType)]
# struct_bed[, featureType := ifelse(grepl("UTR", featureType), "UTR", paste0("CDS", featureType))]
# struct_bed[, c("start", "end") := .(as.numeric(start), as.numeric(end))]
# struct_df_max_length <- struct_bed[, max(end)]
# 
# ref_bed <- pre.process.bed_wt(bed)
# ref_bed_max_length <- ref_bed$end[nrow(ref_bed)]
# ref_bed_df <- data.table(genc_id = paste0(gene_name, "_wt"), start = ref_bed$start, end = ref_bed$end, featureType = ref_bed$featureType)
# new_ref_max <- struct_df_max_length + round((ref_bed_max_length - struct_df_max_length)/4, 0)
# ref_bed_df$end[nrow(ref_bed_df)] <- new_ref_max
# all_df <- rbind(ref_bed_df, struct_bed)
# fwrite(all_df, file = paste0(gene_name, "_structure.bed"), sep = "\t", col.names = FALSE)

# ============================================================================
# Generate an upset plot of the number of samples that share a unique genc_id
# upset_lst <- list()
# 
# new_lst <- lapply(unique(struct_df$sample), function(samp) {
#   sample_list <- c(unique(struct_df$genc_id[struct_df$sample==samp]))
#   sample_list
# })
# names(new_lst) <- unique(struct_df$sample)
# 
# 
# m1 = make_comb_mat(new_lst)
# 
# names(list) <- unique(struct_df$sample)
# 
# upset_df <- data.frame(matrix(ncol=length(unique(struct_df$sample)), nrow=length(unique(struct_df$genc_id))))
# names(upset_df) <- unique(struct_df$sample)
# rownames(upset_df) <- unique(struct_df$genc_id)
# 
# 
# 
# UpSet(m1)
# ============================================================================
# Trouble-shooting