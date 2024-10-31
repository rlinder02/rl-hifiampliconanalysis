#!/usr/bin/env Rscript

##############################################################################

# Project:     Make circos plot of consensus amplicon sequences
# Author:      Robert Linder
# Date:        2024-10-24
# Title:       generate_circos_plots
# Description: Generates circos plots of each consensus amplicon species using the circlize package. 
# Version:     1.0.0

##############################################################################

# ============================================================================
# Parse command line inputs

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("Usage: generate_circos_plots.R <vcfs> <bed> <bounds> <total_reads> <file_name>", call.=FALSE)
}

vcfs <- args[1]
bed <- args[2]
bounds <- args[3]
total_reads <- args[4]
file_name <- args[5]

# ============================================================================
# For trouble-shooting locally

# vcf_list <- dir(pattern = "HU_PCR_SMARCA5_AMPLICON_HPCPS_CTL_cluster.*.vcf.gz")
# bed <- "hSmarca5_cDNA.bed"
# bounds <- "hSmarca5_cDNA.txt"
# total_reads <- "HU_PCR_SMARCA5_AMPLICON_HPCPS_CTL_total_aligned_reads.txt"
# file_name <- "HU_PCR_SMARCA5_AMPLICON_HPCPS_CTL"

# ============================================================================
# Load packages and sourced files
library(circlize)
library(data.table)
library(tidyverse)
library(ComplexHeatmap)
library(gridBase)

# ============================================================================
# Set global options

options(digits = 10)
projectDir <- getwd()
# setwd("/Users/rlinder/Library/CloudStorage/OneDrive-SanfordBurnhamPrebysMedicalDiscoveryInstitute/Chun_lab/Pipelines/HiFi_PCR_analysis/tests/circos_plot_sandbox")


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
  ref_bed_dt[, start := c(0, V5[-length(V5)]+1)]
  ref_bed_dt[, end := V5]
  ref_bed_dt[, feature := gsub("Exon_", "", V1)]
  columns <- c("feature", "start", "end")
  ref_bed_dt <- ref_bed_dt[, ..columns]
  # limit plotting to primer bounds
  if(ref_bounds_dt$V1[1] < ref_bed_dt$end[1]) {
    start_row <- 1
  } else {
    start_row <- max(which(ref_bed_dt$end < ref_bounds_dt$V1[1]))+1
  }
  ref_bed_dt$start[start_row] <- ref_bounds_dt$V1[1]
  if(ref_bounds_dt$V1[2] > ref_bed_dt$start[nrow(ref_bed_dt)]) {
    end_row <- nrow(ref_bed_dt)
  } else {
    end_row <- min(which(ref_bed_dt$start > ref_bounds_dt$V1[2]))-1
  }
  ref_bed_dt$end[end_row] <- ref_bounds_dt$V1[2]
  ref_bed_dt <- ref_bed_dt[start_row:end_row]
  ref_bed_dt
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
  ref_bed_dt[, feature := gsub("Exon_", "", V1)]
  columns <- c("feature", "start", "end")
  ref_bed_dt <- ref_bed_dt[, ..columns]
  ref_bed_dt[, featureType := "exon"]
  ref_bed_dt
}

pre.process.vcf.structure <- function(vcf_file, ref_bed_dt) {
  vcf_dt <- fread(vcf_file)
  vcf_max_depth <- max(as.numeric(str_match(vcf_dt$INFO, 'DP=(\\d+)')[,2]), na.rm = TRUE)
  names(vcf_dt)[10] <- "SAMPLE"
  # subtract one from the vcf file coordinates so is in bed coordinate space (0-based)
  vcf_dt[, POS := POS - 1 ]
  vcf_dt[ref_bed_dt, on=.(POS >= start, POS <= end), feature := i.feature]
  vcf_dt <- vcf_dt[!is.na(feature)]
  vcf_dt[, c("start", "end") := .(min(POS), max(POS)), by = feature]
  #print(vcf_dt)
  struct_columns <- c("feature", "start", "end")
  vcf_dt_structure <- unique(vcf_dt[, ..struct_columns])
  vcf_dt_structure[, maxDepth := vcf_max_depth]
  vcf_dt_structure
}

pre.process.vcf.mutations <- function(vcf_file, ref_bed_dt) {
  vcf_dt <- fread(vcf_file)
  names(vcf_dt)[10] <- "SAMPLE"
  # subtract one from the vcf file coordinates so is in bed coordinate space (0-based)
  vcf_dt[, POS := POS - 1 ]
  vcf_dt[ref_bed_dt, on=.(POS >= start, POS <= end), feature := i.feature]
  vcf_dt <- vcf_dt[!is.na(feature)]
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
  vcf_dt_muts_dt
}

identify_identical_dfs <- function(df_list) {
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

find_common_values <- function(lst1, lst2) {
  common_values <- list()
  # Iterate over each element in the first list
  for (element1 in lst1) {
    # Iterate over each element in the second list
    for (element2 in lst2) {
      find_intersection <- intersect(element1, element2)
      flush.console()
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

# ============================================================================
# Load data

vcf_list <- fread(vcfs, header = F)
total_reads_dt <- fread(total_reads)
total_reads_num <- as.numeric(total_reads_dt$V1[1])
# ============================================================================
# Preprocess bed file
ref_bed_dt <- pre.process.bed(bed, bounds)
base_name <- paste(strsplit(file_name, "_")[[1]][c(3,4)], collapse = "_")
gene_name <- strsplit(file_name, "_")[[1]][3]
# ============================================================================
# Combine consensus sequences if they end up having the same structure and mutation profile after filtering in the callconsensus module
vcf_muts <- lapply(vcf_list$V1, function(vcf) {
  vcf_muts_df <- pre.process.vcf.mutations(vcf, ref_bed_dt)
  muts_columns <- c("feature", "start", "end", "symbol", "REF", "ALT")
  vcf_muts_df <- vcf_muts_df[, ..muts_columns]
  vcf_muts_df 
} )

vcf_muts_values <- lapply(vcf_list$V1, function(vcf) {
  vcf_muts_df <- pre.process.vcf.mutations(vcf, ref_bed_dt)
  muts_columns <- c("feature", "start", "end", "value", "symbol", "REF", "ALT")
  vcf_muts_df <- vcf_muts_df[, ..muts_columns]
  vcf_muts_df 
} )

vcf_structs <- lapply(vcf_list$V1, function(vcf) {
  vcf_struct_df <- pre.process.vcf.structure(vcf, ref_bed_dt)
  struct_columns <- c("feature", "start", "end")
  vcf_struct_df <- vcf_struct_df[, ..struct_columns]
  vcf_struct_df
} )

vcf_structs_depth <- lapply(vcf_list$V1, function(vcf) {
  vcf_struct_df <- pre.process.vcf.structure(vcf, ref_bed_dt)
  struct_columns <- c("feature", "start", "end", "maxDepth")
  vcf_struct_df <- vcf_struct_df[, ..struct_columns]
  vcf_struct_df
} ) 
# Identify groups of identical data frames
identical_mut_groups <- identify_identical_dfs(vcf_muts)
identical_struct_groups <- identify_identical_dfs(vcf_structs)
# Find the intersection of data frames that are identical between both the structural and mutation data frames, outputting a list
common_values <- find_common_values(identical_struct_groups, identical_mut_groups)
# Find the unique data frames that are not duplicated
if(length(common_values) > 0) {
  unique_mut_dfs <- vcf_muts_values[-c(unique(unlist(common_values)))]
  unique_struct_dfs <- vcf_structs_depth[-c(unique(unlist(common_values)))]
  # Merge identical mutation data frames such that the value column is the average of all identical value columns
  new_mut_dfs <- lapply(common_values, function(idx) {
    all_values <- do.call('rbind', vcf_muts_values[c(idx)])
    mean_values <- all_values[, .(value = mean(value)), by = .(feature, start, end, symbol, REF, ALT)]
    setcolorder(mean_values, c("feature", "start", "end", "value", "symbol", "REF", "ALT"))
    mean_values
  })
  new_struct_dfs <- lapply(common_values, function(idx) {
    all_values <- do.call('rbind', vcf_structs_depth[c(idx)])
    combined_depth <- all_values[, .(maxDepth = sum(maxDepth)), by = .(feature, start, end)]
    setcolorder(combined_depth, c("feature", "start", "end", "maxDepth"))
    combined_depth
  })
  # if there are identical data frames, then overwrite the list of data frames stored in vcf_muts_values and vcf_structs_depth
  vcf_muts_values <- append(new_mut_dfs, unique_mut_dfs)
  vcf_structs_depth <- append(new_struct_dfs, unique_struct_dfs)
}

# ============================================================================
# Generate Circos plot

col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))

lgd_muts = Legend(at = c("SNV", "INDEL"), type = "points", pch = c(16,17), title_position = "topleft", title = "Mutation type")
lgd_reads = Legend(col_fun = col_fun, title_position = "topleft", title = "Fraction of reads")
lgd_list_vertical = packLegend(lgd_muts, lgd_reads)

fileName <- paste0(file_name, "_circos_plot.png")
png(fileName, height = 12, width = 8, units = "in", res = 1200)
circos.par("track.height" = 0.05, circle.margin = c(0.1, 0.1, 0.1, 0.1), "start.degree" = 90)
circos.genomicInitialize(ref_bed_dt, plotType = NULL)
# outermost track of wild-type exon structure
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = "gray")
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.075, bg.border = NA)
counter <- 1
cluster_counter <- 0
struct_dfs <- lapply(1:length(vcf_structs_depth), function(idx) {
  vcf_struct_df <- vcf_structs_depth[[idx]]
  vcf_muts_df <- vcf_muts_values[[idx]]
  vcf_max_depth <- vcf_struct_df$maxDepth[1]
  # remove the depth column from vcf_struct_df
  struct_cols <- c("feature", "start", "end")
  vcf_struct_df <- vcf_struct_df[, ..struct_cols]
  vcf_track_col <- vcf_max_depth/total_reads_num
  counter <<- counter + 1
  print(counter)
  flush.console()
  cluster_counter <<- cluster_counter + 1
  circos.genomicTrack(vcf_struct_df, ylim = c(0, 1), track.height = 0.05, bg.border = NA, panel.fun = function(region, value, ...) {
                        i = getI(...)
                        xlim = CELL_META$xlim
                        circos.rect(region$start, 0, region$end, 1, col = add.alpha(col_fun(vcf_track_col), 0.5), border = "black", track.index = counter)
  })
  circos.genomicTrack(vcf_muts_df, numeric.column = 4, ylim = c(0, 1), track.height = 0.05, bg.border = NA, panel.fun = function(region, value, ...) {
                        i = getI(...)
                        xlim = CELL_META$xlim
                        circos.genomicPoints(region, value, pch = value$symbol, cex = 0.7, col = "black", track.index = counter, ...)
  })
  vcf_struct_gsds <- data.table(gene_id = paste0(base_name, "_", cluster_counter), start = vcf_struct_df$start, end = vcf_struct_df$end, featureType = "exon")
  return(vcf_struct_gsds)
})
circos.clear()
draw(lgd_list_vertical, x = unit(0.03, "npc"), y = unit(0.75, "npc"), just = c("left", "top"))
dev.off()

# ============================================================================
# Generate a bed file of wild-type and all amplicons' exon structures

ref_bed <- pre.process.bed_wt(bed)
ref_bed_max_length <- ref_bed$end[nrow(ref_bed)]
ref_bed_df <- data.table(gene_id = paste0(gene_name, "_wt_mRNA"), start = ref_bed$start, end = ref_bed$end, featureType = "exon")
struct_df <- do.call('rbind', struct_dfs)
struct_df_max_length <- struct_df[, max(end)]
new_ref_max <- struct_df_max_length + (ref_bed_max_length - struct_df_max_length)/4
ref_bed_df$end[nrow(ref_bed_df)] <- new_ref_max
all_df <- rbind(ref_bed_df, struct_df)
fwrite(all_df, file = paste0(file_name, "_structure.bed"), sep = "\t", col.names = FALSE)
# ============================================================================
# Trouble-shooting