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

vcfs <- "vcf_fofn.txt"
bed <- "hTARDBP_cDNA_full.bed"
bounds <- "hTARDBP_cDNA.txt"
total_reads <- "total_reads_fofn.txt"
gene_name <- "TARDBP"
orfs <- "orf_fofn.txt"


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

setwd("/Users/rlinder/Library/CloudStorage/OneDrive-SanfordBurnhamPrebysMedicalDiscoveryInstitute/Chun_lab/Projects/gencDNA/PCR_Southerns/Human/TARDBP/2024-12-12_run")


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
  vcf_dt <- vcf_dt[depth > 0]
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
  if(ref_bed_dt$end[1] - orf_dt_cds_longest$V2 == 0) {
    orf_dt_cds_longest[, c("V2", "V3") := .(V2 + 1, V3 + 1)]
  }
  expanded_dt <- data.table(POS = orf_dt_cds_longest$V2:orf_dt_cds_longest$V3, strand = orf_dt_cds_longest$V6)
  expanded_dt[vcf_dt_structure, on=.(POS >= start, POS <= end), c("feature", "runs") := .(i.feature, i.runs)]
  expanded_dt <- expanded_dt[!is.na(feature)]
  expanded_dt[, c("start", "end") := .(min(POS), max(POS)), by = c("runs", "feature")]
  struct_columns <- c("feature", "start", "end", "strand")
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
total_reads_dfs <- lapply(sort(total_reads_list$V1), function(depth) {
  reads <- fread(depth, header = F)
  reads[, sample := gsub("_total.*", "", depth)]
  reads
})
total_reads_dfs <- do.call('rbind', total_reads_dfs)

orf_dfs <- Map(pre.process.orf, sort(orf_list$V1), sort(vcf_list$V1), rep(list(ref_bed_dt), length(vcf_list$V1)))

vcf_structs <- lapply(sort(vcf_list$V1), function(vcf) {
  vcf_struct_df <- pre.process.vcf.structure(vcf, ref_bed_dt)
  vcf_struct_df
} ) 

vcf_muts <- lapply(sort(vcf_list$V1), function(vcf) {
  vcf_muts_df <- pre.process.vcf.mutations(vcf, ref_bed_dt)
  vcf_muts_df 
} )

# Calculate coverage across all aligned sites 
sample_covs <- lapply(sort(vcf_list$V1), function(vcf) {
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

# make an initial dataframe of sample/cluster names and corresponding genc ID ann row index
base_names <- gsub(".vcf.gz", "", sort(vcf_list$V1))
sample_names <- gsub("_cluster.*|_pp.*", "", sort(vcf_list$V1))
cluster_id <- unlist(lapply(gsub("_modified.*", "", sort(vcf_list$V1)), function(clust) {
  splitting <- strsplit(clust, "_")[[1]]
  splitting[length(splitting)]
} ) )
id_dt <- data.table("sample" = sample_names, "sample_cluster" = base_names, "cluster_id" = cluster_id, "genc_id" = paste0(gene_name, "_", seq(1,length(base_names))), "struct_id" = seq(1, length(base_names)))
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
orf_dfs <- do.call('rbind', orf_dfs)
# ============================================================================
# Generate Circos plot

num_sectors <- ref_bed_dt[, uniqueN(feature)]
first_sector <- ref_bed_dt$feature[1]
col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))

lgd_muts = Legend(at = c("SNV", "INDEL"), type = "points", pch = c(16,17), title_position = "topleft", title = "Mutation\ntype", labels_gp = gpar(fontsize = 8), title_gp = gpar(fontsize = 8, fontface = "bold"))
lgd_reads = Legend(col_fun = col_fun, title_position = "topleft", title = "Read\nfraction", labels_gp = gpar(fontsize = 8), title_gp = gpar(fontsize = 8, fontface = "bold"))
lgd_orfs = Legend(labels = "ORF", type = 'points', pch = 26, legend_gp = gpar(col = "slategrey", alpha = 0.5), title_position = "topleft", title = "Longest\nORF", labels_gp = gpar(fontsize = 8), title_gp = gpar(fontsize = 8, fontface = "bold"))

#lgd_orfs = Legend(at = c("ORF"), type = "lines", legend_gp = gpar(col = "purple", lwd = 1, ), title_position = "topleft", title = "Longest ORF")
lgd_list_vertical = packLegend(lgd_orfs, lgd_muts, lgd_reads)

# loop through all samples to create one circos plot per sample that will be combined into a single plot 
# layout(matrix(c(rep(1,15), rep(2, 11), rep(4,4), rep(3, 15)), ncol = 15, byrow = TRUE))
# layout(matrix(1:9, 3, 3))

num_rows <- max(c(floor(length(unique(id_dt$sample))/2), 1))
num_cols <- ceiling(length(unique(id_dt$sample))/num_rows)
png(paste0(gene_name, "_circos_plot_v3.png"), height = min(c(12, num_rows*4)) , width = min(c(12, num_cols*4)), units = "in", res = 1200)
#pdf(paste0(gene_name, "_circos_plot.pdf"), height = min(c(12, num_rows*4)) , width = min(c(12, num_cols*4)))

par(mfrow = c(num_rows, num_cols))

plot_counter <- 0
sample_loop <- lapply(unique(id_dt$sample), function(samp) {
  print(samp)
  flush.console()
  coverage_dt <- sample_coverage_dt[sample == samp]
  vcf_struct_dt <- vcf_structs[sample == samp]
  genc_id_dt <- id_dt[sample == samp]
  vcf_mut_dt <- vcf_muts[sample == samp]
  orf_dt <- orf_dfs[sample == samp]
  total_reads_dt <- total_reads_dfs[sample == samp]
  total_reads_num <- total_reads_dt$V1
  plot_counter <<- plot_counter + 1
  #fileName <- paste0(samp, "_circos_plot.png")
  #png(fileName, height = 12, width = 8, units = "in", res = 1200)
  #par(mar = c(0, 0, 0, 0))
  if(plot_counter %in% seq(1, num_rows*num_cols, num_cols)) {
    par(mar = c(0, 2, 0, 0))
  } else {
    par(mar = c(0, 2, 0, 0))
  }
  circos.par("track.height" = 0.05, track.margin = c(0.001, 0.001), circle.margin = c(0.1, 0.1, 0.1, 0.1), "start.degree" = 90, gap.after = c(rep(1, num_sectors-1), 20), points.overflow.warning = FALSE, cell.padding = c(0.01, 1, 0.01, 1))
  circos.genomicInitialize(ref_bed_dt, plotType = NULL)
  # outermost track of wild-type exon structure
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    #circos.rect(xlim[1], 0, xlim[2], 1, col = "gray")
    circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "black",
                facing = "inside", niceFacing = TRUE)
  }, track.height = 0.05, bg.border = NA)
  
  circos.genomicTrack(coverage_dt, numeric.column = 5, ylim = range(coverage_dt$read_depth, na.rm = TRUE), track.height = 0.075, panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value, col = "blue", lwd = 1, ...)
    if (CELL_META$sector.index == first_sector & CELL_META$track.index == 2) {
      circos.yaxis(at = range(coverage_dt$read_depth), labels = format(range(coverage_dt$read_depth), scientific = TRUE, digits = 2), labels.cex = 0.5)
    }
  } )
  counter <- 2
  genc_dfs <- lapply(unique(genc_id_dt$genc_id), function(genc) {
    genc_dt <- genc_id_dt[genc_id == genc]
    vcf_struct_df <- vcf_struct_dt[cluster %in% unique(genc_dt$cluster_id)]
    vcf_muts_df <- vcf_mut_dt[cluster %in% unique(genc_dt$cluster_id)]
    orf_df <- orf_dt[cluster %in% unique(genc_dt$cluster_id)]
    vcf_max_depth <- vcf_struct_df$maxDepth[1]
    # remove the depth column from vcf_struct_df
    struct_cols <- c("feature", "start", "end")
    vcf_struct_df <- vcf_struct_df[, ..struct_cols]
    orf_struct_df <- orf_df[, ..struct_cols]
    vcf_track_col <- vcf_max_depth/total_reads_num
    counter <<- counter + 1
    print(counter)
    flush.console()
    circos.genomicTrack(vcf_struct_df, ylim = c(0, 1), track.height = 0.05, bg.border = NA, panel.fun = function(region, value, ...) {
      i = getI(...)
      xlim = CELL_META$xlim
      circos.rect(region$start, 0, region$end, 1, col = add.alpha(col_fun(vcf_track_col), 0.4), border = "black", track.index = counter)
      if (CELL_META$sector.index == first_sector) {
        circos.yaxis(at = 0.5, labels = genc, labels.cex = 0.5, track.index = counter)
      }
    })
    circos.genomicTrack(vcf_muts_df, numeric.column = 4, ylim = c(0, 1), track.height = 0.05, bg.border = NA, panel.fun = function(region, value, ...) {
      i = getI(...)
      xlim = CELL_META$xlim
      circos.genomicPoints(region, value, pch = value$symbol, cex = 0.7, col = "darkred", track.index = counter, ...)
    })
    # adding predicted longest orf
    circos.genomicTrack(orf_struct_df, ylim = c(0, 1), track.height = 0.05, bg.border = NA, panel.fun = function(region, value, ...) {
      i = getI(...)
      xlim = CELL_META$xlim
      circos.rect(region$start, 0.25, region$end, 0.75, col = add.alpha("slategrey", 0.5), border = NA, lwd = 1.5, density = 45, angle = 45, track.index = counter)
    })
    sample <- strsplit(samp, "_")[[1]]
    sample <- paste(sample[1:(length(sample)-1)], collapse = "_")
    vcf_struct_gsds <- data.table(sample = sample, genc_id = genc, struct_id = genc_dt$struct_id[1], start = vcf_struct_df$start, end = vcf_struct_df$end, featureType = vcf_struct_df$feature)
    return(vcf_struct_gsds)
  } )
  genc_dfs_combined <- do.call('rbind', genc_dfs)
  title(main = samp, cex.main = 1, line = -2, xpd = NA)
  circos.clear()
  
  #draw(lgd_list_vertical, x = unit(0.03, "npc"), y = unit(0.75, "npc"), just = c("left", "top"))
  #dev.off()
  
  return(genc_dfs_combined)
} )
# plot(1, type = "n", axes=FALSE, xlab="", ylab="", xlim = c(0,1), ylim = c(0,1))
# legend(x = 0.5, y = 1, legend = c("SNV", "INDEL"), pch = c(16,17), cex=1, horiz = TRUE, title="Mutation type", bty="n")
# legend(x = 0.5, y = 0.85, legend = c("ORF"), pch = "/", col = "slategrey", cex=1, horiz = TRUE, title="Longest ORF", bty="n")

#lgd_orfs = Legend(labels = "ORF", type = 'points', pch = 26, legend_gp = gpar(col = "slategrey", alpha = 0.5), title_position = "topleft", title = "Longest ORF")

draw(lgd_list_vertical, x = unit(0.01, "npc"), y = unit(0.99, "npc"), just = c("left", "top")) # need to expand canvas to the left
dev.off()

# ============================================================================
# Generate a bed file of wild-type and all amplicons' exon structures; may want to just keep individual gencDNAs by id, then can have separate plot showing in which samples they were detected in
struct_df <- do.call('rbind', sample_loop)
repeated_features <- struct_df[, rle(featureType), by = genc_id]
repeat_id <- unlist(sapply(repeated_features$lengths, function(x) seq(1, x)))
struct_df[, featureType := paste0(featureType, ".", repeat_id)]
struct_same_dt <- struct_df[, lapply(.SD, function(x) paste(unique(x), collapse = ";")), by = struct_id]
struct_same_df <- struct_same_dt %>% separate_longer_delim(c(start, end, featureType), delim = ";")

struct_bed <- struct_same_df[, c("genc_id", "start", "end", "featureType")]
struct_bed <- setDT(struct_bed)[, featureType := gsub("\\..*", "", featureType)]
struct_bed[, featureType := ifelse(grepl("UTR", featureType), "UTR", paste0("CDS", featureType))]
struct_bed[, c("start", "end") := .(as.numeric(start), as.numeric(end))]
struct_df_max_length <- struct_bed[, max(end)]

ref_bed <- pre.process.bed_wt(bed)
ref_bed_max_length <- ref_bed$end[nrow(ref_bed)]
ref_bed_df <- data.table(genc_id = paste0(gene_name, "_wt"), start = ref_bed$start, end = ref_bed$end, featureType = ref_bed$featureType)
new_ref_max <- struct_df_max_length + round((ref_bed_max_length - struct_df_max_length)/4, 0)
ref_bed_df$end[nrow(ref_bed_df)] <- new_ref_max
all_df <- rbind(ref_bed_df, struct_bed)
fwrite(all_df, file = paste0(gene_name, "_structure.bed"), sep = "\t", col.names = FALSE)

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