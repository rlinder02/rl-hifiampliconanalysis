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

pre.process.vcf.structure <- function(vcf_file) {
  vcf_dt <- fread(vcf_file)
  names(vcf_dt)[10] <- "SAMPLE"
  # subtract one from the vcf file coordinates so is in bed coordinate space (0-based)
  vcf_dt[, POS := POS - 1 ]
  vcf_dt[ref_bed_dt, on=.(POS >= start, POS <= end), feature := i.feature]
  vcf_dt[, c("start", "end") := .(min(POS), max(POS)), by = feature]
  struct_columns <- c("feature", "start", "end")
  vcf_dt_structure <- unique(vcf_dt[, ..struct_columns])
  vcf_dt_structure
}

pre.process.vcf.mutations <- function(vcf_file, ref_bed_dt) {
  vcf_dt <- fread(vcf_file)
  names(vcf_dt)[10] <- "SAMPLE"
  # subtract one from the vcf file coordinates so is in bed coordinate space (0-based)
  vcf_dt[, POS := POS - 1 ]
  vcf_dt[ref_bed_dt, on=.(POS >= start, POS <= end), feature := i.feature]
  vcf_dt[, c("start", "end") := .(min(POS), max(POS)), by = feature]
  vcf_dt[, total_reads := as.numeric(str_match(INFO, 'DP=(\\d+)')[,2])]
  # keep only positions with called mutations 
  vcf_dt_muts <- vcf_dt[grepl("AC=", INFO) & grepl("PASS", FILTER)]
  vcf_dt_muts[, TYPE := ifelse(grepl("INDEL", INFO), "INDEL", "SNV")]
  vcf_dt_muts[, symbol := ifelse(grepl("SNV", TYPE), 16, 17)]
  vcf_dt_muts[, alt_count := as.numeric(gsub(".*,", "", SAMPLE))]
  vcf_dt_muts[, value := alt_count/total_reads]
  vcf_dt_muts[, POS_END := POS +1]
  muts_columns <- c("feature", "POS", "POS_END", "value", "symbol")
  vcf_dt_muts_dt <- vcf_dt_muts[, ..muts_columns]
  setnames(vcf_dt_muts_dt, old = c("POS", "POS_END"), new = c("start", "end"))
  vcf_dt_muts_dt
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

# ============================================================================
# Load data

vcf_list <- fread(vcfs, header = F)
total_reads_dt <- fread(total_reads)
total_reads_num <- as.numeric(total_reads_dt$V1[1])
# ============================================================================
# Preprocess bed file
ref_bed_dt <- pre.process.bed(bed, bounds)

# ============================================================================
# Generate Circos plot

fileName <- paste0(file_name, "_circos_plot.png")
png(fileName, height = 12, width = 8, units = "in", res = 1200)
lgd = Legend(at = c("SNV", "INDEL"), type = "points", pch = c(16,17), title_position = "topleft")
circos.par("track.height" = 0.05, circle.margin = c(0.1, 0.1, 0.1, 0.1), "start.degree" = 90, cell.padding = c(0, 0))
circos.genomicInitialize(ref_bed_dt, plotType = NULL)
# outermost track of wild-type exon structure
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = "gray")
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.05, bg.border = NA)
counter <- 1
lapply(vcf_list$V1[c(1:5)], function(vcf) {
  vcf_struct_df <- pre.process.vcf.structure(vcf)
  vcf_muts_df <- pre.process.vcf.mutations(vcf, ref_bed_dt)
  vcf_max_depth <- vcf.read.depth(vcf)
  vcf_track_height <- vcf_max_depth/total_reads_num
  rescaled_track_height <- rescale(vcf_track_height, 0, 0.02, 1, 0.1)
  print(vcf_track_height)
  print(rescaled_track_height)
  counter <<- counter + 1
  print(counter)
  circos.genomicTrack(vcf_struct_df, ylim = c(0.02, .1), track.height = rescaled_track_height, bg.border = NA, panel.fun = function(region, value, ...) {
                        i = getI(...)
                        xlim = CELL_META$xlim
                        circos.rect(region$start, 0, region$end, 1, col = "white", border = "black", track.index = counter)
  })
  circos.genomicTrack(vcf_muts_df, numeric.column = 4, ylim = c(0.02, 0.1), track.height = rescaled_track_height, bg.border = NA, panel.fun = function(region, value, ...) {
                        i = getI(...)
                        xlim = CELL_META$xlim
                        circos.genomicPoints(region, value, pch = value$symbol, cex = 0.5, col = "red", track.index = counter, ...)
  })
})
circos.clear()
draw(lgd, x = unit(0.03, "npc"), y = unit(0.75, "npc"), just = c("left", "top"))
dev.off()
# ============================================================================
# Trouble-shooting