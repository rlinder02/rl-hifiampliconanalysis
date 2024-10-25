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
if (length(args) < 2) {
  stop("Usage: generate_circos_plots.R <fasta> <vcf>", call.=FALSE)
}

fasta <-  args[1]
vcf <- as.numeric(args[2])

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
#projectDir <- getwd()
setwd("/Users/rlinder/Library/CloudStorage/OneDrive-SanfordBurnhamPrebysMedicalDiscoveryInstitute/Chun_lab/Pipelines/HiFi_PCR_analysis/tests/circos_plot_sandbox")


# ============================================================================
# Custom functions

# ============================================================================
# Load data

ref_fasta <- "hSmarca5_cDNA.fasta"
ref_bed <- "hSmarca5_cDNA.bed"
ref_bounds <- "hSmarca5_cDNA.txt"
vcf1 <- "HU_PCR_SMARCA5_AMPLICON_HPCPS_CTL_cluster7_modified.vcf.gz"
fasta1 <- "HU_PCR_SMARCA5_AMPLICON_HPCPS_CTL_cluster7.fasta"
vcf2 <- "HU_PCR_SMARCA5_AMPLICON_HPCPS_CTL_cluster39_modified.vcf.gz"
fasta2 <- "HU_PCR_SMARCA5_AMPLICON_HPCPS_CTL_cluster39.fasta"

ref_bounds_dt <- fread(ref_bounds)

ref_bed_dt <- fread(ref_bed)
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

vcf_1 <- fread(vcf1)
# subtract one from the vcf file coordinates so is in bed coordinate space (0-based)
vcf_1[, POS := POS - 1 ]
vcf_1[ref_bed_dt, on=.(POS >= start, POS <= end), feature := i.feature]
vcf_1[, c("start", "end") := .(min(POS), max(POS)), by = feature]
struct_columns <- c("feature", "start", "end")
vcf_1_structure <- unique(vcf_1[, ..struct_columns])

# limit plotting to primer bounds 
ref_bed_dt$start[1] <- ref_bounds_dt$V1[1]
ref_bed_dt$end[nrow(ref_bed_dt)] <- ref_bounds_dt$V1[2]

# ============================================================================
# Preprocess data

# ============================================================================
# Generate Circos plot

circos.par("track.height" = 0.1, circle.margin = c(0.1, 0.1, 0.1, 0.1), "start.degree" = 90)
circos.genomicInitialize(ref_bed_dt, plotType = NULL)
# outermost track of wild-type exon structure
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = "gray")
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.1, bg.border = NA)
# track for next species 
circos.genomicTrack(vcf_1_structure, stack = TRUE, track.height = 0.05, bg.border = NA,
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      xlim = CELL_META$xlim
                      ylim = CELL_META$ylim
                      print(region$start)
                      print(region$end)
                      flush.console()
                      circos.rect(region$start, 0, region$end, 1, col = "red")
})
circos.clear()


circos.Track(vcf_1_structure, stack = TRUE, track.height = 0.05, bg.border = NA,
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      xlim = CELL_META$xlim
                      ylim = CELL_META$ylim
                      circos.rect(xlim[1], 0, xlim[2], 1, col = "gray")
                      circos.segments(region$LF_start, y0, region$LF_start, y1, straight=TRUE, col = "red", lwd = 1, lty = 1)
                      ifelse(value$Gene == "Smarca5", cex_val <- 1, cex_val <- 0.75)
                      ifelse(value$Gene == "Smarca5", y_adj <- rep(2.6, nrow(value)), y_adj <- rep(2.3, nrow(value)))
                      by_gene <- lapply(value$Gene, function(gene) {
                        ifelse(gene %in% unique_hit_genes$Gene, font_type <- 2, font_type <- 1)
                        print(value$Gene)
                        print(font_type)
                        circos.genomicText(region, value$Gene, y = y_adj, labels = value$Gene, track.index = 1, facing = "clockwise", niceFacing = TRUE, cex = cex_val, font = font_type)
                      })
                    })
circos.clear()






lgd_points = Legend(at = c("MOPC", "J", "Q", "TSM", "Z"), type = "lines", legend_gp = gpar(col = mycols[1:5], lwd = 3), title_position = "topleft", title = "Cell lines")
lgd_shapes = Legend(at = c("TSD present", "No TSD"), type = "points", pch = c(16,17), title_position = "topleft")
lgd_list_vertical = packLegend(lgd_points, lgd_shapes)
#ins_df_dcast <- dcast(ins_df_clean, LF_chr + LF_start + LF_end ~ Cell_line, value.var = "value")

plot_circle <- circle_plot(ins_dfs, "mm10", lgd_list_vertical)
genome = "mm10"

circle_plot <- function(ins_dfs, genome, legend) {
  fileName <- paste0(getwd(), "/Plots/Genomewide_putatitve_gencDNA/putatitve_mouse_cell_line_gencDNA.png")
  png(fileName, height = 12, width = 8, units = "in", res = 1200)
  circos.par("track.height" = 0.1, circle.margin = c(0.1, 0.1, 0.1, 0.1))
  circos.initializeWithIdeogram(species = genome, plotType = NULL)
  
  
  
  col_counter <- 0
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    col_counter <<- col_counter + 1
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.rect(xlim[1], 0, xlim[2], 1, col = "gray")
    circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
                facing = "inside", niceFacing = TRUE)
  }, track.height = 0.1, bg.border = NA)
  y0 <- -10.5
  y1 <- 3
  counter <- 0
  lapply(ins_dfs, function(df) {
    counter <<- counter + 1
    df <- df[, c("LF_chr", "LF_start", "LF_end", "value", "Gene", "symbol")]
    #df <- rbind(df, empty_df)
    circos.genomicTrack(df, numeric.column = 4, stack = TRUE, track.height = 0.05, bg.border = NA,
                        panel.fun = function(region, value, ...) {
                          i = getI(...)
                          xlim = CELL_META$xlim
                          ylim = CELL_META$ylim
                          circos.genomicPoints(region, value, pch = value$symbol, cex = 0.75, col = "black", ...)
                          circos.segments(region$LF_start, y0, region$LF_start, y1, straight=TRUE, col = "red", lwd = 1, lty = 1)
                          ifelse(value$Gene == "Smarca5", cex_val <- 1, cex_val <- 0.75)
                          ifelse(value$Gene == "Smarca5", y_adj <- rep(2.6, nrow(value)), y_adj <- rep(2.3, nrow(value)))
                          by_gene <- lapply(value$Gene, function(gene) {
                            ifelse(gene %in% unique_hit_genes$Gene, font_type <- 2, font_type <- 1)
                            print(value$Gene)
                            print(font_type)
                            circos.genomicText(region, value$Gene, y = y_adj, labels = value$Gene, track.index = 1, facing = "clockwise", niceFacing = TRUE, cex = cex_val, font = font_type)
                          })
                        })
    y0 <<- y0 + 2.5
    y1 <<- y1 + 2.5
  } )
  for(sect in get.all.sector.index()) {
    for(sn in get.all.track.index()[c(2:6)]) {
      set.current.cell(sector.index = sect, track.index = sn)
      circos.rect(CELL_META$cell.xlim[1], CELL_META$cell.ylim[1], CELL_META$cell.xlim[2], CELL_META$cell.ylim[2], col = mycols[sn-1])
    }
  }
  circos.genomicIdeogram(species = genome, track.height = 0.1)
  lapply(ins_dfs, function(df) {
    df <- df[, c("LF_chr", "LF_start", "Ins_chr", "Ins_start")]
    lapply(1:nrow(df), function (x) {
      circos.link(sector.index1=df$Ins_chr[x], df$Ins_start[x], sector.index2=df$LF_chr[x], df$LF_start[x], directional = 1, arr.width = 0.1, arr.length = 0.1, arr.col = "black", col = "red")
    } )
  } )
  circos.clear()
  draw(legend, x = unit(0.03, "npc"), y = unit(0.75, "npc"), just = c("left", "top"))
  dev.off()
}
# ============================================================================
# Trouble-shooting
