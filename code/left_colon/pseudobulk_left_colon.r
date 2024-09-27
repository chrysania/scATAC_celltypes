# psudobulk_left_colon.r

cell_type <- "left_colon"

setwd('/charonfs/scratch/users/astar/gis/limchr/scATAC_celltypes/code')
library(Signac)
library(Seurat)
library(Matrix)
library(patchwork)
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)
source("utilities.R")
library(tidyr)

obj <- readRDS(paste0("../data/",cell_type,"/",cell_type,"_integrated_filtered.rds"))
obj

ccre <- read.table("../data/combined_cre.bed", sep="\t")
ccre_gr <- GRanges(seqnames = ccre$V1,
                   ranges = IRanges(start = ccre$V2, end = ccre$V3))

# ccre x cells matrix
atac_assay <- obj[["ATAC"]]
frags <- Fragments(atac_assay)

ccre_counts <- FeatureMatrix( # this will take a while
  fragments = frags,
  features = ccre_gr,
  cells = colnames(obj)
)

chrom_assay <- CreateChromatinAssay(ccre_counts, fragments = frags)
obj_ccre <- CreateSeuratObject(chrom_assay, assay = "ATAC")
obj_ccre

# BinaryIdentMatrix
binary_matrix <- Signac:::BinaryIdentMatrix(object = obj_ccre)
# normalize by number of cells
rowsum_bm <- rowSums(binary_matrix)
bm_norm <- (binary_matrix / rowsum_bm) * 1000

# cell type x ccre matrix
ct_ccre <- bm_norm %*% t(ccre_counts)

saveRDS(ct_ccre, paste0("../data/",cell_type,"/",cell_type,"_pseudobulk.rds"))
