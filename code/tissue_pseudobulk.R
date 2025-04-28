library(Signac)
library(Seurat)
library(Matrix)
library(GenomicRanges)
source("code/utilities.R")
library(tidyr)

obj <- readRDS(snakemake@input[['tissue_integrated']])

ccre <- read.table("data/combined_cre.bed", sep="\t")
ccre_gr <- GRanges(seqnames = ccre$V1,
                   ranges = IRanges(start = ccre$V2, end = ccre$V3))

# ccre x cells matrix
atac_assay <- obj[["ATAC"]]
frags <- Fragments(atac_assay)

ccre_counts <- FeatureMatrix(
  fragments = frags,
  features = ccre_gr,
  cells = colnames(obj)
)

chrom_assay <- CreateChromatinAssay(ccre_counts, fragments = frags)
obj_ccre <- CreateSeuratObject(chrom_assay, assay = "ATAC")
obj_ccre

# BinaryIdentMatrix
big_clusters <- names(which(table(Idents(obj_ccre)) >= 100))
binary_matrix <- Signac:::BinaryIdentMatrix(object = obj,
                                           idents = big_clusters)

# normalize by number of cells
rowsum_bm <- rowSums(binary_matrix)
bm_norm <- (binary_matrix / rowsum_bm) * 1000

# cell type x ccre matrix
ct_ccre <- bm_norm %*% t(ccre_counts)

saveRDS(ct_ccre, file = snakemake@output[['pseudobulk']])