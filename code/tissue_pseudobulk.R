library(Signac)
library(Seurat)
library(Matrix)
library(GenomicRanges)
source("code/utilities.R")
library(tidyr)

obj <- readRDS(snakemake@input[['tissue_integrated']])
tissue_name <- snakemake@params[['tissue_name']]
clustering_resolution <- as.numeric(snakemake@params[['clustering_resolution']])

message(clustering_resolution)
message(class(clustering_resolution))

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

saveRDS(obj_ccre, file = paste0("objects/",tissue_name,"_ccre.rds"))
#obj_ccre <- readRDS(paste0("objects/",tissue_name,"_ccre.rds"))

obj_ccre <- RunTFIDF(obj_ccre)
obj_ccre <- FindTopFeatures(obj_ccre)  
obj_ccre <- RunSVD(obj_ccre)
obj_ccre <- RunUMAP(obj_ccre, reduction = 'lsi', dims = 2:20, verbose = FALSE, reduction.name = 'umap.atac')
obj_ccre <- FindNeighbors(obj_ccre, reduction = "lsi", dims = 2:20)
obj_ccre <- FindClusters(obj_ccre, resolution = clustering_resolution) 

# BinaryIdentMatrix
big_clusters <- names(which(table(Idents(obj_ccre)) >= 100))
binary_matrix <- Signac:::BinaryIdentMatrix(object = obj_ccre,
                                           idents = big_clusters)

# normalize by number of cells
rowsum_bm <- rowSums(binary_matrix)
bm_norm <- (binary_matrix / rowsum_bm) * 1000

# get cell type x ccre matrix
counts_matrix <- GetAssayData(obj_ccre, assay = 'ATAC', layer = 'data')
ct_ccre <- bm_norm %*% t(counts_matrix)

# add rownames
matrix_rownames <- paste0(tissue_name, "_cluster", rownames(ct_ccre))
rownames(ct_ccre) <- matrix_rownames

saveRDS(ct_ccre, file = snakemake@output[['pseudobulk']])