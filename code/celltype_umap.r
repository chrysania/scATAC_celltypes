library(Signac)
library(Seurat)
library(Matrix)
library(patchwork)
library(ggplot2)
library(rtracklayer)
source("utilities.R")

args = commandArgs(trailingOnly = TRUE)
sample_id <- as.numeric(args[1])
cell_type <- as.character(args[2])

annotations <- readRDS("../data/annotations.rds")

# counts_dir <- paste0("../data/",cell_type,"/",cell_type,sample_id,"/peak_matrix/matrix.mtx.gz")
# barcodes_dir <- paste0("../data/",cell_type,"/",cell_type,sample_id,"/peak_matrix/barcodes.tsv")
# peaks_dir <- paste0("../data/",cell_type,"/",cell_type,sample_id,"/peak_matrix/features.tsv.gz")

counts <- ReadCounts(dir = paste0("../data/",cell_type,"/",cell_type,sample_id,"/peak_matrix/"))

# counts <- Matrix::readMM(counts_dir)
# barcodes <- readLines(barcodes_dir)
# peaks <- readLines(peaks_dir)
# peaknames <- paste(peaks$V1, peaks$V2, peaks$V3, collapse="-")
# colnames(counts) <- barcodes
# rownames(counts) <- peaknames

obj <- CreateSeuratObject(
  counts = counts,
  assay = "peaks"
)

obj_qc <- obj[, obj$nCount_peaks > 1000 & obj$nCount_peaks < 30000]

# histogram
#pdf(paste0("../data/plots/histograms/hist_nCount_peaks_", cell_type, sample_id,".pdf"))
#hist(obj_qc$nCount_peaks, 
#     breaks = 100, 
#     main = "Histogram of nCount_peaks", 
#     xlab = "nCount_peaks", 
#     col = "lightblue")
#hist(log10(obj_qc$nCount_peaks + 1), 
#     breaks = 100, 
#     main = "Histogram of Log-Transformed nCount_peaks", 
#     xlab = "Log10(nCount_peaks + 1)", 
#     col = "lightblue")
#dev.off()

# UMAP
obj_qc <- RunTFIDF(obj_qc)
obj_qc <- FindTopFeatures(obj_qc)
obj_qc <- RunSVD(obj_qc)
obj_qc <- RunUMAP(obj_qc, reduction = 'lsi', dims = 2:20, verbose = FALSE)
obj_qc <- FindNeighbors(object = obj_qc, reduction = 'lsi', dims = 2:20)
obj_qc <- FindClusters(object = obj_qc, verbose = FALSE, algorithm = 3)

p1 <- DimPlot(obj_qc, label=TRUE) + ggtitle(paste0(cell_type," dataset ", sample_id))
p2 <- FeaturePlot(obj_qc, 'nCount_peaks', max.cutoff = 'q99')
dimplot <- p1 | p2

ggsave(filename = paste0("../data/plots/",cell_type,"_dataset_", sample_id,"_UMAP.pdf"), plot = dimplot, height = 4, width = 8)