# combine_esophagus.r

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

esophagus2 <- readRDS("../data/esophagus/esophagus2/esophagus2_preQC.rds")
esophagus3 <- readRDS("../data/esophagus/esophagus3/esophagus3_preQC.rds")
# esophagus1 is bad

### QC ###
esophagus2_qc <- subset(esophagus2, subset = nCount_ATAC > 3000 & nCount_ATAC < 35000 &
              TSS.enrichment > 4 &
              nucleosome_signal < 1.5)
esophagus2_qc
esophagus3_qc <- subset(esophagus3, subset = nCount_ATAC > 2000 & nCount_ATAC < 30000 &
              TSS.enrichment > 4 &
              nucleosome_signal < 1.5)
esophagus3_qc


### Combine ###
# create common peak set
rownames_obj <- rownames(esophagus2_qc)
split_chr_coords <- strsplit(rownames_obj, "-")
obj_coordinates <- do.call(rbind, lapply(split_chr_coords, function(x) {
  c(x[1], x[2], x[3])
}))
obj_features <- as.data.frame(obj_coordinates, stringsAsFactors = FALSE)
colnames(obj_features) <- c("chr", "start", "end")
obj_features$start <- as.numeric(obj_features$start)
obj_features$end <- as.numeric(obj_features$end)
obj_features_clean <- na.omit(obj_features)
gr.esophagus2 <- makeGRangesFromDataFrame(obj_features_clean)

rownames_obj <- rownames(esophagus3_qc)
split_chr_coords <- strsplit(rownames_obj, "-")
obj_coordinates <- do.call(rbind, lapply(split_chr_coords, function(x) {
  c(x[1], x[2], x[3])
}))
obj_features <- as.data.frame(obj_coordinates, stringsAsFactors = FALSE)
colnames(obj_features) <- c("chr", "start", "end")
obj_features$start <- as.numeric(obj_features$start)
obj_features$end <- as.numeric(obj_features$end)
obj_features_clean <- na.omit(obj_features)
gr.esophagus3 <- makeGRangesFromDataFrame(obj_features_clean)

combined.peaks <- reduce(x = c(gr.esophagus2, gr.esophagus3))
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

## create new seurat obj
# esophagus
#sample_IDs=("ENCSR453TVZ-1" "ENCSR757EGB-1" "ENCSR164GSH-1")
#sample_n=("1" "2" "3")
#cell_type="esophagus"

# create fragment objects
frags.esophagus2 <- CreateFragmentObject(
  path = "../data/esophagus/encode_scatac_dcc_2/results/ENCSR757EGB-1/fragments/fragments.tsv.gz",
  cells = colnames(esophagus2_qc)
)
frags.esophagus3 <- CreateFragmentObject(
  path = "../data/esophagus/encode_scatac_dcc_2/results/ENCSR164GSH-1/fragments/fragments.tsv.gz",
  cells = colnames(esophagus3_qc)
)

# create peak x cell matrix
counts.esophagus2 <- FeatureMatrix(
  fragments = frags.esophagus2,
  features = combined.peaks,
  cells = colnames(esophagus2_qc)
)
counts.esophagus3 <- FeatureMatrix(
  fragments = frags.esophagus3,
  features = combined.peaks,
  cells = colnames(esophagus3_qc)
)

# create new seurat obj with combined.peaks
esophagus2_assay <- CreateChromatinAssay(counts.esophagus2, fragments = frags.esophagus2)
esophagus2_tomerge <- CreateSeuratObject(esophagus2_assay, assay = "ATAC")
esophagus2_tomerge

esophagus3_assay <- CreateChromatinAssay(counts.esophagus3, fragments = frags.esophagus3)
esophagus3_tomerge <- CreateSeuratObject(esophagus3_assay, assay = "ATAC")
esophagus3_tomerge

# add information to identify dataset of origin
esophagus2_tomerge$dataset <- 'esophagus2'
esophagus3_tomerge$dataset <- 'esophagus3'

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = esophagus2_tomerge,
  y = esophagus3_tomerge,
  add.cell.ids = c("2", "3")
)
combined[["ATAC"]]

# dimension reduction
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 10)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:30, reduction = 'lsi')
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:30)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3)

p1 <- DimPlot(object = combined, label = TRUE) + NoLegend()
p2 <- DimPlot(combined, group.by = 'dataset', label = TRUE) + NoLegend()
combined_plot <- p1 | p2
ggsave(filename = paste0("../data/plots/esophagus_combined.pdf"), plot = combined_plot, height = 6, width = 12)

saveRDS(combined, file = "../data/esophagus/esophagus_combined.rds")

### Integration ###
esophagus2_tomerge <- RenameCells(esophagus2_tomerge, add.cell.id = 2)
esophagus2_tomerge <- FindTopFeatures(esophagus2_tomerge, min.cutoff = 0)
esophagus2_tomerge <- RunTFIDF(esophagus2_tomerge)
esophagus2_tomerge <- RunSVD(esophagus2_tomerge)
esophagus2_tomerge <- RunUMAP(esophagus2_tomerge, dims = 2:50, reduction = 'lsi')

esophagus3_tomerge <- RenameCells(esophagus3_tomerge, add.cell.id = 3)
esophagus3_tomerge <- FindTopFeatures(esophagus3_tomerge, min.cutoff = 0)
esophagus3_tomerge <- RunTFIDF(esophagus3_tomerge)
esophagus3_tomerge <- RunSVD(esophagus3_tomerge)
esophagus3_tomerge <- RunUMAP(esophagus3_tomerge, dims = 2:50, reduction = 'lsi')

# find integration anchors 
integration.anchors <- FindIntegrationAnchors(
  object.list = list(esophagus2_tomerge, esophagus3_tomerge),
  anchor.features = rownames(esophagus2_tomerge),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
p2 <- DimPlot(integrated, group.by = 'dataset', pt.size = 0.1)

integrated <- FindNeighbors(object = integrated, reduction = 'integrated_lsi', dims = 2:30)
integrated <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3)
p3 <- DimPlot(object = integrated, label = TRUE) + NoLegend()

combined_plot <- p2 | p3
ggsave(filename = paste0("../data/plots/esophagus_integrated.pdf"), plot = combined_plot, height = 6, width = 12)
saveRDS(integrated, file = "../data/esophagus/esophagus_integrated.rds")

# change clustering resolution
integrated_modif <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3, resolution = 0.05)

# filter out clusters with <100 cells
table(integrated$seurat_clusters)

# Get cluster names that have more than 100 cells
big_clusters <- names(which(table(Idents(integrated)) >= 100))

# Subset the Seurat object to exclude these small clusters
integrated_filtered <- subset(integrated, idents = big_clusters)
integrated_filtered

saveRDS(integrated_filtered, file = paste0("../data/",cell_type,"/",cell_type,"_integrated_filtered.rds"))