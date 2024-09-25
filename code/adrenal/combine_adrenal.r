# combine_adrenal.r

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

adrenal1 <- readRDS("../data/adrenal/adrenal1/adrenal1_preQC.rds")
adrenal2 <- readRDS("../data/adrenal/adrenal2/adrenal2_preQC.rds")
adrenal3 <- readRDS("../data/adrenal/adrenal3/adrenal3_preQC.rds")

### QC ###
adrenal1_qc <- subset(adrenal1, subset = nCount_ATAC > 3000 & nCount_ATAC < 75000 &
              TSS.enrichment > 3.5 &
              nucleosome_signal < 1)
adrenal1_qc
adrenal2_qc <- subset(adrenal2, subset = nCount_ATAC > 3000 & nCount_ATAC < 100000 &
              TSS.enrichment > 3 &
              nucleosome_signal < 1)
adrenal2_qc
adrenal3_qc <- subset(adrenal3, subset = nCount_ATAC > 2000 & nCount_ATAC < 60000 &
              TSS.enrichment > 3 &
              nucleosome_signal < 1)
adrenal3_qc

### Combine ###
# create common peak set
rownames_adrenal <- rownames(adrenal1_qc)
split_chr_coords <- strsplit(rownames_adrenal, "-")
adrenal_coordinates <- do.call(rbind, lapply(split_chr_coords, function(x) {
  c(x[1], x[2], x[3])
}))
adrenal1_features <- as.data.frame(adrenal_coordinates, stringsAsFactors = FALSE)
colnames(adrenal1_features) <- c("chr", "start", "end")
adrenal1_features$start <- as.numeric(adrenal1_features$start)
adrenal1_features$end <- as.numeric(adrenal1_features$end)
adrenal1_features_clean <- na.omit(adrenal1_features)
gr.adrenal1 <- makeGRangesFromDataFrame(adrenal1_features_clean)

rownames_adrenal <- rownames(adrenal2_qc)
split_chr_coords <- strsplit(rownames_adrenal, "-")
adrenal_coordinates <- do.call(rbind, lapply(split_chr_coords, function(x) {
  c(x[1], x[2], x[3])
}))
adrenal2_features <- as.data.frame(adrenal_coordinates, stringsAsFactors = FALSE)
colnames(adrenal2_features) <- c("chr", "start", "end")
adrenal2_features$start <- as.numeric(adrenal2_features$start)
adrenal2_features$end <- as.numeric(adrenal2_features$end)
adrenal2_features_clean <- na.omit(adrenal2_features)
gr.adrenal2 <- makeGRangesFromDataFrame(adrenal2_features_clean)

rownames_adrenal <- rownames(adrenal3_qc)
split_chr_coords <- strsplit(rownames_adrenal, "-")
adrenal_coordinates <- do.call(rbind, lapply(split_chr_coords, function(x) {
  c(x[1], x[2], x[3])
}))
adrenal3_features <- as.data.frame(adrenal_coordinates, stringsAsFactors = FALSE)
colnames(adrenal3_features) <- c("chr", "start", "end")
adrenal3_features$start <- as.numeric(adrenal3_features$start)
adrenal3_features$end <- as.numeric(adrenal3_features$end)
adrenal3_features_clean <- na.omit(adrenal3_features)
gr.adrenal3 <- makeGRangesFromDataFrame(adrenal3_features_clean)

combined.peaks <- reduce(x = c(gr.adrenal1, gr.adrenal2, gr.adrenal3))
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# create new seurat obj
# adrenal gland
#sample_IDs=("ENCSR420EWQ-1" "ENCSR693GAD-1" "ENCSR194KHA-1")
#sample_n=("1" "2" "3")
#cell_type="adrenal"

# create fragment objects
frags.adrenal1 <- CreateFragmentObject(
  path = "../data/adrenal/encode_scatac_dcc_2/results/ENCSR420EWQ-1/fragments/fragments.tsv.gz",
  cells = colnames(adrenal1_qc)
)
frags.adrenal2 <- CreateFragmentObject(
  path = "../data/adrenal/encode_scatac_dcc_2/results/ENCSR693GAD-1/fragments/fragments.tsv.gz",
  cells = colnames(adrenal2_qc)
)
frags.adrenal3 <- CreateFragmentObject(
  path = "../data/adrenal/encode_scatac_dcc_2/results/ENCSR194KHA-1/fragments/fragments.tsv.gz",
  cells = colnames(adrenal3_qc)
)

# create peak x cell matrix
counts.adrenal1 <- FeatureMatrix(
  fragments = frags.adrenal1,
  features = combined.peaks,
  cells = colnames(adrenal1_qc)
)
counts.adrenal2 <- FeatureMatrix(
  fragments = frags.adrenal2,
  features = combined.peaks,
  cells = colnames(adrenal2_qc)
)
counts.adrenal3 <- FeatureMatrix(
  fragments = frags.adrenal3,
  features = combined.peaks,
  cells = colnames(adrenal3_qc)
)

# create new seurat obj with combined.peaks
adrenal1_assay <- CreateChromatinAssay(counts.adrenal1, fragments = frags.adrenal1)
adrenal1_tomerge <- CreateSeuratObject(adrenal1_assay, assay = "ATAC")

adrenal2_assay <- CreateChromatinAssay(counts.adrenal2, fragments = frags.adrenal2)
adrenal2_tomerge <- CreateSeuratObject(adrenal2_assay, assay = "ATAC")

adrenal3_assay <- CreateChromatinAssay(counts.adrenal3, fragments = frags.adrenal3)
adrenal3_tomerge <- CreateSeuratObject(adrenal3_assay, assay = "ATAC")

adrenal1_tomerge
adrenal2_tomerge
adrenal3_tomerge

# add information to identify dataset of origin
adrenal1_tomerge$dataset <- 'adrenal1'
adrenal2_tomerge$dataset <- 'adrenal2'
adrenal3_tomerge$dataset <- 'adrenal3'

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = adrenal1_tomerge,
  y = list(adrenal2_tomerge, adrenal3_tomerge),
  add.cell.ids = c("1", "2", "3")
)
combined[["ATAC"]]

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 10)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:30, reduction = 'lsi')
DimPlot(combined, group.by = 'dataset', pt.size = 0.1)

saveRDS(combined, file = "../data/adrenal/adrenal_combined.rds")

### Integrate ###
adrenal1_tomerge <- RenameCells(adrenal1_tomerge, add.cell.id = 1)
adrenal1_tomerge <- FindTopFeatures(adrenal1_tomerge, min.cutoff = 0)
adrenal1_tomerge <- RunTFIDF(adrenal1_tomerge)
adrenal1_tomerge <- RunSVD(adrenal1_tomerge)
adrenal1_tomerge <- RunUMAP(adrenal1_tomerge, dims = 2:50, reduction = 'lsi')

adrenal2_tomerge <- RenameCells(adrenal2_tomerge, add.cell.id = 2)
adrenal2_tomerge <- FindTopFeatures(adrenal2_tomerge, min.cutoff = 0)
adrenal2_tomerge <- RunTFIDF(adrenal2_tomerge)
adrenal2_tomerge <- RunSVD(adrenal2_tomerge)
adrenal2_tomerge <- RunUMAP(adrenal2_tomerge, dims = 2:50, reduction = 'lsi')

adrenal3_tomerge <- RenameCells(adrenal3_tomerge, add.cell.id = 3)
adrenal3_tomerge <- FindTopFeatures(adrenal3_tomerge, min.cutoff = 0)
adrenal3_tomerge <- RunTFIDF(adrenal3_tomerge)
adrenal3_tomerge <- RunSVD(adrenal3_tomerge)
adrenal3_tomerge <- RunUMAP(adrenal3_tomerge, dims = 2:50, reduction = 'lsi')

# find integration anchors 
integration.anchors <- FindIntegrationAnchors(
  object.list = list(adrenal1_tomerge, adrenal2_tomerge, adrenal3_tomerge),
  anchor.features = rownames(adrenal1_tomerge),
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
DimPlot(integrated, group.by = 'dataset', pt.size = 0.1)

integrated <- FindNeighbors(object = integrated, reduction = 'integrated_lsi', dims = 2:30)
integrated <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3)
DimPlot(object = integrated, label = TRUE) + NoLegend()

saveRDS(integrated, file = "../data/adrenal/adrenal_integrated.rds")