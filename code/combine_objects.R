library(Signac)
library(Seurat)
library(Matrix)
library(patchwork)
library(ggplot2)
library(GenomicRanges)
source("code/utilities.R")
library(tidyr)
options(scipen=999)

tissue_name <- snakemake@params[['tissue_name']]   
clustering_res <- as.numeric(snakemake@params[['clustering_res']])

if (tissue_name == "adrenal") {
    sample_names <- c("adrenal1", "adrenal2")
}
if (tissue_name == "esophagus") {
    sample_names <- c("esophagus2", "esophagus3")
}
if (tissue_name == "heartRV") {
    sample_names <- c("heartRV1", "heartRV2", "heartRV5", "heartRV9")
}
if (tissue_name == "heart_fetal") {
    sample_names <- c("heart_fetal1", "heart_fetal3", "heart_fetal7", "heart_fetal9")
}
if (tissue_name == "left_colon") {
    sample_names <- c("left_colon1", "left_colon3")
}
if (tissue_name == "liver") {
    sample_names <- c("liver3", "liver7")
}
if (tissue_name == "psoas_muscle") {
    sample_names <- c("psoas_muscle1", "psoas_muscle2", "psoas_muscle3")
}

# load objects
obj_list <- list()
features_list <- list() # features
barcodes_list <- list() # cell barcodes
for (i in seq_along(sample_names)) {
    obj <- readRDS(paste0("objects/",sample_names[i],"_peaks.rds"))
    obj_list[[sample_names[i]]] <- obj
    features_list[[sample_names[i]]] <- rownames(obj)
    barcodes_list[[sample_names[i]]] <- colnames(obj)
}

# make common feature set
gr_list <- list()
for (i in seq_along(sample_names)) {
    features_obj <- features_list[[sample_names[i]]]
    split_chr_coords <- strsplit(features_obj, "-")
    obj_coordinates <- do.call(rbind, lapply(split_chr_coords, function(x) {
      c(x[1], x[2], x[3])
    }))
    obj_features_df <- as.data.frame(obj_coordinates, stringsAsFactors = FALSE)
    colnames(obj_features_df) <- c("chr", "start", "end")

    obj_features_df$start <- as.numeric(obj_features_df$start)
    obj_features_df$end <- as.numeric(obj_features_df$end)
    obj_features_clean <- na.omit(obj_features_df)
    gr <- makeGRangesFromDataFrame(obj_features_clean)

    gr_list[[sample_names[i]]] <- gr
}

combined <- do.call(c, unname(gr_list))
combined.peaks <- reduce(combined)
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]

# create new seurat object
merger_list <- list()
for (i in seq_along(sample_names)) {
    fragpath <- paste0("data/",tissue_name,"/",sample_names[i],"/fragments.tsv.gz")
    frags <- CreateFragmentObject(
        path = fragpath,
        cells = barcodes_list[[sample_names[i]]]
    )

    counts <- FeatureMatrix(
        fragments = frags,
        features = combined.peaks,
        cells = barcodes_list[[sample_names[i]]]
    )

    chrom_assay <- CreateChromatinAssay(counts, fragments = frags)

    new_obj <- CreateSeuratObject(chrom_assay, assay = "ATAC")
    new_obj$dataset <- sample_names[i]

    merger_list[sample_names[i]] <- new_obj
}

# merge all datasets
cell_id_vector <- as.integer(sub(".*([0-9])$", "\\1", sample_names))
combined <- merge(
    x = merger_list[[1]],
    y = merger_list[-1],
    add.cell.ids = cell_id_vector
)

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 10)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:30, reduction = 'lsi')
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:30)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3)

saveRDS(combined, file = snakemake@output[['combined']])

# integrate
for (i in seq_along(sample_names)) {
    merger_list[[sample_names[i]]] <- RenameCells(merger_list[[sample_names[i]]], add.cell.id = cell_id_vector[i])
    merger_list[[sample_names[i]]] <- FindTopFeatures(merger_list[[sample_names[i]]], min.cutoff = 0)
    merger_list[[sample_names[i]]] <- RunTFIDF(merger_list[[sample_names[i]]])
    merger_list[[sample_names[i]]] <- RunSVD(merger_list[[sample_names[i]]])
    merger_list[[sample_names[i]]] <- RunUMAP(merger_list[[sample_names[i]]], dims = 2:50, reduction = 'lsi')
}

# find integration anchors 
integration.anchors <- FindIntegrationAnchors(
  object.list = merger_list,
  anchor.features = rownames(merger_list[[1]]),
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
integrated <- FindNeighbors(object = integrated, reduction = 'integrated_lsi', dims = 2:30)
integrated <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3, resolution = clustering_res)

# filter out clusters with <100 cells
# Get cluster names that have more than 100 cells
big_clusters <- names(which(table(Idents(integrated)) >= 100))
integrated_filtered <- subset(integrated, idents = big_clusters)

integrated_filtered

# save object as .rds
saveRDS(integrated_filtered, file = snakemake@output[['integrated']])