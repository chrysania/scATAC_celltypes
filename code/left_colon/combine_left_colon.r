# combine_left_colon.r
# dataset 1 3

cell_type <- "left_colon"
sample_id <- c(1,3)

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

ct_1 <- readRDS(paste0("../data/",cell_type,"/",cell_type,sample_id[1],"/",cell_type,sample_id[1],"_preQC.rds"))
ct_2 <- readRDS(paste0("../data/",cell_type,"/",cell_type,sample_id[2],"/",cell_type,sample_id[2],"_preQC.rds"))


### QC ###
ct_1_qc <- subset(ct_1,
                  subset = nCount_ATAC > 8000 & nCount_ATAC < 100000 &
                  TSS.enrichment > 3 &
                  nucleosome_signal < 1)
ct_1_qc
ct_2_qc <- subset(ct_2,
                  subset = nCount_ATAC > 10000 & nCount_ATAC < 200000 &
                  TSS.enrichment > 3.5 &
                  nucleosome_signal < 1)
ct_2_qc


### Combine ###
# create common peak set
rownames_list <- list()
rownames_list[[1]] <- rownames(ct_1_qc)
rownames_list[[2]] <- rownames(ct_2_qc)

gr_list <- list()
for (i in 1:length(sample_id)) {
    rownames_obj <- rownames_list[[i]]
    split_chr_coords <- strsplit(rownames_obj, "-")
    obj_coordinates <- do.call(rbind, lapply(split_chr_coords, function(x) {
      c(x[1], x[2], x[3])
    }))
    obj_features <- as.data.frame(obj_coordinates, stringsAsFactors = FALSE)
    colnames(obj_features) <- c("chr", "start", "end")
    obj_features$start <- as.numeric(obj_features$start)
    obj_features$end <- as.numeric(obj_features$end)
    obj_features_clean <- na.omit(obj_features)
    gr <- makeGRangesFromDataFrame(obj_features_clean)

    gr_list[[i]] <- gr
}

combined.peaks <- reduce(x = c(gr_list[[1]], gr_list[[2]]))
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# create new seurat obj
# left colon
sample_IDs <- c("ENCSR830FPR-1","ENCSR916RYB-1","ENCSR904WIW-1")
#sample_n=("1" "2" "3")
#cell_type="left_colon"

colnames_list <- list()
colnames_list[[1]] <- colnames(ct_1_qc)
colnames_list[[2]] <- colnames(ct_2_qc)


# create fragment objects
frags.1 <- CreateFragmentObject(
  path = paste0("../data/",cell_type,"/encode_scatac_dcc_2/results/",sample_IDs[sample_id[1]],"/fragments/fragments.tsv.gz"),
  cells = colnames_list[[1]]
)
frags.2 <- CreateFragmentObject(
  path = paste0("../data/",cell_type,"/encode_scatac_dcc_2/results/",sample_IDs[sample_id[2]],"/fragments/fragments.tsv.gz"),
  cells = colnames_list[[2]]
)


# create peak x cell matrix
counts.1 <- FeatureMatrix(
  fragments = frags.1,
  features = combined.peaks,
  cells = colnames_list[[1]]
)
counts.2 <- FeatureMatrix(
  fragments = frags.2,
  features = combined.peaks,
  cells = colnames_list[[2]]
)


# create new seurat obj with combined.peaks
chrom_assay1 <- CreateChromatinAssay(counts.1, fragments = frags.1)
obj1_tomerge <- CreateSeuratObject(chrom_assay1, assay = "ATAC")
chrom_assay2 <- CreateChromatinAssay(counts.2, fragments = frags.2)
obj2_tomerge <- CreateSeuratObject(chrom_assay2, assay = "ATAC")

obj1_tomerge
obj2_tomerge

# add information to identify dataset of origin
obj1_tomerge$dataset <- paste0(cell_type,sample_id[1])
obj2_tomerge$dataset <- paste0(cell_type,sample_id[2])

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = obj1_tomerge,
  y = obj2_tomerge,
  add.cell.ids = c("1", "3")
)
combined[["ATAC"]]

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 10)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:30, reduction = 'lsi')
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:30)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3)

p1 <- DimPlot(object = combined, label = TRUE) + NoLegend()
p2 <- DimPlot(combined, group.by = 'dataset', label = TRUE) + NoLegend()
combined_plot <- p1 | p2
ggsave(filename = paste0("../data/plots/",cell_type,"_combined.pdf"), plot = combined_plot, height = 6, width = 12)

saveRDS(combined, file = paste0("../data/",cell_type,"/",cell_type,"_combined.rds"))

### Integrate ###
obj1_tomerge <- RenameCells(obj1_tomerge, add.cell.id = 1)
obj1_tomerge <- FindTopFeatures(obj1_tomerge, min.cutoff = 0)
obj1_tomerge <- RunTFIDF(obj1_tomerge)
obj1_tomerge <- RunSVD(obj1_tomerge)
obj1_tomerge <- RunUMAP(obj1_tomerge, dims = 2:50, reduction = 'lsi')

obj2_tomerge <- RenameCells(obj2_tomerge, add.cell.id = 3)
obj2_tomerge <- FindTopFeatures(obj2_tomerge, min.cutoff = 0)
obj2_tomerge <- RunTFIDF(obj2_tomerge)
obj2_tomerge <- RunSVD(obj2_tomerge)
obj2_tomerge <- RunUMAP(obj2_tomerge, dims = 2:50, reduction = 'lsi')


# find integration anchors 
integration.anchors <- FindIntegrationAnchors(
  object.list = list(obj1_tomerge, obj2_tomerge),
  anchor.features = rownames(obj1_tomerge),
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
p1 <- DimPlot(integrated, group.by = 'dataset', pt.size = 0.1)

integrated <- FindNeighbors(object = integrated, reduction = 'integrated_lsi', dims = 2:30)
integrated <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3)
p2 <- DimPlot(object = integrated, label = TRUE) + NoLegend()

combined_plot <- p1 | p2
ggsave(filename = paste0("../data/plots/",cell_type,"_integrated.pdf"), plot = combined_plot, height = 6, width = 12)
saveRDS(integrated, file = paste0("../data/",cell_type,"/",cell_type,"_integrated.rds"))
