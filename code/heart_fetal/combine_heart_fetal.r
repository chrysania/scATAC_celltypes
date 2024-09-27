# combine_heart_fetal.r
# dataset 1 3 7 9

cell_type <- "heart_fetal"
sample_id <- c(1,3,7,9)

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

#heart_fetal1 <- readRDS("../data/heart_fetal/heart_fetal1/heart_fetal1_preQC.rds")
heart_fetal1 <- readRDS(paste0("../data/",cell_type,"/",cell_type,sample_id[1],"/",cell_type,sample_id[1],"_preQC.rds"))
heart_fetal3 <- readRDS(paste0("../data/",cell_type,"/",cell_type,sample_id[2],"/",cell_type,sample_id[2],"_preQC.rds"))
heart_fetal7 <- readRDS(paste0("../data/",cell_type,"/",cell_type,sample_id[3],"/",cell_type,sample_id[3],"_preQC.rds"))
heart_fetal9 <- readRDS(paste0("../data/",cell_type,"/",cell_type,sample_id[4],"/",cell_type,sample_id[4],"_preQC.rds"))

### QC ###
heart_fetal1_qc <- subset(heart_fetal1, 
                          subset = nCount_ATAC > 5000 & nCount_ATAC < 100000 &
                          TSS.enrichment > 3.5 &
                          nucleosome_signal < 1)
heart_fetal1_qc
heart_fetal3_qc <- subset(heart_fetal3, 
                          subset = nCount_ATAC > 5000 & nCount_ATAC < 100000 &
                          TSS.enrichment > 3 &
                          nucleosome_signal < 1)
heart_fetal3_qc
heart_fetal7_qc <- subset(heart_fetal7, 
                          subset = nCount_ATAC > 6000 & nCount_ATAC < 100000 &
                          TSS.enrichment > 4 &
                          nucleosome_signal < 1.5)
heart_fetal7_qc
heart_fetal9_qc <- subset(heart_fetal9, 
                          subset = nCount_ATAC > 10000 & nCount_ATAC < 100000 &
                          TSS.enrichment > 3.5 &
                          nucleosome_signal < 1.5)
heart_fetal9_qc

### Combine ###
# create common peak set
rownames_list <- list()
rownames_list[[1]] <- rownames(heart_fetal1_qc)
rownames_list[[2]] <- rownames(heart_fetal3_qc)
rownames_list[[3]] <- rownames(heart_fetal7_qc)
rownames_list[[4]] <- rownames(heart_fetal9_qc)

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

combined.peaks <- reduce(x = c(gr_list[[1]], gr_list[[2]], gr_list[[3]], gr_list[[4]]))
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# create new seurat obj
# heart_fetal 
sample_IDs <- c("ENCSR515SNH-1","ENCSR715JSZ-1","ENCSR282FAK-1","ENCSR890TGR-1","ENCSR024TGD-1","ENCSR004IAY-1","ENCSR376IBI-1","ENCSR306FRQ-1","ENCSR805DID-1")
#sample_n=("1" "2" "3" "4" "5" "6" "7" "8" "9")
#cell_type="heart_fetal"

colnames_list <- list()
colnames_list[[1]] <- colnames(heart_fetal1_qc)
colnames_list[[2]] <- colnames(heart_fetal3_qc)
colnames_list[[3]] <- colnames(heart_fetal7_qc)
colnames_list[[4]] <- colnames(heart_fetal9_qc)

# create fragment objects
frags.1 <- CreateFragmentObject(
  path = paste0("../data/",cell_type,"/encode_scatac_dcc_2/results/",sample_IDs[sample_id[1]],"/fragments/fragments.tsv.gz"),
  cells = colnames_list[[1]]
)
frags.2 <- CreateFragmentObject(
  path = paste0("../data/",cell_type,"/encode_scatac_dcc_2/results/",sample_IDs[sample_id[2]],"/fragments/fragments.tsv.gz"),
  cells = colnames_list[[2]]
)
frags.3 <- CreateFragmentObject(
  path = paste0("../data/",cell_type,"/encode_scatac_dcc_2/results/",sample_IDs[sample_id[3]],"/fragments/fragments.tsv.gz"),
  cells = colnames_list[[3]]
)
frags.4 <- CreateFragmentObject(
  path = paste0("../data/",cell_type,"/encode_scatac_dcc_2/results/",sample_IDs[sample_id[4]],"/fragments/fragments.tsv.gz"),
  cells = colnames_list[[4]]
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
counts.3 <- FeatureMatrix(
  fragments = frags.3,
  features = combined.peaks,
  cells = colnames_list[[3]]
)
counts.4 <- FeatureMatrix(
  fragments = frags.4,
  features = combined.peaks,
  cells = colnames_list[[4]]
)

# create new seurat obj with combined.peaks
chrom_assay1 <- CreateChromatinAssay(counts.1, fragments = frags.1)
obj1_tomerge <- CreateSeuratObject(chrom_assay1, assay = "ATAC")
chrom_assay2 <- CreateChromatinAssay(counts.2, fragments = frags.2)
obj2_tomerge <- CreateSeuratObject(chrom_assay2, assay = "ATAC")
chrom_assay3 <- CreateChromatinAssay(counts.3, fragments = frags.3)
obj3_tomerge <- CreateSeuratObject(chrom_assay3, assay = "ATAC")
chrom_assay4 <- CreateChromatinAssay(counts.4, fragments = frags.4)
obj4_tomerge <- CreateSeuratObject(chrom_assay4, assay = "ATAC")

obj1_tomerge
obj2_tomerge
obj3_tomerge
obj4_tomerge

# add information to identify dataset of origin
obj1_tomerge$dataset <- paste0(cell_type,sample_id[1])
obj2_tomerge$dataset <- paste0(cell_type,sample_id[2])
obj3_tomerge$dataset <- paste0(cell_type,sample_id[3])
obj4_tomerge$dataset <- paste0(cell_type,sample_id[4])

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = obj1_tomerge,
  y = list(obj2_tomerge, obj3_tomerge, obj4_tomerge),
  add.cell.ids = c("1", "3", "7", "9")
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

obj3_tomerge <- RenameCells(obj3_tomerge, add.cell.id = 7)
obj3_tomerge <- FindTopFeatures(obj3_tomerge, min.cutoff = 0)
obj3_tomerge <- RunTFIDF(obj3_tomerge)
obj3_tomerge <- RunSVD(obj3_tomerge)
obj3_tomerge <- RunUMAP(obj3_tomerge, dims = 2:50, reduction = 'lsi')

obj4_tomerge <- RenameCells(obj4_tomerge, add.cell.id = 9)
obj4_tomerge <- FindTopFeatures(obj4_tomerge, min.cutoff = 0)
obj4_tomerge <- RunTFIDF(obj4_tomerge)
obj4_tomerge <- RunSVD(obj4_tomerge)
obj4_tomerge <- RunUMAP(obj4_tomerge, dims = 2:50, reduction = 'lsi')

# find integration anchors 
integration.anchors <- FindIntegrationAnchors(
  object.list = list(obj1_tomerge, obj2_tomerge, obj3_tomerge, obj4_tomerge),
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

# change clustering resolution
integrated_modif <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3, resolution = 0.25)

# filter out clusters with <100 cells
table(integrated$seurat_clusters)

# Get cluster names that have more than 100 cells
big_clusters <- names(which(table(Idents(integrated)) >= 100))

# Subset the Seurat object to exclude these small clusters
integrated_filtered <- subset(integrated, idents = big_clusters)
integrated_filtered

saveRDS(integrated_filtered, file = paste0("../data/",cell_type,"/",cell_type,"_integrated_filtered.rds"))