# combined heartLV heartRV

cell_type1 <- "heartLV"
sample_id1 <- c(4,7,8)

cell_type2 <- "heartRV"
sample_id2 <- c(1,2,5,9)

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

# heartLV
ct1_1 <- readRDS(paste0("../data/",cell_type1,"/",cell_type1,sample_id1[1],"/",cell_type1,sample_id1[1],"_preQC.rds"))
ct1_2 <- readRDS(paste0("../data/",cell_type1,"/",cell_type1,sample_id1[2],"/",cell_type1,sample_id1[2],"_preQC.rds"))
ct1_3 <- readRDS(paste0("../data/",cell_type1,"/",cell_type1,sample_id1[3],"/",cell_type1,sample_id1[3],"_preQC.rds"))

ct1_1_qc <- subset(ct1_1,
                  subset = nCount_ATAC > 3000 & nCount_ATAC < 35000 &
                  TSS.enrichment > 4 &
                  nucleosome_signal < 1)
ct1_1_qc
ct1_2_qc <- subset(ct1_2,
                  subset = nCount_ATAC > 6000 & nCount_ATAC < 100000 &
                  TSS.enrichment > 3 &
                  nucleosome_signal < 1)
ct1_2_qc
ct1_3_qc <- subset(ct1_3, 
                  subset = nCount_ATAC > 5000 & nCount_ATAC < 50000 &
                  TSS.enrichment > 3 &
                  nucleosome_signal < 1)
ct1_3_qc

# heartRV
ct2_1 <- readRDS(paste0("../data/",cell_type2,"/",cell_type2,sample_id2[1],"/",cell_type2,sample_id2[1],"_preQC.rds"))
ct2_2 <- readRDS(paste0("../data/",cell_type2,"/",cell_type2,sample_id2[2],"/",cell_type2,sample_id2[2],"_preQC.rds"))
ct2_3 <- readRDS(paste0("../data/",cell_type2,"/",cell_type2,sample_id2[3],"/",cell_type2,sample_id2[3],"_preQC.rds"))
ct2_4 <- readRDS(paste0("../data/",cell_type2,"/",cell_type2,sample_id2[4],"/",cell_type2,sample_id2[4],"_preQC.rds"))

ct2_1_qc <- subset(ct2_1,
                  subset = nCount_ATAC > 5000 & nCount_ATAC < 60000 &
                  TSS.enrichment > 3 &
                  nucleosome_signal < 1.5)
ct2_1_qc
ct2_2_qc <- subset(ct2_2,
                  subset = nCount_ATAC > 5000 & nCount_ATAC < 50000 &
                  TSS.enrichment > 3.5 &
                  nucleosome_signal < 1.5)
ct2_2_qc
ct2_3_qc <- subset(ct2_3, 
                  subset = nCount_ATAC > 10000 & nCount_ATAC < 150000 &
                  TSS.enrichment > 3 &
                  nucleosome_signal < 1.5)
ct2_3_qc
ct2_4_qc <- subset(ct2_4, 
                  subset = nCount_ATAC > 5000 & nCount_ATAC < 100000 &
                  TSS.enrichment > 3.5 &
                  nucleosome_signal < 1.5)
ct2_4_qc

### Combine ###
# create common peak set
rownames_list <- list()
rownames_list[[1]] <- rownames(ct1_1_qc)
rownames_list[[2]] <- rownames(ct1_2_qc)
rownames_list[[3]] <- rownames(ct1_3_qc)
rownames_list[[4]] <- rownames(ct2_1_qc)
rownames_list[[5]] <- rownames(ct2_2_qc)
rownames_list[[6]] <- rownames(ct2_3_qc)
rownames_list[[7]] <- rownames(ct2_4_qc)

gr_list <- list()
for (i in 1:7) {
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

combined.peaks <- reduce(x = c(gr_list[[1]], gr_list[[2]], gr_list[[3]], gr_list[[4]], gr_list[[5]], gr_list[[6]], gr_list[[7]]))
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# create new seurat obj
# heartLV (left ventricle)
sample_IDs1 <- c("ENCSR321AHR-1","ENCSR769WLL-1","ENCSR701JAT-1","ENCSR270CUT-1","ENCSR960IDI-1","ENCSR506ROZ-1","ENCSR862IWS-1","ENCSR088ZOL-1","ENCSR627IOJ-1","ENCSR020FAW-1")
#sample_n=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10")
#cell_type="heartLV"
# heartRV (right ventricle)
sample_IDs2 <- c("ENCSR588PEE-1","ENCSR681OLJ-1","ENCSR169BCG-1","ENCSR604PDO-1","ENCSR814OLA-1","ENCSR520ZUD-1","ENCSR579TPC-1","ENCSR615TSN-1","ENCSR517QNQ-1","ENCSR454YDZ-1")
#sample_n=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10")
#cell_type="heartRV"

colnames_list <- list()
colnames_list[[1]] <- colnames(ct1_1_qc)
colnames_list[[2]] <- colnames(ct1_2_qc)
colnames_list[[3]] <- colnames(ct1_3_qc)
colnames_list[[4]] <- colnames(ct2_1_qc)
colnames_list[[5]] <- colnames(ct2_2_qc)
colnames_list[[6]] <- colnames(ct2_3_qc)
colnames_list[[7]] <- colnames(ct2_4_qc)

# create fragment objects
frags.1 <- CreateFragmentObject(
  path = paste0("../data/",cell_type1,"/encode_scatac_dcc_2/results/",sample_IDs1[sample_id1[1]],"/fragments/fragments.tsv.gz"),
  cells = colnames_list[[1]]
)
frags.2 <- CreateFragmentObject(
  path = paste0("../data/",cell_type1,"/encode_scatac_dcc_2/results/",sample_IDs1[sample_id1[2]],"/fragments/fragments.tsv.gz"),
  cells = colnames_list[[2]]
)
frags.3 <- CreateFragmentObject(
  path = paste0("../data/",cell_type1,"/encode_scatac_dcc_2/results/",sample_IDs1[sample_id1[3]],"/fragments/fragments.tsv.gz"),
  cells = colnames_list[[3]]
)
frags.4 <- CreateFragmentObject(
  path = paste0("../data/",cell_type2,"/encode_scatac_dcc_2/results/",sample_IDs2[sample_id2[1]],"/fragments/fragments.tsv.gz"),
  cells = colnames_list[[4]]
)
frags.5 <- CreateFragmentObject(
  path = paste0("../data/",cell_type2,"/encode_scatac_dcc_2/results/",sample_IDs2[sample_id2[2]],"/fragments/fragments.tsv.gz"),
  cells = colnames_list[[5]]
)
frags.6 <- CreateFragmentObject(
  path = paste0("../data/",cell_type2,"/encode_scatac_dcc_2/results/",sample_IDs2[sample_id2[3]],"/fragments/fragments.tsv.gz"),
  cells = colnames_list[[6]]
)
frags.7 <- CreateFragmentObject(
  path = paste0("../data/",cell_type2,"/encode_scatac_dcc_2/results/",sample_IDs2[sample_id2[4]],"/fragments/fragments.tsv.gz"),
  cells = colnames_list[[7]]
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
counts.5 <- FeatureMatrix(
  fragments = frags.5,
  features = combined.peaks,
  cells = colnames_list[[5]]
)
counts.6 <- FeatureMatrix(
  fragments = frags.6,
  features = combined.peaks,
  cells = colnames_list[[6]]
)
counts.7 <- FeatureMatrix(
  fragments = frags.7,
  features = combined.peaks,
  cells = colnames_list[[7]]
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
chrom_assay5 <- CreateChromatinAssay(counts.5, fragments = frags.5)
obj5_tomerge <- CreateSeuratObject(chrom_assay5, assay = "ATAC")
chrom_assay6 <- CreateChromatinAssay(counts.6, fragments = frags.6)
obj6_tomerge <- CreateSeuratObject(chrom_assay6, assay = "ATAC")
chrom_assay7 <- CreateChromatinAssay(counts.7, fragments = frags.7)
obj7_tomerge <- CreateSeuratObject(chrom_assay7, assay = "ATAC")

obj1_tomerge
obj2_tomerge
obj3_tomerge
obj4_tomerge
obj5_tomerge
obj6_tomerge
obj7_tomerge

# add information to identify dataset of origin
obj1_tomerge$dataset <- paste0(cell_type1,sample_id1[1])
obj2_tomerge$dataset <- paste0(cell_type1,sample_id1[2])
obj3_tomerge$dataset <- paste0(cell_type1,sample_id1[3])
obj4_tomerge$dataset <- paste0(cell_type2,sample_id2[1])
obj5_tomerge$dataset <- paste0(cell_type2,sample_id2[2])
obj6_tomerge$dataset <- paste0(cell_type2,sample_id2[3])
obj7_tomerge$dataset <- paste0(cell_type2,sample_id2[4])

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = obj1_tomerge,
  y = list(obj2_tomerge, obj3_tomerge, obj4_tomerge, obj5_tomerge, obj6_tomerge, obj7_tomerge),
  add.cell.ids = c("4", "7", "8", "1", "2", "5", "9")
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
ggsave(filename = paste0("../data/plots/hearts_combined.pdf"), plot = combined_plot, height = 6, width = 12)

saveRDS(combined, file = paste0("../data/hearts_combined.rds"))

### Integrate ###
obj1_tomerge <- RenameCells(obj1_tomerge, add.cell.id = 4)
obj1_tomerge <- FindTopFeatures(obj1_tomerge, min.cutoff = 0)
obj1_tomerge <- RunTFIDF(obj1_tomerge)
obj1_tomerge <- RunSVD(obj1_tomerge)
obj1_tomerge <- RunUMAP(obj1_tomerge, dims = 2:50, reduction = 'lsi')

obj2_tomerge <- RenameCells(obj2_tomerge, add.cell.id = 7)
obj2_tomerge <- FindTopFeatures(obj2_tomerge, min.cutoff = 0)
obj2_tomerge <- RunTFIDF(obj2_tomerge)
obj2_tomerge <- RunSVD(obj2_tomerge)
obj2_tomerge <- RunUMAP(obj2_tomerge, dims = 2:50, reduction = 'lsi')

obj3_tomerge <- RenameCells(obj3_tomerge, add.cell.id = 8)
obj3_tomerge <- FindTopFeatures(obj3_tomerge, min.cutoff = 0)
obj3_tomerge <- RunTFIDF(obj3_tomerge)
obj3_tomerge <- RunSVD(obj3_tomerge)
obj3_tomerge <- RunUMAP(obj3_tomerge, dims = 2:50, reduction = 'lsi')

obj4_tomerge <- RenameCells(obj4_tomerge, add.cell.id = 1)
obj4_tomerge <- FindTopFeatures(obj4_tomerge, min.cutoff = 0)
obj4_tomerge <- RunTFIDF(obj4_tomerge)
obj4_tomerge <- RunSVD(obj4_tomerge)
obj4_tomerge <- RunUMAP(obj4_tomerge, dims = 2:50, reduction = 'lsi')

obj5_tomerge <- RenameCells(obj5_tomerge, add.cell.id = 2)
obj5_tomerge <- FindTopFeatures(obj5_tomerge, min.cutoff = 0)
obj5_tomerge <- RunTFIDF(obj5_tomerge)
obj5_tomerge <- RunSVD(obj5_tomerge)
obj5_tomerge <- RunUMAP(obj5_tomerge, dims = 2:50, reduction = 'lsi')

obj6_tomerge <- RenameCells(obj6_tomerge, add.cell.id = 5)
obj6_tomerge <- FindTopFeatures(obj6_tomerge, min.cutoff = 0)
obj6_tomerge <- RunTFIDF(obj6_tomerge)
obj6_tomerge <- RunSVD(obj6_tomerge)
obj6_tomerge <- RunUMAP(obj6_tomerge, dims = 2:50, reduction = 'lsi')

obj7_tomerge <- RenameCells(obj7_tomerge, add.cell.id = 9)
obj7_tomerge <- FindTopFeatures(obj7_tomerge, min.cutoff = 0)
obj7_tomerge <- RunTFIDF(obj7_tomerge)
obj7_tomerge <- RunSVD(obj7_tomerge)
obj7_tomerge <- RunUMAP(obj7_tomerge, dims = 2:50, reduction = 'lsi')

# find integration anchors 
integration.anchors <- FindIntegrationAnchors(
  object.list = list(obj1_tomerge, obj2_tomerge, obj3_tomerge, obj4_tomerge, obj5_tomerge, obj6_tomerge, obj7_tomerge),
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
ggsave(filename = paste0("../data/plots/hearts_integrated.pdf"), plot = combined_plot, height = 6, width = 12)
saveRDS(integrated, file = paste0("../data/hearts_integrated.rds"))