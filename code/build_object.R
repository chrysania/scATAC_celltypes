# build seurat object

library(Signac)
library(Seurat)
library(Matrix)
source("code/utilities.R")

frag_path <- snakemake@input[['frags']]
barcodes <- snakemake@input[['barcodes']]
peak_counts_dir <- paste0(snakemake@input[['peak_counts']], "/")
annotations <- readRDS(snakemake@input[['annotations']])

# QC params
nCount_ATAC_above <- as.numeric(snakemake@params[["nCount_ATAC_above"]])
nCount_ATAC_below <- as.numeric(snakemake@params[["nCount_ATAC_below"]])
TSS_above <- as.numeric(snakemake@params[["TSS_above"]])
nucleosome_signal_below <- as.numeric(snakemake@params[["nucleosome_signal_below"]])

# read counts matrix
counts <- ReadCounts(peak_counts_dir)
ca <- CreateChromatinAssay(
    counts = counts,
    annotation = annotations
)

# make fragment object
bc <- read.table(barcodes, sep = "\t")
frag_obj <- CreateFragmentObject(frag_path, bc$V1)

# make seurat object
obj <- CreateSeuratObject(counts = ca, assay = 'ATAC')
Annotation(obj[["ATAC"]]) <- annotations
Fragments(obj) <- frag_obj

obj

# ATAC dimension reduc
DefaultAssay(obj) <- 'ATAC'

obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj)  
obj <- RunSVD(obj)
obj <- RunUMAP(obj, reduction = 'lsi', dims = 2:20, verbose = FALSE, reduction.name = 'umap.atac')
obj <- FindNeighbors(obj, reduction = "lsi", dims = 2:20)
obj <- FindClusters(obj)

obj <- TSSEnrichment(obj)
obj <- NucleosomeSignal(obj)

# ATAC QC
obj_subset <- subset(obj, subset = nCount_ATAC > nCount_ATAC_above & nCount_ATAC < nCount_ATAC_below & TSS.enrichment > TSS_above & nucleosome_signal < nucleosome_signal_below)

obj_subset

# save object
saveRDS(obj_subset, file = snakemake@output[['object']])