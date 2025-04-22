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
nucleosome_signal <- as.numeric(snakemake@params[["nucleosome_signal"]])

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

### debug
message(class(nCount_ATAC_above))
message(class(nCount_ATAC_below))
message(class(TSS_above))
message(class(nucleosome_signal))

message(print(nCount_ATAC_above))
message(print(nCount_ATAC_below))
message(print(TSS_above))
message(print(nucleosome_signal))

p1 <- DimPlot(obj, reduction = 'umap.atac')
p2 <- VlnPlot(obj, c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), pt.size=0)
p3 <- DensityScatter(obj, "nCount_ATAC", "TSS.enrichment", log_x = TRUE)

pdf(
  paste0("~/scratch/scATAC_celltypes/debug_QC.pdf"), 
  width = 6,      # Width of the PDF in inches
  height = 6       # Height of the PDF in inches
)
p1
p2
p3
dev.off()
# debug

# ATAC QC
obj_subset <- subset(obj, subset = nCount_ATAC > nCount_ATAC_above & nCount_ATAC < nCount_ATAC_below &
                     TSS.enrichment > TSS_above &
                     nucleosome_signal < nucleosome_signal)

#obj_subset <- subset(obj, subset = nCount_ATAC > 2000 & nCount_ATAC < 30000 &
#                     TSS.enrichment > 4 &
#                     nucleosome_signal < 1.5)

# save object
saveRDS(obj_subset, file = snakemake@output[['object']])