# QC
args = commandArgs(trailingOnly = TRUE)
sample_id <- as.numeric(args[1])
cell_type <- as.character(args[2])

library(Signac)
library(Seurat)
library(Matrix)
library(patchwork)
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)
source("utilities.R")
library(tidyr)
options(scipen = 999)

annotations <- readRDS("../data/annotations.rds")

### datasets ###
if (cell_type == "adrenal") { # adrenal gland
    sample_IDs <- c("ENCSR420EWQ-1","ENCSR693GAD-1","ENCSR194KHA-1")
    sample_n <- 1:3
}
if (cell_type == "esophagus") { # esophagus
    sample_IDs <- c("ENCSR453TVZ-1","ENCSR757EGB-1","ENCSR164GSH-1")
    sample_n <- 1:3
}
if (cell_type == "heart_fetal") { # heart_fetal
    sample_IDs <- c("ENCSR515SNH-1","ENCSR715JSZ-1","ENCSR282FAK-1","ENCSR890TGR-1","ENCSR024TGD-1","ENCSR004IAY-1","ENCSR376IBI-1","ENCSR306FRQ-1","ENCSR805DID-1")
    sample_n <- 1:9
}
if (cell_type == "heartLV") { # heartLV
    sample_IDs <- c("ENCSR321AHR-1","ENCSR769WLL-1","ENCSR701JAT-1","ENCSR270CUT-1","ENCSR960IDI-1","ENCSR506ROZ-1","ENCSR862IWS-1","ENCSR088ZOL-1","ENCSR627IOJ-1","ENCSR020FAW-1")
    sample_n <- 1:10
}
if (cell_type == "heartRV") { # heartRV
    sample_IDs <- c("ENCSR588PEE-1","ENCSR681OLJ-1","ENCSR169BCG-1","ENCSR604PDO-1","ENCSR814OLA-1","ENCSR520ZUD-1","ENCSR579TPC-1","ENCSR615TSN-1","ENCSR517QNQ-1","ENCSR454YDZ-1")
    sample_n <- 1:10
}
if (cell_type == "left_colon") { # left colon
    sample_IDs <- c("ENCSR830FPR-1","ENCSR916RYB-1","ENCSR904WIW-1")
    sample_n <- 1:3
}
if (cell_type == "leg_muscle") { # leg_muscle (gastrocnemius medialis)
    sample_IDs <- c("ENCSR696YOC-1","ENCSR819EGE-1","ENCSR244GZL-1","ENCSR023FME-1","ENCSR139TIQ-1")
    sample_n <- 1:5
}
if (cell_type == "leg_skin") { # leg_skin
    sample_IDs <- c("ENCSR733SZL-1","ENCSR397ODX-1","ENCSR474TGL-1","ENCSR513HZN-1")
    sample_n <- 1:4
}
if (cell_type == "liver") { # liver
    sample_IDs <- c("ENCSR594BJF-1","ENCSR659ANG-1","ENCSR074DOR-1","ENCSR497ZST-1","ENCSR825WYQ-1",
                    "ENCSR118PUH-1","ENCSR650BBI-1","ENCSR552QVU-1","ENCSR630YEA-1","ENCSR139ATR-1")
    sample_n <- 1:10
}
if (cell_type == "lung") { # lung
    sample_IDs <- c("ENCSR816NWE-1","ENCSR391BWM-1","ENCSR824OCY-1")
    sample_n <- 1:3
}
if (cell_type == "omental_fat") { # omental_fat
    sample_IDs <- c("ENCSR181XXQ-1","ENCSR492GGN-1","ENCSR644SCP-1","ENCSR274HQD-1")
    sample_n <- 1:4
}
if (cell_type == "pancreas") { # pancreas
    sample_IDs <- c("ENCSR868CRK-1","ENCSR690NZI-1","ENCSR229VVY-1","ENCSR496XXB-1")
    sample_n <- 1:4
}
if (cell_type == "peyer") { # peyer's patch
    sample_IDs <- c("ENCSR052DKH-1","ENCSR652WJE-1","ENCSR101JHK-1")
    sample_n <- 1:3
}
if (cell_type == "psoas_muscle") { # psoas muscle
    sample_IDs <- c("ENCSR000XQD-1","ENCSR869LEV-1","ENCSR916EDP-1","ENCSR332XEW-1")
    sample_n <- 1:4
}
if (cell_type == "stomach") { # stomach
    sample_IDs <- c("ENCSR848EID-1","ENCSR549JEQ-1","ENCSR266ZPZ-1","ENCSR940QRM-1")
    sample_n <- 1:4
}
if (cell_type == "thyroid") { # thyroid
    sample_IDs <- c("ENCSR796RXX-1","ENCSR909OXO-1","ENCSR817VFO-1")
    sample_n <- 1:3
}
if (cell_type == "tibial_nerve") { # tibial_nerve
    sample_IDs <- c("ENCSR205TUH-1","ENCSR453BVR-1","ENCSR726QTF-1")
    sample_n <- 1:3
}
if (cell_type == "transverse_colon") { # transverse_colon
    sample_IDs <- c("ENCSR434SXE-1","ENCSR997YNO-1","ENCSR506YMX-1","ENCSR007QIO-1","ENCSR349XKD-1")
    sample_n <- 1:5
}
if (cell_type == "uterus") { # uterus
    sample_IDs <- c("ENCSR828HVB-1","ENCSR455CVZ-1","ENCSR028RFK-1")
    sample_n <- 1:3
}
#if (cell_type == "ovary") { # ovary
#    sample_IDs <- c("ENCSR322LEY-1","ENCSR533CSY-1","ENCSR305VES-1","ENCSR105VKN-1","ENCSR197GLE-1")
#    sample_n <- 1:5
#}

# read counts matrix
counts <- ReadCounts(paste0("../data/",cell_type,"/",cell_type,sample_id,"/peak_matrix/"))
ca <- CreateChromatinAssay(
    counts = counts,
    annotation = annotations
)

# make fragment object
bc <- read.table(paste0("../data/",cell_type,"/",cell_type,sample_id,"/10000_barcodes_only.tsv"), sep = "\t")
frag_path <- paste0("../data/",cell_type,"/encode_scatac_dcc_2/results/",sample_IDs[sample_id],"/fragments/fragments.tsv.gz")
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
obj <- obj[, obj$nCount_ATAC > 1000]
obj

p1 <- DimPlot(obj, reduction = 'umap.atac')
p2 <- VlnPlot(obj, c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), pt.size=0)
p3 <- DensityScatter(obj, "nCount_ATAC", "TSS.enrichment", log_x = TRUE)

pdf(
  paste0("~/scratch/scATAC_celltypes/data/plots/check_QC/", cell_type, sample_id, "_QC.pdf"), 
  width = 6,      # Width of the PDF in inches
  height = 6       # Height of the PDF in inches
)
p1
p2
p3
dev.off()

# save object
saveRDS(obj, paste0("../data/",cell_type,"/",cell_type,sample_id,"/",cell_type,sample_id,"_preQC.rds"))