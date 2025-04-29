library(Signac)
library(Seurat)
library(Matrix)
library(GenomicRanges)
source("code/utilities.R")
library(tidyr)

# load pseudobulk matrix
adrenal <- readRDS(snakemake@input[['adrenal']])
esophagus <- readRDS(snakemake@input[['esophagus']])
heartRV <- readRDS(snakemake@input[['heartRV']])
heart_fetal <- readRDS(snakemake@input[['heart_fetal']])
left_colon <- readRDS(snakemake@input[['left_colon']])
liver <- readRDS(snakemake@input[['liver']])
psoas_muscle <- readRDS(snakemake@input[['psoas_muscle']])

# combine matrices
combined_pseudobulk <- rbind(adrenal, esophagus, heart_fetal, heartRV, 
                             left_colon, liver, psoas_muscle)

# split pseudobulk per chromosome
ccre <- read.table("data/combined_cre.bed", sep="\t")
ccre_gr <- GRanges(seqnames = ccre$V1,
                   ranges = IRanges(start = ccre$V2, end = ccre$V3))
ccre_id <- ccre$V4

colnames(combined_pseudobulk) <- ccre_id

saveRDS(combined_pseudobulk, snakemake@output[['tissues_pseudobulk']])

# split matrix per chromosome
tissue_names <- c("adrenal", "esophagus", "heartRV", "heart_fetal", "left_colon", "liver", "psoas_muscle")
chromosome_names <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                     "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
                     "chr21", "chr22", "chrX", "chrY")

for (l in 1:length(tissue_names)) {
    tissue <- tissue_names[l]
    message(tissue)
    pseudobulk_matrix <- readRDS(paste0("data/pseudobulk/",tissue,"_pseudobulk.rds"))
    colnames(pseudobulk_matrix) <- ccre_id
    for (i in 1:length(chromosome_names)) {
        chr <- chromosome_names[i]
        message(chr)
        chr_start <- head(which(ccre$V1 == chr), 1)
        chr_end <- tail(which(ccre$V1 == chr), 1)

        matrix_split <- pseudobulk_matrix[,chr_start:chr_end]
        saveRDS(matrix_split, file = paste0("data/",tissue,"/",tissue,"_",chr,"_pseudobulk.rds"))
    }
}