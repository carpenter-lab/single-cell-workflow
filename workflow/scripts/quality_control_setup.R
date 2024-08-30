snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(Seurat)
library(scRepertoire)
library(Trex)

seurat <- readRDS(snakemake@input[["seurat"]])

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
seurat[["percent.tcr"]] <- PercentageFeatureSet(seurat, pattern = "^TR[BD][VDJ]*|^TR[AG][VJ]*")
seurat[["clonal_expansion_type"]] <- data.table::fcase(
    seurat[["clonalFrequency"]] == 0, "Absent",
    seurat[["clonalFrequency"]] < 2, "Not Expanded",
    seurat[["clonalFrequency"]] == 2, "Undetermined",
    seurat[["clonalFrequency"]] >= 3, "Expanded"
)

saveRDS(seurat, file = snakemake@output[["seurat"]])