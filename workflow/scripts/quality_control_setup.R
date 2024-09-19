snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(Seurat)
library(scRepertoire)

seurat <- readRDS(snakemake@input[["seurat"]])

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
seurat[["percent.tcr"]] <- PercentageFeatureSet(seurat, pattern = "^TR[BD][VDJ]*|^TR[AG][VJ]*")
seurat[["clonal_expansion_type"]] <- data.table::fcase(
    seurat[["Frequency"]] == 0, "Absent",
    seurat[["Frequency"]] < 2, "Not Expanded",
    seurat[["Frequency"]] == 2, "Undetermined",
    seurat[["Frequency"]] >= 3, "Expanded"
)

saveRDS(seurat, file = snakemake@output[["seurat"]])