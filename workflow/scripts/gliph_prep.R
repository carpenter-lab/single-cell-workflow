
snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(tidyverse)
library(glue)
library(SeuratObject)

seurat <- readRDS(snakemake@input[["seurat"]])

data <- seurat[[c("CTgene", "CTaa", "clonalFrequency", "patient_id", "condition")]]

data |>
    mutate(
        TRBV = str_extract(CTgene, "TRBV[0-9|-]+"),
        TRBJ = str_extract(CTgene, "TRBJ[0-9|-]+"),
        "subject:condition" = glue("{patient_id}:{condition}")
    ) |>
    separate(CTaa, into = c("CDR3a", "CDR3b")) |>
    select(CDR3b, TRBV, TRBJ, CDR3a, "subject:condition", count = clonalFrequency) |>
    filter(!is.na(CDR3b)) |>
    write_tsv(snakemake@output[["tsv"]])
