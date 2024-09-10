snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(Seurat)
library(stringr)
library(dplyr)
library(purrr)
library(tidyr)
library(data.table)

CollapseTCRPatternString <- function(v) unique(v) |> paste(collapse = "-") |> str_remove("-?NA-?")

DetectTCRMotif <- function(motif, pattern, data, ...) if_else(str_detect(data, pattern), motif, NA_character_)

seurat <- readRDS(snakemake@input[["seurat"]])

tcr_patterns <- snakemake@params[["tcr_patterns"]] |>
    as_tibble() |>
    pivot_longer(cols = everything(), names_to = "name", values_to = "motif") |>
    distinct() |>
    mutate(pattern = str_replace(motif, "%", "."))

cols_to_keep <- c("tcr_motif", "tcr_name")

seurat[["cell_id"]] <- rownames(seurat[[]])
tcr <- copy(seurat[[]])

setDT(tcr)


tcr_patterns[["motif"]] -> motif_use

tcr[ , c(motif_use) := pmap(tcr_patterns, DetectTCRMotif, data = CTaa)
    ][ , tcr_motif := .(apply(.SD, 1, unique)), .SDcols = tcr_patterns[["motif"]]
    ][ , tcr_motif := lapply(tcr_motif, discard, is.na)
    ][ , tcr_name := lapply(tcr_motif, function(x) tcr_patterns[["name"]][tcr_patterns[["motif"]] %chin% unlist(x)])
    ][ , c("tcr_motif", "tcr_name") := lapply(.SD, function(x) lapply(x, CollapseTCRPatternString)), .SDcols = c("tcr_motif", "tcr_name")
    ][ , tcr_patterns[["motif"]] := NULL
]

setDF(tcr)

rownames(tcr) <- tcr[["cell_id"]]

tcr <- tcr[ , c("tcr_motif", "tcr_name")]

seurat <- AddMetaData(seurat, tcr)

saveRDS(seurat, snakemake@output[["seurat"]])