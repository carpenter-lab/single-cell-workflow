snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(Seurat)
library(dplyr)
library(glue)
library(rlang)

seurat <- readRDS(snakemake@input[["seurat"]])

labels <- snakemake@params[["cluster_labs"]]
clustering <- snakemake@params[["group_by"]]
new_labels <- snakemake@params[["new_group_by"]]

labels2 <- glue("'{from}' ~ '{from}|{to}'", from = names(labels), to = labels)

seurat[[new_labels]] <- seurat[[clustering]] |>
    mutate(
        "{new_labels}" := case_match(.data[[clustering]], !!!parse_exprs(labels2)),
        "{new_labels}" := forcats::fct(.data[[new_labels]], levels = glue("{from}|{to}", from = names(labels), to = labels))
    ) |>
    pull(all_of(new_labels))

saveRDS(seurat, snakemake@output[["seurat"]])