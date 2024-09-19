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

test_regex_matches <- function(vec, regex_list) {
    matches <- sapply(vec, function(element) {
        matched_regex <- sapply(regex_list, function(regex) {
            if (grepl(regex, element)) {
                return(regex)
            } else {
                return(NA)
            }
        })
    matched_regex <- matched_regex[!is.na(matched_regex)]
    if (length(matched_regex) > 0) {
        return(matched_regex[1])
    } else {
        return(NA)
    }
    })
    names(matches) <- names(vec)
    return(matches)
}

seurat[["tcr_motif"]] <- test_regex_matches(seurat$CTaa, tcr_patterns$pattern)

saveRDS(seurat, snakemake@output[["seurat"]])