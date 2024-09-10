
library(Seurat)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

seurat <- readRDS(snakemake@input[["seurat"]])

data <- seurat[[c("CTgene", "CTaa", "clonalFrequency", "patient_id", "condition")]]

data |>
    as_tibble() |>
    mutate(
        TRBV = str_extract(CTgene, "TRBV[0-9|-]+"),
        TRBJ = str_extract(CTgene, "TRBJ[0-9|-]+"),
        TRAJ = str_extract(CTgene, "TRAJ[0-9|-]+"),
        TRAV = str_extract(CTgene, "TRAV[0-9|-]+"),
        across(.cols = c(TRBV, TRBJ, TRAV, TRAJ), .fns = ~ paste0(.x, "*01"))
    ) |>
    separate(CTaa, into = c("CDR3a", "CDR3b"))|>
    select(
        cdr3_b_aa = CDR3b, cdr3_a_aa = CDR3a,
        v_b_gene=TRBV, j_b_gene=TRBJ,
        v_a_gene = TRAV, j_a_gene = TRAJ,
        subject = patient_id,
        condition,
        count = clonalFrequency
    ) |>
    filter(!is.na(cdr3_b_aa)) |>
    write_csv(snakemake@output[["csv"]])
