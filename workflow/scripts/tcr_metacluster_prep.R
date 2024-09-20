snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(tidyverse)
library(glue)

cell_ranger_tcr_spec <- cols_only(
    barcode = col_character(),
    chain = col_character(),
    v_gene = col_character(),
    d_gene = col_character(),
    j_gene = col_character(),
    c_gene = col_character(),
    full_length = col_logical(),
    productive = col_logical(),
    cdr1 = col_character(),
    cdr2 = col_character(),
    cdr3 = col_character(),
    cdr3_nt = col_character()
)

input_matrix <- list(
    file = snakemake@input[["tcr"]],
    patient_id = snakemake@params[["patient_id"]],
    condition = snakemake@params[["condition"]]
)

tcr_data <- pmap(
    input_matrix,
    function(file, patient_id, condition) {
        read_csv(file, col_types = cell_ranger_tcr_spec) |>
            mutate(patient_id = patient_id, condition = condition)
    }
) |>
    bind_rows() |>
    filter(chain != "Multi", productive, cdr3 != "None") |>
    mutate(vdj_cdr3 = paste0(v_gene, "|", d_gene,"|", j_gene, "|", cdr3)) |>
    pivot_wider(id_cols = c(barcode, patient_id, condition), names_from = chain, values_from = vdj_cdr3) |>
    unnest(c(TRA, TRB)) |>
    separate_wider_delim(cols = c(TRA, TRB), delim = "|", names = c("V", "D", "J", "CDR3"), names_sep = "") |>
    mutate(across(starts_with("TR"), ~ if_else(.x == "NA", NA_character_, .x))) |>
    group_by(pick(starts_with("TR"))) |>
    mutate(count = n()) |>
    ungroup()

data_clean <- if (snakemake@params[["tool"]] == "gliph") {
    tcr_data |>
        mutate("subject:condition" = glue("{patient_id}:{condition}")) |>
        select(CDR3b = TRBCDR3, TRBV, TRBJ, CDR3a = TRACDR3, "subject:condition", count) |>
        filter(!is.na(CDR3b))
} else if (snakemake@params[["tool"]] == "tcrdist") {
    tcr_data |>
        mutate(across(.cols = starts_with("TR"), .fns = ~ paste0(.x, "*01"))) |>
        select(
            cdr3_b_aa = TRBCDR3, cdr3_a_aa = TRACDR3,
            v_b_gene = TRBV, j_b_gene = TRBJ,
            v_a_gene = TRAV, j_a_gene = TRAJ,
            subject = patient_id,
            condition,
            count
        )
}

outputs <- snakemake@output
outputs[names(outputs) == ""] <- NULL

names(outputs) <- case_match(
    names(outputs),
    "csv" ~ "write_csv",
    "tsv" ~ "write_tsv"
)

iwalk(outputs, function(path, fun) {
    fun <- utils::getFromNamespace(fun, ns = "readr")
    fun(data_clean, path)
})