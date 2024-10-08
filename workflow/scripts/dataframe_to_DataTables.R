snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(DT)
library(data.table)

data <- fread(snakemake@input[["table"]])

dt <- data |>
    datatable(
        rownames = FALSE,
        filter = "top",
        options = list(pageLength = 15, lengthMenu = c(5, 10, 15, 20))
    )

dt <- tryCatch(
    dt |>
        formatPercentage(c("pct.1", "pct.2"), 1) |>
        formatSignif(c("p_val", "p_val_adj")) |>
        formatRound("avg_log2FC", digits = 3),
    error = function(e) dt
)

dt |>
    saveWidget(
        snakemake@output[["html"]],
        libdir = fs::path_rel(
            path = snakemake@output[["html_library"]],
            start = dirname(snakemake@output[["html"]])
        )
    )

