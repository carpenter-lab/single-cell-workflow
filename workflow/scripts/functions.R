SinkAllOutput <- function(log, .envir = parent.frame()) {
    if (is.null(log) || length(log) < 1) return(invisible(NULL))

    if (is.list(log)) log <- log[[1]]
    assign("sink_connection", file(log), envir = .envir)
    sink(sink_connection, append = TRUE)
    sink(sink_connection, append = TRUE, type = "message")
    return(invisible(NULL))
}