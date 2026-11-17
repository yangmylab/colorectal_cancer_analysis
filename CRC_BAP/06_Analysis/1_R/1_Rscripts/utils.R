# Shared helpers for the CRC Track A pipeline.
# Worker count is set by run_pipeline.R from config$params$n_workers (one level of parallelism).

read_matrix_any <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
    if (is.matrix(obj)) return(obj)
    if (is.data.frame(obj)) {
      m <- as.matrix(obj[, -1, drop = FALSE])
      rownames(m) <- obj[[1]]
      storage.mode(m) <- "numeric"
      return(m)
    }
    if (requireNamespace("SummarizedExperiment", quietly = TRUE) && methods::is(obj, "SummarizedExperiment")) {
      return(SummarizedExperiment::assay(obj))
    }
    stop("Unsupported RDS object in ", path)
  }
  dt <- data.table::fread(path, data.table = FALSE)
  if (is.character(dt[[1]]) || is.factor(dt[[1]])) {
    m <- as.matrix(dt[, -1, drop = FALSE])
    rownames(m) <- dt[[1]]
  } else {
    m <- as.matrix(dt)
  }
  storage.mode(m) <- "numeric"
  m
}

save_plot <- function(plot_fun, file_base, width = 7, height = 5) {
  grDevices::pdf(paste0(file_base, ".pdf"), width = width, height = height)
  plot_fun()
  grDevices::dev.off()
  grDevices::jpeg(paste0(file_base, ".jpg"), width = width, height = height, units = "in", res = 300, quality = 95)
  plot_fun()
  grDevices::dev.off()
}

zscore_rows <- function(m) {
  t(scale(t(m)))
}

top_variable <- function(m, n = 1000) {
  v <- apply(m, 1, stats::var, na.rm = TRUE)
  v[is.na(v)] <- 0
  keep <- order(v, decreasing = TRUE)[seq_len(min(n, length(v)))]
  m[keep, , drop = FALSE]
}

ari <- function(a, b) {
  mclust::adjustedRandIndex(a, b)
}

log_decision <- function(step, decision, rationale) {
  path <- get0("DECISION_LOG_PATH", envir = .GlobalEnv, ifnotfound = file.path("Logs", "DECISION_LOG.tsv"))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  row <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    step = step,
    decision = decision,
    rationale = rationale,
    stringsAsFactors = FALSE
  )
  if (file.exists(path)) {
    old <- data.table::fread(path, data.table = FALSE)
    row <- rbind(old, row)
  }
  data.table::fwrite(row, path, sep = "\t")
}

require_file <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path, call. = FALSE)
  invisible(TRUE)
}
