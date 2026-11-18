# Shared helpers for the CRC Track A pipeline.
# Worker count is set by run_pipeline.R from config$params$n_workers (one level of parallelism).

# One palette for the two ecotypes in every figure (Okabe-Ito blue / vermillion: colour-blind safe).
ECOTYPE_COL <- c(E1 = "#0072B2", E2 = "#D55E00")

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
  # The composites are laid out at a fixed printed width, so a panel's authored size decides its printed
  # type size: a panel authored at 3.4 in prints its 11 pt text at >=7 pt in a two-column figure.
  # CRC_PANEL_MAX_W caps the authored width (inches; 0 = off) and CRC_PANEL_SCALE multiplies both
  # dimensions, so a re-authored panel set needs no edits to the individual plot calls.
  scale <- as.numeric(Sys.getenv("CRC_PANEL_SCALE", "1"))
  max_w <- as.numeric(Sys.getenv("CRC_PANEL_MAX_W", "0"))
  width <- width * scale
  height <- height * scale
  if (max_w > 0 && width > max_w) {
    shrink <- max_w / width
    width <- width * shrink
    height <- height * shrink
  }
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

# Hallmark gene sets, offline-safe.
# msigdbr 26.x needs the `msigdbdf` data package or an online release download; msigdbdf is on no
# reachable mirror here and the zenodo release times out, so fall back to the frozen GMT exported from
# the discovery run (`03_Raw_Data/hallmark_frozen.gmt`, 50 sets, Ensembl IDs).
read_gmt_sets <- function(path) {
  lines <- readLines(path)
  sets <- list()
  for (l in lines) {
    f <- strsplit(l, "\t")[[1]]
    if (length(f) >= 3) sets[[f[1]]] <- unique(f[-c(1, 2)])
  }
  sets
}

map_ensembl_to_symbol <- function(ens) {
  map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = unique(ens),
                               keytype = "ENSEMBL", columns = "SYMBOL")
  map <- map[!is.na(map$SYMBOL) & map$SYMBOL != "", , drop = FALSE]
  map <- map[!duplicated(map$ENSEMBL), , drop = FALSE]
  stats::setNames(map$SYMBOL, map$ENSEMBL)
}

map_symbol_to_ensembl <- function(sym) {
  map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = unique(sym),
                               keytype = "SYMBOL", columns = "ENSEMBL")
  map <- map[!is.na(map$ENSEMBL) & map$ENSEMBL != "", , drop = FALSE]
  map <- map[!duplicated(map$SYMBOL), , drop = FALSE]
  stats::setNames(map$ENSEMBL, map$SYMBOL)
}

hallmark_gene_sets <- function(id_type = c("ensembl", "symbol"),
                               gmt = "Environment/hallmark_frozen.gmt") {
  id_type <- match.arg(id_type)
  # Only call msigdbr when its data are actually present locally (msigdbdf package or a populated cache);
  # otherwise every call would sit in a 10 s network timeout before falling back.
  cache_dir <- tryCatch(tools::R_user_dir("msigdbr", "cache"), error = function(e) NA_character_)
  have_data <- requireNamespace("msigdbdf", quietly = TRUE) ||
    (!is.na(cache_dir) && dir.exists(cache_dir) && length(list.files(cache_dir, recursive = TRUE)) > 0)
  from_msigdbr <- if (have_data) tryCatch({
    msig <- msigdbr::msigdbr(db_species = "HS", collection = "H")
    key <- if (id_type == "ensembl" && "ensembl_gene" %in% colnames(msig)) "ensembl_gene" else "gene_symbol"
    list(sets = split(msig[[key]], msig$gs_name),
         source = paste0("MSigDB Hallmark (msigdbr; ", key, ")"), fallback = FALSE)
  }, error = function(e) NULL) else NULL
  if (!is.null(from_msigdbr)) return(from_msigdbr)
  if (!file.exists(gmt)) {
    stop("msigdbr is unavailable (no msigdbdf / no network) and no frozen GMT at ", gmt, call. = FALSE)
  }
  sets <- read_gmt_sets(gmt)
  is_ens <- grepl("^ENSG", unique(unlist(sets))[1])
  if (id_type == "symbol" && is_ens) {
    m <- map_ensembl_to_symbol(unlist(sets))
    sets <- lapply(sets, function(g) unique(unname(m[g[g %in% names(m)]])))
  } else if (id_type == "ensembl" && !is_ens) {
    stop("Frozen GMT holds symbols but Ensembl IDs are required", call. = FALSE)
  }
  list(sets = sets, source = paste0("frozen GMT (", gmt, "; ", id_type, ")"), fallback = TRUE)
}
