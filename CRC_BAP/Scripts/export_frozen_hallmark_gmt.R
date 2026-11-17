# export_frozen_hallmark_gmt.R
# msigdbr 26.1.1 needs the `msigdbdf` data package or an online release download; `msigdbdf` is not on
# any reachable mirror here (tuna CRAN 404, CRAN/r-universe unreachable) and the zenodo release times out.
# So the 50 Hallmark sets used by the discovery run are exported once to a plain GMT (Ensembl IDs, the
# same id space as 03_Raw_Data/expression.tsv). Pipeline steps can then read the GMT offline and the
# analysis stays reproducible without network access.

suppressPackageStartupMessages(library(data.table))
sig <- readRDS("04_Processed_Data/signature_scores.rds")
gs <- sig$gene_sets
stopifnot(length(gs) > 0)
cat("sets:", length(gs), "| id type:",
    if (all(grepl("^ENSG", unlist(gs)[1:50]))) "ensembl" else "symbol", "\n")

lines <- vapply(names(gs), function(nm) {
  paste(c(nm, "frozen_discovery_run", unique(gs[[nm]])), collapse = "\t")
}, character(1))
out <- "Environment/hallmark_frozen.gmt"
writeLines(unname(lines), out)
cat("wrote", out, "|", length(lines), "sets |",
    sum(lengths(gs)), "gene entries |", round(file.size(out) / 1024, 1), "KB\n")
cat("sizes: min", min(lengths(gs)), "median", median(lengths(gs)), "max", max(lengths(gs)), "\n")
