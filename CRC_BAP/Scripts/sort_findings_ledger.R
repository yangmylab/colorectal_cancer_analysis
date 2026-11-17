# Keep the findings ledger in F01, F02, ... order regardless of the order in which layer scripts ran.
suppressPackageStartupMessages(library(data.table))
p <- "13_Reports/Findings_Ledger.csv"
led <- fread(p, data.table = FALSE)
led <- led[order(as.integer(sub("^F", "", led$finding_id))), , drop = FALSE]
fwrite(led, p, sep = ",")
cat("ids in order:", paste(led$finding_id, collapse = ","), "\n")
cat("rows:", nrow(led), "| cols:", ncol(led), "| duplicate ids:", any(duplicated(led$finding_id)), "\n")
