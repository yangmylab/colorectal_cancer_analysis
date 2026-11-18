#!/usr/bin/env Rscript
# compose_figures.R
# Assemble the 44 published panels into 6 main + 2 supplementary composite figures, driven by
# 12_Manuscript/figures_plan.tsv (the single source of truth for panel order).
#
# Typography policy (rewritten 2026-11-18 after measuring the panels):
#   * panels are placed as VECTOR PDF pages through pdflatex + TikZ, so panel text stays selectable
#     and editable; the previous version rescaled every panel to a fixed 900 px width and shipped a
#     raster montage inside a PDF.
#   * one typographic scale per figure: every panel is multiplied by the same factor, so an 11 pt
#     label keeps its size relative to every other panel (`--mode=uniform`, the default). The old
#     fixed-pixel rescaling spread panel type size by 1.14x-1.67x within a single figure.
#   * `--mode=equal` keeps equal-width cells instead (fills the page better, but panel type then
#     varies with each panel's authored width). It exists for comparison, not as the default.
#   * panel letters are set at a fixed 8 pt Helvetica bold in a reserved strip above each panel, so
#     they are identical across all eight figures and never overlay panel content.
#   * every figure is written at a fixed printed width (180 mm) under a height budget (240 mm); the
#     resulting printed type size of every panel is audited against a floor (7 pt) in
#     12_Manuscript/compose_figures_typography.{tsv,md}.
#
# Honest limitation: the panels were authored at 6-10 in with ~11 pt type. At 180 mm that type prints
# at 3-5.5 pt, i.e. below the floor, and no montage can repair it. The audit names the panels that
# must be re-authored at print width (~3.4 in, base_size ~8) for a compliant figure set.
#
# Parallelism: 6 workers (6 physical cores; nproc reports 12 SMT threads), one level only, with
# OMP/OPENBLAS/MKL pinned to 1 thread per worker.

library(data.table)
library(pdftools)

PLAN <- "12_Manuscript/figures_plan.tsv"
SRC <- "09_Figures/plots"
OUT <- "12_Manuscript/Figures"
AUDIT_TSV <- "12_Manuscript/compose_figures_typography.tsv"
AUDIT_MD <- "12_Manuscript/compose_figures_typography.md"
REF_DOC <- Sys.getenv("CRC_REF_DOC", "/run/media/meng/Seqdata/MN/Paper_1/reference/tnr_fangsong_redcode.docx")
WORKERS <- 6L
ASSUMED_BODY_PT <- 11
ASSUMED_TITLE_PT <- 13.2
MM_PER_IN <- 25.4

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}
MODE <- arg("mode", "uniform")
TARGET_W <- as.numeric(arg("width_mm", 180)) / MM_PER_IN
MAX_H <- as.numeric(arg("height_mm", 240)) / MM_PER_IN
FLOOR_PT <- as.numeric(arg("floor_pt", 7))
GUTTER <- as.numeric(arg("gutter_mm", 4)) / MM_PER_IN
LETTER_STRIP <- as.numeric(arg("letter_strip_mm", 4)) / MM_PER_IN
LETTER_PT <- as.numeric(arg("letter_pt", 8))
stopifnot(MODE %in% c("uniform", "equal"), TARGET_W > 0, MAX_H > 0, GUTTER >= 0, LETTER_STRIP > 0)

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
for (tool in c("pdflatex", "pdftoppm")) {
  if (!nzchar(Sys.which(tool))) {
    stop("compose_figures.R needs ", tool,
         " for the vector montage; install texlive-latex-recommended (TikZ) and poppler-utils",
         call. = FALSE)
  }
}

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
OUT_ABS <- normalizePath(OUT, mustWork = TRUE)

plan <- fread(PLAN, data.table = FALSE)
stopifnot(all(c("figure", "panel", "source") %in% names(plan)))
stopifnot(!any(duplicated(paste(plan$figure, plan$panel))))
panels <- data.frame(
  figure = plan$figure,
  panel = plan$panel,
  plot = plan$source,
  src = file.path(SRC, paste0(plan$source, ".pdf")),
  stringsAsFactors = FALSE
)
missing <- unique(panels$src[!file.exists(panels$src)])
if (length(missing)) stop("missing panel PDFs: ", paste(missing, collapse = ", "), call. = FALSE)
sizes <- t(vapply(panels$src, function(p) {
  s <- pdftools::pdf_pagesize(p)
  c(w = s$width[1] / 72, h = s$height[1] / 72)
}, numeric(2)))
panels$w <- sizes[, "w"]
panels$h <- sizes[, "h"]

# One candidate layout per column count. In uniform mode the column widths follow the authored panel
# widths, so a single factor scales every panel and the printed type size is identical across panels.
layout_one <- function(w, h, ncol) {
  n <- length(w)
  nrow <- ceiling(n / ncol)
  col <- ((seq_len(n) - 1L) %% ncol) + 1L
  row <- ((seq_len(n) - 1L) %/% ncol) + 1L
  if (MODE == "uniform") {
    col_w_auth <- vapply(seq_len(ncol), function(c) max(w[col == c]), numeric(1))
    row_h_auth <- vapply(seq_len(nrow), function(r) max(h[row == r]), numeric(1))
    scale <- rep((TARGET_W - GUTTER * (ncol - 1)) / sum(col_w_auth), n)
    col_w <- scale[1] * col_w_auth
    row_h <- scale[1] * row_h_auth + LETTER_STRIP
  } else {
    col_w <- rep((TARGET_W - GUTTER * (ncol - 1)) / ncol, ncol)
    scale <- col_w[1] / w
    row_h <- vapply(seq_len(nrow), function(r) max(h[row == r] * scale[row == r]) + LETTER_STRIP,
                    numeric(1))
  }
  list(ncol = ncol, nrow = nrow, col = col, row = row, scale = scale,
       col_x = c(0, cumsum(col_w + GUTTER))[seq_len(ncol)],
       row_y = c(0, cumsum(row_h + GUTTER))[seq_len(nrow)],
       panel_w = w * scale, panel_h = h * scale,
       canvas_w = TARGET_W, canvas_h = sum(row_h) + GUTTER * (nrow - 1),
       printed_pt = ASSUMED_BODY_PT * scale)
}

# Fewest columns that still fit the height budget; otherwise the tallest-but-smallest layout.
choose_layout <- function(w, h) {
  cand <- lapply(seq_len(min(3L, length(w))), function(k) layout_one(w, h, k))
  ok <- Filter(function(L) L$canvas_h <= MAX_H + 1e-9, cand)
  if (!length(ok)) ok <- cand
  best <- ok[[which.max(vapply(ok, function(L) min(L$printed_pt), numeric(1)))]]
  best$within_height <- best$canvas_h <= MAX_H + 1e-9
  best
}

tex_for <- function(L, files, letters) {
  lines <- c(
    "\\documentclass{article}",
    sprintf("\\usepackage[paperwidth=%.4fin,paperheight=%.4fin,margin=0pt]{geometry}",
            L$canvas_w, L$canvas_h),
    "\\usepackage{graphicx}",
    "\\usepackage{tikz}",
    "\\usepackage{helvet}",
    "\\renewcommand{\\familydefault}{\\sfdefault}",
    "\\pagestyle{empty}",
    "\\setlength{\\parindent}{0pt}",
    "\\setlength{\\topskip}{0pt}",
    "\\begin{document}\\null",
    "\\begin{tikzpicture}[remember picture,overlay,x=1in,y=1in]"
  )
  for (i in seq_along(files)) {
    x <- L$col_x[L$col[i]]
    y <- L$row_y[L$row[i]] + LETTER_STRIP
    lines <- c(
      lines,
      sprintf(paste0("\\node[anchor=north west,inner sep=0] at ",
                     "([xshift=%.4fin,yshift=-%.4fin]current page.north west) ",
                     "{\\includegraphics[width=%.4fin,height=%.4fin]{%s}};"),
              x, y, L$panel_w[i], L$panel_h[i], files[i]),
      sprintf(paste0("\\node[anchor=south west,inner sep=0] at ",
                     "([xshift=%.4fin,yshift=-%.4fin]current page.north west) ",
                     "{\\fontsize{%.1f}{%.1f}\\selectfont\\bfseries %s};"),
              x, y, LETTER_PT, LETTER_PT * 1.2, letters[i])
    )
  }
  c(lines, "\\end{tikzpicture}", "\\end{document}")
}

build_figure <- function(fig) {
  sub <- panels[panels$figure == fig, , drop = FALSE]
  sub <- sub[order(match(sub$panel, LETTERS)), , drop = FALSE]
  if (!identical(sub$panel, LETTERS[seq_len(nrow(sub))])) {
    stop("panels for ", fig, " are not a contiguous A.. sequence", call. = FALSE)
  }
  L <- choose_layout(sub$w, sub$h)
  dir <- file.path(tempdir(), paste0("crcfig_", gsub("[^A-Za-z0-9]", "_", fig)))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  staged <- paste0("panel_", sub$panel, ".pdf")
  if (!all(file.copy(sub$src, file.path(dir, staged), overwrite = TRUE))) {
    stop("could not stage panels for ", fig, call. = FALSE)
  }
  writeLines(tex_for(L, sub(".pdf$", "", staged), sub$panel), file.path(dir, "composite.tex"))

  old_wd <- setwd(dir)
  on.exit(setwd(old_wd), add = TRUE)
  for (pass in seq_len(2L)) {
    log <- suppressWarnings(system2("pdflatex",
                                    c("-interaction=nonstopmode", "-halt-on-error", "composite.tex"),
                                    stdout = TRUE, stderr = TRUE))
    status <- attr(log, "status")
    if (!is.null(status) && status != 0L) {
      stop("pdflatex failed for ", fig, " (pass ", pass, "):\n",
           paste(tail(log, 20), collapse = "\n"), call. = FALSE)
    }
  }
  if (!file.exists("composite.pdf")) stop("pdflatex produced no PDF for ", fig, call. = FALSE)

  base <- file.path(OUT_ABS, gsub(" ", "_", fig))
  file.copy("composite.pdf", paste0(base, ".pdf"), overwrite = TRUE)
  suppressWarnings(system2("pdftoppm",
                           c("-jpeg", "-singlefile", "-r", "300", "-jpegopt", "quality=95",
                             "composite.pdf", "render"),
                           stdout = TRUE, stderr = TRUE))
  if (!file.exists("render.jpg")) stop("pdftoppm produced no JPEG for ", fig, call. = FALSE)
  file.copy("render.jpg", paste0(base, ".jpg"), overwrite = TRUE)

  got <- pdftools::pdf_pagesize(paste0(base, ".pdf"))
  if (abs(got$width[1] / 72 - L$canvas_w) > 0.02 || abs(got$height[1] / 72 - L$canvas_h) > 0.02) {
    stop("composite size check failed for ", fig, call. = FALSE)
  }
  cat(sprintf("%-24s %d panels, %d cols, canvas %.0f x %.0f mm, body type %.2f pt%s\n",
              fig, nrow(sub), L$ncol, L$canvas_w * MM_PER_IN, L$canvas_h * MM_PER_IN,
              min(L$printed_pt), if (L$within_height) "" else "  [OVER HEIGHT BUDGET]"))
  data.frame(
    figure = fig,
    panel = sub$panel,
    plot = sub$plot,
    mode = MODE,
    cols = L$ncol,
    authored_in = sprintf("%.2f x %.2f", sub$w, sub$h),
    printed_in = sprintf("%.2f x %.2f", L$panel_w, L$panel_h),
    scale = round(L$scale, 3),
    printed_body_pt = round(L$printed_pt, 2),
    printed_title_pt = round(L$printed_pt * (ASSUMED_TITLE_PT / ASSUMED_BODY_PT), 2),
    floor_pt = FLOOR_PT,
    meets_floor = L$printed_pt >= FLOOR_PT - 1e-9,
    canvas_mm = sprintf("%.0f x %.0f", L$canvas_w * MM_PER_IN, L$canvas_h * MM_PER_IN),
    height_ok = L$within_height,
    stringsAsFactors = FALSE
  )
}

figures <- unique(panels$figure)
cat("panels in plan:", nrow(panels), "| figures:", length(figures),
    "| mode:", MODE, "| width_mm:", TARGET_W * MM_PER_IN, "| floor_pt:", FLOOR_PT, "\n")
res <- parallel::mclapply(figures, function(f) try(build_figure(f), silent = TRUE),
                          mc.cores = min(WORKERS, length(figures)), mc.preschedule = FALSE)
failed <- vapply(res, inherits, logical(1), what = "try-error")
if (any(failed)) {
  for (i in which(failed)) cat(figures[i], ":\n", as.character(res[[i]]), "\n")
  stop("composite build failed for: ", paste(figures[failed], collapse = ", "), call. = FALSE)
}
audit <- do.call(rbind, res)

write.table(audit, AUDIT_TSV, sep = "\t", quote = FALSE, row.names = FALSE)

below <- audit[!audit$meets_floor, , drop = FALSE]
fig_rows <- unique(audit[, c("figure", "cols", "canvas_mm", "height_ok")])
fig_rows$min_body_pt <- vapply(fig_rows$figure, function(f) min(audit$printed_body_pt[audit$figure == f]), numeric(1))
fig_rows$max_body_pt <- vapply(fig_rows$figure, function(f) max(audit$printed_body_pt[audit$figure == f]), numeric(1))
fig_rows$panels_below <- vapply(fig_rows$figure, function(f) sum(!audit$meets_floor[audit$figure == f]), integer(1))

md <- c(
  "# Composite figure typography audit",
  "",
  sprintf("Generated by `Scripts/compose_figures.R` (%s) on %s.", MODE, format(Sys.Date())),
  "",
  sprintf("- Printed width: **%.0f mm**; height budget: **%.0f mm**; legibility floor: **%.1f pt**.",
          TARGET_W * MM_PER_IN, MAX_H * MM_PER_IN, FLOOR_PT),
  sprintf("- Panel panels are placed as vector PDF pages; type sizes below are derived from each panel's authored width (%.1f pt body, %.1f pt title assumed).",
          ASSUMED_BODY_PT, ASSUMED_TITLE_PT),
  sprintf("- Panels below the floor: **%d of %d**.", nrow(below), nrow(audit)),
  "",
  "## Per figure",
  "",
  "| Figure | Panels | Columns | Canvas (mm) | Body type (pt) | Below floor | Height OK |",
  "|---|---|---|---|---|---|---|"
)
for (i in seq_len(nrow(fig_rows))) {
  r <- fig_rows[i, ]
  md <- c(md, sprintf("| %s | %d | %d | %s | %.2f-%.2f | %d | %s |",
                      r$figure, sum(audit$figure == r$figure), r$cols, r$canvas_mm,
                      r$min_body_pt, r$max_body_pt, r$panels_below,
                      if (r$height_ok) "yes" else "**no**"))
}
md <- c(md, "", "## Panels below the legibility floor", "")
if (nrow(below)) {
  md <- c(md, "| Figure | Panel | Plot | Authored (in) | Printed (in) | Body (pt) | Title (pt) |",
          "|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(below))) {
    r <- below[i, ]
    md <- c(md, sprintf("| %s | %s | %s | %s | %s | %.2f | %.2f |",
                        r$figure, r$panel, r$plot, r$authored_in, r$printed_in,
                        r$printed_body_pt, r$printed_title_pt))
  }
} else {
  md <- c(md, "None.")
}
md <- c(md, "",
  "## What this means",
  "",
  "The panels were authored at 6-10 in with ~11 pt type, so at 180 mm every panel prints its text",
  "between roughly 3 and 5.5 pt. `--mode=uniform` makes that size identical across panels (the old",
  "fixed-pixel montage let it vary by 1.14x-1.67x within a figure); it cannot make it larger.",
  "Reaching the floor requires re-authoring the panels at print width (~3.4 in with `base_size = 8`),",
  "after which the same montage yields 7-8 pt. Regenerating panels means re-running the layer scripts,",
  "so it is tracked separately from this composer.")
writeLines(md, AUDIT_MD)
if (nzchar(Sys.which("pandoc")) && file.exists(REF_DOC)) {
  invisible(suppressWarnings(system2("pandoc",
                                     c(shQuote(AUDIT_MD), paste0("--reference-doc=", shQuote(REF_DOC)),
                                       "-o", shQuote(sub("\\.md$", ".docx", AUDIT_MD))),
                                     stdout = TRUE, stderr = TRUE)))
}

source("06_Analysis/1_R/1_Rscripts/utils.R")
log_decision("compose_figures", "PASS",
             sprintf("mode=%s; width=%.0fmm; floor=%.1fpt; figures=%d; panels=%d; below_floor=%d",
                     MODE, TARGET_W * MM_PER_IN, FLOOR_PT, length(figures), nrow(audit), nrow(below)))
cat(sprintf("done: %d figures, %d panels, %d below the %.1f pt floor (see %s)\n",
            length(figures), nrow(audit), nrow(below), FLOOR_PT, AUDIT_TSV))
