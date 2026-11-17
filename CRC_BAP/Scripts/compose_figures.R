#!/usr/bin/env Rscript
# compose_figures.R
# Assemble the 44 published panels into 6 main + 2 supplementary composite figures, driven by
# 12_Manuscript/figures_plan.tsv (the single source of truth for panel order).
#
# Output: 09_Figures/composite/Figure_1.{pdf,jpg} ... Figure_6, SuppFigure_1, SuppFigure_2
# Each panel is lettered A, B, C ... in the order given in the plan, so the manuscript can cite
# Figure 1A -> 1B -> ... and the numbering is reproducible.
#
# Parallelism: none needed (8 figures x a handful of images); ImageMagick handles the resampling.

suppressPackageStartupMessages({
  library(data.table)
  library(magick)
})

PLAN <- "12_Manuscript/figures_plan.tsv"
SRC <- "09_Figures/plots"
OUT <- "12_Manuscript/Figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

plan <- fread(PLAN, data.table = FALSE)
stopifnot(all(c("figure", "panel", "source") %in% names(plan)))
stopifnot(!any(duplicated(paste(plan$figure, plan$panel))))
cat("panels in plan:", nrow(plan), "| figures:", length(unique(plan$figure)), "\n")

figures <- unique(plan$figure)
for (fig in figures) {
  p <- plan[plan$figure == fig, , drop = FALSE]
  p <- p[order(match(p$panel, LETTERS)), , drop = FALSE]
  files <- file.path(SRC, paste0(p$source, ".jpg"))
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("missing panels for ", fig, ": ", paste(missing, collapse = ", "), call. = FALSE)

  imgs <- lapply(seq_along(files), function(i) {
    im <- magick::image_read(files[i])
    im <- magick::image_scale(im, "900")                      # common width, keeps aspect
    magick::image_annotate(im, p$panel[i], size = 46, weight = 700, color = "black",
                           boxcolor = "white", gravity = "northwest", location = "+12+8")
  })

  n <- length(imgs)
  ncol <- if (n <= 4) 2 else if (n == 5) 2 else 3
  nrow <- ceiling(n / ncol)
  tile <- paste0(ncol, "x", nrow)
  mont <- magick::image_montage(magick::image_join(imgs), geometry = "900x900+14+14", tile = tile,
                                bg = "white", shadow = FALSE)
  mont <- magick::image_convert(mont, "jpg")
  mont <- magick::image_trim(mont)
  mont <- magick::image_border(mont, "white", "20x20")

  base <- file.path(OUT, gsub(" ", "_", fig))
  magick::image_write(mont, paste0(base, ".jpg"), format = "jpeg", quality = 95, density = 300)
  magick::image_write(mont, paste0(base, ".pdf"), format = "pdf")
  cat(sprintf("%-22s %d panels (%s) -> %s.jpg\n", fig, n, paste(p$panel, collapse = ""), base))
}
cat("done\n")
