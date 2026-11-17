#!/usr/bin/env python3
"""Build a compact cell-type signature matrix from GSE132465 for TCGA deconvolution.

Streams the dense UMI matrix and extracts only marker genes, so the full 63k x 34k
matrix is never loaded into memory.
"""

import gzip
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "03_Raw_Data" / "GSE132465"
OUT = ROOT / "06_Analysis" / "09_SingleCell"
FIG = ROOT / "09_Figures" / "plots"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

ann = pd.read_csv(RAW / "GSE132465_cell_annotation.txt.gz", sep="\t")
type_order = ["Epithelial cells", "T cells", "B cells", "Myeloids", "Stromal cells", "Mast cells"]
type_idx = {t: i for i, t in enumerate(type_order)}

markers = {
    "Epithelial cells": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "CDH1", "VIL1", "MUC2", "TFF3", "LGR5", "OLFM4", "CEACAM5", "CEACAM6", "KRT13", "KRT17"],
    "T cells": ["CD3D", "CD3E", "CD2", "CD8A", "CD8B", "CD4", "IL7R", "GZMB", "NKG7", "GNLY", "KLRD1", "TRAC", "TRBC2", "CD52", "IL32"],
    "B cells": ["MS4A1", "CD79A", "CD79B", "CD19", "MZB1", "JCHAIN", "IGHG1", "IGKC", "DERL3", "TNFRSF17", "SDC1", "CD37", "BANK1"],
    "Myeloids": ["LYZ", "CD68", "CD14", "FCGR3A", "ITGAX", "C1QA", "C1QB", "CSF1R", "S100A8", "S100A9", "IL1B", "NLRP3", "FCER1G", "TYROBP", "AIF1"],
    "Stromal cells": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "THY1", "ACTA2", "TAGLN", "FAP", "POSTN", "MMP2", "PECAM1", "VWF", "CDH5", "CLDN5", "ENG", "RAMP2", "PLVAP"],
    "Mast cells": ["TPSAB1", "TPSB2", "CPA3", "MS4A2", "KIT", "HDC", "CMA1"],
}
all_markers = {g for gs in markers.values() for g in gs}

with gzip.open(RAW / "GSE132465_natural_log_TPM_matrix.txt.gz", "rt") as f:
    header = f.readline().rstrip("\n").split("\t")
    cell_ids = header[1:]
    ann_idx = ann.set_index("Index")
    cell_types = ann_idx.loc[cell_ids, "Cell_type"].values
    codes = np.array([type_idx.get(t, -1) for t in cell_types])
    n_by_type = np.array([(codes == i).sum() for i in range(len(type_order))])

    sums = {g: np.zeros(len(type_order), dtype=np.float64) for g in all_markers}
    frac = {g: np.zeros(len(type_order), dtype=np.float64) for g in all_markers}
    found = set()
    for line in f:
        parts = line.rstrip("\n").split("\t")
        gene = parts[0]
        if gene not in all_markers:
            continue
        found.add(gene)
        vals = np.asarray(parts[1:], dtype=np.float64)
        for i in range(len(type_order)):
            m = codes == i
            if m.sum() == 0:
                continue
            # values are natural-log TPM; convert to log2(mean TPM + 1)
            sums[gene][i] = np.log2(np.exp(vals[m]).mean() + 1.0)
            frac[gene][i] = (vals[m] > 0).mean()

sig = pd.DataFrame({t: [sums[g][i] for g in sorted(found)] for i, t in enumerate(type_order)}, index=sorted(found))
sig.index.name = "gene_symbol"
sig.to_csv(OUT / "celltype_signature_symbols.tsv", sep="\t")

comp = pd.DataFrame({"cell_type": type_order, "n_cells": n_by_type})
comp["fraction"] = comp["n_cells"] / comp["n_cells"].sum()
comp.to_csv(OUT / "celltype_composition.tsv", sep="\t", index=False)

# Marker heatmap: log1p mean counts, gene-wise z-score
mat = sig.values.copy()
z = (mat - mat.mean(axis=1, keepdims=True)) / (mat.std(axis=1, keepdims=True) + 1e-8)
fig, ax = plt.subplots(figsize=(7, 10))
im = ax.imshow(z, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
ax.set_xticks(range(len(type_order)))
ax.set_xticklabels([t.replace(" cells", "") for t in type_order], rotation=45, ha="right")
ax.set_yticks(range(len(sig.index)))
ax.set_yticklabels(sig.index, fontsize=6)
ax.set_title("GSE132465 marker signatures")
fig.colorbar(im, ax=ax, shrink=0.5, label="z-score")
fig.tight_layout()
fig.savefig(FIG / "plot_1.pdf")
fig.savefig(FIG / "plot_1.jpg", dpi=300, facecolor="white")
plt.close(fig)

print("markers found:", len(found), "of", len(all_markers))
print("cell types:", dict(zip(type_order, n_by_type)))
print("wrote:", OUT / "celltype_signature_symbols.tsv")
print("wrote:", FIG / "plot_1.pdf")
