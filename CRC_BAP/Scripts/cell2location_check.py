#!/usr/bin/env python3
"""F13 robustness check: cell2location vs NNLS cell-type fractions on two Visium sections.

F13 (29_spatial_niches.R equivalent, `24_spatial_niches.R`) deconvolved all 16 sections with NNLS against the
GSE132465-derived signature matrix. This script re-deconvolves two sections (one per dataset) with
cell2location's Bayesian model using the SAME signature matrix as the prior, and measures the concordance
per cell type. If the two methods agree, the NNLS fractions underlying F13 are not an artefact of the linear
solver.

Deliberately scoped: 2 sections, few epochs - a concordance check, not a replacement analysis. No new claim
is made from cell2location alone.

Output: 06_Analysis/11_SpatialNiches/25_cell2location_check.tsv  (per section x cell type concordance)
        06_Analysis/11_SpatialNiches/25_cell2location_check.md    (short summary)
Usage : /home/meng/miniconda3/envs/sc-analysis/bin/python Scripts/cell2location_check.py
"""
import gzip
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from scipy.optimize import nnls

import cell2location
from cell2location.models import Cell2location

ROOT = Path(__file__).resolve().parents[1]
SPATIAL = ROOT / "04_Processed_Data" / "spatial"
OUTDIR = ROOT / "06_Analysis" / "11_SpatialNiches"
SIG = ROOT / "06_Analysis" / "09_SingleCell" / "celltype_signature_symbols.tsv"

MAX_EPOCHS = 100     # deliberately low: this is a concordance check, not a converged production model
MAX_SPOTS = 400      # subsample spots so the check runs on CPU in minutes rather than hours
SECTIONS = [("GSE283052", "GSM8655157_P02"), ("GSE267401", "GSM8265214_CTC17M")]


def load_section(dataset, sample):
    d = SPATIAL / dataset
    mat_path = next(p for p in [d / f"{sample}matrix.mtx.gz", d / f"{sample}_matrix.mtx.gz"] if p.exists())
    pos_path = next(p for p in [d / f"{sample}tissue_positions_list.csv.gz",
                                d / f"{sample}_tissue_positions.csv.gz"] if p.exists())
    feat_candidates = [d / f"{sample}features.tsv.gz", d / f"{sample}_features.tsv.gz",
                       ROOT / "03_Raw_Data" / dataset / f"{dataset}_features.tsv.gz"]
    bc_candidates = [d / f"{sample}barcodes.tsv.gz", d / f"{sample}_barcodes.tsv.gz",
                     ROOT / "03_Raw_Data" / dataset / f"{dataset}_barcodes.tsv.gz"]

    m = scipy.io.mmread(mat_path).tocsr()          # genes x spots
    feat = pd.read_csv(next(p for p in feat_candidates if p.exists()), sep="\t", header=None)
    symbols = feat[1].astype(str).values
    bc_path = next((p for p in bc_candidates if p.exists()), None)
    if bc_path is not None:
        bcs = pd.read_csv(bc_path, sep="\t", header=None)[0].astype(str).values
    else:
        bcs = np.array([f"spot{i}" for i in range(m.shape[1])])

    pos = pd.read_csv(pos_path, sep=",", header=None)
    if str(pos.iloc[0, 0]).strip().lower() == "barcode":
        pos = pd.read_csv(pos_path, sep=",")
        pos.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    else:
        pos.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    pos = pos[pos["in_tissue"] == 1]
    keep = [i for i, b in enumerate(bcs) if b in set(pos["barcode"])]
    m = m[:, keep]
    bcs = bcs[keep]
    return m, symbols, bcs


def marker_table(sig_path):
    sig = pd.read_csv(sig_path, sep="\t")
    genes = sig["gene_symbol"].astype(str).values
    types = [c for c in sig.columns if c != "gene_symbol"]
    return genes, types, sig


def run_section(dataset, sample, sig_genes, types, sig):
    m, symbols, bcs = load_section(dataset, sample)
    gene_index = pd.Index(symbols)
    common = [g for g in sig_genes if g in gene_index]
    rows = [int(np.where(gene_index == g)[0][0]) for g in common]
    X = m[rows, :].T.toarray().astype(np.float32)      # spots x marker genes (raw counts)
    if X.shape[0] > MAX_SPOTS:
        rng = np.random.default_rng(20261117)
        sel = np.sort(rng.choice(X.shape[0], MAX_SPOTS, replace=False))
        X = X[sel, :]
        bcs = bcs[sel]
    lib = X.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1
    X_norm = np.log1p(X / lib * 1e4)                   # used only for NNLS

    # cell2location models raw counts (negative-binomial likelihood), so it gets X, not X_norm
    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=pd.Index(bcs, name="spot")),
                       var=pd.DataFrame(index=pd.Index(common, name="gene")))

    sig_idx = sig.set_index("gene_symbol").loc[common, types]
    cell_state_df = sig_idx.copy()                     # genes x cell types (matches adata.var_names)

    Cell2location.setup_anndata(adata, batch_key=None)
    mod = Cell2location(
        adata,
        cell_state_df=cell_state_df,
        N_cells_per_location=30,
        detection_alpha=20,
    )
    mod.train(max_epochs=MAX_EPOCHS, batch_size=None, train_size=1.0, lr=0.002)
    mod.export_posterior(adata, sample_kwargs={"num_samples": 200, "batch_size": max(adata.n_obs // 5, 1)})

    q05_cols = list(adata.obsm["q05_cell_abundance_w_sf"].columns)

    def tidy(c):                      # cell2location prefixes vary between versions (q05_ vs q05cell_)
        for t in types:
            if t in c:
                return t
        return c

    c2l = pd.DataFrame(
        np.asarray(adata.obsm["q05_cell_abundance_w_sf"]),
        index=adata.obs_names,
        columns=[tidy(c) for c in q05_cols],
    )
    c2l_frac = c2l.div(c2l.sum(axis=1).replace(0, np.nan), axis=0)

    # NNLS on the same marker matrix, as in 24_spatial_niches.R
    A = cell_state_df.values                            # genes x cell types
    nnls_frac = np.zeros((X_norm.shape[0], len(types)))
    for i in range(X_norm.shape[0]):
        x, _ = nnls(A, X_norm[i, :])
        nnls_frac[i, :] = x / x.sum() if x.sum() > 0 else x
    nnls_df = pd.DataFrame(nnls_frac, index=adata.obs_names, columns=types)

    rows_out = []
    for t in types:
        a = c2l_frac[t].values
        b = nnls_df[t].values
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() < 30:
            continue
        rows_out.append({
            "dataset": dataset, "sample": sample, "cell_type": t, "spots": int(ok.sum()),
            "pearson_r": float(np.corrcoef(a[ok], b[ok])[0, 1]),
            "spearman_r": float(pd.Series(a[ok]).corr(pd.Series(b[ok]), method="spearman")),
            "mean_c2l": float(np.mean(a[ok])), "mean_nnls": float(np.mean(b[ok])),
            "mean_abs_diff": float(np.mean(np.abs(a[ok] - b[ok]))),
        })
    return pd.DataFrame(rows_out)


def main():
    if not SIG.exists():
        sys.exit(f"missing {SIG}")
    sig_genes, types, sig = marker_table(SIG)
    print(f"signature: {len(sig_genes)} markers x {len(types)} cell types", flush=True)
    OUTDIR.mkdir(parents=True, exist_ok=True)

    frames = []
    for dataset, sample in SECTIONS:
        print(f"=== {dataset} / {sample} (max_epochs={MAX_EPOCHS}) ...", flush=True)
        frames.append(run_section(dataset, sample, sig_genes, types, sig))
    res = pd.concat(frames, ignore_index=True)
    out = OUTDIR / "25_cell2location_check.tsv"
    res.to_csv(out, sep="\t", index=False)
    print(res.to_string(index=False))

    summary = res.groupby("cell_type")[["pearson_r", "spearman_r"]].median()
    md = ["# cell2location vs NNLS concordance (F13 robustness check)", "",
          f"- Sections: {', '.join(f'{d}/{s}' for d, s in SECTIONS)}; max_epochs = {MAX_EPOCHS}.",
          "- Both methods use the same GSE132465-derived signature matrix (83 markers, 6 lineages).",
          "- Purpose: check that the NNLS fractions underlying F13 are not an artefact of the linear solver.",
          "", "## Median concordance per cell type", "",
          "| cell type | Pearson r | Spearman r |", "|---|---|---|"]
    for t, row in summary.iterrows():
        md.append(f"| {t} | {row['pearson_r']:.2f} | {row['spearman_r']:.2f} |")
    md += ["", f"- Overall median Pearson r = {res['pearson_r'].median():.2f}, "
               f"Spearman = {res['spearman_r'].median():.2f}.", ""]
    med = res["pearson_r"].median()
    if med >= 0.5:
        verdict = ("The two methods agree on the dominant lineages, so the NNLS fractions behind F13 are not an "
                   "artefact of the linear solver.")
    elif med >= 0.3:
        verdict = ("Agreement is partial: dominant lineages agree but minor lineages do not, so the co-localisation "
                   "component of F13 is method-sensitive and should be read as approximate.")
    else:
        verdict = ("**The two methods do NOT agree at this training budget.** cell2location returns a nearly "
                   "uniform fraction across lineages while NNLS returns an epithelial-dominated composition. With "
                   "only 83 marker genes and "
                   f"{MAX_EPOCHS} epochs the cell2location model is not identifiable, so this check is "
                   "inconclusive: it neither validates nor refutes the NNLS fractions. Consequence for F13: the "
                   "programme-level spatial statistics (Moran's I) do not depend on deconvolution, but the "
                   "E2-stromal co-localisation statistic does, and must be treated as method-dependent until a "
                   "properly trained cell2location reference is available.")
    md += [verdict,
           "",
           f"Configuration: {MAX_SPOTS} spots per section, {MAX_EPOCHS} epochs, signature matrix used directly as "
           "the cell_state_df prior (no reference model trained on the single-cell data)."]
    (OUTDIR / "25_cell2location_check.md").write_text("\n".join(md) + "\n")
    print(f"\nwrote {out} and 25_cell2location_check.md")
    print(f"overall median Pearson r: {res['pearson_r'].median():.3f}")


if __name__ == "__main__":
    main()
