#!/usr/bin/env python3
"""Train a cell2location reference model on GSE132465 and export a per-lineage expression prior.

Why this exists: `cell2location_check.py` (2026-11-17) used the 83-marker signature matrix directly as the
`cell_state_df` prior, with 400 spots and 100 epochs. That model is under-determined by construction and returned
near-uniform fractions, so it neither supported nor refuted the NNLS fractions behind F13. cell2location's own
workflow trains a reference on labelled single-cell data first; this script does exactly that and exports the
per-cell-type mean expression profile that `cell2location_check.py --signature` then uses as its prior.

Reference data (both on disk, no network):
  03_Raw_Data/GSE132465/GSE132465_cell_annotation.txt.gz       63,689 cells x 6 lineages, 23 patients
  03_Raw_Data/GSE132465/GSE132465_raw_UMI_count_matrix.txt.gz   genes x cells, raw UMI counts

Raw counts (not the 1.07 GB log-TPM matrix) are used because cell2location models counts with a negative-binomial
likelihood.

Gene panel: genes detected in >=5% of cells, excluding mitochondrial/ribosomal/haemoglobin genes, restricted to the
36,591 symbols measured on the Visium features, ranked by variance, capped at --genes; the 83 F13 markers are always
retained. Cells are capped per lineage (--cells-per-type) because the reference only estimates a per-lineage mean -
a few thousand cells per type is ample and keeps CPU training tractable.

Outputs:
  04_Processed_Data/singlecell/c2l_reference.h5ad            cached single-cell reference (gitignored dir)
  06_Analysis/09_SingleCell/cell2location_reference_signature.tsv   genes x 6 lineages, trained prior
  06_Analysis/09_SingleCell/cell2location_reference_summary.tsv     training diagnostics

Usage:
  /home/meng/miniconda3/envs/sc-analysis/bin/python Scripts/cell2location_reference.py
  ... --cells-per-type 2000 --genes 2000 --epochs 250 [--force]

Parallelism: one process, torch threads pinned to 6 (6 physical cores; nproc reports 12 SMT threads).
"""
import argparse
import gzip
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch

ROOT = Path(__file__).resolve().parents[1]
ANNOT = ROOT / "03_Raw_Data" / "GSE132465" / "GSE132465_cell_annotation.txt.gz"
UMI = ROOT / "03_Raw_Data" / "GSE132465" / "GSE132465_raw_UMI_count_matrix.txt.gz"
FEATURES = ROOT / "03_Raw_Data" / "GSE283052" / "GSE283052_features.tsv.gz"
MARKERS = ROOT / "06_Analysis" / "09_SingleCell" / "celltype_signature_symbols.tsv"
CACHE = ROOT / "04_Processed_Data" / "singlecell" / "c2l_reference.h5ad"
OUT_SIG = ROOT / "06_Analysis" / "09_SingleCell" / "cell2location_reference_signature.tsv"
OUT_SUM = ROOT / "06_Analysis" / "09_SingleCell" / "cell2location_reference_summary.tsv"
EXCLUDE_PREFIXES = ("MT-", "RPL", "RPS", "HB")


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def read_annotation(path):
    ann = pd.read_csv(path, sep="\t")
    ann = ann.set_index("Index")
    if ann.index.duplicated().any():
        sys.exit("duplicated cell ids in the annotation")
    return ann


def load_matrix(umi_path):
    """Single streaming pass: per-gene stats plus the whole matrix as sparse CSR (genes x cells)."""
    stats = {"sum": [], "sumsq": [], "detected": []}
    gene_names = []
    blocks = []
    n_chunk = 0
    t0 = time.time()
    for chunk in pd.read_csv(umi_path, sep="\t", index_col=0, chunksize=500):
        n_chunk += 1
        values = chunk.to_numpy(dtype=np.float64)
        gene_names.extend(chunk.index.astype(str).tolist())
        stats["sum"].append(values.sum(axis=1))
        stats["sumsq"].append((values ** 2).sum(axis=1))
        stats["detected"].append((values > 0).sum(axis=1))
        blocks.append(sp.csr_matrix(values.astype(np.float32)))
        if n_chunk % 10 == 0:
            log(f"  scanned {n_chunk * 500} genes ({time.time() - t0:.0f}s, "
                f"{sum(b.nnz for b in blocks) / 1e6:.0f}M non-zeros)")
        del values, chunk
    stats = pd.DataFrame({k: np.concatenate(v) for k, v in stats.items()})
    stats["gene"] = gene_names
    X = sp.vstack(blocks, format="csr")
    return stats, X, np.asarray(gene_names)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells-per-type", type=int, default=2000)
    ap.add_argument("--genes", type=int, default=2000)
    ap.add_argument("--epochs", type=int, default=250)
    ap.add_argument("--batch-size", type=int, default=2500)
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--force", action="store_true", help="retrain even if the cached reference exists")
    ap.add_argument("--umi", type=Path, default=UMI, help="raw count matrix (genes x cells)")
    ap.add_argument("--annot", type=Path, default=ANNOT, help="cell annotation table")
    ap.add_argument("--features", type=Path, default=FEATURES, help="Visium feature table")
    ap.add_argument("--markers", type=Path, default=MARKERS, help="F13 marker signature")
    ap.add_argument("--cache", type=Path, default=CACHE, help="cached reference h5ad")
    ap.add_argument("--out-sig", type=Path, default=OUT_SIG, help="exported lineage prior")
    ap.add_argument("--out-sum", type=Path, default=OUT_SUM, help="training diagnostics")
    args = ap.parse_args()
    torch.set_num_threads(args.threads)

    from cell2location.models import RegressionModel

    args.cache.parent.mkdir(parents=True, exist_ok=True)
    ann = read_annotation(args.annot)
    cells = list(ann.index)
    log(f"annotation: {len(cells)} cells, lineages: {ann['Cell_type'].value_counts().to_dict()}")

    if args.cache.exists() and not args.force:
        log(f"loading cached reference {args.cache}")
        ref = ad.read_h5ad(args.cache)
    else:
        markers = pd.read_csv(args.markers, sep="\t")["gene_symbol"].astype(str)
        feat = pd.read_csv(args.features, sep="\t", header=None)
        visium = set(feat[1].astype(str))
        header = pd.read_csv(args.umi, sep="\t", nrows=0).columns[1:].astype(str)
        log(f"markers: {len(markers)} | visium symbols: {len(visium)} | matrix cells: {len(header)}")
        if len(header) != len(ann):
            sys.exit(f"matrix has {len(header)} cells but annotation has {len(ann)}")

        t0 = time.time()
        stats, X_all, genes = load_matrix(args.umi)
        log(f"read {len(genes)} genes x {X_all.shape[1]} cells in {time.time() - t0:.0f}s "
            f"({X_all.nnz / 1e6:.0f}M non-zeros)")
        stats["mean"] = stats["sum"] / len(cells)
        stats["var"] = stats["sumsq"] / len(cells) - stats["mean"] ** 2
        stats["det_frac"] = stats["detected"] / len(cells)
        stats = stats.drop_duplicates("gene").set_index("gene")

        eligible = stats[
            (stats["det_frac"] >= 0.05)
            & (stats["var"] > 0)
            & stats.index.isin(visium)
            & ~stats.index.str.upper().str.startswith(EXCLUDE_PREFIXES)
        ]
        ranked = eligible.sort_values("var", ascending=False).index[: args.genes].tolist()
        panel = sorted(set(ranked) | (set(markers) & set(eligible.index)))
        log(f"gene panel: {len(panel)} genes ({len(set(markers) & set(panel))} markers retained)")

        pos = {g: i for i, g in enumerate(genes)}
        missing = [g for g in panel if g not in pos]
        if missing:
            log(f"WARNING: {len(missing)} panel genes absent from the matrix: {missing[:5]}")
        panel = [g for g in panel if g in pos]
        X = X_all[[pos[g] for g in panel], :]
        del X_all
        col_pos = pd.Index(header).get_indexer(ann.index)
        if (col_pos < 0).any():
            sys.exit("cell ids in the matrix header and the annotation do not match")
        X = X[:, col_pos]                     # columns now follow ann.index exactly
        log(f"selected {X.shape[0]} genes x {X.shape[1]} cells")

        rng = np.random.default_rng(20261118)
        cell_pos = pd.Series(np.arange(len(ann)), index=ann.index)
        keep_cells = []
        for lineage, labels in ann.groupby("Cell_type").groups.items():
            idx = cell_pos.loc[labels].to_numpy()
            if len(idx) > args.cells_per_type:
                idx = rng.choice(idx, args.cells_per_type, replace=False)
            keep_cells.extend(idx)
        keep_cells = np.sort(np.asarray(keep_cells))
        obs = ann.iloc[keep_cells].copy()
        counts = X[:, keep_cells].T.toarray()
        ref = ad.AnnData(
            X=counts,
            obs=obs[["Patient", "Sample", "Class", "Cell_type"]].copy(),
            var=pd.DataFrame(index=pd.Index(panel, name="gene")),
        )
        ref.obs["cell_type"] = ref.obs["Cell_type"].astype("category")
        ref.obs["sample"] = ref.obs["Sample"].astype("category")
        ref.write_h5ad(args.cache)
        log(f"cached reference: {ref.n_obs} cells x {ref.n_vars} genes -> {args.cache}")

    log("setting up the regression model")
    RegressionModel.setup_anndata(ref, batch_key="sample", labels_key="cell_type")
    mod = RegressionModel(ref)
    log(f"training: {args.epochs} epochs, batch {args.batch_size}")
    t0 = time.time()
    mod.train(max_epochs=args.epochs, batch_size=args.batch_size, train_size=1.0, lr=0.002)
    elapsed = time.time() - t0
    log(f"trained in {elapsed / 60:.1f} min")

    ref = mod.export_posterior(ref, sample_kwargs={"num_samples": 300,
                                                   "batch_size": args.batch_size})
    key = next((k for k in ("means_cell_abundance_w_sf", "means_per_cluster_mu_fg",
                            "q05_per_cluster_mu_fg") if k in ref.varm), None)
    if key is None:
        sys.exit(f"no per-lineage expression in varm: {list(ref.varm.keys())}")
    aver = np.asarray(ref.varm[key])
    names = [str(c) for c in mod.adata.obs["cell_type"].cat.categories]
    if aver.shape[1] != len(names):
        sys.exit(f"prior has {aver.shape[1]} columns but {len(names)} lineage labels")
    log(f"exported prior from varm['{key}']: {aver.shape[0]} genes x {aver.shape[1]} lineages")
    if not np.isfinite(aver).all():
        bad = int((~np.isfinite(aver)).sum())
        sys.exit(f"prior contains {bad} non-finite values - refusing to write it")
    sig = pd.DataFrame(aver, index=ref.var_names, columns=names)
    sig.insert(0, "gene_symbol", sig.index)
    sig.to_csv(args.out_sig, sep="\t", index=False)
    log(f"wrote {args.out_sig} ({sig.shape[0]} genes x {len(names)} lineages)")

    summary = pd.DataFrame([{
        "cells": ref.n_obs,
        "genes": ref.n_vars,
        "lineages": "|".join(names),
        "epochs": args.epochs,
        "batch_size": args.batch_size,
        "threads": args.threads,
        "train_minutes": round(elapsed / 60, 1),
        "mean_prior": float(np.mean(aver)),
        "per_lineage_mean": "|".join(f"{n}:{float(np.mean(sig[n])):.2f}" for n in names),
    }])
    summary.to_csv(args.out_sum, sep="\t", index=False)
    log(summary.to_string(index=False))


if __name__ == "__main__":
    main()
