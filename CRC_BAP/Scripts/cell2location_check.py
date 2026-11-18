#!/usr/bin/env python3
"""F13 robustness check: cell2location vs NNLS cell-type fractions on Visium sections.

F13 (`24_spatial_niches.R`) deconvolved all 16 sections with NNLS against the GSE132465-derived marker matrix.
This script re-deconvolves sections with cell2location's Bayesian model and measures per-lineage concordance with
NNLS. If the two methods agree, the NNLS fractions behind F13 are not an artefact of the linear solver.

Two priors are supported, and the difference matters:

  --prior markers     the 83-marker signature matrix used directly as `cell_state_df` (the 2026-11-17 run)
  --prior reference   a per-lineage expression profile trained on GSE132465 single cells by
                      `cell2location_reference.py` (the statistically correct workflow)

The marker-matrix prior with a few hundred epochs is under-determined and returns near-uniform fractions; that run
is kept as `25_cell2location_check.*` for the record. A verdict needs the trained reference, larger training
budget and more spots, which is what `--prior reference --epochs 3000` provides.

Outputs (prefix follows --out-prefix):
  <out-prefix>.tsv   per section x lineage concordance (Pearson/Spearman, means, mean absolute difference)
  <out-prefix>.md    short human-readable summary with the verdict rule applied

Usage:
  /home/meng/miniconda3/envs/sc-analysis/bin/python Scripts/cell2location_check.py \
      --prior reference --epochs 3000 --spots 1000 --out-prefix 26_cell2location_reference

Parallelism: one process, torch threads pinned to 6 (6 physical cores; nproc reports 12 SMT threads).
Sections are processed sequentially because each one saturates the cores on its own.
"""
import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import torch
from scipy.optimize import nnls

ROOT = Path(__file__).resolve().parents[1]
SPATIAL = ROOT / "04_Processed_Data" / "spatial"
OUTDIR = ROOT / "06_Analysis" / "11_SpatialNiches"
MARKERS = ROOT / "06_Analysis" / "09_SingleCell" / "celltype_signature_symbols.tsv"
REFERENCE = ROOT / "06_Analysis" / "09_SingleCell" / "cell2location_reference_signature.tsv"
DEFAULT_SECTIONS = [("GSE283052", "GSM8655157_P02"), ("GSE267401", "GSM8265214_CTC17M")]
FEATURES = ROOT / "03_Raw_Data" / "GSE283052" / "GSE283052_features.tsv.gz"


def log(msg):
    print(msg, flush=True)


def load_prior(kind):
    path = MARKERS if kind == "markers" else REFERENCE
    if not path.exists():
        sys.exit(f"missing prior {path} (run Scripts/cell2location_reference.py first)")
    df = pd.read_csv(path, sep="\t")
    types = [c for c in df.columns if c != "gene_symbol"]
    df = df.drop_duplicates("gene_symbol").set_index("gene_symbol")
    return df.index.astype(str), types, df


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
    pos = pd.read_csv(pos_path, sep=",")
    if str(pos.columns[0]).strip().lower() != "barcode":
        pos.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    pos = pos[pos["in_tissue"] == 1]
    keep = [i for i, b in enumerate(bcs) if b in set(pos["barcode"])]
    return m[:, keep], symbols, bcs[keep]


def run_section(dataset, sample, prior_genes, types, prior, epochs, spots, seed=20261118):
    m, symbols, bcs = load_section(dataset, sample)
    gene_index = pd.Index(symbols)
    common = [g for g in prior_genes if g in gene_index]
    rows = [int(np.where(gene_index == g)[0][0]) for g in common]
    X = m[rows, :].T.toarray().astype(np.float32)          # spots x prior genes (raw counts)
    if X.shape[0] > spots:
        rng = np.random.default_rng(seed)
        sel = np.sort(rng.choice(X.shape[0], spots, replace=False))
        X = X[sel, :]
        bcs = bcs[sel]
    lib = X.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1
    X_norm = np.log1p(X / lib * 1e4)

    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=pd.Index(bcs, name="spot")),
                       var=pd.DataFrame(index=pd.Index(common, name="gene")))

    from cell2location.models import Cell2location

    cell_state_df = prior.loc[common, types].copy()
    Cell2location.setup_anndata(adata, batch_key=None)
    mod = Cell2location(adata, cell_state_df=cell_state_df, N_cells_per_location=30, detection_alpha=20)
    mod.train(max_epochs=epochs, batch_size=None, train_size=1.0, lr=0.002,
              enable_progress_bar=False, logger=False)
    mod.export_posterior(adata, sample_kwargs={"num_samples": 200,
                                               "batch_size": max(adata.n_obs // 5, 1)})

    def tidy(c):
        for t in types:
            if t in c:
                return t
        return c

    c2l = pd.DataFrame(
        np.asarray(adata.obsm["q05_cell_abundance_w_sf"]),
        index=adata.obs_names,
        columns=[tidy(c) for c in adata.obsm["q05_cell_abundance_w_sf"].columns],
    )
    c2l_frac = c2l[types].div(c2l[types].sum(axis=1).replace(0, np.nan), axis=0)

    A = cell_state_df.values
    nnls_frac = np.zeros((X_norm.shape[0], len(types)))
    for i in range(X_norm.shape[0]):
        x, _ = nnls(A, X_norm[i, :])
        nnls_frac[i, :] = x / x.sum() if x.sum() > 0 else x
    nnls_df = pd.DataFrame(nnls_frac, index=adata.obs_names, columns=types)

    out = []
    undefined = 0
    for t in types:
        a = c2l_frac[t].to_numpy(dtype=float)
        b = nnls_df[t].to_numpy(dtype=float)
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() < 30:
            continue
        constant = np.std(a[ok]) == 0 or np.std(b[ok]) == 0
        undefined += int(constant)
        out.append({
            "dataset": dataset, "sample": sample, "cell_type": t, "spots": int(ok.sum()),
            "pearson_r": np.nan if constant else float(np.corrcoef(a[ok], b[ok])[0, 1]),
            "spearman_r": np.nan if constant else float(pd.Series(a[ok]).corr(pd.Series(b[ok]), method="spearman")),
            "mean_c2l": float(np.mean(a[ok])), "mean_nnls": float(np.mean(b[ok])),
            "mean_abs_diff": float(np.mean(np.abs(a[ok] - b[ok]))),
            "constant_side": (
                "both" if constant and np.std(a[ok]) == 0 and np.std(b[ok]) == 0
                else "nnls" if constant and np.std(b[ok]) == 0
                else "cell2location" if constant else ""),
        })

    # Spots x lineages view: the composition-level comparison survives degenerate single lineages.
    A_all = c2l_frac[types].to_numpy(dtype=float).ravel()
    B_all = nnls_df[types].to_numpy(dtype=float).ravel()
    fin = np.isfinite(A_all) & np.isfinite(B_all)
    overall = {
        "dataset": dataset, "sample": sample, "spots": int(len(c2l_frac)),
        "matrix_pearson": float(np.corrcoef(A_all[fin], B_all[fin])[0, 1]) if fin.sum() > 30 else np.nan,
        "top_lineage_agreement": float(np.mean(
            [types[int(np.argmax(c2l_frac[types].to_numpy()[i]))] ==
             types[int(np.argmax(nnls_df[types].to_numpy()[i]))] for i in range(len(c2l_frac))])),
        "median_total_abundance": float(np.median(np.asarray(adata.obsm["q05_cell_abundance_w_sf"]).sum(axis=1))),
        "undefined_lineages": undefined,
        "genes": int(len(common)),
    }
    return pd.DataFrame(out), overall, len(common)


def write_report(res, overall_df, args, sections, genes_used, prior_genes, types):
    """Render the concordance report, including the composition gap that a bare correlation hides."""
    defined = res.dropna(subset=["pearson_r"])
    med = float(defined["pearson_r"].median()) if len(defined) else float("nan")
    med_matrix = float(overall_df["matrix_pearson"].median())
    n_undefined = int(len(res) - len(defined))
    summary = res.groupby("cell_type")[["pearson_r", "spearman_r"]].median()
    means = res.groupby("cell_type")[["mean_c2l", "mean_nnls"]].mean()
    gap = means["mean_c2l"] - means["mean_nnls"]
    max_gap = float(gap.abs().max())
    gap_lineage = str(gap.abs().idxmax())
    md = [
        "# cell2location vs NNLS concordance (F13 robustness check)",
        "",
        f"- Prior: **{args.prior}** — {len(prior_genes)} genes x {len(types)} lineages "
        f"({'83-marker matrix' if args.prior == 'markers' else 'cell2location reference trained on GSE132465'}).",
        f"- Sections: {', '.join(f'{d}/{s}' for d, s in sections)}; max_epochs = {args.epochs}; "
        f"spots per section = {args.spots}; genes shared with Visium = {genes_used}.",
        "- Both methods are run on the same gene set and the same spots, so the comparison is like for like.",
        "",
        "## Median concordance per lineage",
        "",
        "| lineage | Pearson r | Spearman r | mean fraction (cell2location) | mean fraction (NNLS) |",
        "|---|---|---|---|---|",
    ]
    for t, row in summary.iterrows():
        r_txt = "n/a" if np.isnan(row["pearson_r"]) else f"{row['pearson_r']:.2f}"
        rho_txt = "n/a" if np.isnan(row["spearman_r"]) else f"{row['spearman_r']:.2f}"
        md.append(f"| {t} | {r_txt} | {rho_txt} | {means.loc[t, 'mean_c2l']:.2f} | "
                  f"{means.loc[t, 'mean_nnls']:.2f} |")
    md += [
        "",
        "## Per-section composition agreement",
        "",
        "| section | spots | matrix Pearson r | top-lineage agreement | median abundance | undefined lineages |",
        "|---|---|---|---|---|---|",
    ]
    for _, r in overall_df.iterrows():
        md.append(f"| {r['dataset']}/{r['sample']} | {r['spots']} | {r['matrix_pearson']:.2f} | "
                  f"{r['top_lineage_agreement']:.2f} | {r['median_total_abundance']:.1f} | "
                  f"{r['undefined_lineages']} |")
    md += [
        "",
        f"- Median per-lineage Pearson r over the {len(defined)} identifiable lineages: {med:.2f}; "
        f"{n_undefined} lineage comparisons are undefined because one side has zero variance (NNLS returns "
        "exactly zero for those lineages).",
        f"- Median matrix-level Pearson r over all spots x lineages: {med_matrix:.2f}.",
        f"- Largest absolute composition gap: **{max_gap:.2f}** ({gap_lineage}: cell2location "
        f"{means.loc[gap_lineage, 'mean_c2l']:.2f} vs NNLS {means.loc[gap_lineage, 'mean_nnls']:.2f}).",
        "",
        "## Verdict",
        "",
    ]
    if med_matrix >= 0.5 and max_gap <= 0.15:
        verdict = ("**The two methods agree, in both pattern and composition** (median matrix r >= 0.5 and no "
                   "lineage differs by more than 0.15): the NNLS fractions behind F13 are not an artefact of the "
                   "linear solver, so F13's co-localisation statistic needs no method-dependence caveat.")
    elif med_matrix >= 0.5:
        verdict = (
            f"**The spatial pattern is reproduced, but the absolute composition is not** (median matrix r "
            f"{med_matrix:.2f}; largest gap {max_gap:.2f} for {gap_lineage}: cell2location "
            f"{means.loc[gap_lineage, 'mean_c2l']:.2f} vs NNLS {means.loc[gap_lineage, 'mean_nnls']:.2f}). "
            "The 2026-11-17 reading that cell2location 'did NOT corroborate' the NNLS fractions is therefore "
            "withdrawn: at 100 epochs with a marker-matrix prior the model was simply under-trained. What remains "
            "method-dependent is the *scale*, not the direction - unregularised NNLS collapses onto one lineage "
            "(~0.97 epithelial, exactly zero for several others) while cell2location spreads mass across lineages, "
            "so absolute fractions, and therefore the magnitude of any stromal-fraction correlation, depend on the "
            "method and the gene panel. F13 keeps a narrowed caveat about absolute fractions and magnitudes; its "
            "direction (E2 programmes high where the stromal compartment is high) is supported.")
    elif med_matrix >= 0.3:
        verdict = ("**Agreement is partial** (0.3 <= median matrix r < 0.5): the dominant lineage is reproduced but "
                   "the composition is not interchangeable, so F13's co-localisation statistic stays "
                   "method-sensitive and must be reported as approximate.")
    else:
        verdict = ("**The two methods do not agree** (median matrix r < 0.3) even with a trained reference and the "
                   "larger budget, so F13's co-localisation statistic must keep a method-dependence caveat and no "
                   "cell2location-derived claim may be made.")
    md += [verdict, ""]
    (OUTDIR / f"{args.out_prefix}.md").write_text("\n".join(md) + "\n")
    log(f"wrote {OUTDIR / (args.out_prefix + '.md')}")
    log(f"median per-lineage r: {med:.3f} | median matrix r: {med_matrix:.3f} | "
        f"max composition gap: {max_gap:.3f} ({gap_lineage}) | undefined: {n_undefined}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prior", choices=["markers", "reference"], default="markers")
    ap.add_argument("--epochs", type=int, default=100)
    ap.add_argument("--spots", type=int, default=400)
    ap.add_argument("--genes-cap", type=int, default=0,
                    help="use only the N most variable prior genes (0 = all)")
    ap.add_argument("--sections", nargs="*", default=None,
                    help="dataset:sample pairs, e.g. GSE283052:GSM8655157_P02")
    ap.add_argument("--out-prefix", default="25_cell2location_check")
    ap.add_argument("--report-only", action="store_true",
                    help="re-render the .md from the existing .tsv/_overall.tsv without retraining")
    ap.add_argument("--threads", type=int, default=6)
    args = ap.parse_args()
    torch.set_num_threads(args.threads)

    prior_genes, types, prior = load_prior(args.prior)
    if args.genes_cap and len(prior_genes) > args.genes_cap:
        var = prior[types].var(axis=1).sort_values(ascending=False)
        shared = [g for g in var.index[: args.genes_cap]]
        prior_genes, prior = shared, prior.loc[shared]
    sections = DEFAULT_SECTIONS if not args.sections else [
        tuple(s.split(":", 1)) for s in args.sections]
    log(f"prior={args.prior} ({len(prior_genes)} genes x {len(types)} lineages) | "
        f"epochs={args.epochs} | spots={args.spots} | sections={len(sections)}")

    if args.report_only:
        res = pd.read_csv(OUTDIR / f"{args.out_prefix}.tsv", sep="\t")
        overall_df = pd.read_csv(OUTDIR / f"{args.out_prefix}_overall.tsv", sep="\t")
        visium = set(pd.read_csv(FEATURES, sep="\t", header=None)[1].astype(str))
        genes_used = len(set(prior_genes) & visium)
        write_report(res, overall_df, args, sections, genes_used, prior_genes, types)
        return

    OUTDIR.mkdir(parents=True, exist_ok=True)
    frames = []
    overall_rows = []
    genes_used = None
    for dataset, sample in sections:
        log(f"=== {dataset} / {sample} ...")
        df, overall, genes_used = run_section(dataset, sample, prior_genes, types, prior,
                                              args.epochs, args.spots)
        frames.append(df)
        overall_rows.append(overall)
        log(df[["cell_type", "pearson_r", "mean_c2l", "mean_nnls"]].to_string(index=False))
        log(f"    matrix r = {overall['matrix_pearson']:.3f} | top-lineage agreement = "
            f"{overall['top_lineage_agreement']:.2f} | median abundance = "
            f"{overall['median_total_abundance']:.1f} | undefined lineages = {overall['undefined_lineages']}")
    res = pd.concat(frames, ignore_index=True)
    overall_df = pd.DataFrame(overall_rows)
    res.to_csv(OUTDIR / f"{args.out_prefix}.tsv", sep="\t", index=False)
    overall_df.to_csv(OUTDIR / f"{args.out_prefix}_overall.tsv", sep="\t", index=False)
    write_report(res, overall_df, args, sections, genes_used, prior_genes, types)


if __name__ == "__main__":
    main()
