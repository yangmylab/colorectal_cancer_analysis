#!/usr/bin/env python3
"""Build gene-level promoter methylation from the full TCGA 450K beta matrix.

Uses 03_Raw_Data/TCGA/methylation/probe2gene_promoter.tsv.gz (promoter probes with a
unique gene, exported from the 450K annotation) so that gene selection is annotation
based, not variance based (avoids circular differential-methylation testing).

Writes:
  03_Raw_Data/TCGA/methylation/gene_beta.tsv.gz    (genes x samples, mean beta)
  03_Raw_Data/TCGA/methylation/global_beta.tsv     (mean beta over all probes per sample)
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
METH = ROOT / "03_Raw_Data" / "TCGA" / "methylation"

matrix = np.load(METH / "beta_matrix.npy", mmap_mode="r")
probes = pd.read_csv(METH / "beta_matrix_probes.txt", header=None)[0].values
samples = pd.read_csv(METH / "beta_matrix_samples.txt", header=None)[0].values
print(f"beta matrix: {matrix.shape[0]} probes x {matrix.shape[1]} samples", flush=True)

global_beta = pd.DataFrame(
    {"sample_id": samples, "global_beta": np.nanmean(matrix, axis=0)}
)
global_beta.to_csv(METH / "global_beta.tsv", sep="\t", index=False, float_format="%.4f")
print(f"global beta range: {global_beta.global_beta.min():.3f}-{global_beta.global_beta.max():.3f}", flush=True)

mapping = pd.read_csv(METH / "probe2gene_promoter.tsv.gz", sep="\t")
pos = pd.Series(np.arange(len(probes)), index=probes)
rows = pos.reindex(mapping["probe_id"].values).values
ok = ~pd.isna(rows)
rows = rows[ok].astype(int)
genes = mapping["gene"].values[ok]
print(f"promoter probes matched in matrix: {len(rows)} of {len(mapping)}", flush=True)

sub = np.asarray(matrix[rows], dtype="float32")
order = np.argsort(genes, kind="stable")
sub = sub[order]
genes_sorted = genes[order]
uniq, start = np.unique(genes_sorted, return_index=True)
gene_beta = np.add.reduceat(sub, start, axis=0) / np.diff(np.append(start, len(genes_sorted)))[:, None]
gene_beta = np.round(gene_beta, 4)

out = pd.DataFrame(gene_beta, index=uniq, columns=samples)
out.index.name = "gene"
out.to_csv(METH / "gene_beta.tsv.gz", sep="\t", compression="gzip", float_format="%.4f")
print(f"wrote gene_beta: {out.shape[0]} genes x {out.shape[1]} samples", flush=True)
