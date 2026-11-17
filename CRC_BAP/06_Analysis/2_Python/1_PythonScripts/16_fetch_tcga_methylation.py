#!/usr/bin/env python3
"""Download TCGA-COAD/READ 450K methylation beta values from GDC and build a beta matrix.

Open-access files only (SeSAMe 'Methylation Beta Value', Illumina Human Methylation 450).
Writes:
  03_Raw_Data/TCGA/methylation/raw/<file_id>.txt
  03_Raw_Data/TCGA/methylation/methylation_files.tsv   (file <-> case/sample manifest)
  03_Raw_Data/TCGA/methylation/beta_matrix.npy         (probes x samples, float32)
  03_Raw_Data/TCGA/methylation/beta_matrix_probes.txt  (probe ids, one per line)
  03_Raw_Data/TCGA/methylation/beta_matrix_samples.txt (sample ids, one per line)
"""

import argparse
import json
import subprocess
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
METH = ROOT / "03_Raw_Data" / "TCGA" / "methylation"
RAW = METH / "raw"
RAW.mkdir(parents=True, exist_ok=True)

GDC_FILES = "https://api.gdc.cancer.gov/files"
GDC_DATA = "https://api.gdc.cancer.gov/data/{}"

FIELDS = ",".join([
    "file_id", "file_name", "file_size", "md5sum",
    "cases.submitter_id", "cases.samples.submitter_id", "cases.samples.sample_type",
])


def query_files() -> list:
    project = {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-COAD", "TCGA-READ"]}}
    filters = {"op": "and", "content": [
        project,
        {"op": "in", "content": {"field": "data_category", "value": ["DNA Methylation"]}},
        {"op": "in", "content": {"field": "data_type", "value": ["Methylation Beta Value"]}},
        {"op": "in", "content": {"field": "platform", "value": ["Illumina Human Methylation 450"]}},
        {"op": "in", "content": {"field": "access", "value": ["open"]}},
    ]}
    payload = {"filters": json.dumps(filters), "fields": FIELDS, "format": "JSON", "size": 5000}
    req = urllib.request.Request(
        GDC_FILES,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=180) as fh:
        return json.load(fh)["data"]["hits"]


def download_one(hit: dict, retries: int = 4) -> tuple:
    file_id = hit["file_id"]
    expected = int(hit.get("file_size", 0))
    dest = RAW / f"{file_id}.txt"
    if dest.exists() and dest.stat().st_size == expected:
        return file_id, "cached"
    url = GDC_DATA.format(file_id)
    for attempt in range(1, retries + 1):
        try:
            tmp = dest.with_suffix(".txt.part")
            subprocess.run(
                ["curl", "-L", "--retry", "3", "--silent", "--show-error", "--max-time", "600", "-o", str(tmp), url],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
            )
            if expected and tmp.stat().st_size != expected:
                raise IOError(f"size mismatch {tmp.stat().st_size} != {expected}")
            tmp.replace(dest)
            return file_id, "downloaded"
        except Exception as exc:  # noqa: BLE001
            if attempt == retries:
                return file_id, f"FAILED: {exc}"
            time.sleep(2 * attempt)
    return file_id, "FAILED"


def read_beta(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep="\t", header=0, dtype=str)
    col = None
    for candidate in df.columns:
        if "beta" in candidate.lower():
            col = candidate
            break
    if col is None:
        if df.shape[1] < 2:
            raise ValueError(f"unexpected methylation file layout: {path}")
        col = df.columns[1]
    probe_col = df.columns[0]
    out = pd.Series(pd.to_numeric(df[col], errors="coerce").values, index=df[probe_col].values, dtype="float32")
    return out[out.index.str.startswith("cg", na=False)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--download-only", action="store_true")
    ap.add_argument("--threads", type=int, default=6)
    args = ap.parse_args()

    hits = query_files()
    total_gb = sum(int(h.get("file_size", 0)) for h in hits) / 1e9
    print(f"{len(hits)} methylation beta files, {total_gb:.2f} GB total", flush=True)

    rows, primary = [], []
    for h in hits:
        case = (h.get("cases") or [{}])[0]
        samples = case.get("samples") or [{}]
        sample = samples[0]
        rows.append({
            "file_id": h["file_id"],
            "file_name": h.get("file_name", ""),
            "file_size": int(h.get("file_size", 0)),
            "patient_id": case.get("submitter_id", ""),
            "sample_id": sample.get("submitter_id", ""),
            "sample_type": sample.get("sample_type", ""),
        })
        if sample.get("sample_type") == "Primary Tumor":
            primary.append(h["file_id"])
    manifest = pd.DataFrame(rows).sort_values(["sample_type", "patient_id"])
    manifest.to_csv(METH / "methylation_files.tsv", sep="\t", index=False)
    print(f"primary-tumour files: {len(primary)} of {len(manifest)}", flush=True)

    if not args.download_only:
        with ThreadPoolExecutor(max_workers=args.threads) as pool:
            futures = [pool.submit(download_one, h) for h in hits]
            done = 0
            for fut in as_completed(futures):
                file_id, status = fut.result()
                done += 1
                if status.startswith("FAILED"):
                    print(f"[{done}/{len(hits)}] {status}", flush=True)
                elif done % 25 == 0:
                    print(f"[{done}/{len(hits)}] {status}", flush=True)

    missing = [p for p in primary if not (RAW / f"{p}.txt").exists()]
    if missing:
        print(f"WARNING: {len(missing)} primary files missing; building the matrix from the rest", flush=True)
    primary = [p for p in primary if (RAW / f"{p}.txt").exists()]
    if not primary:
        print("ERROR: no primary-tumour files available", flush=True)
        return 1

    # Build the matrix on primary tumours, one patient per sample id.
    meta = manifest[manifest["file_id"].isin(primary)].drop_duplicates("patient_id")
    sample_ids = meta["sample_id"].tolist()
    probe_index = None
    matrix = None
    for i, (file_id, sample_id) in enumerate(zip(meta["file_id"], sample_ids), start=1):
        series = read_beta(RAW / f"{file_id}.txt")
        if probe_index is None:
            probe_index = series.index.values
            matrix = np.full((len(probe_index), len(sample_ids)), np.nan, dtype="float32")
            pos = pd.Series(np.arange(len(probe_index)), index=probe_index)
        aligned = pos.reindex(series.index).values
        ok = ~pd.isna(aligned)
        matrix[aligned[ok].astype(int), i - 1] = series.values[ok]
        if i % 50 == 0:
            print(f"assembled {i}/{len(sample_ids)} samples", flush=True)

    keep = np.isfinite(matrix).mean(axis=1) >= 0.9
    matrix = matrix[keep]
    probes = probe_index[keep]
    np.save(METH / "beta_matrix.npy", matrix)
    pd.Series(probes).to_csv(METH / "beta_matrix_probes.txt", index=False, header=False)
    pd.Series(sample_ids).to_csv(METH / "beta_matrix_samples.txt", index=False, header=False)
    print(f"beta matrix: {matrix.shape[0]} probes x {matrix.shape[1]} samples", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
