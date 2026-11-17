#!/usr/bin/env python3
"""Download TCGA-COAD/READ RNA-seq from GDC and build CRC_BAP input tables.

This script downloads open-access gene-expression quantification files, parses them,
builds a sample-level clinical table, collapses duplicate patients, filters genes and
writes:
  03_Raw_Data/expression.tsv
  03_Raw_Data/clinical.tsv
  03_Raw_Data/TCGA/expression_all.tsv
  03_Raw_Data/TCGA/clinical_all.tsv

No survival values are invented: OS is derived from GDC vital_status and
days_to_death / days_to_last_follow_up. CSS is left missing because cause of death
is not available in the open clinical table.
"""

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "03_Raw_Data"
TCGA = RAW / "TCGA"
RNA = TCGA / "rna"
RNA.mkdir(parents=True, exist_ok=True)

GDC_DATA = "https://api.gdc.cancer.gov/data/{}"


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download_one(hit: dict, retries: int = 3) -> str:
    file_id = hit["file_id"]
    expected = int(hit.get("file_size", 0))
    dest = RNA / f"{file_id}.tsv"
    if dest.exists() and dest.stat().st_size == expected:
        return "cached"
    url = GDC_DATA.format(file_id)
    for attempt in range(1, retries + 1):
        try:
            tmp = dest.with_suffix(".tsv.part")
            subprocess.run(
                ["curl", "-L", "--retry", "2", "--silent", "--show-error", "--max-time", "300", "-o", str(tmp), url],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
            )
            if expected and tmp.stat().st_size != expected:
                raise IOError(f"size mismatch {tmp.stat().st_size} != {expected}")
            if hit.get("md5sum") and md5(tmp) != hit["md5sum"]:
                raise IOError("md5 mismatch")
            tmp.replace(dest)
            return "downloaded"
        except Exception as e:
            if attempt == retries:
                return f"FAILED: {file_id}: {e}"
            time.sleep(2 * attempt)
    return f"FAILED: {file_id}"


def parse_tpm(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep="\t", comment="#", dtype={"gene_id": str, "gene_name": str, "gene_type": str})
    df = df[df["gene_id"].astype(str).str.startswith("ENSG")].copy()
    df["gene_id"] = df["gene_id"].str.split(".").str[0]
    return df.groupby("gene_id")["tpm_unstranded"].sum()


def map_stage(x):
    if not isinstance(x, str):
        return np.nan
    s = x.replace("Stage", "").strip().upper()
    for level in ["IV", "III", "II", "I"]:
        if s.startswith(level):
            return level
    return np.nan


def map_side(x):
    if not isinstance(x, str):
        return np.nan
    s = x.lower()
    if any(k in s for k in ["cecum", "ascending", "hepatic flexure", "transverse"]):
        return "right"
    if any(k in s for k in ["splenic flexure", "descending", "sigmoid"]):
        return "left"
    if "rect" in s:
        return "rectum"
    return np.nan


def map_sex(x):
    if not isinstance(x, str):
        return np.nan
    s = x.lower()
    if s.startswith("f"):
        return "F"
    if s.startswith("m"):
        return "M"
    return np.nan


def build_clinical(hits, cases):
    rows = []
    for h in hits:
        case = cases.get(h["cases"][0]["case_id"], {})
        sample = (h.get("cases") or [{}])[0].get("samples", [{}])[0]
        demo = case.get("demographic") or {}
        diag = (case.get("diagnoses") or [{}])[0]
        project = (case.get("project") or {}).get("project_id", "")
        sample_type = sample.get("sample_type", "")
        vital = demo.get("vital_status")
        death = demo.get("days_to_death")
        follow = diag.get("days_to_last_follow_up")
        if vital == "Dead" and death is not None:
            os_time, os_event = float(death), 1
        elif vital == "Alive" and follow is not None:
            os_time, os_event = float(follow), 0
        else:
            os_time, os_event = np.nan, np.nan
        rfs_time = diag.get("days_to_recurrence")
        rfs_event = 1 if str(diag.get("progression_or_recurrence", "")).lower().startswith("yes") else 0
        rows.append({
            "sample_id": sample.get("submitter_id", h["file_id"]),
            "patient_id": case.get("submitter_id", ""),
            "group": "tumour" if "Primary Tumor" in sample_type else ("control" if "Normal" in sample_type else "other"),
            "group_source": "GDC",
            "disease_status": project,
            "treatment_status": "na",
            "timepoint": "baseline",
            "tissue": "colon" if project == "TCGA-COAD" else "rectum",
            "cell_type": "epithelial",
            "technology": "RNA-seq",
            "platform": "Illumina",
            "library_prep": "unknown",
            "repository_id": h["file_id"],
            "stage": map_stage(diag.get("ajcc_pathologic_stage")),
            "msi_status": np.nan,
            "sidedness": map_side(diag.get("site_of_resection_or_biopsy")),
            "age": demo.get("age_at_index"),
            "sex": map_sex(demo.get("sex_at_birth") or demo.get("gender")),
            "batch": project,
            "purity": np.nan,
            "os_time_days": os_time,
            "os_event": os_event,
            "css_time_days": np.nan,
            "css_event": np.nan,
            "rfs_time_days": float(rfs_time) if rfs_time is not None else np.nan,
            "rfs_event": rfs_event,
            "primary_diagnosis": diag.get("primary_diagnosis"),
            "site_of_resection": diag.get("site_of_resection_or_biopsy"),
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="download/process only first N files")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--write-ready", action="store_true")
    args = ap.parse_args()

    hits = json.load(open(TCGA / "file_manifest.json"))
    if args.limit:
        hits = hits[: args.limit]
    cases = {c["case_id"]: c for c in json.load(open(TCGA / "clinical_cases.json"))}
    print(f"files: {len(hits)}  workers: {args.workers}", flush=True)

    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(download_one, h): h for h in hits}
        for i, fut in enumerate(as_completed(futs), 1):
            res = fut.result()
            results.append(res)
            if i % 25 == 0 or i == len(hits):
                print(f"  downloaded/checked {i}/{len(hits)}", flush=True)

    failures = [r for r in results if str(r).startswith("FAILED")]
    if failures:
        print("\n".join(failures), flush=True)
        print(f"{len(failures)} download failures; stop.", flush=True)
        sys.exit(1)

    print("parsing TPM...", flush=True)
    series = {}
    for i, h in enumerate(hits, 1):
        series[h["file_id"]] = parse_tpm(RNA / f"{h['file_id']}.tsv")
        if i % 50 == 0 or i == len(hits):
            print(f"  parsed {i}/{len(hits)}", flush=True)
    expr = pd.DataFrame(series)
    expr.index.name = "gene_id"

    clinical = build_clinical(hits, cases)
    clinical = clinical.drop_duplicates("sample_id")

    total_tpm = expr.sum(axis=0)
    clinical["total_tpm"] = clinical["repository_id"].map(total_tpm)
    primary = clinical[clinical["group"] == "tumour"].copy()
    primary = primary.sort_values("total_tpm", ascending=False).drop_duplicates("patient_id", keep="first")
    sample_to_file = dict(zip(clinical["sample_id"], clinical["repository_id"]))
    selected_ids = [s for s in primary["sample_id"] if s in sample_to_file]
    selected_files = [sample_to_file[s] for s in selected_ids]
    discovery = expr[selected_files].copy()
    discovery.columns = selected_ids
    clinical_discovery = primary.set_index("sample_id").loc[selected_ids].reset_index()

    expressed = (discovery > 1).mean(axis=1) >= 0.2
    discovery = discovery.loc[expressed]
    if discovery.shape[0] > 15000:
        top = discovery.var(axis=1).sort_values(ascending=False).head(15000).index
        discovery = discovery.loc[top]
    discovery = discovery.sort_index()

    expr_out = discovery.reset_index()
    expr_out.columns = ["gene_id"] + list(discovery.columns)
    expr_out.to_csv(RAW / "expression.tsv", sep="\t", index=False)

    cols = ["sample_id", "patient_id", "group", "group_source", "disease_status", "treatment_status", "timepoint", "tissue", "cell_type", "technology", "platform", "library_prep", "repository_id", "stage", "msi_status", "sidedness", "age", "sex", "batch", "purity", "os_time_days", "os_event", "css_time_days", "css_event", "rfs_time_days", "rfs_event"]
    clinical_discovery[cols].to_csv(RAW / "clinical.tsv", sep="\t", index=False)

    file_to_sample = dict(zip(clinical["repository_id"], clinical["sample_id"]))
    expr_all = expr.copy()
    expr_all.columns = [file_to_sample.get(c, c) for c in expr_all.columns]
    all_expr = expr_all.loc[discovery.index].reset_index()
    all_expr.to_csv(TCGA / "expression_all.tsv", sep="\t", index=False)
    clinical.to_csv(TCGA / "clinical_all.tsv", sep="\t", index=False)

    n_events = int(clinical_discovery["os_event"].sum()) if "os_event" in clinical_discovery else 0
    print("\n=== BUILD SUMMARY ===")
    print(f"samples selected: {len(selected_ids)}")
    print(f"patients: {clinical_discovery['patient_id'].nunique()}")
    print(f"genes: {discovery.shape[0]}")
    print(f"OS events: {n_events}")
    print(f"stage missing: {clinical_discovery['stage'].isna().sum()}")
    print(f"OS time missing: {clinical_discovery['os_time_days'].isna().sum()}")
    print(f"wrote: {RAW/'expression.tsv'} ({discovery.shape[0]} x {discovery.shape[1]})")
    print(f"wrote: {RAW/'clinical.tsv'} ({len(clinical_discovery)} rows)")

    if args.write_ready:
        (RAW / "READY").write_text("TCGA-COAD/READ built\n")
        print("wrote READY")


if __name__ == "__main__":
    main()
