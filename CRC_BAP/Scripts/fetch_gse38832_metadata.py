#!/usr/bin/env python3
"""Fetch GSE38832 sample metadata from GEO's SOFT endpoint and write a series-matrix-compatible file.

Why this exists: `ftp.ncbi.nlm.nih.gov` is blocked from this host, so the usual
`GSE38832_series_matrix.txt.gz` download (ftp/https /geo/series/.../matrix/) fails, and the
`/geo/download/?format=file` servlet only serves `suppl` files (404 for the matrix). The GEO SOFT
endpoint on `www.ncbi.nlm.nih.gov` *is* reachable, and each per-sample brief record carries the
survival fields we need (ajcc_stage, dfs_event, dfs_time, dss_event, dss_time).

Output: 03_Raw_Data/GSE38832/GSE38832_series_matrix.txt.gz containing the header rows that
`22c_build_gse38832_dfs.R` expects (!Sample_geo_accession, !Sample_title, !Sample_characteristics_ch1),
synthesised from the 122 sample records. No expression data is written: 22c reads expression from the
RMA matrix it builds from the CEL files.

Usage: python3 Scripts/fetch_gse38832_metadata.py
"""
import gzip
import re
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "03_Raw_Data" / "GSE38832" / "GSE38832_series_matrix.txt.gz"
BASE = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=brief"
UA = {"User-Agent": "Mozilla/5.0 (compatible; CRC-TrackA metadata fetch)"}


def get(acc, tries=4):
    url = BASE.format(acc=acc)
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001 - retry any transport failure
            if attempt == tries - 1:
                raise
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError("unreachable")


def parse_sample(text):
    acc = re.search(r"^!Sample_geo_accession = (\S+)", text, re.M)
    title = re.search(r"^!Sample_title = (.+)$", text, re.M)
    chars = re.findall(r"^!Sample_characteristics_ch1 = (.+)$", text, re.M)
    return {
        "acc": acc.group(1) if acc else None,
        "title": title.group(1).strip() if title else "",
        "chars": [c.strip() for c in chars],
    }


def main():
    series = get("GSE38832")
    ids = re.findall(r"^!Series_sample_id = (GSM\d+)", series, re.M)
    if not ids:
        sys.exit("could not read sample ids from the series record")
    print(f"series record lists {len(ids)} samples; fetching brief records ...", flush=True)

    with ThreadPoolExecutor(max_workers=4) as pool:
        samples = list(pool.map(lambda a: (time.sleep(0.1), parse_sample(get(a)))[1], ids))
    samples = [s for s in samples if s["acc"]]
    print(f"parsed {len(samples)} sample records", flush=True)

    order = {a: i for i, a in enumerate(ids)}
    samples.sort(key=lambda s: order.get(s["acc"], 10**9))

    acc_row = "\t".join(['!Sample_geo_accession'] + [f'"{s["acc"]}"' for s in samples])
    title_row = "\t".join(['!Sample_title'] + [f'"{s["title"]}"' for s in samples])
    char_row = "\t".join(
        ['!Sample_characteristics_ch1'] + ['"' + "; ".join(s["chars"]) + '"' for s in samples]
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(OUT, "wt") as fh:
        fh.write("^SERIES = GSE38832\n")
        fh.write('!Series_title = "NFAT transcription factors in colorectal cancer"\n')
        fh.write(acc_row + "\n")
        fh.write(title_row + "\n")
        fh.write(char_row + "\n")

    n_dfs = sum(1 for s in samples if any(c.lower().startswith("dfs_time") and not c.endswith("NA") for c in s["chars"]))
    print(f"wrote {OUT} ({OUT.stat().st_size/1024:.0f} KB)")
    print(f"samples with a non-missing dfs_time: {n_dfs} of {len(samples)}")
    print("fields seen:", sorted({c.split(":")[0] for s in samples for c in s["chars"]}))


if __name__ == "__main__":
    main()
