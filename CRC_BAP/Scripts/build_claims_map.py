#!/usr/bin/env python3
"""Build 12_Manuscript/claims_to_evidence.tsv: every manuscript claim mapped to its graded finding.

The ledger is the authority; this script only adds the manuscript section in which each finding is used and a
short quotable statement. Edit SECTION_MAP when the manuscript structure changes, then re-run.

Usage: python3 Scripts/build_claims_map.py
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "13_Reports" / "Findings_Ledger.csv"
OUT = ROOT / "12_Manuscript" / "claims_to_evidence.tsv"

SECTION_MAP = {
    "F01": "Results 4 / Abstract / Discussion",
    "F02": "Results 1 (post hoc, explicitly labelled)",
    "F03": "Results 4 / Abstract",
    "F04": "Results 1 (not replicated)",
    "F05": "Results 1 (context only)",
    "F06": "Results 2 / Abstract",
    "F07": "Results 3 / Abstract",
    "F08": "Results 3 / Abstract / What we do not claim",
    "F09": "Results 2 / Abstract",
    "F11": "Results 5 / Abstract",
    "F12": "Results 2 (exploratory)",
    "F13": "Results 3 (descriptive; caveats stated)",
    "F14": "Results 5 / Abstract / Discussion",
    "F15": "Results 4 / Abstract",
}

WORDING = {
    "F01": "no incremental association with overall survival after adjustment",
    "F02": "K=4 is exploratory only (post hoc, did not replicate)",
    "F03": "overall-survival association did not replicate (meta p=0.081)",
    "F04": "the K=4 association did not replicate either",
    "F05": "single-cell deconvolution is descriptive context, not a finding about the axis",
    "F06": "protein-level reproduction of the ecotype assignment",
    "F07": "the axis does not reproduce in promoter methylation",
    "F08": "the bulk classifier does not transfer to spatial spots",
    "F09": "phosphoproteome-level reproduction, permutation-controlled",
    "F11": "disease-free association replicates across three cohorts (E2)",
    "F12": "kinase activity differences are inferred and exploratory (E1)",
    "F13": "spatial autocorrelation and stromal co-localisation, reported descriptively",
    "F14": "no incremental value beyond stage, MSI and CMS",
    "F15": "no association with cancer-specific death in a third cohort",
}

FORBIDDEN = {
    "F01": "must not be cited as 'no association at all' - it is a null for the adjusted OS model",
    "F02": "must never be presented as confirmatory",
    "F03": "must not be described as replication",
    "F07": "must not be described as 'methylation-independent biology' - only promoter 450K was tested",
    "F08": "must not be used to claim spots lack ecotype biology",
    "F11": "must always be reported together with F14 (non-incremental over CMS)",
    "F12": "must not be phrased as measured kinase activity",
    "F13": "must carry the permutation-ceiling and cell2location caveats",
    "F14": "must not be over-generalised: one cohort with CMS",
}


def main():
    with open(LEDGER, newline="", encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh)]
    fields = ["finding_id", "grade", "manuscript_section", "statement", "n", "datasets", "figures",
              "ledger_statistic", "citation_caveat"]
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            fid = r["finding_id"]
            w.writerow({
                "finding_id": fid,
                "grade": r["grade"],
                "manuscript_section": SECTION_MAP.get(fid, "unmapped"),
                "statement": WORDING.get(fid, r["claim"][:120]),
                "n": r["n"],
                "datasets": r["dataset_s"],
                "figures": r["figure_table"],
                "ledger_statistic": r["statistic"],
                "citation_caveat": FORBIDDEN.get(fid, ""),
            })
    print(f"wrote {OUT} with {len(rows)} graded claims")
    unmapped = [r["finding_id"] for r in rows if r["finding_id"] not in SECTION_MAP]
    if unmapped:
        print("WARNING unmapped findings:", ", ".join(unmapped))


if __name__ == "__main__":
    main()
