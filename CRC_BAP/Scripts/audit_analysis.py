#!/usr/bin/env python3
"""Audit the analysis artefacts: scripts, tables, plots, ledger, figures and endpoint discipline.

Produces 13_Reports/Consistency_Audit.md with pass/fail checks so the state of the project can be verified in one
place. Read-only: it never rewrites an analysis artefact.

Usage: python3 Scripts/audit_analysis.py
"""
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "13_Reports" / "Consistency_Audit.md"
checks = []          # (name, ok, detail)


def check(name, ok, detail=""):
    checks.append((name, bool(ok), detail))


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------- scripts
scripts = sorted((ROOT / "06_Analysis" / "1_R" / "1_Rscripts").glob("*.R"))
helpers = sorted((ROOT / "Scripts").glob("*.*"))
check("analysis R scripts present", len(scripts) >= 20, f"{len(scripts)} scripts: " +
      ", ".join(s.stem for s in scripts))
check("helper scripts present", len(helpers) >= 6, ", ".join(h.name for h in helpers))

# ---------------------------------------------------------------- plots
manifest_path = ROOT / "09_Figures" / "plots" / "plot_manifest.tsv"
manifest = read_tsv(manifest_path)
plot_ids = [m["plot"] for m in manifest]
missing_files = [p for p in plot_ids
                 if not all((ROOT / "09_Figures" / "plots" / f"{p}.{e}").exists() for e in ("pdf", "jpg"))]
check("every manifest plot has pdf+jpg", not missing_files, f"{len(plot_ids)} plots; missing: {missing_files or 'none'}")
check("manifest ids unique", len(plot_ids) == len(set(plot_ids)), f"{len(plot_ids)} rows")

# ---------------------------------------------------------------- figure plan
plan = read_tsv(ROOT / "12_Manuscript" / "figures_plan.tsv")
plan_sources = [p["source"] for p in plan]
check("plan sources all exist in manifest", set(plan_sources) <= set(plot_ids),
      f"unknown: {sorted(set(plan_sources) - set(plot_ids)) or 'none'}")
check("every plot assigned to exactly one panel", sorted(plan_sources) == sorted(plot_ids),
      f"unassigned: {sorted(set(plot_ids) - set(plan_sources)) or 'none'}; duplicated: "
      f"{sorted({s for s in plan_sources if plan_sources.count(s) > 1}) or 'none'}")

from collections import Counter
per_fig = Counter(p["figure"] for p in plan)
main = {k: v for k, v in per_fig.items() if not k.lower().startswith("supplementary")}
supp = {k: v for k, v in per_fig.items() if k.lower().startswith("supplementary")}
check("6 main figures with 5-7 panels", len(main) == 6 and all(5 <= v <= 7 for v in main.values()),
      str(main))
check("2 supplementary figures with 4-6 panels", len(supp) == 2 and all(4 <= v <= 6 for v in supp.values()),
      str(supp))

composite = ROOT / "12_Manuscript" / "Figures"
comp_files = sorted(p.stem for p in composite.glob("*.jpg"))
check("composite figures rendered", len(comp_files) == 8, ", ".join(comp_files))

# ---------------------------------------------------------------- ledger
ledger = read_csv(ROOT / "13_Reports" / "Findings_Ledger.csv")
ids = [r["finding_id"] for r in ledger]
check("ledger has 12 fields", all(len(r) == 12 for r in ledger), f"{len(ledger)} rows")
check("ledger ids unique", len(ids) == len(set(ids)), ", ".join(ids))
check("ledger in id order", ids == sorted(ids, key=lambda x: int(x[1:])), ", ".join(ids))
grades = Counter(r["grade"] for r in ledger)
check("every finding graded", all(r["grade"].strip() for r in ledger), str(dict(grades)))
allowed = {"E1", "E2", "E3", "descriptive", "not_supported"}
check("grades from the allowed vocabulary", set(grades) <= allowed, str(dict(grades)))

# ledger figures must exist
ledger_figs = set()
for r in ledger:
    ledger_figs |= {f.strip() for f in r["figure_table"].split(";") if f.strip().startswith("plot_")}
check("ledger-cited plots all exist", ledger_figs <= set(plot_ids),
      f"missing: {sorted(ledger_figs - set(plot_ids)) or 'none'}")

# ---------------------------------------------------------------- number propagation
def contains(path, needle):
    p = ROOT / path
    return p.exists() and needle in p.read_text(encoding="utf-8", errors="replace")


prop = [
    ("F06 83.3% in Proteomics_Report", "07_Validation/Proteomics_Report.md", "83.3"),
    ("F09 87.5% in Phospho_Report", "07_Validation/21_Phospho_Report.md", "87.5"),
    ("F11 HR 1.50 in DFS replication report", "07_Validation/22d_DFS_Replication_Report.md", "1.50"),
    ("F14 HR 1.42 in incremental report", "07_Validation/25_Incremental_Value_Report.md", "1.42"),
    ("F15 HR 1.00 in GSE38832 report", "07_Validation/22e_GSE38832_Report.md", "1.00"),
]
for name, path, needle in prop:
    check(name, contains(path, needle), path)

# ---------------------------------------------------------------- endpoint discipline
gse38832_clin = ROOT / "03_Raw_Data" / "GSE38832" / "22c_GSE38832_clinical_dfs.tsv"
if gse38832_clin.exists():
    header = gse38832_clin.read_text(encoding="utf-8").splitlines()[0].split("\t")
    forbidden = [c for c in header if c.startswith(("os_", "css_", "rfs_"))]
    check("GSE38832 DFS file carries no OS/CSS/RFS columns", not forbidden, f"forbidden: {forbidden or 'none'}")

os_script = (ROOT / "06_Analysis" / "1_R" / "1_Rscripts" / "04_survival_validation.R").read_text(encoding="utf-8")
check("OS validation reads only 03_Raw_Data/validation", "validation_dir" in os_script or "validation" in os_script,
      "cohort scan confined to the validation directory")

# ---------------------------------------------------------------- report
lines = ["# Consistency audit — CRC Track A", "",
         f"Generated by `Scripts/audit_analysis.py`. {len(checks)} checks.", "",
         "| Check | Result | Detail |", "|---|---|---|"]
for name, ok, detail in checks:
    lines.append(f"| {name} | {'PASS' if ok else '**FAIL**'} | {str(detail)[:180]} |")
n_fail = sum(1 for _, ok, _ in checks if not ok)
lines += ["", f"**{len(checks) - n_fail} passed, {n_fail} failed.**"]
OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"wrote {OUT}: {len(checks) - n_fail} passed, {n_fail} failed")
for name, ok, detail in checks:
    if not ok:
        print("  FAIL:", name, "-", detail)
