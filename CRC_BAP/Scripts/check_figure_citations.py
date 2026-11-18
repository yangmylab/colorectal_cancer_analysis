#!/usr/bin/env python3
"""Check that the manuscript cites every figure panel, in an order that matches the argument.

Rule (revised 2026-11-18 at the author's request): supplementary panels are no longer cited in a block after the
main figures. Each one is cited inside the Results passage for the main figure it supports, so the rule is now:

  1. every panel in figures_plan.tsv is cited at least once in Results;
  2. main-figure panels are *first* cited in strict order (1A, 1B, ... 6E);
  3. each supplementary panel is first cited after the main panel it supports
     (the `supports_main` column of figures_plan.tsv, e.g. "Figure 1A-1D" -> anchor 1A);
  4. no panel is cited that is absent from the plan.

Later references back to a panel are normal prose and are not constrained.

Usage: python3 Scripts/check_figure_citations.py
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MS = ROOT / "12_Manuscript" / "manuscript.md"
PLAN = ROOT / "12_Manuscript" / "figures_plan.tsv"

text = MS.read_text(encoding="utf-8")
results = text.split("## Results", 1)[1].split("## Discussion", 1)[0]

pattern = re.compile(r"(Supplementary\s+)?Figure\s+(\d+)([A-Z])")
cites = []
for m in pattern.finditer(results):
    kind = "S" if m.group(1) else "M"
    cites.append((kind, int(m.group(2)), m.group(3)))

expected = []
anchors = {}
for line in PLAN.read_text(encoding="utf-8").splitlines()[1:]:
    fig, panel, source, desc, supports, supports_main = (line.split("\t") + [""] * 6)[:6]
    kind = "S" if fig.lower().startswith("supplementary") else "M"
    num = int(re.findall(r"\d+", fig)[0])
    expected.append((kind, num, panel))
    if kind == "S" and supports_main.strip() not in ("", "--"):
        m = re.search(r"Figure\s+(\d+)([A-Z])", supports_main)
        if m:
            anchors[(kind, num, panel)] = ("M", int(m.group(1)), m.group(2))

print(f"panels expected in plan: {len(expected)}")
print(f"citations found in Results: {len(cites)}")

problems = []
missing = [e for e in expected if e not in cites]
if missing:
    problems.append("NOT CITED: " + ", ".join(f"{k}{n}{p}" for k, n, p in missing))

dupes = {c for c in cites if cites.count(c) > 1}
if dupes:
    print("note: repeated citations (allowed, prose back-references):",
          ", ".join(f"{k}{n}{p}" for k, n, p in sorted(dupes)))

extra = [c for c in cites if c not in expected]
if extra:
    problems.append("CITED BUT NOT IN PLAN: " + ", ".join(f"{k}{n}{p}" for k, n, p in extra))

first = {}
for i, c in enumerate(cites):
    first.setdefault(c, i)

mains = [(k, n, p) for k, n, p in expected if k == "M"]
main_order = [first[c] for c in mains if c in first]
if main_order != sorted(main_order):
    problems.append("MAIN PANELS OUT OF ORDER: " + " ".join(f"{k}{n}{p}" for k, n, p in mains))

for supp in [c for c in expected if c[0] == "S"]:
    if supp not in first:
        continue
    anchor = anchors.get(supp)
    if anchor is None:
        problems.append(f"{supp[0]}{supp[1]}{supp[2]}: no parsable supports_main anchor")
    elif anchor not in first:
        problems.append(f"{supp[0]}{supp[1]}{supp[2]}: anchor {anchor[0]}{anchor[1]}{anchor[2]} never cited")
    elif first[supp] < first[anchor]:
        problems.append(
            f"{supp[0]}{supp[1]}{supp[2]} cited before the main panel it supports "
            f"({anchor[0]}{anchor[1]}{anchor[2]})")

if problems:
    print("\n".join(problems))
    sys.exit(1)

print("PASS: every panel cited; main panels first-cited in order; every supplementary panel cited after the")
print("      main panel it supports (supports_main column).")
print("first-citation order:", " ".join(f"{k}{n}{p}" for k, n, p in sorted(first, key=first.get)))
