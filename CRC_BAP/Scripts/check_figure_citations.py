#!/usr/bin/env python3
"""Check that the manuscript cites every figure panel exactly once, in strict sequential order.

Rule (user requirement): panels must be cited as Figure 1A, 1B, ... 6E, then Supplementary Figure 1A ... 2D,
with no panel skipped, none cited twice, and no panel cited out of order.

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

# normalise "Supplementary Figure 2C" and "Figure 1A"
pattern = re.compile(r"(Supplementary\s+)?Figure\s+(\d+)([A-Z])")
cites = []
for m in pattern.finditer(results):
    kind = "S" if m.group(1) else "M"
    cites.append((kind, int(m.group(2)), m.group(3)))

expected = []
for line in PLAN.read_text(encoding="utf-8").splitlines()[1:]:
    fig, panel, *_ = line.split("\t")
    kind = "S" if fig.lower().startswith("supplementary") else "M"
    num = int(re.findall(r"\d+", fig)[0])
    expected.append((kind, num, panel))

print(f"panels expected in plan: {len(expected)}")
print(f"citations found in Results: {len(cites)}")

problems = []
missing = [e for e in expected if e not in cites]
if missing:
    problems.append("NOT CITED: " + ", ".join(f"{k}{n}{p}" for k, n, p in missing))

dupes = {c for c in cites if cites.count(c) > 1}
if dupes:
    problems.append("CITED MORE THAN ONCE: " + ", ".join(f"{k}{n}{p}" for k, n, p in sorted(dupes)))

extra = [c for c in cites if c not in expected]
if extra:
    problems.append("CITED BUT NOT IN PLAN: " + ", ".join(f"{k}{n}{p}" for k, n, p in extra))

# order: main figures 1..N with panels A.., then supplementary 1..N
def sort_key(c):
    kind, num, panel = c
    return (0 if kind == "M" else 1, num, panel)

first_seen = []
for c in cites:
    if c not in first_seen:
        first_seen.append(c)
if first_seen != sorted(first_seen, key=sort_key):
    problems.append("OUT OF ORDER: " + " -> ".join(f"{k}{n}{p}" for k, n, p in first_seen))

if problems:
    print("\nFAIL")
    for p in problems:
        print("  -", p)
    sys.exit(1)

print("\nPASS: every panel cited exactly once, in strict sequential order")
print("order:", " ".join(f"{k}{n}{p}" for k, n, p in first_seen))
