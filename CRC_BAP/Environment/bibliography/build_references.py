#!/usr/bin/env python3
"""Fetch metadata for the selected DOIs and print the numbered list (dry run; nothing is written to the docx)."""
import json, time, urllib.parse, urllib.request

# assignment: paragraph index -> DOIs cited in that paragraph, in reading order
PLAN = [
 (29, ["10.1038/nm.3967", "10.1038/ncomms15107", "10.1158/1078-0432.ccr-23-3964"]),
 (30, ["10.1186/s13073-023-01176-5"]),
 (34, ["10.1038/nature11252", "10.1016/j.cell.2019.03.030"]),
 (40, ["10.1016/j.cels.2015.12.004", "10.1186/1471-2105-14-7"]),
 (42, ["10.1038/nmeth.2810"]),
 (48, ["10.1093/nar/gkv007"]),
 (56, ["10.1038/s41587-021-01139-4"]),
 (73, ["10.1158/0008-5472.can-22-2794"]),
 (76, ["10.1038/s41416-024-02650-6"]),
 (78, ["10.1038/s41588-025-02193-3"]),
 (80, ["10.1038/s41586-024-07769-3"]),
 (84, ["10.1200/jco.22.02582"]),
 (88, ["10.1016/j.ccell.2023.02.016", "10.1038/s41467-023-39953-w"]),
 (89, ["10.1038/s41591-023-02324-5", "10.1038/s41586-024-07747-9"]),
]

def fetch(doi):
    q = f'DOI:"{doi}"'
    url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search?" +
           urllib.parse.urlencode({"query": q, "format": "json", "pageSize": 1, "resultType": "lite"}))
    with urllib.request.urlopen(url, timeout=30) as fh:
        res = json.load(fh).get("resultList", {}).get("result", [])
    return res[0] if res else None

def vancouver(r):
    authors = [a.strip() for a in (r.get("authorString") or "").rstrip(".").split(",") if a.strip()]
    if len(authors) > 6:
        au = ", ".join(authors[:6]) + ", et al."
    else:
        au = ", ".join(authors) + "."
    title = (r.get("title") or "").rstrip(".")
    journal = r.get("journalTitle", "")
    year = r.get("pubYear", "")
    vol = r.get("journalVolume", "")
    pages = r.get("pageInfo", "")
    loc = f"{year}"
    if vol:
        loc += f";{vol}"
    if pages:
        loc += f":{pages}"
    return f"{au} {title}. {journal}. {loc}. doi:{r.get('doi','')}"

order, seen, entries = [], {}, {}
for pidx, dois in PLAN:
    nums = []
    for doi in dois:
        if doi not in seen:
            rec = fetch(doi)
            if not rec:
                print(f"  MISSING: {doi}")
                continue
            seen[doi] = len(order) + 1
            order.append(doi)
            entries[doi] = vancouver(rec)
            time.sleep(0.2)
        nums.append(seen[doi])
    print(f"paragraph {pidx}: cites {nums}")

print(f"\n{len(order)} references")
with open("Environment/bibliography/plan.json", "w", encoding="utf-8") as fh:
    json.dump({"plan": PLAN, "numbers": seen, "entries": entries}, fh, indent=1)
for doi in order:
    print(f"{seen[doi]:2}. {entries[doi][:170]}")
