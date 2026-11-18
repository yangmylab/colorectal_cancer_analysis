#!/usr/bin/env python3
"""Query Europe PMC for citation candidates on the manuscript's claims, and write a compact TSV.

Records only what the API returns (title, journal, year, volume, pages, DOI, PMID, citation count); nothing is
written from memory. Journals are restricted to a high-impact whitelist and years to 2023-2026, except where a
query targets a framework or method original, which the author agreed may be older.
"""
import json, time, urllib.parse, urllib.request
from pathlib import Path

BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
JOURNALS = ['"Nature"', '"Nature medicine"', '"Nature cancer"', '"Nature genetics"', '"Nature communications"',
            '"Cancer cell"', '"Cell"', '"Gut"', '"Journal of clinical oncology"', '"Annals of oncology"',
            '"Gastroenterology"', '"Lancet oncology"', '"JAMA oncology"', '"Cancer discovery"',
            '"Clinical cancer research"', '"Nature reviews. Cancer"', '"Nature reviews. Clinical oncology"']

QUERIES = {
 "taxonomy": '(TITLE_ABS:"colorectal cancer") AND (TITLE_ABS:subtype OR TITLE_ABS:classification OR TITLE_ABS:CMS)',
 "reproducibility": '(TITLE_ABS:colorectal) AND (TITLE_ABS:reproducibility OR TITLE_ABS:reproducible OR TITLE_ABS:classifier OR TITLE_ABS:consensus)',
 "microenvironment": '(TITLE_ABS:colorectal) AND (TITLE_ABS:stroma OR TITLE_ABS:stromal OR TITLE_ABS:microenvironment OR TITLE_ABS:fibroblast)',
 "singlecell": '(TITLE_ABS:colorectal) AND (TITLE_ABS:"single-cell" OR TITLE_ABS:"single cell")',
 "spatial": '(TITLE_ABS:colorectal) AND (TITLE_ABS:"spatial transcriptomics" OR TITLE_ABS:Visium OR TITLE_ABS:"spatially resolved")',
 "proteogenomics": '(TITLE_ABS:colorectal) AND (TITLE_ABS:proteogenomic OR TITLE_ABS:proteomic OR TITLE_ABS:phosphoproteomic)',
 "methylation": '(TITLE_ABS:colorectal) AND (TITLE_ABS:methylation OR TITLE_ABS:CIMP)',
 "immunemsi": '(TITLE_ABS:colorectal) AND (TITLE_ABS:"microsatellite" OR TITLE_ABS:MSI OR TITLE_ABS:immune)',
 "prognosis": '(TITLE_ABS:"colorectal cancer") AND (TITLE_ABS:prognosis OR TITLE_ABS:survival OR TITLE_ABS:recurrence)',
 "cms": '(TITLE_ABS:"consensus molecular subtypes")',
 "subtypeutility": '(TITLE_ABS:"molecular subtype" OR TITLE_ABS:"cancer subtype") AND (TITLE_ABS:"clinical utility" OR TITLE_ABS:"prognostic value" OR TITLE_ABS:reproducibility)',
}

def search(query, page_size=25, years="[2023 TO 2026]"):
    q = f"({query})"
    if years:
        q += f" AND ({' OR '.join('JOURNAL:' + j for j in JOURNALS)}) AND PUB_YEAR:{years}"
    url = BASE + "?" + urllib.parse.urlencode(
        {"query": q, "format": "json", "pageSize": page_size, "sort": "CITED desc", "resultType": "lite"})
    with urllib.request.urlopen(url, timeout=40) as fh:
        return json.load(fh).get("resultList", {}).get("result", [])

rows = []
for name, q in QUERIES.items():
    try:
        for r in search(q):
            rows.append({
                "query": name, "title": r.get("title", "").rstrip("."),
                "journal": r.get("journalTitle", ""), "year": r.get("pubYear", ""),
                "volume": r.get("journalVolume", ""), "pages": r.get("pageInfo", ""),
                "doi": r.get("doi", ""), "pmid": r.get("pmid", ""),
                "cites": r.get("citedByCount", 0),
                "authors": r.get("authorString", ""),
            })
        print(f"  {name}: {len(rows)} cumulative", flush=True)
        time.sleep(0.4)
    except Exception as exc:
        print(f"  {name}: FAILED {exc}", flush=True)

seen, out = set(), []
for r in sorted(rows, key=lambda x: -int(x["cites"] or 0)):
    if r["doi"].lower() in seen:
        continue
    seen.add(r["doi"].lower())
    out.append(r)
Path("Environment/bibliography/candidates.tsv").write_text(
    "\t".join(["query", "year", "journal", "title", "volume", "pages", "doi", "pmid", "cites", "authors"]) + "\n" +
    "\n".join("\t".join(str(r[k]) for k in ["query","year","journal","title","volume","pages","doi","pmid","cites","authors"]) for r in out),
    encoding="utf-8")
print("unique candidates:", len(out))
