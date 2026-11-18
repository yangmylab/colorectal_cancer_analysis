#!/usr/bin/env python3
"""Second pass: claim-specific Europe PMC queries plus exact-title lookups for framework/method originals."""
import json, time, urllib.parse, urllib.request, csv
from pathlib import Path

BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
HI = ['"Nature"','"Nature medicine"','"Nature cancer"','"Nature genetics"','"Nature communications"','"Cancer cell"','"Cell"','"Gut"','"Journal of clinical oncology"','"Annals of oncology"','"Gastroenterology"','"Lancet oncology"','"JAMA oncology"','"Cancer discovery"','"Clinical cancer research"','"Nature reviews. Cancer"','"Nature reviews. Clinical oncology"','"Journal of pathology"','"British journal of cancer"','"Molecular oncology"','"Cell reports medicine"','"Cancer research"','"Genome biology"','"Genome medicine"','"Nature biotechnology"','"Nature methods"','"Nucleic acids research"','"BMC bioinformatics"','"Cell systems"']

def run(q, years="[2023 TO 2026]", hi=True, size=20):
    full = q
    if hi:
        full += " AND (" + " OR ".join("JOURNAL:" + j for j in HI) + ")"
    if years:
        full += " AND PUB_YEAR:" + years
    url = BASE + "?" + urllib.parse.urlencode({"query": full, "format":"json","pageSize":size,
                                               "sort":"CITED desc","resultType":"lite"})
    with urllib.request.urlopen(url, timeout=40) as fh:
        return json.load(fh).get("resultList", {}).get("result", [])

TOPIC = {
 "cms_mesenchymal": '(TITLE_ABS:colorectal) AND (TITLE_ABS:CMS4 OR TITLE_ABS:"mesenchymal subtype" OR TITLE_ABS:"epithelial-mesenchymal")',
 "intrinsic_subtype": '(TITLE_ABS:colorectal) AND (TITLE_ABS:intrinsic AND TITLE_ABS:subtype)',
 "multiomics": '(TITLE_ABS:colorectal) AND (TITLE_ABS:"multi-omics" OR TITLE_ABS:integrative OR TITLE_ABS:"multiomic")',
 "stage_beyond": '(TITLE_ABS:colorectal) AND (TITLE_ABS:prognostic AND (TITLE_ABS:stage OR TITLE_ABS:"beyond"))',
 "subtype_repro": '(TITLE_ABS:cancer) AND (TITLE_ABS:subtype AND (TITLE_ABS:reproducibility OR TITLE_ABS:robustness OR TITLE_ABS:transferability OR TITLE_ABS:benchmark))',
 "stroma_emt": '(TITLE_ABS:colorectal) AND (TITLE_ABS:stroma OR TITLE_ABS:fibroblast OR TITLE_ABS:EMT)',
 "spatial_crc": '(TITLE_ABS:colorectal) AND (TITLE_ABS:"spatial transcriptomics" OR TITLE_ABS:Visium)',
 "singlecell_crc": '(TITLE_ABS:colorectal) AND (TITLE_ABS:"single-cell" OR TITLE_ABS:scRNA)',
 "phospho_crc": '(TITLE_ABS:colorectal) AND (TITLE_ABS:phosphoproteomic OR TITLE_ABS:kinase)',
 "meth_cimp": '(TITLE_ABS:colorectal) AND (TITLE_ABS:methylation OR TITLE_ABS:CIMP)',
 "cms_general": '(TITLE_ABS:"consensus molecular subtypes")',
 "prognosis_crc": '(TITLE_ABS:"colorectal cancer") AND (TITLE_ABS:prognosis OR TITLE_ABS:"prognostic")',
}
EXACT = {
 "guinney2015": 'TITLE:"consensus molecular subtypes of colorectal cancer"',
 "isella2017": 'TITLE:"cancer-cell intrinsic transcriptional traits" AND TITLE:colorectal',
 "phipps2015": 'TITLE:"prognostic subtypes" AND TITLE:colorectal',
 "gsva": 'TITLE:"GSVA" AND AUTH:"Hänzelmann"',
 "limma": 'TITLE:"limma powers differential expression"',
 "snf": 'TITLE:"Similarity network fusion"',
 "cell2location": 'TITLE:"cell2location"',
 "msigdb": 'TITLE:"The Molecular Signatures Database"',
 "cptac_crc": 'TITLE:"Proteogenomic Analysis of Human Colon Cancer"',
 "tcga_crc": 'TITLE:"Comprehensive molecular characterization of human colon and rectal cancer"',
 "ksea": 'TITLE:"Kinase Substrate Enrichment Analysis" AND AUTH:"Wiredja"',
}

rows=[]
def add(tag, r):
    rows.append({"query":tag,"title":(r.get("title") or "").rstrip("."),"journal":r.get("journalTitle",""),
                 "year":r.get("pubYear",""),"volume":r.get("journalVolume",""),"pages":r.get("pageInfo",""),
                 "doi":r.get("doi",""),"pmid":r.get("pmid",""),"cites":r.get("citedByCount",0),
                 "authors":r.get("authorString","")})

for tag,q in TOPIC.items():
    try:
        for r in run(q): add(tag,r)
    except Exception as e: print(tag,"FAIL",e)
    time.sleep(0.3)
for tag,q in EXACT.items():
    try:
        res = run(q, years=None, hi=False, size=3)
        if res: add(tag,res[0])
        else: print(tag,"no hit")
    except Exception as e: print(tag,"FAIL",e)
    time.sleep(0.3)

seen=set(); out=[]
for r in rows:
    k=r["doi"].lower()
    if not k or k in seen: continue
    seen.add(k); out.append(r)
Path("Environment/bibliography/candidates2.tsv").write_text(
  "\t".join(["query","year","journal","title","volume","pages","doi","pmid","cites","authors"])+"\n"+
  "\n".join("\t".join(str(r[k]) for k in ["query","year","journal","title","volume","pages","doi","pmid","cites","authors"]) for r in out), encoding="utf-8")
print("rows:", len(out))
