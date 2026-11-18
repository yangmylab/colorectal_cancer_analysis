#!/usr/bin/env python3
"""Insert numbered superscript citations and a References section into the docx by direct XML surgery.

Body paragraphs are addressed by their index among <w:p> elements. Superscripts are appended as a new run at the
end of each target paragraph, and the References section is inserted immediately before the final <w:sectPr>.
"""
import json, re, shutil, zipfile
from pathlib import Path

DOCX = Path("12_Manuscript/manuscript.docx")
plan = json.loads(Path("Environment/bibliography/plan.json").read_text())

def clean(s):
    s = re.sub(r"<[^>]+>", "", s)
    return (s.replace("&lt;", "<").replace("&gt;", ">").replace("&amp;", "&").rstrip("."))

xml = zipfile.ZipFile(DOCX).read("word/document.xml").decode("utf-8")
paras = list(re.finditer(r"<w:p\b[^>]*>.*?</w:p>", xml, re.S))
print("paragraphs found:", len(paras))

def sup_run(nums):
    txt = ",".join(str(n) for n in nums)
    return ('<w:r><w:rPr><w:vertAlign w:val="superscript"/></w:rPr>'
            f'<w:t xml:space="preserve">{txt}</w:t></w:r>')

# 1. insert superscript runs, last paragraph first so offsets stay valid
edits = []
for pidx, dois in plan["plan"]:
    nums = [plan["numbers"][d] for d in dois if d in plan["numbers"]]
    if not nums:
        continue
    edits.append((paras[pidx].end() - len("</w:p>"), sup_run(nums)))
for pos, run in sorted(edits, key=lambda x: -x[0]):
    xml = xml[:pos] + run + xml[pos:]
print("superscripts inserted in", len(edits), "paragraphs")

# 2. build the References section
NS = 'xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"'
def para(style, text):
    return (f'<w:p><w:pPr><w:pStyle w:val="{style}"/></w:pPr><w:r><w:t xml:space="preserve">'
            f'{text}</w:t></w:r></w:p>')

entries = []
for doi, num in sorted(plan["numbers"].items(), key=lambda kv: kv[1]):
    entries.append(para("3", f'{num}. {clean(plan["entries"][doi])}'))
section = para("4", "References") + "".join(entries)

m = re.search(r"<w:sectPr\b", xml)
if m:
    xml = xml[:m.start()] + section + xml[m.start():]
else:
    xml = xml.replace("</w:body>", section + "</w:body>")
print("references appended:", len(entries))

shutil.copy(DOCX, str(DOCX) + ".bak")
z = zipfile.ZipFile(str(DOCX) + ".bak")
items = [(i, z.read(i.filename)) for i in z.infolist()]
z.close()
with zipfile.ZipFile(DOCX, "w", zipfile.ZIP_DEFLATED) as out:
    for info, data in items:
        if info.filename == "word/document.xml":
            data = xml.encode("utf-8")
        out.writestr(info, data)
print("written")
