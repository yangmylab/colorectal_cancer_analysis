#!/usr/bin/env python3
"""List the docx body paragraphs with their style and opening text, so citations can be anchored by paragraph."""
import sys, zipfile
import xml.etree.ElementTree as ET
W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
NS = {"w": W}
root = ET.fromstring(zipfile.ZipFile(sys.argv[1]).read("word/document.xml"))
body = root.find("w:body", NS)
for i, p in enumerate(body.findall("w:p", NS)):
    style = p.find("./w:pPr/w:pStyle", NS)
    s = style.get("{%s}val" % W, "") if style is not None else ""
    text = "".join(t.text or "" for t in p.iter("{%s}t" % W))
    if not text.strip():
        continue
    print(f"{i:3} [{s:>6}] {text[:110]}")
