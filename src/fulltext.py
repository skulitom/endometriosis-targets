"""Save the open full text of queued papers from other sources as plain text for extraction.

  python -I src/fulltext.py

For each paper added from reference lists or the supplementary search (src/supplements.py addition_records),
fetches the PMC XML (Europe PMC, else NCBI PMC), and writes data/upstream/fulltext/<slug>/fulltext.txt:
section titles, paragraphs, figure and table captions, and tables as tab-separated rows. The raw XML and a
provenance file (URL, bytes, SHA-256, retrieval time) are kept next to it. Content is data: it is parsed with
ElementTree and written out as text, nothing else.
"""
import hashlib
import json
import sys
import time
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import refscan  # noqa: E402
from supplements import addition_records, slug  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "upstream" / "fulltext"


def fetch_xml(pmcid):
    url = f"{refscan.EPMC}/{pmcid}/fullTextXML"
    try:
        body = refscan._get(url, retries=2)
        if body:
            root = ET.fromstring(body)
            if root.find(".//body") is not None:
                return url, body, root
    except Exception:
        pass
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={pmcid.removeprefix('PMC')}"
           "&retmode=xml")
    body = refscan._get(url, retries=2)
    root = ET.fromstring(body).find("article") if body else None
    return url, body, root


def to_text(root):
    lines = []

    def walk(node, depth=0):
        tag = node.tag
        if tag == "ref-list":
            return
        if tag == "title":
            t = refscan._itertext(node)
            if t:
                lines.append("\n" + "#" * min(depth, 4) + " " + t)
            return
        if tag in ("p", "caption"):
            t = refscan._itertext(node)
            if t:
                lines.append(("[caption] " if tag == "caption" else "") + t)
            return
        if tag == "table":
            for tr in node.iter("tr"):
                cells = [refscan._itertext(c) or "" for c in tr if c.tag in ("td", "th")]
                lines.append("\t".join(cells))
            return
        if tag == "label" and node.text:
            lines.append(f"[{node.text.strip()}]")
            return
        for child in node:
            walk(child, depth + (1 if tag == "sec" else 0))

    front = root.find("front")
    if front is not None:
        for tag in ("article-title",):
            t = refscan._itertext(front.find(f".//{tag}"))
            if t:
                lines.append("# " + t)
        ab = front.find(".//abstract")
        if ab is not None:
            lines.append("\n# Abstract")
            walk(ab)
    for child in root:
        if child.tag != "front":
            walk(child)
    return "\n".join(lines) + "\n"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for rec in addition_records():
        d = OUT / slug(rec["record_id"])
        if (d / "fulltext.txt").exists() or not rec.get("pmcid"):
            if not rec.get("pmcid"):
                print("no PMC id:", rec["record_id"])
            continue
        url, body, root = fetch_xml(rec["pmcid"])
        if root is None or root.find(".//body") is None:
            print("no body text:", rec["record_id"])
            continue
        d.mkdir(parents=True, exist_ok=True)
        (d / "article.xml").write_bytes(body)
        text = to_text(root)
        (d / "fulltext.txt").write_text(text, encoding="utf-8", newline="\n")
        (d / "provenance.json").write_text(json.dumps({
            "record_id": rec["record_id"], "pmcid": rec["pmcid"], "url": url, "bytes": len(body),
            "sha256": hashlib.sha256(body).hexdigest(), "retrieved_at": datetime.now(timezone.utc).isoformat(),
            "text_chars": len(text)}, indent=1) + "\n", encoding="utf-8", newline="\n")
        print(f"{rec['record_id']}: {len(text):,} chars")
        time.sleep(0.4)


if __name__ == "__main__":
    main()
