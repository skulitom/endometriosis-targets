"""Run the stored literature searches and merge results into one de-duplicated record list.

  python src/literature.py search   # query PubMed and Europe PMC, save raw results with date and counts
  python src/literature.py merge    # de-duplicate saved results into data/literature/records.json
  python src/literature.py supplementary   # one-off Europe PMC full-text search (PROTOCOL.md section 1, Search)

Standard library only.
"""
import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LIT = ROOT / "data" / "literature"
UA = {"User-Agent": "EndometriosisTargetCheck/0.1 (research literature audit)"}

# Broad on purpose: every endometriosis MR study is screened, so gene-level claims are not missed.
QUERIES = {
    "pubmed_mr": {
        "source": "pubmed",
        "query": 'endometriosis[tiab] AND (mendelian[tiab] OR "summary-data-based"[tiab] '
                 'OR "drug target"[tiab] OR "drug-target"[tiab] OR druggable[tiab])',
    },
    "europepmc_mr": {
        "source": "europepmc",
        "query": '(TITLE:"endometriosis" OR ABSTRACT:"endometriosis") AND (TITLE:"mendelian" '
                 'OR ABSTRACT:"mendelian" OR ABSTRACT:"summary-data-based" OR ABSTRACT:"druggable" '
                 'OR ABSTRACT:"drug target" OR ABSTRACT:"drug-target")',
    },
}

# Supplementary search approved 2026-10-31: matches open full texts as well as titles and abstracts.
SUPPLEMENTARY = {
    "europepmc_fulltext": {
        "source": "europepmc",
        "query": '"endometriosis" AND ("summary-data-based" OR "summary data-based" OR HEIDI '
                 'OR "mendelian randomization" OR "mendelian randomisation")',
    },
}

# Terms suggesting a gene-level exposure. A screening hint only; every record is screened.
GENE_LEVEL_HINTS = re.compile(
    r"\b(pqtl|eqtl|sqtl|cis-|druggable|drug[- ]target|proteom|plasma protein|"
    r"summary-data-based|smr\b|coloc|gene expression|transcriptom|inhibitor|"
    r"lipid-lowering|antihypertensive|hmgcr|pcsk9|il-?6r)", re.IGNORECASE)


def _get(url, retries=3):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                return r.read()
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(2 * (attempt + 1))


def _text(node):
    return "".join(node.itertext()).strip() if node is not None else ""


def search_pubmed(query):
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    params = urllib.parse.urlencode({"db": "pubmed", "term": query, "retmax": 5000, "retmode": "json"})
    ids = json.loads(_get(f"{base}esearch.fcgi?{params}"))["esearchresult"]["idlist"]
    records = []
    for start in range(0, len(ids), 200):
        batch = ids[start:start + 200]
        params = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(batch), "retmode": "xml"})
        root = ET.fromstring(_get(f"{base}efetch.fcgi?{params}"))
        for article in root.findall(".//PubmedArticle"):
            ids_node = {i.get("IdType"): _text(i) for i in article.findall(".//ArticleIdList/ArticleId")}
            records.append({
                "pmid": _text(article.find(".//PMID")),
                "doi": ids_node.get("doi", "").lower() or None,
                "pmcid": ids_node.get("pmc") or None,
                "title": _text(article.find(".//ArticleTitle")),
                "abstract": " ".join(_text(a) for a in article.findall(".//Abstract/AbstractText")),
                "journal": _text(article.find(".//Journal/Title")),
                "year": _text(article.find(".//JournalIssue/PubDate/Year")) or _text(article.find(".//ArticleDate/Year")),
                "publication_types": [_text(p) for p in article.findall(".//PublicationType")],
                "preprint": False,
            })
        time.sleep(0.4)
    return records


def search_europepmc(query):
    base = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    records, cursor = [], "*"
    while True:
        params = urllib.parse.urlencode({"query": query, "format": "json", "resultType": "core",
                                         "pageSize": 1000, "cursorMark": cursor})
        data = json.loads(_get(f"{base}?{params}"))
        for item in data.get("resultList", {}).get("result", []):
            records.append({
                "pmid": item.get("pmid"),
                "doi": (item.get("doi") or "").lower() or None,
                "pmcid": item.get("pmcid"),
                "europepmc_id": f"{item.get('source')}:{item.get('id')}",
                "title": item.get("title", ""),
                "abstract": re.sub(r"<[^>]+>", " ", item.get("abstractText", "") or "").strip(),
                "journal": (item.get("journalInfo") or {}).get("journal", {}).get("title", "") or item.get("bookOrReportDetails", {}).get("publisher", ""),
                "year": item.get("pubYear"),
                "publication_types": (item.get("pubTypeList") or {}).get("pubType", []),
                "preprint": item.get("source") == "PPR",
            })
        nxt = data.get("nextCursorMark")
        if not nxt or nxt == cursor or not data.get("resultList", {}).get("result"):
            break
        cursor = nxt
    return records


def record_key(record):
    if record.get("doi"):
        return "doi:" + record["doi"].strip().lower()
    if record.get("pmid"):
        return "pmid:" + str(record["pmid"]).strip()
    if record.get("europepmc_id"):
        return "epmc:" + record["europepmc_id"]
    return "title:" + re.sub(r"\W+", " ", record.get("title", "")).strip().lower()


def merge(result_sets):
    """De-duplicate by DOI, then PMID. Keeps the longest abstract and every source that found it."""
    merged, pmid_index = {}, {}
    for source_name, records in result_sets:
        for rec in records:
            key = record_key(rec)
            if key not in merged and rec.get("pmid") and str(rec["pmid"]) in pmid_index:
                key = pmid_index[str(rec["pmid"])]
            if key in merged:
                existing = merged[key]
                existing["found_by"] = sorted(set(existing["found_by"]) | {source_name})
                for field, value in rec.items():
                    if value and not existing.get(field):
                        existing[field] = value
                if len(rec.get("abstract") or "") > len(existing.get("abstract") or ""):
                    existing["abstract"] = rec["abstract"]
            else:
                merged[key] = {**rec, "found_by": [source_name]}
            if rec.get("pmid"):
                pmid_index[str(rec["pmid"])] = key
    out = []
    for key, rec in sorted(merged.items()):
        text = f"{rec.get('title', '')} {rec.get('abstract', '')}"
        out.append({"record_id": key, **rec, "gene_level_hint": bool(GENE_LEVEL_HINTS.search(text))})
    return out


def cmd_search():
    LIT.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    for name, spec in QUERIES.items():
        runner = search_pubmed if spec["source"] == "pubmed" else search_europepmc
        records = runner(spec["query"])
        payload = {"name": name, "source": spec["source"], "query": spec["query"],
                   "run_at": datetime.now(timezone.utc).isoformat(), "count": len(records), "records": records}
        path = LIT / f"search_{name}_{stamp}.json"
        path.write_text(json.dumps(payload, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
        print(f"{name}: {len(records)} records -> {path.name}")


def cmd_supplementary():
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    for name, spec in SUPPLEMENTARY.items():
        path = LIT / f"search_{name}_{stamp}.json"
        if path.exists():
            raise FileExistsError(f"{path.name} exists; the supplementary search is run once")
        records = search_europepmc(spec["query"])
        payload = {"name": name, "source": spec["source"], "query": spec["query"],
                   "run_at": datetime.now(timezone.utc).isoformat(), "count": len(records), "records": records}
        path.write_text(json.dumps(payload, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
        print(f"{name}: {len(records)} records -> {path.name}")


def _norm_title(text):
    return re.sub(r"[^a-z0-9]+", " ", re.sub(r"<[^>]+>", " ", text or "").lower()).strip()


def supplementary_only(main, supplementary):
    """Supplementary records not already among the main records (by DOI, PMID or normalised title)."""
    seen = set()
    for r in main:
        seen |= {("doi", (r.get("doi") or "").lower()), ("pmid", str(r.get("pmid") or "")),
                 ("title", _norm_title(r.get("title")))}
    seen -= {("doi", ""), ("pmid", ""), ("title", "")}
    return [r for r in supplementary
            if not {("doi", (r.get("doi") or "").lower()), ("pmid", str(r.get("pmid") or "")),
                    ("title", _norm_title(r.get("title")))} & seen]


def cmd_merge():
    """records.json from the database searches only; other sources (reference lists, the supplementary
    full-text search) go to records_supplementary.json so the original records never change."""
    main_sets, other_sets = [], []
    for path in sorted(LIT.glob("search_*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        (main_sets if payload["name"] in QUERIES else other_sets).append((payload["name"], payload["records"]))
    records = merge(main_sets)
    (LIT / "records.json").write_text(json.dumps(records, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    hinted = sum(r["gene_level_hint"] for r in records)
    print(f"{len(records)} unique records ({hinted} with gene-level hint terms) -> records.json")
    if other_sets:
        extra = supplementary_only(records, merge(other_sets))
        (LIT / "records_supplementary.json").write_text(
            json.dumps(extra, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
        print(f"{len(extra)} records from other sources not in records.json -> records_supplementary.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["search", "merge", "supplementary"])
    args = parser.parse_args()
    try:
        {"search": cmd_search, "merge": cmd_merge, "supplementary": cmd_supplementary}[args.command]()
    except Exception as exc:
        print(f"Literature step failed: {exc}", file=sys.stderr)
        sys.exit(1)
