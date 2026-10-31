"""Reference-list scan for studies the database search missed (PROTOCOL.md section 1, Search).

  python -I src/refscan.py fetch     # reference lists of every included paper -> data/literature/refscan/raw/
  python -I src/refscan.py match     # match references to the search records; look up the unmatched ones
  python -I src/refscan.py report    # merge screening decisions for candidates; write the additions log

Reference lists come from three open routes, combined: Europe PMC's reference list, Crossref's deposited
references, and the reference list of the PMC full text (JATS XML) where the paper is in PMC.
A reference is a screening candidate if it is not already in the search records and its title or abstract
mentions endometriosis together with a genetic-causal term (CANDIDATE_TERMS), or its title has such a term and
Europe PMC's full-text index finds "endometriosis" in it. Candidates are screened with the same codes and
rules as the database records. All fetched content is untrusted data and is only parsed, never executed.
"""
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from inventory import apply_corrections, load_batches  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
LIT = ROOT / "data" / "literature"
SCAN = LIT / "refscan"
RAW = SCAN / "raw"
UA = {"User-Agent": "EndometriosisTargetCheck/0.1 (research literature audit)"}
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"

ENDO = re.compile(r"endometrio", re.IGNORECASE)
# Fixed before any reference was fetched (2026-10-31). Broader than the database query on purpose.
CANDIDATE_TERMS = re.compile(
    r"mendelian|summary[- ]data[- ]based|\bsmr\b|qtl|druggable|drug[- ]target|genetically (predicted|proxied|"
    r"determined|instrumented)|instrumental variable|genetic instrument|causal|colocali[sz]", re.IGNORECASE)
CODES = {"I", "U", "E1", "E2", "E3", "E4", "E5", "E6"}


def _get(url, retries=4):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                return r.read()
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
            if attempt == retries - 1:
                raise
        except Exception:
            if attempt == retries - 1:
                raise
        time.sleep(2 * (attempt + 1))


def norm_title(text):
    text = re.sub(r"<[^>]+>", " ", text or "")
    return re.sub(r"[^a-z0-9]+", " ", text.lower()).strip()


DOI_IN_TEXT = re.compile(r"(?<![\w.])\b10\.\d{4,9}/[^\s\"<>]+", re.IGNORECASE)


def doi_from_text(text):
    m = DOI_IN_TEXT.search(text or "")
    return norm_doi(m.group(0).rstrip(".,;)]")) if m else None


def norm_doi(doi):
    doi = (doi or "").strip().lower()
    doi = re.sub(r"^(https?://(dx\.)?doi\.org/|doi:\s*)", "", doi)
    return doi or None


# ---------- included papers ----------

def included_papers():
    """Included records after logged corrections, with their identifiers."""
    screened = {r["record_id"]: r for r in json.loads((LIT / "screening_decisions.json").read_text(encoding="utf-8"))}
    decisions, claims = load_batches("fulltext_decisions"), load_batches("claims")
    corrections = []
    for name in ("corrections.json", "corrections_second.json", "corrections_lead.json"):
        path = LIT / name
        if path.exists():
            corrections += json.loads(path.read_text(encoding="utf-8"))
    decisions, _ = apply_corrections(decisions, claims, corrections)
    out = []
    for d in decisions:
        if d["final_decision"] == "include":
            r = screened[d["record_id"]]
            out.append({k: r.get(k) for k in ("record_id", "pmid", "doi", "pmcid", "europepmc_id", "title", "year")})
    return sorted(out, key=lambda r: r["record_id"])


# ---------- fetching ----------

def refs_europepmc(paper):
    src_id = paper.get("europepmc_id") or (f"MED:{paper['pmid']}" if paper.get("pmid") else None)
    if not src_id:
        return None
    src, pid = src_id.split(":", 1)
    refs, page = [], 1
    while True:
        body = _get(f"{EPMC}/{src}/{pid}/references?page={page}&pageSize=1000&format=json")
        if body is None:
            return None if page == 1 else refs
        data = json.loads(body)
        batch = (data.get("referenceList") or {}).get("reference", [])
        for x in batch:
            refs.append({"title": x.get("title"), "doi": norm_doi(x.get("doi")),
                         "pmid": x.get("id") if x.get("source") == "MED" else None,
                         "year": x.get("pubYear"), "journal": x.get("journalAbbreviation"),
                         "authors": x.get("authorString"), "text": None})
        if len(refs) >= int(data.get("hitCount") or 0) or not batch:
            return refs
        page += 1


def refs_crossref(paper):
    if not paper.get("doi"):
        return None
    body = _get(f"https://api.crossref.org/works/{urllib.parse.quote(paper['doi'], safe='/')}")
    if body is None:
        return None
    msg = json.loads(body).get("message", {})
    if "reference" not in msg:
        return None
    return [{"title": x.get("article-title"), "doi": norm_doi(x.get("DOI")), "pmid": None,
             "year": x.get("year"), "journal": x.get("journal-title"), "authors": x.get("author"),
             "text": x.get("unstructured")} for x in msg["reference"]]


def _itertext(node):
    return re.sub(r"\s+", " ", "".join(node.itertext())).strip() if node is not None else None


def refs_jats(paper):
    if not paper.get("pmcid"):
        return None
    body = _get(f"{EPMC}/{paper['pmcid']}/fullTextXML")
    if body is None:
        return None
    root = ET.fromstring(body)
    refs = []
    for ref in root.iter("ref"):
        cit = next((c for c in ref if c.tag in ("element-citation", "mixed-citation", "citation")), ref)
        ids = {p.get("pub-id-type"): (p.text or "").strip() for p in cit.iter("pub-id")}
        refs.append({"title": _itertext(cit.find(".//article-title")), "doi": norm_doi(ids.get("doi")),
                     "pmid": ids.get("pmid") or None, "year": _itertext(cit.find(".//year")),
                     "journal": _itertext(cit.find(".//source")), "authors": None, "text": _itertext(cit)})
    return refs


def cmd_fetch():
    RAW.mkdir(parents=True, exist_ok=True)
    papers = included_papers()
    for p in papers:
        out = RAW / (re.sub(r"[^A-Za-z0-9._-]+", "_", p["record_id"]) + ".json")
        if out.exists():
            continue
        routes = {}
        for name, fn in (("europepmc", refs_europepmc), ("crossref", refs_crossref), ("pmc_jats", refs_jats)):
            try:
                refs = fn(p)
                routes[name] = {"status": "ok" if refs else "none", "references": refs or []}
            except Exception as exc:  # recorded, not hidden
                routes[name] = {"status": f"error: {type(exc).__name__}: {exc}"[:200], "references": []}
            time.sleep(0.3)
        out.write_text(json.dumps({"paper": p, "retrieved_at": datetime.now(timezone.utc).isoformat(),
                                   "routes": routes}, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
        print(p["record_id"], {k: len(v["references"]) for k, v in routes.items()})
    print(f"{len(papers)} included papers")


# ---------- matching ----------

def record_index(records):
    idx = {}
    for r in records:
        for key in (("doi", norm_doi(r.get("doi"))), ("pmid", str(r["pmid"]) if r.get("pmid") else None),
                    ("title", norm_title(r.get("title")) or None)):
            if key[1]:
                idx[key] = r["record_id"]
    return idx


def unique_references(raw_files):
    """Combine routes and papers into unique references (by DOI, then PMID, then normalised title)."""
    refs, alias = {}, {}
    for path in raw_files:
        data = json.loads(path.read_text(encoding="utf-8"))
        citing = data["paper"]["record_id"]
        for route, block in data["routes"].items():
            for x in block["references"]:
                if not x.get("doi") and x.get("text"):
                    x = {**x, "doi": doi_from_text(x["text"])}
                keys = [k for k in (("doi", x.get("doi")), ("pmid", x.get("pmid")),
                                    ("title", norm_title(x.get("title")) or None)) if k[1]]
                if not keys:
                    keys = [("text", norm_title(x.get("text")))]
                    if not keys[0][1]:
                        continue
                uid = next((alias[k] for k in keys if k in alias), None)
                if uid is None:
                    uid = f"{keys[0][0]}:{keys[0][1]}"
                    refs[uid] = {"ref_id": uid, "doi": None, "pmid": None, "title": None, "text": None,
                                 "year": None, "cited_by": set(), "routes": set()}
                r = refs[uid]
                for f in ("doi", "pmid", "title", "text", "year"):
                    if x.get(f) and not r[f]:
                        r[f] = x[f]
                r["cited_by"].add(citing)
                r["routes"].add(route)
                for k in keys:
                    alias.setdefault(k, uid)
                for k in (("doi", r["doi"]), ("pmid", r["pmid"]), ("title", norm_title(r["title"]) or None)):
                    if k[1]:
                        alias.setdefault(k, uid)
    return refs


def match_record(ref, idx):
    for key in (("doi", ref.get("doi")), ("pmid", ref.get("pmid")), ("title", norm_title(ref.get("title")) or None)):
        if key[1] and key in idx:
            return idx[key]
    return None


CACHE_PATH = SCAN / "europepmc_lookups.json"
_cache = None


def epmc_lookup(query, page_size=100):
    """Europe PMC search; results kept in a dated cache so reruns repeat no queries and the lookups are on file."""
    global _cache
    if _cache is None:
        _cache = json.loads(CACHE_PATH.read_text(encoding="utf-8")) if CACHE_PATH.exists() else {}
    key = f"{page_size}|{query}"
    if key not in _cache:
        params = urllib.parse.urlencode({"query": query, "format": "json", "resultType": "core", "pageSize": page_size})
        _cache[key] = {"at": datetime.now(timezone.utc).isoformat(),
                       "result": json.loads(_get(f"{EPMC}/search?{params}")).get("resultList", {}).get("result", [])}
        time.sleep(0.2)
    return _cache[key]["result"]


def save_cache():
    if _cache is not None:
        SCAN.mkdir(parents=True, exist_ok=True)
        CACHE_PATH.write_text(json.dumps(_cache, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")


def _clean(item):
    return {"europepmc_id": f"{item.get('source')}:{item.get('id')}", "pmid": item.get("pmid"),
            "doi": norm_doi(item.get("doi")), "pmcid": item.get("pmcid"), "title": item.get("title", ""),
            "abstract": re.sub(r"<[^>]+>", " ", item.get("abstractText", "") or "").strip(),
            "year": item.get("pubYear"), "journal": (item.get("journalInfo") or {}).get("journal", {}).get("title"),
            "publication_types": (item.get("pubTypeList") or {}).get("pubType", []),
            "preprint": item.get("source") == "PPR"}


def resolve(refs):
    """Look up each unmatched reference in Europe PMC (by DOI or PMID in batches; by title otherwise)."""
    found = {}
    by_doi = {r["doi"]: uid for uid, r in refs.items() if r.get("doi")}
    by_pmid = {r["pmid"]: uid for uid, r in refs.items() if r.get("pmid") and uid not in by_doi.values()}
    for field, mapping in (("DOI", by_doi), ("EXT_ID", by_pmid)):
        keys = sorted(mapping)
        for i in range(0, len(keys), 40):
            chunk = keys[i:i + 40]
            q = " OR ".join(f'{field}:"{k}"' for k in chunk)
            if field == "EXT_ID":
                q = f"({q}) AND SRC:MED"
            for item in epmc_lookup(q):
                c = _clean(item)
                key = c["doi"] if field == "DOI" else c["pmid"]
                if key in mapping:
                    found.setdefault(mapping[key], c)
            time.sleep(0.3)
    for uid, r in refs.items():
        if uid in found or not r.get("title") or len(norm_title(r["title"])) < 20:
            continue
        title = re.sub(r'["\\]', " ", r["title"])[:250]
        for item in epmc_lookup(f'TITLE:"{title}"', page_size=5):
            c = _clean(item)
            if norm_title(c["title"]) == norm_title(r["title"]):
                found[uid] = c
                break
        time.sleep(0.2)
    return found


# Added 2026-10-31 after the first match run, before screening: Rahmioglu 2023 reports gene-level SMR only in
# its full text, so endometriosis references are also checked for MR or SMR anywhere in their open full text.
FULLTEXT_MR = ('("mendelian randomization" OR "mendelian randomisation" OR "summary-data-based" OR SMR '
               'OR "drug target" OR "drug-target")')


def fulltext_mentions_endometriosis(ids, clause="endometriosis"):
    """Europe PMC ids whose indexed text (full text where open) matches the clause."""
    hits = set()
    ids = sorted(ids)
    for i in range(0, len(ids), 40):
        chunk = ids[i:i + 40]
        q = "(" + " OR ".join(f'(SRC:{s} AND EXT_ID:{e})' for s, e in (x.split(":", 1) for x in chunk)) + f") AND {clause}"
        for item in epmc_lookup(q):
            hits.add(f"{item.get('source')}:{item.get('id')}")
        time.sleep(0.3)
    return hits


def cmd_match():
    records = json.loads((LIT / "records.json").read_text(encoding="utf-8"))
    idx = record_index(records)
    raw_files = sorted(RAW.glob("*.json"))
    refs = unique_references(raw_files)
    matched = {uid: match_record(r, idx) for uid, r in refs.items()}
    unmatched = {uid: r for uid, r in refs.items() if not matched[uid]}
    found = resolve(unmatched)
    # A resolved reference may turn out to be a search record after all (identifier missing in the citation).
    for uid, c in found.items():
        rid = match_record(c, idx)
        if rid:
            matched[uid] = rid
    unmatched = {uid: r for uid, r in refs.items() if not matched[uid]}
    title_terms = {uid for uid in unmatched
                   if CANDIDATE_TERMS.search((found.get(uid) or {}).get("title") or unmatched[uid].get("title")
                                             or unmatched[uid].get("text") or "")}
    ft_hits = fulltext_mentions_endometriosis({found[u]["europepmc_id"] for u in title_terms if u in found})
    endo_only = {found[u]["europepmc_id"] for u in unmatched if u in found and u not in title_terms
                 and ENDO.search(" ".join(filter(None, [found[u].get("title"), found[u].get("abstract")])))}
    ft_mr_hits = fulltext_mentions_endometriosis(endo_only, FULLTEXT_MR)
    save_cache()
    # References that resolve to the same Europe PMC record are one study (cited in different styles).
    groups = defaultdict(list)
    for uid in sorted(unmatched):
        groups[found[uid]["europepmc_id"] if uid in found else uid].append(uid)
    candidates, rest = [], Counter()
    for gid, uids in sorted(groups.items()):
        uid, r = uids[0], unmatched[uids[0]]
        r = {**r, "cited_by": set().union(*(unmatched[u]["cited_by"] for u in uids)),
             "routes": set().union(*(unmatched[u]["routes"] for u in uids))}
        meta = found.get(uid)
        text = " ".join(filter(None, [(meta or {}).get("title"), (meta or {}).get("abstract")]
                               + [unmatched[u].get(f) for u in uids for f in ("title", "text")]))
        endo, terms = bool(ENDO.search(text)), bool(CANDIDATE_TERMS.search(text))
        ft = bool(meta and meta["europepmc_id"] in ft_hits)
        ft_mr = bool(meta and meta["europepmc_id"] in ft_mr_hits)
        if (endo and terms) or (any(u in title_terms for u in uids) and ft) or (endo and ft_mr):
            candidates.append({"ref_id": uid, "same_study_as": uids[1:], "cited_by": sorted(r["cited_by"]), "routes": sorted(r["routes"]),
                               "why_candidate": "title/abstract: endometriosis + genetic-causal term" if endo and terms
                               else "title has genetic-causal term; Europe PMC full text mentions endometriosis" if ft
                               else "endometriosis in title/abstract; Europe PMC full text mentions MR, SMR or drug target",
                               "citation": {k: r.get(k) for k in ("doi", "pmid", "title", "text", "year")},
                               "europepmc": meta})
        else:
            rest["endometriosis, no genetic-causal term" if endo else
                 "genetic-causal term, endometriosis not found" if terms else "neither"] += 1
    SCAN.mkdir(parents=True, exist_ok=True)
    (SCAN / "candidates.json").write_text(json.dumps(candidates, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    per_paper = {}
    for path in raw_files:
        data = json.loads(path.read_text(encoding="utf-8"))
        per_paper[data["paper"]["record_id"]] = {k: (v["status"], len(v["references"])) for k, v in data["routes"].items()}
    summary = {
        "run_at": datetime.now(timezone.utc).isoformat(),
        "included_papers": len(raw_files),
        "papers_without_any_reference_list": sorted(p for p, v in per_paper.items() if not any(n for _, n in v.values())),
        "unique_references": len(refs),
        "already_in_search_records": sum(1 for v in matched.values() if v),
        "not_in_search_records": len(unmatched),
        "distinct_studies_not_in_search_records": len(groups),
        "resolved_in_europepmc": sum(1 for u in unmatched if u in found),
        "not_candidates": dict(rest),
        "candidates": len(candidates),
        "per_paper_routes": per_paper,
    }
    (SCAN / "match_summary.json").write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "per_paper_routes"}, indent=1))


SNIPPET_TERMS = re.compile(r"mendelian|summary[- ]data[- ]based|\bsmr\b|heidi|drug[- ]target|instrumental variable|"
                           r"qtl|colocali[sz]", re.IGNORECASE)


def cmd_snippets():
    """Passages and table captions naming MR, SMR or QTL terms, from each candidate's open PMC full text.

    Screening evidence only: saved to refscan/fulltext_snippets.json so each decision can be checked."""
    candidates = json.loads((SCAN / "candidates.json").read_text(encoding="utf-8"))
    out = {}
    for c in candidates:
        pmcid = (c.get("europepmc") or {}).get("pmcid")
        if not pmcid:
            out[c["ref_id"]] = {"pmcid": None, "snippets": []}
            continue
        try:
            body = _get(f"{EPMC}/{pmcid}/fullTextXML")
        except urllib.error.HTTPError as exc:
            out[c["ref_id"]] = {"pmcid": pmcid, "full_text": f"error HTTP {exc.code}", "snippets": []}
            continue
        time.sleep(0.3)
        if body is None:
            out[c["ref_id"]] = {"pmcid": pmcid, "full_text": "not available", "snippets": []}
            continue
        root = ET.fromstring(body)
        snippets = []
        for node in root.iter():
            if node.tag in ("p", "caption", "title", "label"):
                text = _itertext(node)
                if text and SNIPPET_TERMS.search(text) and text not in snippets:
                    snippets.append(text[:1500])
        out[c["ref_id"]] = {"pmcid": pmcid, "full_text": "europepmc_xml", "snippets": snippets}
    (SCAN / "fulltext_snippets.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print({k: len(v["snippets"]) for k, v in out.items()})


MR_PASSAGE = re.compile(r"mendelian|summary[- ]data[- ]based|\bheidi\b|\bsmr\b", re.IGNORECASE)
# Own-analysis wording, used only to pick exclusions for a lead re-check (PROTOCOL.md section 1, Search).
OWN_ANALYSIS = re.compile(
    r"\b(we|this study|the present study|our study|our analysis)\b[^.]{0,150}?\b(perform|conduct|appl|us|carri|undert|implement)"
    r"\w*[^.]{0,150}?(mendelian|\bsmr\b|summary[- ]data)|(mendelian|\bsmr\b|summary[- ]data)[^.]{0,120}?\b(was|were)\s+"
    r"(performed|conducted|applied|used|carried out|undertaken)", re.IGNORECASE)


OUTCOME_WORDS = re.compile(r"outcome|gwas|finngen|n14|cases|odds ratio|\bor\b|associat|causal|mendelian|\bsmr\b|"
                           r"phenome|ukb|biobank", re.IGNORECASE)


def _ncbi_pmc_xml(pmcid):
    """Open-access PMC XML from NCBI E-utilities; used when Europe PMC has no XML for the article."""
    body = _get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={pmcid.removeprefix('PMC')}"
                "&retmode=xml", retries=2)
    if not body:
        return None
    article = ET.fromstring(body).find("article")
    return article


def mr_passages(pmcid, route="europepmc"):
    """(status, MR passages, endometriosis sentences) for one PMC full text.

    MR passages: paragraphs and captions naming MR, SMR or HEIDI. Endometriosis sentences: up to 8 sentences
    (body text, captions and table cells) that mention endometriosis, to show what role it plays."""
    try:
        if route == "europepmc":
            body = _get(f"{EPMC}/{pmcid}/fullTextXML")
            root = ET.fromstring(body) if body else None
        else:
            root = _ncbi_pmc_xml(pmcid)
            time.sleep(0.4)
    except urllib.error.HTTPError as exc:
        return f"error HTTP {exc.code}", [], []
    except (urllib.error.URLError, OSError) as exc:
        return f"error network {type(exc).__name__}", [], []
    except ET.ParseError:
        return "unparseable", [], []
    if root is None:
        return "not available", [], []
    # Body, back matter (minus references) and floating tables; the abstract is already in the record.
    parts = []
    for child in root:
        if child.tag == "front":
            continue
        if child.tag == "back":
            parts.extend(c for c in child if c.tag != "ref-list")
        else:
            parts.append(child)
    out, endo = [], []
    for part in parts:
        for node in part.iter():
            if node.tag not in ("p", "caption", "title", "td", "th"):
                continue
            text = _itertext(node)
            if not text:
                continue
            if node.tag not in ("td", "th") and MR_PASSAGE.search(text) and text[:1500] not in out:
                out.append(text[:1500])
            for sent in re.split(r"(?<=[.;])\s+", text):
                if ENDO.search(sent):
                    rank = 0 if node.tag in ("td", "th") else 1 if OUTCOME_WORDS.search(sent) else 2
                    endo.append((rank, len(endo), sent[:350]))
    seen, picked = set(), []
    for _, _, sent in sorted(endo):
        if sent not in seen:
            seen.add(sent)
            picked.append(sent)
        if len(picked) == 8:
            break
    return ("europepmc_xml" if route == "europepmc" else "ncbi_pmc_xml"), out, picked


def _passage_row(r, route="europepmc"):
    status, passages, endo = mr_passages(r["pmcid"], route) if r.get("pmcid") else ("no PMC full text", [], [])
    text = " ".join(passages)
    return r["record_id"], {"pmcid": r.get("pmcid"), "status": status, "passages": passages,
                            "endometriosis_sentences": endo,
                            "own_analysis_with_endometriosis": bool(
                                any(OWN_ANALYSIS.search(p) for p in passages) and ENDO.search(text))}


def cmd_passages():
    """Full-text MR passages for supplementary-search records that still need screening (4 parallel requests)."""
    from concurrent.futures import ThreadPoolExecutor
    path = SCAN / "supplementary_passages.json"
    done = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    recs = json.loads((LIT / "records_supplementary.json").read_text(encoding="utf-8"))
    todo = [r for r in recs if r["found_by"] == ["europepmc_fulltext"] and r["record_id"] not in done]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for i, (rid, row) in enumerate(pool.map(_passage_row, todo)):
            done[rid] = row
            if i % 100 == 99:
                path.write_text(json.dumps(done, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
                print(f"{i + 1}/{len(todo)}", flush=True)
    # Europe PMC has no XML for some PMC articles (HTTP 500); NCBI serves the same open-access XML.
    by_id = {r["record_id"]: r for r in recs}
    for rid, row in list(done.items()):
        if row["status"].startswith("error") and by_id[rid].get("pmcid"):
            done[rid] = _passage_row(by_id[rid], route="ncbi")[1]
    path.write_text(json.dumps(done, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"records": len(done), "status": dict(Counter(v["status"] for v in done.values())),
                      "own_analysis_with_endometriosis": sum(v["own_analysis_with_endometriosis"] for v in done.values())}))


# ---------- screening merge ----------

def merge_decisions(candidates, decisions):
    by_id = {}
    for d in decisions:
        if d["ref_id"] in by_id:
            raise ValueError(f"Duplicate decision for {d['ref_id']}")
        expected = {"include": {"I"}, "unsure": {"U"}, "exclude": CODES - {"I", "U"}}.get(d.get("decision"))
        if not expected or d.get("reason_code") not in expected or not d.get("reason"):
            raise ValueError(f"Invalid decision for {d['ref_id']}")
        by_id[d["ref_id"]] = d
    ids = {c["ref_id"] for c in candidates}
    if ids != set(by_id):
        raise ValueError(f"Decisions do not match candidates: {len(ids - set(by_id))} missing, {len(set(by_id) - ids)} unknown")
    return [{**c, "screening": by_id[c["ref_id"]]} for c in candidates]


def cmd_report():
    candidates = json.loads((SCAN / "candidates.json").read_text(encoding="utf-8"))
    decisions = json.loads((SCAN / "decisions.json").read_text(encoding="utf-8"))
    merged = merge_decisions(candidates, decisions)
    additions = [{"ref_id": m["ref_id"], "title": (m["europepmc"] or {}).get("title") or m["citation"]["title"],
                  "doi": (m["europepmc"] or {}).get("doi") or m["citation"]["doi"],
                  "pmid": (m["europepmc"] or {}).get("pmid") or m["citation"]["pmid"],
                  "found_via_reference_list_of": m["cited_by"], "decision": m["screening"]["decision"],
                  "reason": m["screening"]["reason"]}
                 for m in merged if m["screening"]["decision"] in ("include", "unsure")]
    out = {"note": "Reference-list scan (PROTOCOL.md section 1). Studies not found by the database search; each "
                   "addition names the included paper(s) whose reference list it came from.",
           "candidates_screened": len(merged),
           "decisions": dict(Counter(m["screening"]["decision"] for m in merged)),
           "reason_codes": dict(Counter(m["screening"]["reason_code"] for m in merged)),
           "additions_for_full_text": additions}
    # Every screened candidate becomes a record from "other sources" (merged by src/literature.py merge).
    summary = json.loads((SCAN / "match_summary.json").read_text(encoding="utf-8"))
    recs = []
    for m in merged:
        meta = m["europepmc"] or {"europepmc_id": None, "pmid": m["citation"]["pmid"], "doi": m["citation"]["doi"],
                                  "pmcid": None, "title": m["citation"]["title"] or m["citation"]["text"],
                                  "abstract": "", "year": m["citation"]["year"], "journal": None,
                                  "publication_types": [], "preprint": False}
        recs.append({**meta, "cited_by_included": m["cited_by"], "refscan_id": m["ref_id"]})
    (LIT / f"search_reference_lists_{summary['run_at'][:10]}.json").write_text(json.dumps(
        {"name": "reference_lists", "source": "reference lists of included papers (src/refscan.py)",
         "query": "candidate rules (a)-(c), PROTOCOL.md section 1 (Search)", "run_at": summary["run_at"],
         "count": len(recs), "records": recs}, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    # The same decisions keyed by record_id, for src/screening.py (screened by the lead, so lead_checked).
    from literature import record_key  # noqa: E402
    rows = [{"record_id": record_key(rec), "decision": m["screening"]["decision"],
             "reason_code": m["screening"]["reason_code"], "reason": m["screening"]["reason"],
             "evidence": m["screening"].get("evidence"), "source": "reference_lists", "lead_checked": True}
            for rec, m in zip(recs, merged)]
    (LIT / "screening" / "decisions_supp_refscan.json").write_text(
        json.dumps(rows, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    (SCAN / "screened.json").write_text(json.dumps(merged, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    (SCAN / "additions.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({k: v for k, v in out.items() if k != "note"}, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    cmds = {"fetch": cmd_fetch, "match": cmd_match, "snippets": cmd_snippets, "report": cmd_report,
            "passages": cmd_passages}
    if len(sys.argv) != 2 or sys.argv[1] not in cmds:
        sys.exit(__doc__)
    cmds[sys.argv[1]]()
