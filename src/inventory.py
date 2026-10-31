"""Merge full-text decisions and extracted claims into the claim inventory.

  python src/inventory.py

Reads data/literature/screening_decisions.json and data/literature/extraction/{fulltext_decisions,claims}_*.json,
then applies the logged edits in data/literature/corrections.json (raw extraction files are never edited).
Writes:
  data/literature/claims.json         validated, de-duplicated claim table (input to the frozen analysis)
  results/inventory_flow.json         counts at each step (PRISMA-style) with exclusion reasons
  results/inventory_summary.json      claim summary, genes claimed by several papers
  results/claims_review.csv           flat table for human review
"""
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from claims import summary, validate_table, write_csv  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
LIT = ROOT / "data" / "literature"
RESULTS = ROOT / "results"
DATABASE_SEARCHES = {"pubmed_mr", "europepmc_mr"}


def load_batches(prefix):
    rows = []
    for path in sorted((LIT / "extraction").glob(f"{prefix}_*.json")):
        rows.extend(json.loads(path.read_text(encoding="utf-8")))
    return rows


def apply_corrections(decisions, claims, corrections):
    """Apply logged corrections without touching the raw extraction files.

    Each correction names its target row, the values it expects to find, and the values to set.
    A correction whose expected values no longer match fails loudly instead of silently re-applying.
    """
    tables = {"decision": ({d["record_id"]: d for d in decisions}),
              "claim": ({c["claim_id"]: c for c in claims})}
    for corr in corrections:
        row = tables[corr["target"]].get(corr["id"])
        if row is None:
            raise ValueError(f"Correction target not found: {corr['target']} {corr['id']}")
        for field, expected in corr.get("expect", {}).items():
            if row.get(field) != expected:
                raise ValueError(f"Correction for {corr['id']} expected {field}={expected!r}, found {row.get(field)!r}")
        if not corr.get("why"):
            raise ValueError(f"Correction for {corr['id']} has no reason")
        row.update(corr["set"])
    return decisions, claims


def merge_instruments(claims, path):
    """Attach variants extracted by src/instruments.py (copied from supplement tables, never retyped)."""
    if not path.exists():
        return claims
    extracted = json.loads(path.read_text(encoding="utf-8"))
    by_id = {c["claim_id"]: c for c in claims}
    for cid, r in extracted.items():
        claim = by_id.get(cid)
        if claim is None:
            raise ValueError(f"instruments.json names unknown claim {cid}")
        claim["instruments_status"] = r["status"]
        claim["instruments_check"] = {"recompute": r.get("recompute"), "flags": r.get("flags") or [],
                                      "problems": r.get("problems") or [], "transcribed": r.get("transcribed", False),
                                      "reason": r.get("reason")}
        if r["status"] == "extracted":
            claim["instruments"] = r["variants"]
    return claims


def merge_rederived(claims, path):
    """Attach instruments rebuilt from public exposure files by src/rederive.py (kept apart from extracted ones)."""
    if not path.exists():
        return claims
    rebuilt = json.loads(path.read_text(encoding="utf-8"))["claims"]
    by_id = {c["claim_id"]: c for c in claims}
    for r in rebuilt:
        claim = by_id.get(r["claim_id"])
        if claim is None:
            raise ValueError(f"rederived.json names unknown claim {r['claim_id']}")
        if claim.get("instruments_status") == "extracted":
            raise ValueError(f"{r['claim_id']} already has extracted instruments; re-derivation is only for the rest")
        record = {k: r.get(k) for k in ("route", "dataset", "dataset_file", "dataset_sha256", "n_reported",
                                        "n_rebuilt", "analysis_role", "reasons", "reading_used")}
        check = dict(claim.get("instruments_check") or {})
        check["rederived"] = record
        claim["instruments_check"] = check
        if not r["variants"]:
            continue
        claim["instruments_status"] = "rederived"
        claim["instruments"] = [{"rsid": v["rsid"], "effect_allele": v["ea"], "other_allele": v["oa"],
                                 "beta_exposure": v["beta"], "se_exposure": v["se"], "pval_exposure": v["p"],
                                 "eaf_exposure": v.get("eaf"), "chrom": v.get("chrom"),
                                 "pos": v.get("pos", v.get("pos_GRCh38")),
                                 "build": "GRCh38" if "pos_GRCh38" in v else "GRCh37",
                                 **({"proxy_for": v["proxy_for"], "proxy_r2": v["proxy_r2"]} if v.get("proxy_for") else {}),
                                 "source": {"file": r["dataset_file"]}} for v in r["variants"]]
    return claims


def build(screened, decisions, claims):
    queue_ids = [r["record_id"] for r in screened if r["screening"]["decision"] in {"include", "unsure"}]
    by_id = {}
    for d in decisions:
        if d["record_id"] in by_id:
            raise ValueError(f"Duplicate full-text decision for {d['record_id']}")
        by_id[d["record_id"]] = d
    # Database-search records must all be decided. Records from other sources (reference lists, the
    # supplementary full-text search) may still be awaiting full-text review; they are counted, not dropped.
    from_db = {r["record_id"] for r in screened if set(r.get("found_by", DATABASE_SEARCHES)) & DATABASE_SEARCHES}
    missing = set(queue_ids) - set(by_id)
    if missing & from_db:
        raise ValueError(f"{len(missing & from_db)} queued records have no full-text decision: {sorted(missing & from_db)[:5]}")
    awaiting = sorted(missing)
    queue_ids = [r for r in queue_ids if r not in missing]
    included = {rid for rid, d in by_id.items() if d["final_decision"] == "include"}
    for rid, d in by_id.items():
        if d.get("duplicate_of") and d["final_decision"] == "include":
            raise ValueError(f"{rid} is marked duplicate but included")
        if d["final_decision"] == "include" and d["reason_code"] != "I":
            raise ValueError(f"{rid} included with reason {d['reason_code']}")

    # Claims are linked to records through DOI/PMID in paper_id.
    def record_for(claim):
        candidates = {claim["paper_id"]}
        if claim.get("doi"):
            candidates.add("doi:" + claim["doi"].lower())
        if claim.get("pmid"):
            candidates.add("pmid:" + str(claim["pmid"]))
        hits = candidates & set(by_id)
        return hits.pop() if len(hits) == 1 else None

    unlinked = [c["claim_id"] for c in claims if record_for(c) is None]
    if unlinked:
        raise ValueError(f"{len(unlinked)} claims do not link to exactly one record: {unlinked[:5]}")
    from_excluded = [c["claim_id"] for c in claims if record_for(c) not in included]
    if from_excluded:
        raise ValueError(f"Claims from excluded records: {from_excluded[:5]}")
    counts = Counter(record_for(c) for c in claims)
    mismatched = [rid for rid in included if counts.get(rid, 0) != by_id[rid]["n_claims"]]
    if mismatched:
        raise ValueError(f"n_claims mismatch for {mismatched[:5]}")
    problems = validate_table(claims)
    if problems:
        raise ValueError(f"{len(problems)} claims fail validation, e.g. {next(iter(problems.items()))}")

    def counts(subset):
        ids = {r["record_id"] for r in subset}
        q = [r for r in queue_ids if r in ids]
        inc = included & ids
        return {
            "records_identified": len(subset),
            "title_abstract_excluded": dict(Counter(r["screening"]["reason_code"] for r in subset
                                                    if r["screening"]["decision"] == "exclude")),
            "full_text_assessed": len(q),
            "full_text_accessed": sum(by_id[r]["full_text_accessed"] for r in q),
            "full_text_excluded": dict(Counter(by_id[r]["reason_code"] for r in q
                                               if by_id[r]["final_decision"] == "exclude")),
            "papers_included": len(inc),
            "papers_included_with_zero_claims": sum(1 for r in inc if by_id[r]["n_claims"] == 0),
        }

    # Database searches first (PRISMA "databases"); reference lists and the supplementary full-text search
    # are "other methods". A record found by both counts as a database record.
    database = [r for r in screened if set(r.get("found_by", DATABASE_SEARCHES)) & DATABASE_SEARCHES]
    other = [r for r in screened if not set(r.get("found_by", DATABASE_SEARCHES)) & DATABASE_SEARCHES]
    flow = {**counts(database), "claims": len(claims)}
    if other:
        flow["other_sources"] = {**counts(other), "awaiting_full_text_review": len(awaiting),
                                 "records_by_source": dict(Counter(" + ".join(sorted(r["found_by"])) for r in other))}
        flow["papers_included_all_sources"] = len(included)
    papers_per_gene = defaultdict(set)
    for c in claims:
        papers_per_gene[c["gene_symbol"]].add(c["paper_id"])
    multi = {g: sorted(p) for g, p in sorted(papers_per_gene.items(), key=lambda kv: -len(kv[1])) if len(p) > 1}
    return flow, {**summary(claims), "genes_claimed_by_several_papers": multi}


def main():
    screened = json.loads((LIT / "screening_decisions.json").read_text(encoding="utf-8"))
    decisions, claims = load_batches("fulltext_decisions"), load_batches("claims")
    corrections = []
    # Order matters: each stage expects the values left by the previous one.
    # corrections.json = first lead review; corrections_second.json = generated by src/reconcile.py;
    # corrections_lead.json = lead fixes made after reconciliation.
    for name in ("corrections.json", "corrections_second.json", "corrections_lead.json"):
        path = LIT / name
        if path.exists():
            corrections += json.loads(path.read_text(encoding="utf-8"))
    decisions, claims = apply_corrections(decisions, claims, corrections)
    merge_instruments(claims, LIT / "instruments" / "instruments.json")
    merge_rederived(claims, LIT / "instruments" / "rederive" / "rederived.json")
    flow, summ = build(screened, decisions, claims)
    flow["corrections_applied"] = len(corrections)
    RESULTS.mkdir(exist_ok=True)
    (LIT / "claims.json").write_text(json.dumps(claims, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    (RESULTS / "inventory_flow.json").write_text(json.dumps(flow, indent=1) + "\n", encoding="utf-8", newline="\n")
    (RESULTS / "inventory_summary.json").write_text(json.dumps(summ, indent=1) + "\n", encoding="utf-8", newline="\n")
    write_csv(claims, RESULTS / "claims_review.csv")
    print(json.dumps(flow, indent=1))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Inventory build failed: {exc}", file=sys.stderr)
        sys.exit(1)
