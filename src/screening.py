"""Merge title/abstract screening decisions and build the full-text queue.

  python src/screening.py

Reads data/literature/records.json and data/literature/screening/decisions_batch_*.json.
Writes:
  data/literature/screening_decisions.json   every record with its decision
  data/literature/fulltext_queue.json        include + unsure records for full-text review
  data/literature/exclusion_review_sample.json  random exclusions for human review (seed fixed)
"""
import json
import random
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LIT = ROOT / "data" / "literature"
DECISIONS = {"include", "exclude", "unsure"}
CODES = {"I", "U", "E1", "E2", "E3", "E4", "E5", "E6"}
SEED = 20261029


def merge_decisions(records, decision_rows):
    by_id = {}
    for row in decision_rows:
        rid = row.get("record_id")
        if rid in by_id:
            raise ValueError(f"Duplicate decision for {rid}")
        if row.get("decision") not in DECISIONS or row.get("reason_code") not in CODES:
            raise ValueError(f"Invalid decision or code for {rid}")
        expected = {"include": {"I"}, "unsure": {"U"}, "exclude": CODES - {"I", "U"}}[row["decision"]]
        if row["reason_code"] not in expected:
            raise ValueError(f"Reason code {row['reason_code']} inconsistent with {row['decision']} for {rid}")
        by_id[rid] = row
    record_ids = {r["record_id"] for r in records}
    missing, extra = record_ids - set(by_id), set(by_id) - record_ids
    if missing or extra:
        raise ValueError(f"Decisions do not match records: {len(missing)} missing, {len(extra)} unknown")
    return [{**r, "screening": by_id[r["record_id"]]} for r in records]


def exclusion_sample(merged, fraction=0.2, minimum=10, seed=SEED):
    excluded = [m for m in merged if m["screening"]["decision"] == "exclude"]
    k = min(len(excluded), max(minimum, round(fraction * len(excluded))))
    return random.Random(seed).sample(excluded, k)


def write_sample(path, sample):
    """Write a human-review sample, keeping any answers already given (human_agrees, notes)."""
    kept = {}
    if path.exists():
        kept = {r["record_id"]: r for r in json.loads(path.read_text(encoding="utf-8"))}
    rows = []
    for s in sample:
        row = {"record_id": s["record_id"], "title": s["title"], "reason_code": s["screening"]["reason_code"],
               "reason": s["screening"]["reason"], "human_agrees": None}
        old = kept.get(s["record_id"], {})
        row.update({k: v for k, v in old.items() if k not in row or (k == "human_agrees" and v is not None)})
        rows.append(row)
    path.write_text(json.dumps(rows, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")


def apply_lead_checks(rows, checks):
    """Lead re-checks override screener rows; each names the screener decision it expects to replace."""
    by_id = {r["record_id"]: r for r in rows}
    for c in checks:
        row = by_id.get(c["record_id"])
        if row is None:
            raise ValueError(f"Lead check for unknown record {c['record_id']}")
        if (row["decision"], row["reason_code"]) != (c["screener_decision"], c["screener_reason_code"]):
            raise ValueError(f"Lead check for {c['record_id']} expects {c['screener_decision']}/{c['screener_reason_code']}, "
                             f"found {row['decision']}/{row['reason_code']}")
        row.update({"screener_decision": row["decision"], "screener_reason_code": row["reason_code"],
                    "screener_reason": row["reason"], "decision": c["decision"], "reason_code": c["reason_code"],
                    "reason": c["reason"], "lead_checked": True})
    return rows


def _rows(pattern):
    rows = []
    for path in sorted((LIT / "screening").glob(pattern)):
        rows.extend(json.loads(path.read_text(encoding="utf-8")))
    return rows


def main():
    records = json.loads((LIT / "records.json").read_text(encoding="utf-8"))
    merged = merge_decisions(records, _rows("decisions_batch_*.json"))
    # Other sources (reference lists, supplementary full-text search): own decision files and own review sample,
    # so the database-search sample offered to Artem never changes.
    supp_path = LIT / "records_supplementary.json"
    merged_supp = []
    if supp_path.exists():
        supp_rows = _rows("decisions_supp_*.json")
        checks_path = LIT / "screening" / "lead_checks_supp.json"
        if checks_path.exists():
            supp_rows = apply_lead_checks(supp_rows, json.loads(checks_path.read_text(encoding="utf-8")))
        merged_supp = merge_decisions(json.loads(supp_path.read_text(encoding="utf-8")), supp_rows)
        unchecked = [m for m in merged_supp if not m["screening"].get("lead_checked")]
        supp_sample = exclusion_sample(unchecked)
        write_sample(LIT / "exclusion_review_sample_supplementary.json", supp_sample)
    sample = exclusion_sample(merged)
    merged = merged + merged_supp
    queue = [m for m in merged if m["screening"]["decision"] in {"include", "unsure"}]
    (LIT / "screening_decisions.json").write_text(json.dumps(merged, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    (LIT / "fulltext_queue.json").write_text(json.dumps(queue, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    write_sample(LIT / "exclusion_review_sample.json", sample)
    print(json.dumps({"records": len(merged),
                      "decisions": dict(Counter(m["screening"]["decision"] for m in merged)),
                      "reason_codes": dict(Counter(m["screening"]["reason_code"] for m in merged).most_common()),
                      "fulltext_queue": len(queue), "exclusion_review_sample": len(sample)}, indent=1))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Screening merge failed: {exc}", file=sys.stderr)
        sys.exit(1)
