"""Compare the blind second extraction with the first and list disagreements.

  python src/second.py

Reads the first extraction with first-stage corrections only (so re-running after reconciliation
gives the same comparison) and data/literature/second/second_*.json. Writes results/second_extraction_comparison.json.
Agreement rules are fixed here so they cannot drift after looking at the results.
"""
import json
import math
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LIT = ROOT / "data" / "literature"
RESULTS = ROOT / "results"

NUMERIC = ["estimate", "ci_low", "ci_high"]
EXACT = ["estimate_scale", "direction", "instrument_scope", "coloc_reported", "headline"]
SOFT = ["outcome_release", "method"]  # reported, not counted as disagreement
COUNTS = ["outcome_cases", "outcome_controls", "n_instruments"]


def same_number(a, b):
    """Equal within rounding: absolute 0.006 or 1% relative."""
    if a is None or b is None:
        return a is None and b is None
    return abs(a - b) <= max(0.006, 0.01 * abs(a))


def same_p(a, b):
    """Equal if both missing, or within a factor of 1.15 (rounding of mantissa)."""
    if a is None or b is None:
        return a is None and b is None
    if a == b:
        return True
    if a <= 0 or b <= 0:
        return False
    return abs(math.log(a) - math.log(b)) <= math.log(1.15)


def compare(first, second):
    diffs = {}
    for f in NUMERIC:
        if not same_number(first.get(f), second.get(f)):
            diffs[f] = [first.get(f), second.get(f)]
    if not same_p(first.get("p_value"), second.get("p_value")):
        diffs["p_value"] = [first.get("p_value"), second.get("p_value")]
    for f in EXACT:
        if f in second and first.get(f) != second.get(f):
            diffs[f] = [first.get(f), second.get(f)]
    if "outcome_cohorts" in second and sorted(first.get("outcome_cohorts") or []) != sorted(second["outcome_cohorts"] or []):
        diffs["outcome_cohorts"] = [first.get("outcome_cohorts"), second["outcome_cohorts"]]
    for f in COUNTS:
        if f in second and first.get(f) != second.get(f):
            diffs[f] = [first.get(f), second.get(f)]
    soft = {f: [first.get(f), second.get(f)] for f in SOFT if f in second and first.get(f) != second.get(f)}
    return diffs, soft


def main():
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from reconcile import stage1_claims  # first extraction + first-stage corrections, never later stages
    first = stage1_claims()
    rows = []
    for path in sorted((LIT / "second").glob("second_*.json")):
        rows += json.loads(path.read_text(encoding="utf-8"))
    out, seen = [], set()
    for s in rows:
        cid = s["claim_id"]
        if cid in seen:
            raise ValueError(f"Duplicate second extraction for {cid}")
        seen.add(cid)
        if cid not in first:
            raise ValueError(f"Unknown claim_id {cid}")
        if s.get("found") is False:
            out.append({"claim_id": cid, "status": "not_found", "notes": s.get("notes")})
            continue
        diffs, soft = compare(first[cid], s)
        out.append({"claim_id": cid, "status": "agree" if not diffs else "disagree", "differences": diffs,
                    "soft_differences": soft, "second_source_location": s.get("source_location"),
                    "first_source_location": first[cid].get("source_location"), "notes": s.get("notes")})
    headline = {cid for cid, c in first.items() if c["headline"]}
    missing = sorted(headline - seen)
    report = {"compared": len(out), "headline_claims": len(headline), "headline_not_second_extracted": missing,
              "status": dict(Counter(r["status"] for r in out)),
              "fields_disagreeing": dict(Counter(f for r in out for f in r.get("differences", {})).most_common()),
              "claims": out}
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "second_extraction_comparison.json").write_text(json.dumps(report, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({k: v for k, v in report.items() if k != "claims"}, indent=1))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Comparison failed: {exc}", file=sys.stderr)
        sys.exit(1)
