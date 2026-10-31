"""Turn the second extraction and adjudication into logged corrections.

  python src/reconcile.py

Reads the claim state after first-stage corrections (raw extraction + data/literature/corrections.json),
results/second_extraction_comparison.json, data/literature/second/second_*.json and
data/literature/adjudication/adjudicated_*.json. Writes data/literature/corrections_second.json,
which src/inventory.py applies after corrections.json. Re-running regenerates the file from scratch.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from inventory import LIT, apply_corrections, load_batches  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BY = "claude-opus-5-5 (lead reconciliation), 2026-10-29"


def stage1_claims():
    decisions, claims = load_batches("fulltext_decisions"), load_batches("claims")
    corrections = json.loads((LIT / "corrections.json").read_text(encoding="utf-8"))
    _, claims = apply_corrections(decisions, claims, corrections)
    return {c["claim_id"]: c for c in claims}


def build(claims, comparison, second, adjudicated):
    out = []
    for row in comparison["claims"]:
        cid, c = row["claim_id"], claims[row["claim_id"]]
        s = second.get(cid, {})
        base = {"target": "claim", "id": cid, "by": BY}
        fill_loc = {}
        if not c.get("instruments_location") and s.get("instruments_location"):
            fill_loc = {"instruments_location": s["instruments_location"]}
        if row["status"] == "agree":
            out.append({**base, "expect": {"second_extraction": c.get("second_extraction")},
                        "set": {"second_extraction": "agree (blind second extraction)", **fill_loc},
                        "why": "Blind second extraction agreed on every compared field."})
        elif row["status"] == "not_found":
            out.append({**base, "expect": {"second_extraction": c.get("second_extraction")},
                        "set": {"second_extraction": f"not found by blind second extractor: {(row.get('notes') or '')[:150]}"},
                        "why": "Second extractor could not locate the estimate; first extraction stays single-extracted."})
        else:
            adj = adjudicated.get(cid)
            if adj is None:
                raise ValueError(f"No adjudication for disagreeing claim {cid}")
            decided = {f["field"]: f for f in adj["fields"]}
            missing = set(row["differences"]) - set(decided)
            if missing:
                raise ValueError(f"Adjudication for {cid} misses fields {sorted(missing)}")
            expect, setv, summary, issues = {}, {}, [], []
            for field, (first_value, _second_value) in row["differences"].items():
                d = decided[field]
                expect[field] = first_value
                if d["value"] != first_value:
                    setv[field] = d["value"]
                summary.append(f"{field}={d['decision']}")
                if d.get("source_inconsistent"):
                    issues.append(f"{field}: {d.get('note') or d.get('evidence')}")
            setv["second_extraction"] = "adjudicated: " + ", ".join(summary)
            if issues:
                setv["source_issues"] = "; ".join(issues)[:400]
            if adj.get("overall_note"):
                setv["notes"] = ((c.get("notes") or "") + f" [Adjudication: {adj['overall_note']}]").strip()
            setv.update(fill_loc)
            out.append({**base, "expect": expect, "set": setv,
                        "why": "Second extraction disagreed; adjudicator checked each field against the source."})
    return out


def main():
    claims = stage1_claims()
    comparison = json.loads((ROOT / "results" / "second_extraction_comparison.json").read_text(encoding="utf-8"))
    second = {}
    for path in sorted((LIT / "second").glob("second_*.json")):
        for s in json.loads(path.read_text(encoding="utf-8")):
            second[s["claim_id"]] = s
    adjudicated = {}
    for path in sorted((LIT / "adjudication").glob("adjudicated_*.json")):
        for a in json.loads(path.read_text(encoding="utf-8")):
            adjudicated[a["claim_id"]] = a
    corrections = build(claims, comparison, second, adjudicated)
    (LIT / "corrections_second.json").write_text(json.dumps(corrections, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"{len(corrections)} second-stage corrections written")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Reconciliation failed: {exc}", file=sys.stderr)
        sys.exit(1)
