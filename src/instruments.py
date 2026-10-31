"""Extract per-variant instruments for claims from supplement tables, driven by per-paper mapping specs.

  python src/instruments.py                    # all specs -> shared outputs
  python src/instruments.py --check <spec>     # one spec, printed only (safe while others work)

A spec (data/literature/instruments/specs/<slug>.json) says which converted table holds the variants,
which columns mean what, and which rows belong to each claim. This script copies values from the
tables (never retyped), checks them, and recomputes the Wald ratio / IVW estimate wherever outcome
effects are present, to verify that the right rows were taken.

Writes data/literature/instruments/instruments.json and results/instrument_verification.json.

Spec format:
{
  "record_id": "doi:...",
  "tables": {"T1": {"file": "text/<converted>.csv", "header_row": 1,
                    "columns": {"rsid": "SNP", "effect_allele": "...", "other_allele": "...",
                                "beta_exposure": "...", "se_exposure": "...", "pval_exposure": "...",
                                "eaf_exposure": "...", "beta_outcome": "...", "se_outcome": "...",
                                "pval_outcome": "...", "effect_allele_outcome": "...", "other_allele_outcome": "..."}}},
  "claims": {"<claim_id>": {"status": "mapped", "table": "T1",
                            "filters": [{"column": "exposure", "op": "equals|contains|regex", "value": "RSPO3"}],
                            "expected_n": 3, "notes": null}
             | {"status": "not_available", "reason": "..."}}
}
Only rsid, effect_allele and beta_exposure columns are required; others may be omitted.
"""
import csv
import json
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LIT = ROOT / "data" / "literature"
SPECS = LIT / "instruments" / "specs"
SUPP = ROOT / "data" / "upstream" / "supplements"
RESULTS = ROOT / "results"

FIELDS = ["rsid", "effect_allele", "other_allele", "beta_exposure", "se_exposure", "pval_exposure",
          "eaf_exposure", "beta_outcome", "se_outcome", "pval_outcome", "effect_allele_outcome",
          "other_allele_outcome"]
NUMERIC = {"beta_exposure", "se_exposure", "pval_exposure", "eaf_exposure", "beta_outcome", "se_outcome", "pval_outcome"}
REQUIRED = {"rsid", "effect_allele", "beta_exposure"}
VARIANT_RE = re.compile(r"^(rs\d+|(chr)?[0-9XYM]{1,2}[:_]\d+([:_][ACGT]+[:_][ACGT]+)?)$", re.IGNORECASE)
ALLELE_RE = re.compile(r"^[ACGT]+$")
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def parse_number(text):
    if text is None:
        return None
    t = str(text).strip().replace("−", "-").replace("–", "-").replace("×", "x")
    if t in {"", "NA", "N/A", "NaN", "nan", "-", "."}:
        return None
    m = re.fullmatch(r"([-+]?\d*\.?\d+)\s*x\s*10\^?\s*([-+]?\d+)", t)  # 1.2 x 10^-5
    if m:
        return float(m.group(1)) * 10 ** int(m.group(2))
    try:
        value = float(t)
    except ValueError:
        raise ValueError(f"not a number: {text!r}")
    if math.isnan(value) or math.isinf(value):
        return None
    return value


def read_table(path, header_row):
    with open(path, newline="", encoding="utf-8", errors="replace") as handle:
        rows = list(csv.reader(handle))
    if header_row < 1 or header_row > len(rows):
        raise ValueError(f"header_row {header_row} outside {path.name}")
    header = [h.strip() for h in rows[header_row - 1]]
    body = []
    for offset, row in enumerate(rows[header_row:], start=header_row + 1):
        if any(cell.strip() for cell in row):
            body.append((offset, dict(zip(header, [c.strip() for c in row]))))
    return header, body


def row_matches(row, filters):
    for f in filters:
        cell = row.get(f["column"], "")
        value = f["value"]
        if f["op"] == "equals" and cell != value:
            return False
        if f["op"] == "contains" and value.lower() not in cell.lower():
            return False
        if f["op"] == "regex" and not re.search(value, cell):
            return False
    return True


def to_variant(row, columns, file_label, line):
    v = {}
    for field in FIELDS:
        col = columns.get(field)
        raw = row.get(col) if col else None
        if field in NUMERIC:
            v[field] = parse_number(raw)
        else:
            v[field] = (raw or "").strip().upper() or None if field != "rsid" else (raw or "").strip() or None
    problems = []
    if not v["rsid"] or not VARIANT_RE.match(v["rsid"]):
        problems.append(f"bad variant id {v['rsid']!r}")
    for allele_field in ("effect_allele", "other_allele", "effect_allele_outcome", "other_allele_outcome"):
        if v[allele_field] is not None and not ALLELE_RE.match(v[allele_field]):
            problems.append(f"bad {allele_field} {v[allele_field]!r}")
    if v["beta_exposure"] is None:
        problems.append("missing beta_exposure")
    if v["se_exposure"] is not None and v["se_exposure"] <= 0:
        problems.append("non-positive se_exposure")
    if v["se_outcome"] is not None and v["se_outcome"] <= 0:
        problems.append("non-positive se_outcome")
    v["source"] = {"file": file_label, "line": line}
    return v, problems


def align_outcome(v):
    """Return beta_outcome aligned to the exposure effect allele, or None if it cannot be aligned safely."""
    by = v.get("beta_outcome")
    if by is None:
        return None
    ea, oa = v.get("effect_allele"), v.get("other_allele")
    eo, oo = v.get("effect_allele_outcome"), v.get("other_allele_outcome")
    if eo is None:
        return by  # harmonised table: outcome effect already on the exposure allele
    if eo == ea and (oo is None or oo == oa):
        return by
    if eo == oa and (oo is None or oo == ea):
        return -by
    if ea and eo == ea.translate(COMPLEMENT) and ea.translate(COMPLEMENT) != oa:
        return by  # strand flip, non-palindromic
    return None


def recompute(variants):
    """Wald ratio (one variant) or fixed-effect IVW. Returns (beta, se) on the outcome log scale, or None."""
    pairs = []
    for v in variants:
        by, bx, sy = align_outcome(v), v.get("beta_exposure"), v.get("se_outcome")
        if by is None or bx in (None, 0) or sy in (None, 0):
            return None
        pairs.append((bx, by, sy))
    if not pairs:
        return None
    if len(pairs) == 1:
        bx, by, sy = pairs[0]
        return by / bx, sy / abs(bx)
    weight = sum(bx * bx / (sy * sy) for bx, _, sy in pairs)
    beta = sum(bx * by / (sy * sy) for bx, by, sy in pairs) / weight
    return beta, 1 / math.sqrt(weight)


def compare_estimate(claim, recomputed):
    """Agreement of a recomputed log-scale estimate with the claim's reported estimate."""
    est = claim.get("estimate")
    if recomputed is None or est is None:
        return None
    reported = math.log(est) if claim["estimate_scale"] in {"OR", "HR"} else est
    if claim["estimate_scale"] in {"OR", "HR"} and est <= 0:
        return None
    beta = recomputed[0]
    tolerance = max(0.01, 0.05 * abs(reported))
    return {"reported_log": reported, "recomputed_log": beta, "recomputed_se": recomputed[1],
            "difference": beta - reported, "agree": abs(beta - reported) <= tolerance,
            "same_sign": (beta > 0) == (reported > 0) or abs(reported) < 1e-9}


def apply_spec(spec, supp_dir, claims):
    tables, out = {}, {}
    for name, t in spec.get("tables", {}).items():
        # Transcribed tables (typed from a PDF or figure) live apart and are always sent for verification.
        base = (LIT / "instruments" / "transcribed") if t.get("transcribed") else supp_dir
        path = (base / t["file"]).resolve()
        if base.resolve() not in path.parents:
            raise ValueError(f"table path escapes its folder: {t['file']}")
        header, body = read_table(path, t.get("header_row", 1))
        missing_cols = [c for f, c in t["columns"].items() if c and c not in header]
        if missing_cols:
            raise ValueError(f"{t['file']}: columns not in header: {missing_cols}")
        if not REQUIRED <= {f for f, c in t["columns"].items() if c}:
            raise ValueError(f"{t['file']}: spec lacks required columns {sorted(REQUIRED)}")
        tables[name] = (t, body)
    for cid, entry in spec["claims"].items():
        if cid not in claims:
            raise ValueError(f"spec names unknown claim {cid}")
        if entry["status"] != "mapped":
            out[cid] = {"status": "not_available", "reason": entry.get("reason")}
            continue
        t, body = tables[entry["table"]]
        variants, problems, seen = [], [], set()
        for line, row in body:
            if row_matches(row, entry["filters"]):
                v, p = to_variant(row, t["columns"], t["file"], line)
                if v["rsid"] in seen:
                    problems.append(f"duplicate variant {v['rsid']} (line {line})")
                    continue
                seen.add(v["rsid"])
                variants.append(v)
                problems += [f"line {line}: {x}" for x in p]
        flags = []
        if not variants:
            problems.append("filters matched no rows")
        if entry.get("expected_n") is not None and len(variants) != entry["expected_n"]:
            flags.append(f"spec expected {entry['expected_n']} variants, found {len(variants)}")
        reported_n = claims[cid].get("n_instruments")
        if reported_n is not None and variants and len(variants) != reported_n:
            flags.append(f"claim reports {reported_n} instruments, table has {len(variants)}")
        check = compare_estimate(claims[cid], recompute(variants)) if variants else None
        out[cid] = {"status": "error" if problems else "extracted", "n": len(variants), "variants": variants,
                    "problems": problems, "flags": flags, "recompute": check, "notes": entry.get("notes"),
                    "transcribed": bool(t.get("transcribed")), "record_id": spec["record_id"]}
    return out


def check_one(spec_path):
    """Run one spec and print a compact report without writing shared outputs (for mapping work)."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from supplements import slug
    claims = {c["claim_id"]: c for c in json.loads((LIT / "claims.json").read_text(encoding="utf-8"))}
    spec = json.loads(Path(spec_path).read_text(encoding="utf-8"))
    result = apply_spec(spec, SUPP / slug(spec["record_id"]), claims)
    for cid, r in result.items():
        rc = r.get("recompute") or {}
        print(json.dumps({"claim_id": cid, "status": r["status"], "n": r.get("n"), "flags": r.get("flags"),
                          "problems": (r.get("problems") or [])[:3], "reason": r.get("reason"),
                          "recompute_agree": rc.get("agree"), "reported_log": rc.get("reported_log"),
                          "recomputed_log": rc.get("recomputed_log")}))


def main():
    claims = {c["claim_id"]: c for c in json.loads((LIT / "claims.json").read_text(encoding="utf-8"))}
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from supplements import slug
    results, spec_errors = {}, {}
    for spec_path in sorted(SPECS.glob("*.json")):
        spec = json.loads(spec_path.read_text(encoding="utf-8"))
        try:
            results.update(apply_spec(spec, SUPP / slug(spec["record_id"]), claims))
        except Exception as exc:
            spec_errors[spec_path.name] = str(exc)[:200]
    (LIT / "instruments" / "instruments.json").write_text(json.dumps(results, indent=1) + "\n", encoding="utf-8")
    primary = [cid for cid, c in claims.items() if c["exposure_type"] != "mQTL"]
    checked = [r["recompute"] for r in results.values() if r.get("recompute")]
    summary = {
        "primary_claims": len(primary),
        "claims_in_specs": len(results),
        "primary_without_spec": sorted(set(primary) - set(results)),
        "status": {s: sum(r["status"] == s for r in results.values()) for s in ("extracted", "error", "not_available")},
        "variants_extracted": sum(r.get("n", 0) for r in results.values() if r["status"] == "extracted"),
        "recomputed": len(checked),
        "recompute_agree": sum(c["agree"] for c in checked),
        "recompute_sign_mismatch": sum(not c["same_sign"] for c in checked),
        "claims_with_count_flags": sum(bool(r.get("flags")) for r in results.values()),
        "spec_errors": spec_errors,
        "per_claim": {cid: {k: r.get(k) for k in ("status", "n", "flags", "problems", "recompute", "reason")}
                      for cid, r in sorted(results.items())},
    }
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "instrument_verification.json").write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k not in ("per_claim", "primary_without_spec")}, indent=1))
    print(f"primary claims without a spec: {len(summary['primary_without_spec'])}")


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--check":
        check_one(sys.argv[2])
    else:
        main()
