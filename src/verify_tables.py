"""Cross-check extracted variants against the original .xlsx files with an independent reader.

  .venv\\Scripts\\python src\\verify_tables.py

src/tables.py converts spreadsheets with its own standard-library parser. This script re-reads the
original workbook with openpyxl (project venv only) and compares every extracted value for claims whose
table came from an .xlsx file. Writes results/table_crosscheck.json.
"""
import json
import math
import re
import sys
from pathlib import Path

import openpyxl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from instruments import LIT, SPECS, SUPP, parse_number, row_matches  # noqa: E402
from supplements import slug  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CHECK_FIELDS = ["effect_allele", "other_allele", "beta_exposure", "se_exposure", "beta_outcome", "se_outcome"]


def sheet_key(name):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)[:60]


def cell_text(value):
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer() and abs(value) < 1e15:
        return str(int(value))
    return str(value).strip()


def same(a, b, numeric):
    if numeric:
        x, y = parse_number(a), parse_number(b)
        if x is None or y is None:
            return x is None and y is None
        return math.isclose(x, y, rel_tol=1e-9, abs_tol=1e-12)
    return (a or "").strip().upper() == (b or "").strip().upper()


def load_sheet(xlsx, wanted_key, header_names):
    wb = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    try:
        sheet = next((wb[s] for s in wb.sheetnames if sheet_key(s) == wanted_key), None)
        if sheet is None:
            return None
        rows = [[cell_text(v) for v in r] for r in sheet.iter_rows(values_only=True)]
    finally:
        wb.close()
    for i, row in enumerate(rows):
        stripped = [c.strip() for c in row]
        if all(h in stripped for h in header_names):
            return [dict(zip(stripped, r)) for r in rows[i + 1:]]
    return None


def main():
    extracted = json.loads((LIT / "instruments" / "instruments.json").read_text(encoding="utf-8"))
    report, cache = {}, {}
    for spec_path in sorted(SPECS.glob("*.json")):
        spec = json.loads(spec_path.read_text(encoding="utf-8"))
        supp = SUPP / slug(spec["record_id"])
        for cid, entry in spec["claims"].items():
            r = extracted.get(cid)
            if entry.get("status") != "mapped" or not r or r["status"] != "extracted":
                continue
            t = spec["tables"][entry["table"]]
            if t.get("transcribed"):
                report[cid] = {"checked": False, "reason": "transcribed table; verified by hand sample"}
                continue
            m = re.match(r"text/(.+?)__(.+)\.csv$", t["file"])
            xlsx = supp / "files" / f"{m.group(1)}.xlsx" if m else None
            if not m or not xlsx.exists():
                report[cid] = {"checked": False, "reason": "table not from an .xlsx file (csv/txt/xls/docx source)"}
                continue
            key = (str(xlsx), m.group(2))
            header_names = [c for c in t["columns"].values() if c] + [f["column"] for f in entry["filters"]]
            if key not in cache:
                cache[key] = load_sheet(xlsx, m.group(2), sorted(set(header_names)))
            rows = cache[key]
            if rows is None:
                report[cid] = {"checked": False, "reason": "sheet or header not found by openpyxl"}
                continue
            candidates = [row for row in rows if row_matches(row, entry["filters"])]
            by_rsid = {}
            for row in candidates:
                by_rsid.setdefault(row.get(t["columns"]["rsid"], "").strip(), []).append(row)
            mismatches, checked = [], 0
            mapped_cols = [c for c in t["columns"].values() if c]
            for v in r["variants"]:
                hits = by_rsid.get(v["rsid"], [])
                # Rows repeated with identical mapped values (src/instruments.py ignores them) count as one.
                if len(hits) > 1 and all([h.get(c) for c in mapped_cols] == [hits[0].get(c) for c in mapped_cols] for h in hits):
                    hits = hits[:1]
                if len(hits) != 1:
                    mismatches.append(f"{v['rsid']}: {len(hits)} matching rows in original")
                    continue
                checked += 1
                for field in CHECK_FIELDS:
                    col = t["columns"].get(field)
                    if not col:
                        continue
                    ours = v.get(field)
                    theirs = hits[0].get(col)
                    numeric = field.startswith(("beta", "se"))
                    if not same(None if ours is None else str(ours), theirs, numeric):
                        mismatches.append(f"{v['rsid']} {field}: extracted {ours!r} vs original {theirs!r}")
            if len(candidates) != r["n"] + len(r.get("identical_repeats_ignored") or []):
                mismatches.append(f"filters match {len(candidates)} rows in original vs {r['n']} extracted")
            report[cid] = {"checked": True, "variants_checked": checked, "mismatches": mismatches}
    summary = {"claims_crosschecked": sum(v["checked"] for v in report.values()),
               "claims_with_mismatches": sum(bool(v.get("mismatches")) for v in report.values()),
               "variants_checked": sum(v.get("variants_checked", 0) for v in report.values()),
               "not_checked_reasons": {}, "claims": report}
    for v in report.values():
        if not v["checked"]:
            summary["not_checked_reasons"][v["reason"]] = summary["not_checked_reasons"].get(v["reason"], 0) + 1
    (ROOT / "results" / "table_crosscheck.json").write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "claims"}, indent=1))


if __name__ == "__main__":
    main()
