"""Write the public versions of the literature records: identifiers, titles, decisions and reasons, no copied text.

  python -I src/public_export.py

The working files (search results, merged records, screening inputs) hold abstracts and full-text excerpts
copied from papers, so they stay out of the public repository (.gitignore). This writes, to
data/literature/public/:
  searches.json    every search: name, source, query, run date, hit count and the identifiers it returned
  records.json     every record screened, from all sources: identifiers, title, year, journal, found_by
  screening.json   every title/abstract decision with its reason (and the lead's re-check, where made)
  fulltext.json    every full-text decision with its reason
Anyone can re-fetch the abstracts from PubMed or Europe PMC with the identifiers.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from inventory import LIT, apply_corrections, load_batches  # noqa: E402

OUT = LIT / "public"
ID_FIELDS = ("record_id", "pmid", "doi", "pmcid", "europepmc_id", "title", "year", "journal", "preprint", "found_by")
SCREEN_FIELDS = ("decision", "reason_code", "reason", "lead_checked", "screener_decision", "screener_reason_code",
                 "screener_reason", "evidence", "source")


def _write(name, obj):
    (OUT / name).write_text(json.dumps(obj, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    searches = []
    for path in sorted(LIT.glob("search_*.json")):
        p = json.loads(path.read_text(encoding="utf-8"))
        searches.append({"file": path.name, "name": p["name"], "source": p["source"], "query": p["query"],
                         "run_at": p["run_at"], "count": p["count"],
                         "hits": [{k: r.get(k) for k in ("pmid", "doi", "pmcid", "europepmc_id", "title")}
                                  for r in p["records"]]})
    _write("searches.json", searches)
    screened = json.loads((LIT / "screening_decisions.json").read_text(encoding="utf-8"))
    _write("records.json", [{k: r.get(k) for k in ID_FIELDS} for r in screened])
    _write("screening.json", [{"record_id": r["record_id"], "title": r.get("title"),
                               **{k: r["screening"][k] for k in SCREEN_FIELDS if k in r["screening"]}}
                              for r in screened])
    corrections = []
    for name in ("corrections.json", "corrections_second.json", "corrections_lead.json"):
        path = LIT / name
        if path.exists():
            corrections += json.loads(path.read_text(encoding="utf-8"))
    decisions, _ = apply_corrections(load_batches("fulltext_decisions"), load_batches("claims"), corrections)
    _write("fulltext.json", sorted(decisions, key=lambda d: d["record_id"]))
    print({"searches": len(searches), "records": len(screened), "fulltext_decisions": len(decisions)})


if __name__ == "__main__":
    main()
