"""Build screening input batches for the supplementary full-text search (PROTOCOL.md section 1, Search).

  python -I src/supp_batches.py [n_batches]

Reads data/literature/records_supplementary.json and data/literature/refscan/supplementary_passages.json;
writes data/literature/screening/input_supp_<k>.json. Excerpts are trimmed so each batch stays readable:
MR passages that mention endometriosis or describe the paper's own analysis come first.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from refscan import ENDO, OWN_ANALYSIS  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
LIT = ROOT / "data" / "literature"
MAX_PASSAGES, PASSAGE_CHARS, ABSTRACT_CHARS = 5, 700, 1600


def pick_passages(passages):
    ranked = sorted(range(len(passages)), key=lambda i: (
        not (OWN_ANALYSIS.search(passages[i]) and ENDO.search(passages[i])),
        not OWN_ANALYSIS.search(passages[i]), not ENDO.search(passages[i]), i))
    return [passages[i][:PASSAGE_CHARS] for i in ranked[:MAX_PASSAGES]]


def build(records, passages, n_batches):
    todo = [r for r in records if r["found_by"] == ["europepmc_fulltext"]]
    missing = [r["record_id"] for r in todo if r["record_id"] not in passages]
    if missing:
        raise ValueError(f"{len(missing)} records have no passage entry, e.g. {missing[:3]}")
    rows = []
    for r in todo:
        p = passages[r["record_id"]]
        rows.append({"record_id": r["record_id"], "title": r.get("title"), "year": r.get("year"),
                     "journal": r.get("journal"), "publication_types": r.get("publication_types"),
                     "preprint": r.get("preprint"), "abstract": (r.get("abstract") or "")[:ABSTRACT_CHARS],
                     "full_text_status": p["status"], "n_mr_passages": len(p["passages"]),
                     "mr_passages": pick_passages(p["passages"]),
                     "endometriosis_sentences": p.get("endometriosis_sentences", [])})
    size = -(-len(rows) // n_batches)
    return [rows[i:i + size] for i in range(0, len(rows), size)]


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    records = json.loads((LIT / "records_supplementary.json").read_text(encoding="utf-8"))
    passages = json.loads((LIT / "refscan" / "supplementary_passages.json").read_text(encoding="utf-8"))
    batches = build(records, passages, n)
    for k, batch in enumerate(batches, 1):
        path = LIT / "screening" / f"input_supp_{k}.json"
        path.write_text(json.dumps(batch, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
        print(path.name, len(batch), sum(len(json.dumps(x)) for x in batch))


if __name__ == "__main__":
    main()
