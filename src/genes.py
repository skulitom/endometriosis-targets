"""Look up genomic positions (GRCh37 and GRCh38) of every gene named in the claim table.

  python src/genes.py          # full lookup
  python src/genes.py --fill   # only resolve what is still unresolved (assembly map, then gene_overrides.json)

Uses the public Ensembl REST API (symbol lookup, then alias lookup for symbols that fail).
Writes data/reference/claimed_genes.json. src/blinding.py refuses to run while any gene is unresolved,
so the masking around claimed genes can never silently miss one.
"""
import json
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "reference" / "claimed_genes.json"
SERVERS = {"GRCh38": "https://rest.ensembl.org", "GRCh37": "https://grch37.rest.ensembl.org"}
HEADERS = {"Content-Type": "application/json", "Accept": "application/json",
           "User-Agent": "EndometriosisTargetCheck/0.1"}
PRIMARY_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}


def _request(url, payload=None):
    data = json.dumps(payload).encode() if payload is not None else None
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, data=data, headers=HEADERS, method="POST" if data else "GET")
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as exc:
            if exc.code in (400, 404):
                return None
            time.sleep(2 * (attempt + 1))
        except Exception:
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"Ensembl request failed: {url}")


def lookup_symbols(server, symbols):
    found = {}
    for i in range(0, len(symbols), 500):
        batch = symbols[i:i + 500]
        result = _request(f"{server}/lookup/symbol/homo_sapiens", {"symbols": batch}) or {}
        for sym, rec in result.items():
            if rec and str(rec.get("seq_region_name")) in PRIMARY_CHROMS:
                found[sym] = rec
        time.sleep(0.5)
    return found


def lookup_alias(server, symbol):
    """Resolve an alias or previous symbol to one gene on a primary chromosome, or None if ambiguous."""
    xrefs = _request(f"{server}/xrefs/symbol/homo_sapiens/{urllib.parse.quote(symbol)}?object_type=gene") or []
    hits = []
    for x in xrefs:
        rec = _request(f"{server}/lookup/id/{x['id']}")
        if rec and str(rec.get("seq_region_name")) in PRIMARY_CHROMS:
            hits.append(rec)
        time.sleep(0.2)
    unique = {(h["id"]) : h for h in hits}
    return next(iter(unique.values())) if len(unique) == 1 else None


def main():
    import urllib.parse  # noqa: F401  (used by lookup_alias)
    claims = json.loads((ROOT / "data" / "literature" / "claims.json").read_text(encoding="utf-8"))
    symbols = sorted({c["gene_symbol"] for c in claims})
    genes, unresolved = {}, []
    for build, server in SERVERS.items():
        found = lookup_symbols(server, symbols)
        for sym in symbols:
            rec = found.get(sym)
            method = "symbol"
            if rec is None:
                rec, method = lookup_alias(server, sym), "alias"
            if rec is None:
                unresolved.append(f"{sym} ({build})")
                continue
            genes.setdefault(sym, {})[build] = {"chrom": str(rec["seq_region_name"]), "start": int(rec["start"]),
                                                "end": int(rec["end"]), "ensembl_id": rec["id"],
                                                "ensembl_name": rec.get("display_name"), "method": method}
    # Newer HGNC symbols are unknown to the GRCh37 archive; look those up by their stable gene ID instead.
    still = []
    for item in unresolved:
        sym, build = item.rsplit(" (", 1)
        build = build.rstrip(")")
        other = "GRCh38" if build == "GRCh37" else "GRCh37"
        gid = genes.get(sym, {}).get(other, {}).get("ensembl_id")
        rec = _request(f"{SERVERS[build]}/lookup/id/{gid}") if gid else None
        if rec and str(rec.get("seq_region_name")) in PRIMARY_CHROMS:
            genes[sym][build] = {"chrom": str(rec["seq_region_name"]), "start": int(rec["start"]), "end": int(rec["end"]),
                                 "ensembl_id": rec["id"], "ensembl_name": rec.get("display_name"),
                                 "method": f"stable id from {other}"}
        else:
            still.append(item)
    genes, unresolved = fill_unresolved(genes, still)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"retrieved_at": datetime.now(timezone.utc).isoformat(),
                               "source": "Ensembl REST (rest.ensembl.org, grch37.rest.ensembl.org)",
                               "genes": genes, "unresolved": unresolved}, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"{len(symbols)} symbols; resolved in both builds: "
          f"{sum(1 for g in genes.values() if len(g) == 2)}; unresolved: {unresolved}")


OVERRIDES = ROOT / "data" / "reference" / "gene_overrides.json"


def map_assembly(rec, src, dst):
    """Ensembl assembly mapping of a gene span; returns the merged mapped span on one chromosome, or None."""
    region = f"{rec['chrom']}:{rec['start']}..{rec['end']}:1"
    data = _request(f"https://rest.ensembl.org/map/human/{src}/{region}/{dst}")
    spans = [m["mapped"] for m in (data or {}).get("mappings", [])]
    chroms = {str(s["seq_region_name"]) for s in spans}
    if len(chroms) != 1 or not chroms <= PRIMARY_CHROMS:
        return None
    return {"chrom": chroms.pop(), "start": min(s["start"] for s in spans), "end": max(s["end"] for s in spans)}


def fill_unresolved(genes, unresolved):
    """Second pass for symbols one build lacks: Ensembl assembly mapping from the other build, then the
    logged overrides in data/reference/gene_overrides.json."""
    overrides = json.loads(OVERRIDES.read_text(encoding="utf-8"))["genes"] if OVERRIDES.exists() else {}
    still = []
    for item in unresolved:
        sym, build = item.rsplit(" (", 1)
        build = build.rstrip(")")
        other = "GRCh38" if build == "GRCh37" else "GRCh37"
        src = genes.get(sym, {}).get(other)
        mapped = map_assembly(src, other, build) if src and "override" not in src.get("method", "") else None
        if mapped:
            genes[sym][build] = {**mapped, "ensembl_id": src.get("ensembl_id"), "ensembl_name": src.get("ensembl_name"),
                                 "method": f"Ensembl assembly map from {other}"}
        elif sym in overrides and build in overrides[sym]:
            genes.setdefault(sym, {})[build] = {**overrides[sym][build], "ensembl_id": None, "ensembl_name": None,
                                                "method": "override (data/reference/gene_overrides.json)"}
        else:
            still.append(item)
    return genes, still


def fill_only():
    """Apply the second pass to the existing file without repeating the Ensembl symbol lookups."""
    data = json.loads(OUT.read_text(encoding="utf-8"))
    genes, unresolved = fill_unresolved(data["genes"], data["unresolved"])
    data.update(genes=genes, unresolved=unresolved, filled_at=datetime.now(timezone.utc).isoformat())
    OUT.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"unresolved after fill: {unresolved}")


if __name__ == "__main__":
    try:
        fill_only() if sys.argv[1:] == ["--fill"] else main()
    except Exception as exc:
        print(f"Gene lookup failed: {exc}", file=sys.stderr)
        sys.exit(1)
