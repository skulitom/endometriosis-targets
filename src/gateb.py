"""Gate B positive controls: do the replication files behave as expected at known endometriosis loci?

  python src/gateb.py

Follows PROTOCOL.md section 2 ("Positive-control and reproduction rules"), fixed before any file was read.
All reads go through src/blinding.py, so nothing within 1 Mb of a claimed gene is ever seen here.
Writes results/gateb_positive_controls.json.
"""
import json
import math
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from blinding import norm_chrom, read_blinded  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
UP = ROOT / "data" / "upstream"
OUT = ROOT / "results" / "gateb_positive_controls.json"
GW = 5e-8
CLUMP = 1_000_000
COMP = str.maketrans("ACGT", "TGCA")

SOURCES = {
    "pujolgualdo2025_estbb_finngen": {"file": "GCST90454213.tsv.gz", "build": "GRCh37", "chrom": "chromosome",
                                      "pos": "base_pair_location", "ea": "effect_allele", "oa": "other_allele",
                                      "beta": "beta", "se": "standard_error", "af": "effect_allele_frequency",
                                      "match": "position"},
    "panukbb_phecode615_females": {"file": "phecode-615-females.tsv.bgz", "build": "GRCh37", "chrom": "chr",
                                   "pos": "pos", "ea": "alt", "oa": "ref", "beta": "beta_EUR", "se": "se_EUR",
                                   "af": "af_controls_EUR", "exclude": ("low_confidence_EUR", "true"),
                                   "match": "position"},
    "finngen_r12_endo": {"file": "finngen_R12_N14_ENDOMETRIOSIS.gz", "build": "GRCh38", "chrom": "chrom", "pos": "pos",
                         "ea": "alt", "oa": "ref", "beta": "beta", "se": "sebeta", "af": "af_alt", "match": "rsid"},
    "finngen_r10_endo": {"file": "finngen_R10_N14_ENDOMETRIOSIS.gz", "build": "GRCh38", "chrom": "chrom", "pos": "pos",
                         "ea": "alt", "oa": "ref", "beta": "beta", "se": "sebeta", "af": "af_alt", "match": "rsid"},
}


def palindromic(a, b):
    return {a, b} in ({"A", "T"}, {"C", "G"})


def align(lead, ea, oa, beta, af):
    """Return beta aligned to the lead's effect allele, or (None, reason)."""
    le, lo = lead["ea"], lead["oa"]
    ea, oa = ea.upper(), oa.upper()
    if palindromic(le, lo):
        freq = lead.get("eaf")
        if freq is None or 0.42 <= freq <= 0.58:
            return None, "palindromic with ambiguous frequency"
    if (ea, oa) == (le, lo):
        return beta, None
    if (ea, oa) == (lo, le):
        return -beta, None
    if not palindromic(le, lo):
        if (ea.translate(COMP), oa.translate(COMP)) == (le, lo):
            return beta, None
        if (ea.translate(COMP), oa.translate(COMP)) == (lo, le):
            return -beta, None
    return None, "alleles do not match"


def binom_upper(k, n):
    """One-sided P(X >= k) for X ~ Binomial(n, 0.5)."""
    return sum(math.comb(n, i) for i in range(k, n + 1)) / 2 ** n


def num(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except (TypeError, ValueError):
        return None


def lead_variants():
    hits = []
    for r in read_blinded(UP / "GCST90205183_buildGRCh37.tsv.gz", "GRCh37", "chromosome", "base_pair_location",
                          "positive controls: genome-wide significant variants outside masked regions",
                          "rahmioglu2023_meta", row_filter=lambda p, i: (num(p[i["p_value"]]) or 1) < GW):
        hits.append({"chrom": norm_chrom(r["chromosome"]), "pos": int(r["base_pair_location"]),
                     "ea": r["effect_allele"].upper(), "oa": r["other_allele"].upper(),
                     "eaf": num(r["effect_allele_frequency"]), "beta": float(r["beta"]),
                     "se": float(r["standard_error"]), "p": float(r["p_value"])})
    leads = []
    for h in sorted(hits, key=lambda x: x["p"]):
        if all(not (h["chrom"] == L["chrom"] and abs(h["pos"] - L["pos"]) <= CLUMP) for L in leads):
            leads.append(h)
    return leads, len(hits)


def ensembl_rsid(lead):
    url = (f"https://grch37.rest.ensembl.org/overlap/region/human/{lead['chrom']}:{lead['pos']}-{lead['pos']}"
           "?feature=variation;content-type=application/json")
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "EndometriosisTargetCheck/0.1"}),
                                        timeout=60) as r:
                rows = json.loads(r.read())
            ids = [x["id"] for x in rows if x["id"].startswith("rs") and x["start"] == lead["pos"]]
            return ids[0] if len(ids) == 1 else None
        except Exception:
            time.sleep(2 * (attempt + 1))
    return None


def scan(key, spec, leads, rsid_of):
    by_pos = {(L["chrom"], L["pos"]): L for L in leads}
    wanted_rsids = {rsid_of[(L["chrom"], L["pos"])]: L for L in leads if rsid_of.get((L["chrom"], L["pos"]))}
    found = {}

    def keep(parts, idx):
        if spec["match"] == "position":
            return (norm_chrom(parts[idx[spec["chrom"]]]), int(parts[idx[spec["pos"]]])) in by_pos
        return any(r in wanted_rsids for r in parts[idx["rsids"]].split(","))

    for r in read_blinded(UP / spec["file"], spec["build"], spec["chrom"], spec["pos"],
                          "positive controls: lookup at lead variants", key, row_filter=keep):
        if spec.get("exclude") and r.get(spec["exclude"][0], "").lower() == spec["exclude"][1]:
            continue
        if spec["match"] == "position":
            lead = by_pos[(norm_chrom(r[spec["chrom"]]), int(r[spec["pos"]]))]
        else:
            lead = next(wanted_rsids[x] for x in r["rsids"].split(",") if x in wanted_rsids)
        beta = num(r[spec["beta"]])
        if beta is None:
            continue
        aligned, why = align(lead, r[spec["ea"]], r[spec["oa"]], beta, num(r.get(spec["af"])))
        key_ = (lead["chrom"], lead["pos"])
        if aligned is None:
            found.setdefault(key_, {"excluded": why})
            continue
        found[key_] = {"beta": aligned, "se": num(r[spec["se"]]), "source_pos": int(r[spec["pos"]])}
    matched = {k: v for k, v in found.items() if "beta" in v}
    concordant = sum((v["beta"] > 0) == (by_pos[k]["beta"] > 0) for k, v in matched.items())
    n = len(matched)
    p = binom_upper(concordant, n) if n else None
    nominal = sum(1 for k, v in matched.items() if v["se"] and abs(v["beta"] / v["se"]) > 1.96
                  and (v["beta"] > 0) == (by_pos[k]["beta"] > 0))
    return {"leads": len(leads), "matched": n, "excluded": {k: v["excluded"] for k, v in
                                                            ((f"{a}:{b}", v) for (a, b), v in found.items()) if "excluded" in v},
            "not_found": len(leads) - len(found), "concordant": concordant,
            "concordance": round(concordant / n, 3) if n else None, "binomial_p": p,
            "concordant_and_p_below_0.05": nominal,
            "passes": bool(n and concordant / n >= 0.70 and p < 0.01),
            "per_lead": {f"{k[0]}:{k[1]}": v for k, v in found.items()}}


def main():
    leads, n_sig = lead_variants()
    print(f"{n_sig} genome-wide significant variants outside masked regions; {len(leads)} independent leads")
    rsid_of = {}
    est = SOURCES["pujolgualdo2025_estbb_finngen"]
    by_pos = {(L["chrom"], L["pos"]) for L in leads}
    for r in read_blinded(UP / est["file"], "GRCh37", est["chrom"], est["pos"],
                          "positive controls: rsIDs at lead positions", "pujolgualdo2025_estbb_finngen",
                          row_filter=lambda p, i: (norm_chrom(p[i["chromosome"]]), int(p[i["base_pair_location"]])) in by_pos):
        if r.get("rs_id", "").startswith("rs"):
            rsid_of[(norm_chrom(r["chromosome"]), int(r["base_pair_location"]))] = r["rs_id"]
    for L in leads:
        k = (L["chrom"], L["pos"])
        if k not in rsid_of:
            rs = ensembl_rsid(L)
            if rs:
                rsid_of[k] = rs
            time.sleep(0.3)
    report = {"rule": "PROTOCOL.md section 2, positive-control and reproduction rules",
              "genome_wide_significant_unmasked": n_sig,
              "leads": [{**L, "rsid": rsid_of.get((L["chrom"], L["pos"]))} for L in leads], "sources": {}}
    for key, spec in SOURCES.items():
        t = time.time()
        report["sources"][key] = scan(key, spec, leads, rsid_of)
        s = report["sources"][key]
        print(f"{key}: matched {s['matched']}/{s['leads']}, concordant {s['concordant']} "
              f"({s['concordance']}), binomial p {s['binomial_p']}, passes={s['passes']} [{time.time() - t:.0f}s]")
    OUT.write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
