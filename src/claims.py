"""Claim table schema, validation, independence check and summary.

  python src/claims.py validate data/literature/claims.json
  python src/claims.py summary  data/literature/claims.json --csv results/claims_review.csv

A claim is one paper x gene x exposure dataset x outcome GWAS (see PROTOCOL.md section 1).
Validation checks internal consistency only; it does not show that an extraction matches its source.
"""
import argparse
import csv
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path

EXPOSURE_TYPES = {"pQTL", "eQTL", "sQTL", "SMR_eQTL", "SMR_pQTL", "drug_target_proxy", "mQTL", "other"}
SECONDARY_ONLY_TYPES = {"mQTL"}  # recorded, but not in the primary replication analysis (PROTOCOL.md section 1)
INSTRUMENT_SCOPES = {"cis", "cis_and_trans", "genome_wide", "unclear"}
# Scale of the reported effect: per higher level of the gene product (or, for drug-target proxies, of the proxied
# biomarker through variants at the gene); per lower level or inhibition; or not stated. None means per_higher.
EFFECT_SCALES = {None, "per_higher", "per_lower", "unclear"}
SCALES = {"OR", "beta", "HR"}
DIRECTIONS = {"risk_increasing", "protective", "null_or_unclear"}
KNOWN_COHORTS = {"FinnGen", "UK Biobank", "23andMe", "BioBank Japan", "Rahmioglu2023_meta",
                 "deCODE", "Estonian Biobank", "Million Veteran Program", "All of Us", "other", "unknown"}

# field: (required, allowed python types)
FIELDS = {
    "claim_id": (True, (str,)),
    "paper_id": (True, (str,)),
    "doi": (False, (str, type(None))),
    "pmid": (False, (str, type(None))),
    "first_author": (True, (str,)),
    "year": (True, (int,)),
    "journal": (False, (str, type(None))),
    "preprint": (True, (bool,)),
    "gene_symbol": (True, (str,)),
    "exposure_type": (True, (str,)),
    "exposure_source": (True, (str,)),
    "exposure_tissue": (False, (str, type(None))),
    "outcome_source": (True, (str,)),
    "outcome_release": (False, (str, type(None))),
    "outcome_cases": (False, (int, type(None))),
    "outcome_controls": (False, (int, type(None))),
    "outcome_cohorts": (True, (list,)),
    "instrument_scope": (True, (str,)),
    "n_instruments": (False, (int, type(None))),
    "instruments": (False, (list, type(None))),
    "instruments_location": (False, (str, type(None))),
    "method": (True, (str,)),
    "estimate_scale": (True, (str,)),
    "estimate": (False, (int, float, type(None))),
    "ci_low": (False, (int, float, type(None))),
    "ci_high": (False, (int, float, type(None))),
    "p_value": (False, (int, float, type(None))),
    "direction": (True, (str,)),
    "coloc_reported": (True, (bool,)),
    "coloc_pph4": (False, (int, float, type(None))),
    "replication_reported": (True, (bool,)),
    "replication_details": (False, (str, type(None))),
    "experimental_validation": (False, (str, type(None))),
    "headline": (True, (bool,)),
    "source_location": (True, (str,)),
    "full_text_accessed": (True, (bool,)),
    "extraction_by": (True, (str,)),
    "second_extraction": (False, (str, type(None))),
    "human_verified": (True, (bool,)),
    "source_issues": (False, (str, type(None))),
    "instruments_status": (False, (str, type(None))),   # extracted | not_available | error (src/instruments.py); rederived (src/rederive.py)
    "instruments_check": (False, (dict, type(None))),   # recompute agreement and count flags
    "effect_scale": (False, (str, type(None))),   # per_higher (default when absent) | per_lower | unclear; PROTOCOL.md section 5
    "notes": (False, (str, type(None))),
}

GENE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9\-\.@]*$")


def validate_claim(claim):
    """Return a list of problems; an empty list means internally consistent."""
    problems = []
    for name, (required, types) in FIELDS.items():
        if name not in claim:
            if required:
                problems.append(f"missing required field {name}")
            continue
        if not isinstance(claim[name], types) or (isinstance(claim[name], bool) and bool not in types):
            problems.append(f"{name} has type {type(claim[name]).__name__}")
    unknown = set(claim) - set(FIELDS)
    if unknown:
        problems.append(f"unknown fields {sorted(unknown)}")
    if problems:
        return problems
    if not GENE_RE.match(claim["gene_symbol"]):
        problems.append("gene_symbol is not a plain gene symbol")
    if claim["exposure_type"] not in EXPOSURE_TYPES:
        problems.append(f"exposure_type {claim['exposure_type']!r} not allowed")
    if claim["instrument_scope"] not in INSTRUMENT_SCOPES:
        problems.append(f"instrument_scope {claim['instrument_scope']!r} not allowed")
    if claim.get("effect_scale") not in EFFECT_SCALES:
        problems.append(f"effect_scale {claim['effect_scale']!r} not allowed")
    if claim["estimate_scale"] not in SCALES:
        problems.append(f"estimate_scale {claim['estimate_scale']!r} not allowed")
    if claim["direction"] not in DIRECTIONS:
        problems.append(f"direction {claim['direction']!r} not allowed")
    if not claim["outcome_cohorts"]:
        problems.append("outcome_cohorts is empty; use ['unknown'] if not reported")
    bad = [c for c in claim["outcome_cohorts"] if c not in KNOWN_COHORTS]
    if bad:
        problems.append(f"unrecognised cohorts {bad}; map to a known name or 'other'")
    if not claim["source_location"].strip():
        problems.append("source_location is blank")
    est, lo, hi, p = claim.get("estimate"), claim.get("ci_low"), claim.get("ci_high"), claim.get("p_value")
    ratio = claim["estimate_scale"] in {"OR", "HR"}
    if est is not None:
        if ratio and est <= 0:
            problems.append("ratio estimate must be positive")
        if lo is not None and hi is not None and not (lo <= est <= hi):
            problems.append("estimate lies outside its confidence interval")
        null = 1.0 if ratio else 0.0
        if claim["direction"] == "risk_increasing" and est < null:
            problems.append("direction says risk_increasing but estimate is below null")
        if claim["direction"] == "protective" and est > null:
            problems.append("direction says protective but estimate is above null")
    if p is not None and not (0 < p <= 1):
        problems.append("p_value outside (0, 1]")
    if claim.get("coloc_pph4") is not None and not (0 <= claim["coloc_pph4"] <= 1):
        problems.append("coloc_pph4 outside [0, 1]")
    if claim.get("coloc_pph4") is not None and not claim["coloc_reported"]:
        problems.append("coloc_pph4 given but coloc_reported is false")
    if claim.get("n_instruments") is not None and claim["n_instruments"] < 1:
        problems.append("n_instruments must be at least 1")
    # Extracted lists must be complete; lists kept from the first extraction on claims without a usable
    # table may be variant IDs only (useful for re-deriving instruments from open exposure data later).
    full = claim.get("instruments_status") in {"extracted", "rederived"}
    for i, v in enumerate(claim.get("instruments") or []):
        if not isinstance(v, dict) or not v.get("rsid"):
            problems.append(f"instrument {i} lacks rsid")
        elif full and (not isinstance(v.get("beta_exposure"), (int, float)) or isinstance(v.get("beta_exposure"), bool)):
            problems.append(f"instrument {v['rsid']} lacks numeric beta_exposure")
    if full and not claim.get("instruments"):
        problems.append(f"instruments_status is {claim['instruments_status']} but no instruments listed")
    if claim.get("instruments_status") == "rederived" and not (claim.get("instruments_check") or {}).get("rederived"):
        problems.append("instruments_status is rederived but the re-derivation record is missing")
    return problems


def validate_table(claims):
    report, seen = {}, set()
    for claim in claims:
        cid = claim.get("claim_id", "<no id>")
        problems = validate_claim(claim)
        if cid in seen:
            problems.append("duplicate claim_id")
        seen.add(cid)
        if problems:
            report[cid] = problems
    return report


# Cohort contents used for overlap checks. Rahmioglu2023_meta (public file without 23andMe) contains UK Biobank,
# Estonian Biobank (EGCUT), deCODE, BioBank Japan, the Northern Finland Birth Cohort (NFBC) and 17 other smaller
# cohorts, per its Supplementary Table 1 (data/reference/rahmioglu2023_cohorts.json). FinnGen is not in it, but
# NFBC1966/1986 are FinnGen legacy cohorts, so FinnGen contains NFBC too.
EXPAND = {"Rahmioglu2023_meta": {"Rahmioglu2023_meta", "UK Biobank", "BioBank Japan", "Estonian Biobank",
                                 "deCODE", "NFBC", "other"},
          "FinnGen": {"FinnGen", "NFBC"}}
# Shared participants accepted as near-independent (Artem, 2026-10-31, before any replication lookup):
# NFBC is 0.8% of Rahmioglu's effective sample; estimated correlation of estimates about 0.009.
ACCEPTED_OVERLAP = {"NFBC": "near-independent (NFBC overlap, 167 cases / 5,235 controls)"}


def shared_cohorts(discovery_cohorts, replication_cohorts):
    def norm(cohorts):
        out = set()
        for c in cohorts:
            out |= EXPAND.get(c, {c})
        return out
    return norm(discovery_cohorts) & norm(replication_cohorts)


def independent(discovery_cohorts, replication_cohorts):
    """True when no cohort is shared, apart from accepted small overlaps, and neither side is unknown.

    Use independence_note() to label results that rest on an accepted overlap.
    """
    if "unknown" in discovery_cohorts or "unknown" in replication_cohorts:
        return False
    return not (shared_cohorts(discovery_cohorts, replication_cohorts) - set(ACCEPTED_OVERLAP))


def independence_note(discovery_cohorts, replication_cohorts):
    """'independent', a near-independent label naming the accepted overlap, or None if not independent."""
    if not independent(discovery_cohorts, replication_cohorts):
        return None
    shared = shared_cohorts(discovery_cohorts, replication_cohorts)
    return "; ".join(ACCEPTED_OVERLAP[c] for c in sorted(shared)) if shared else "independent"


def common_direction(claim):
    """+1 (higher level raises risk), -1 (lowers risk) or None, on the shared per-higher-level scale.

    Claims reported per lower level or inhibition are flipped; claims with an unclear scale or no clear
    direction give None and are left out of conflict checks and drug-direction matching.
    """
    sign = {"risk_increasing": 1, "protective": -1}.get(claim.get("direction"))
    scale = claim.get("effect_scale") or "per_higher"
    if sign is None or scale == "unclear":
        return None
    return -sign if scale == "per_lower" else sign


def log_effect(claim):
    est = claim.get("estimate")
    if est is None:
        return None
    return math.log(est) if claim["estimate_scale"] in {"OR", "HR"} else est


def summary(claims):
    return {
        "claims": len(claims),
        "papers": len({c["paper_id"] for c in claims}),
        "genes": len({c["gene_symbol"] for c in claims}),
        "headline_claims": sum(c["headline"] for c in claims),
        "claims_per_gene": dict(Counter(c["gene_symbol"] for c in claims).most_common()),
        "outcome_cohort_sets": dict(Counter("+".join(sorted(c["outcome_cohorts"])) for c in claims).most_common()),
        "exposure_types": dict(Counter(c["exposure_type"] for c in claims).most_common()),
        "instrument_scopes": dict(Counter(c["instrument_scope"] for c in claims).most_common()),
        "secondary_only_claims": sum(c["exposure_type"] in SECONDARY_ONLY_TYPES for c in claims),
        "claims_without_estimate": sum(c.get("estimate") is None for c in claims),
        "with_reported_instruments": sum(bool(c.get("instruments")) for c in claims),
        "instruments_status": dict(Counter(c.get("instruments_status") or "none" for c in claims).most_common()),
        "instruments_recompute_agree": sum(bool(((c.get("instruments_check") or {}).get("recompute") or {}).get("agree"))
                                           for c in claims),
        "independent_of_UKB_possible": sum(independent(c["outcome_cohorts"], ["UK Biobank"]) for c in claims),
        "full_text_accessed": sum(c["full_text_accessed"] for c in claims),
        "human_verified": sum(c["human_verified"] for c in claims),
        "second_extraction": dict(Counter((c.get("second_extraction") or "none").split(":")[0].split(" (")[0]
                                          for c in claims).most_common()),
        "claims_with_source_issues": sum(bool(c.get("source_issues")) for c in claims),
    }


def write_csv(claims, path):
    cols = [f for f in FIELDS if f != "instruments"]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=cols, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for c in claims:
            writer.writerow({k: ("; ".join(map(str, c[k])) if isinstance(c.get(k), list) else c.get(k))
                             for k in cols})


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["validate", "summary"])
    parser.add_argument("path")
    parser.add_argument("--csv")
    args = parser.parse_args()
    claims = json.loads(Path(args.path).read_text(encoding="utf-8"))
    report = validate_table(claims)
    if args.command == "validate":
        print(json.dumps({"claims": len(claims), "claims_with_problems": len(report), "problems": report}, indent=1))
        sys.exit(1 if report else 0)
    if report:
        print(f"{len(report)} claims have problems; run validate first", file=sys.stderr)
        sys.exit(1)
    print(json.dumps(summary(claims), indent=1))
    if args.csv:
        write_csv(claims, args.csv)


if __name__ == "__main__":
    main()
