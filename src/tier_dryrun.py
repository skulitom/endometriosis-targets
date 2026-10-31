"""Dry run of the replication-source and tier rules on the claim table, before the freeze.

  python src/tier_dryrun.py

Uses no outcome data. For each primary claim it lists the eligible replication sources (PROTOCOL.md section 1,
independence map, and section 3, choice of replication source) and the best tier the claim could reach given its
instruments (section 5, tier rules). Writes results/tier_ceiling_dryrun.json.
"""
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from claims import common_direction, independence_note, independent  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CLAIMS = ROOT / "data" / "literature" / "claims.json"
GENES = ROOT / "data" / "reference" / "claimed_genes.json"
GATEB = ROOT / "results" / "gateb_positive_controls.json"
GATEB_LABEL = "not blind: Gate B positive-control locus"
WINDOW = 1_000_000
OUT = ROOT / "results" / "tier_ceiling_dryrun.json"

# Replication sources that passed the positive controls (results/gateb_positive_controls.json), with the
# cohorts they contain and their endometriosis case counts from each source's own documentation.
SOURCES = {
    "Rahmioglu 2023 (public, no 23andMe)": {"cohorts": ["Rahmioglu2023_meta"], "cases": 23492},
    "FinnGen R12": {"cohorts": ["FinnGen"], "cases": 20190},
    "Estonia + FinnGen": {"cohorts": ["Estonian Biobank", "FinnGen"], "cases": 19588},
    "Pan-UKBB females": {"cohorts": ["UK Biobank"], "cases": 4354},
}
MHC_GRCH37 = ("6", 28_477_797, 33_448_354)


def eligible_sources(claim):
    out = []
    for name, src in SOURCES.items():
        note = independence_note(claim["outcome_cohorts"], src["cohorts"])
        if note:
            out.append({"source": name, "independence": note, "cases": src["cases"]})
    return out


def primary_sources(eligible):
    """Most cases first; add any other eligible source that shares no participants with those already chosen."""
    chosen = []
    for e in sorted(eligible, key=lambda e: -e["cases"]):
        if all(independent(SOURCES[e["source"]]["cohorts"], SOURCES[c["source"]]["cohorts"]) for c in chosen):
            chosen.append(e)
    return chosen


def in_mhc(claim, genes):
    g = genes.get(claim["gene_symbol"], {}).get("GRCh37")
    chrom, lo, hi = MHC_GRCH37
    return bool(g) and g["chrom"] == chrom and g["start"] <= hi and g["end"] >= lo


def gate_b_lead_near(claim, genes, leads):
    """Gate B positive-control leads (GRCh37) within 1 Mb of the claim's gene: replication there is not blind.

    PROTOCOL.md section 5, qualifying condition 7. A gene without GRCh37 coordinates is treated as near a lead,
    so a missing position can never make a claim look blind."""
    g = genes.get(claim["gene_symbol"], {}).get("GRCh37")
    if not g:
        return ["position unknown"]
    return [lead["rsid"] for lead in leads
            if lead["chrom"] == g["chrom"] and g["start"] - WINDOW <= lead["pos"] <= g["end"] + WINDOW]


def ceiling(claim, eligible, genes, leads=()):
    status = claim.get("instruments_status")
    check = claim.get("instruments_check") or {}
    if status not in ("extracted", "rederived"):
        return "not testable", ["no usable instruments"]
    if not eligible:
        return "not testable", ["no independent source"]
    reasons = []
    if gate_b_lead_near(claim, genes, leads):
        reasons.append(GATEB_LABEL)
    if status == "rederived" and check.get("rederived", {}).get("analysis_role") != "primary_rederived":
        reasons.append("sensitivity-only instruments (re-derived)")
    if status == "extracted" and check.get("flags"):
        reasons.append("sensitivity-only instruments (published list problems)")
    if claim["instrument_scope"] != "cis":
        reasons.append(f"instrument scope {claim['instrument_scope']}")
    if in_mhc(claim, genes):
        reasons.append("MHC region")
    return ("Tier 1", []) if not reasons else ("Tier 2", reasons)


def exposure_family(claim):
    """Exposure type for the conflict rule: SMR variants count with their QTL type."""
    return claim["exposure_type"].removeprefix("SMR_")


def discovery_direction_conflicts(claims):
    """Genes whose published claims of one exposure type point both ways on the shared scale (discovery data only)."""
    groups = defaultdict(lambda: {1: [], -1: []})
    for c in claims:
        if c["exposure_type"] == "mQTL":
            continue
        d = common_direction(c)
        if d is not None:
            groups[(c["gene_symbol"], exposure_family(c))][d].append(c["claim_id"])
    return [{"gene": g, "exposure_type": t, "higher_level_raises_risk": v[1], "higher_level_lowers_risk": v[-1]}
            for (g, t), v in sorted(groups.items()) if v[1] and v[-1]]


def main():
    claims = json.loads(CLAIMS.read_text(encoding="utf-8"))
    genes = json.loads(GENES.read_text(encoding="utf-8"))["genes"]
    leads = json.loads(GATEB.read_text(encoding="utf-8"))["leads"]
    gate_b = [{"claim_id": c["claim_id"], "gene": c["gene_symbol"], "exposure_type": c["exposure_type"],
               "leads_within_1Mb": near}
              for c in claims if (near := gate_b_lead_near(c, genes, leads))]
    rows, per_gene = [], defaultdict(list)
    for c in claims:
        if c["exposure_type"] == "mQTL":
            continue
        eligible = eligible_sources(c)
        best, reasons = ceiling(c, eligible, genes, leads)
        primary = primary_sources(eligible) if best != "not testable" else []
        rows.append({"claim_id": c["claim_id"], "gene": c["gene_symbol"], "exposure_type": c["exposure_type"],
                     "discovery_cohorts": c["outcome_cohorts"], "eligible_sources": eligible,
                     "primary_sources": [p["source"] for p in primary],
                     "near_independent": any(p["independence"] != "independent" for p in primary),
                     "ceiling": best, "reasons": reasons})
        per_gene[c["gene_symbol"]].append(best)
    order = {"Tier 1": 0, "Tier 2": 1, "not testable": 2}
    gene_ceiling = {g: min(v, key=order.get) for g, v in sorted(per_gene.items())}
    testable = [g for g, t in gene_ceiling.items() if t != "not testable"]
    tested = [r for r in rows if r["ceiling"] != "not testable"]
    out = {
        "note": "Dry run of PROTOCOL.md source and tier rules before the freeze. No outcome data used: shows eligible "
                "sources and the best tier each claim could reach given its instruments.",
        "distinct_testable_genes": len(testable),
        "alpha_strict": 0.05 / len(testable) if testable else None,
        "claims_by_ceiling": dict(Counter(r["ceiling"] for r in rows)),
        "genes_by_ceiling": dict(Counter(gene_ceiling.values())),
        "primary_source_of_testable_claims": dict(Counter(" + ".join(r["primary_sources"]) for r in tested)),
        "near_independent_testable_claims": sum(r["near_independent"] for r in tested),
        "tier2_reasons": dict(Counter("; ".join(r["reasons"]) for r in rows if r["ceiling"] == "Tier 2")),
        "not_testable_reasons": dict(Counter("; ".join(r["reasons"]) for r in rows if r["ceiling"] == "not testable")),
        "genes_that_could_reach_tier1": [g for g, t in gene_ceiling.items() if t == "Tier 1"],
        "discovery_direction_conflicts": discovery_direction_conflicts(claims),
        "effect_scale_flags": {c["claim_id"]: c["effect_scale"] for c in claims if c.get("effect_scale")},
        "gate_b_locus_label": GATEB_LABEL,
        "gate_b_locus_claims": gate_b,
        "gene_ceiling": gene_ceiling,
        "claims": rows,
    }
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({k: v for k, v in out.items() if k not in ("gene_ceiling", "claims")}, indent=1))


if __name__ == "__main__":
    main()
