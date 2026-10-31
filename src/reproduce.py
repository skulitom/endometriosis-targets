"""Gate B reproduction: recompute published MR estimates from FinnGen R10 and the papers' own instruments.

  python src/reproduce.py

Targets follow PROTOCOL.md section 2 (rule plus logged clarification). Each target gene is unmasked in
the FinnGen R10 file only, through an entry in data/reference/blinding_exemptions.json with a reason.
Writes results/gateb_reproduction.json.
"""
import json
import math
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from blinding import EXEMPTIONS, read_blinded  # noqa: E402
from instruments import recompute  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
R10 = ROOT / "data" / "upstream" / "finngen_R10_N14_ENDOMETRIOSIS.gz"
TARGETS = [("Tao2024-1", "primary"), ("Xie2026FI-1", "secondary: endpoint not named in paper")]


def add_exemption(gene, claim_id, role):
    entries = json.loads(EXEMPTIONS.read_text(encoding="utf-8")) if EXEMPTIONS.exists() else []
    if not any(e["file_key"] == "finngen_r10_endo" and e["gene"] == gene for e in entries):
        entries.append({"file_key": "finngen_r10_endo", "gene": gene, "claim_id": claim_id,
                        "reason": f"Gate B reproduction ({role}): FinnGen R10 is this claim's discovery outcome, and every "
                                  f"claim of {gene} includes FinnGen, so FinnGen near {gene} is never an independent "
                                  "replication source.",
                        "added": datetime.now(timezone.utc).isoformat()})
        EXEMPTIONS.parent.mkdir(parents=True, exist_ok=True)
        EXEMPTIONS.write_text(json.dumps(entries, indent=1) + "\n", encoding="utf-8", newline="\n")


def grch38_position(rsid):
    url = f"https://rest.ensembl.org/variation/human/{rsid}?content-type=application/json"
    with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "EndometriosisTargetCheck/0.1"}),
                                timeout=60) as r:
        data = json.loads(r.read())
    maps = [m for m in data.get("mappings", []) if m.get("assembly_name") == "GRCh38"
            and m["seq_region_name"] in {str(i) for i in range(1, 23)} | {"X"}]
    return [(m["seq_region_name"], m["start"]) for m in maps]


def neighbour_exemptions(claim, claims):
    """Claimed genes within 1 Mb of the claim's instruments, exempted only if all their claims include FinnGen."""
    genes = json.loads((ROOT / "data" / "reference" / "claimed_genes.json").read_text(encoding="utf-8"))["genes"]
    by_gene = {}
    for c in claims.values():
        by_gene.setdefault(c["gene_symbol"], []).append(c)
    safe, blocked = set(), set()
    for v in claim["instruments"]:
        for chrom, pos in grch38_position(v["rsid"]):
            for sym, g in genes.items():
                b = g["GRCh38"]
                if b["chrom"] == chrom and b["start"] - 1_000_000 <= pos <= b["end"] + 1_000_000:
                    ok = all("FinnGen" in c["outcome_cohorts"] for c in by_gene.get(sym, []))
                    (safe if ok else blocked).add(sym)
    return sorted(safe - blocked), sorted(blocked)


def main():
    claims = {c["claim_id"]: c for c in json.loads((ROOT / "data" / "literature" / "claims.json").read_text(encoding="utf-8"))}
    blocked_by = {}
    for cid, role in TARGETS:
        safe, blocked = neighbour_exemptions(claims[cid], claims)
        blocked_by[cid] = blocked
        for gene in sorted(set(safe) | {claims[cid]["gene_symbol"]} - set(blocked)):
            add_exemption(gene, cid, role)
    wanted = {v["rsid"]: cid for cid, _ in TARGETS for v in claims[cid]["instruments"]}
    rows = {}
    for r in read_blinded(R10, "GRCh38", "chrom", "pos", "Gate B reproduction of published estimates",
                          "finngen_r10_endo", row_filter=lambda p, i: any(x in wanted for x in p[i["rsids"]].split(","))):
        for x in r["rsids"].split(","):
            if x in wanted:
                rows[x] = r
    report = {"file": R10.name, "targets": {}}
    for cid, role in TARGETS:
        c = claims[cid]
        variants, missing = [], []
        for v in c["instruments"]:
            r = rows.get(v["rsid"])
            if r is None:
                missing.append(v["rsid"])
                continue
            variants.append({**v, "beta_outcome": float(r["beta"]), "se_outcome": float(r["sebeta"]),
                             "effect_allele_outcome": r["alt"].upper(), "other_allele_outcome": r["ref"].upper(),
                             "finngen_pos": f"{r['chrom']}:{r['pos']}"})
        est = recompute(variants) if variants and not missing else None
        reported = math.log(c["estimate"]) if c["estimate_scale"] == "OR" else c["estimate"]
        out = {"role": role, "gene": c["gene_symbol"], "reported_estimate": c["estimate"], "reported_log": reported,
               "n_instruments": len(c["instruments"]), "found_in_r10": len(variants), "missing": missing,
               "kept_masked_because_of": blocked_by[cid]}
        if est:
            diff = est[0] - reported
            out.update({"recomputed_log": est[0], "recomputed_se": est[1], "recomputed_or": math.exp(est[0]),
                        "difference": diff, "agree": abs(diff) <= max(0.02, 0.05 * abs(reported)),
                        "variants": [{k: v[k] for k in ("rsid", "effect_allele", "other_allele", "beta_exposure",
                                                        "beta_outcome", "effect_allele_outcome", "finngen_pos")}
                                     for v in variants]})
        report["targets"][cid] = out
        print(json.dumps({k: out.get(k) for k in ("role", "gene", "reported_estimate", "recomputed_or", "agree",
                                                   "found_in_r10", "missing", "kept_masked_because_of")}))
    (ROOT / "results" / "gateb_reproduction.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
