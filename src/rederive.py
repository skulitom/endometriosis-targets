"""Re-derive instruments for claims without published per-variant effects (PROTOCOL.md section 3).

Run with the project venv (needs numpy):
  .venv\\Scripts\\python -I src\\rederive.py unpack      # copy the needed members out of the 1000 Genomes and GTEx archives
  .venv\\Scripts\\python -I src\\rederive.py spotcheck   # 1000 Genomes structure and Ensembl allele-frequency check
  .venv\\Scripts\\python -I src\\rederive.py build       # rebuild every claim

Writes data/literature/instruments/rederive/rederived.json and results/rederive_report.json.
Exposure files hold no endometriosis data; nothing here reads an outcome file.
All parameters and readings follow PROTOCOL.md "Fixed details", written before any exposure file was opened.
"""
import gzip
import json
import math
import random
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / "data" / "upstream" / "exposure"
REF_DIR = EXP / "1kg_v3"
GTEX_DIR = EXP / "gtex_v8"
RD = ROOT / "data" / "literature" / "instruments" / "rederive"
OUT = RD / "rederived.json"
REPORT = ROOT / "results" / "rederive_report.json"
SPOT = ROOT / "results" / "ld_reference_check.json"
GENES = ROOT / "data" / "reference" / "claimed_genes.json"
ENSEMBL = "https://rest.ensembl.org"
UA = {"User-Agent": "EndometriosisTargetCheck/0.1"}
SEED = 20261031
COMP = str.maketrans("ACGT", "TGCA")

REF_MEMBERS = {"EUR.bed", "EUR.bim", "EUR.fam"}
GTEX_MEMBERS = {f"{t}.v8.{k}.txt.gz" for t in ("Whole_Blood", "Uterus") for k in ("egenes", "signif_variant_gene_pairs")}

# Column layouts, checked against each file's header on reading.
DATASETS = {
    "zhao2023": {"sep": "\t", "cols": {"chrom": "chromosome", "pos": "base_pair_location", "ea": "effect_allele",
                                       "oa": "other_allele", "beta": "beta", "se": "standard_error",
                                       "eaf": "effect_allele_frequency", "p": "p_value", "rsid": "rsid"}},
    "sun2018": {"sep": "\t", "cols": {"chrom": "chromosome", "pos": "base_pair_location", "ea": "effect_allele",
                                      "oa": "other_allele", "beta": "beta", "se": "standard_error", "p": "p_value",
                                      "log10p": "log(P)"}},
    "aholaolli2017": {"sep": None, "cols": {"rsid": "MarkerName", "chrom": "Chromosome", "pos": "Position",
                                            "oa": "OtherAllele", "ea": "EffectAllele", "beta": "Effect",
                                            "se": "StdErr", "p": "P.value"}},
    "glgc2021": {"sep": "\t", "cols": {"rsid": "rsID", "chrom": "CHROM", "pos": "POS_b37", "oa": "REF", "ea": "ALT",
                                       "eaf": "POOLED_ALT_AF", "beta": "EFFECT_SIZE", "se": "SE", "p": "pvalue",
                                       "neglog10p": "pvalue_neg_log10"}},
}

FILES = {  # download key -> (dataset layout, filename)
    **{f"zhao2023_{g}": ("zhao2023", f"GCST902{acc}.tsv.gz") for g, acc in [
        ("CCL4", "74770"), ("CD274", "74832"), ("CDCP1", "74775"), ("CXCL11", "74781"), ("CXCL5", "74782"),
        ("CXCL6", "74783"), ("FGF21", "74788"), ("IL15RA", "74800"), ("IL17A", "74801"), ("NGF", "74762"),
        ("SLAMF1", "74835"), ("TNFSF10", "74843"), ("TNFSF11", "74844")]},
    "sun2018_CTSH_8465_52_3": ("sun2018", "GCST90240624_GRCh37.tsv.gz"),
    "sun2018_CTSH_3737_6_3": ("sun2018", "GCST90240625_GRCh37.tsv.gz"),
    "aholaolli2017_TRAIL": ("aholaolli2017", "TRAIL.data.gz"),
    "glgc2021_logTG_EUR": ("glgc2021", "logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
}

TWO_SAMPLE_MR_KB = 10_000  # TwoSampleMR clump_data default distance, used where only r2 is stated
PROXY_R2 = 0.8
PROXY_KB = 500
SMR_P = 5e-8


def norm_chrom(c):
    c = str(c).strip().upper().removeprefix("CHR")
    return {"23": "X", "24": "Y"}.get(c, c)


def palindromic(a, b):
    return {a, b} in ({"A", "T"}, {"C", "G"})


def f_stat(beta, se):
    return (beta / se) ** 2


# ---------- archives ----------

def unpack_members(archive, wanted, dest):
    """Copy wanted regular-file members (by basename) out of a tar archive; nothing else is written."""
    dest.mkdir(parents=True, exist_ok=True)
    done = {}
    with tarfile.open(archive) as tar:
        for member in tar:
            base = Path(member.name).name
            if base not in wanted:
                continue
            if not member.isfile():
                raise ValueError(f"{member.name} is not a regular file")
            if base in done:
                raise ValueError(f"{base} appears twice in {archive.name}")
            target = dest / base
            if target.exists() and target.stat().st_size == member.size:
                done[base] = {"bytes": member.size, "status": "existing"}
                continue
            partial = target.with_name(base + ".partial")
            with tar.extractfile(member) as src, partial.open("xb") as out:
                for chunk in iter(lambda: src.read(1 << 20), b""):
                    out.write(chunk)
            if partial.stat().st_size != member.size:
                raise ValueError(f"{base}: size differs from archive entry")
            partial.rename(target)
            done[base] = {"bytes": member.size, "status": "extracted", "archive_member": member.name}
    missing = wanted - set(done)
    if missing:
        raise ValueError(f"{archive.name} lacks {sorted(missing)}")
    return done


def unpack():
    out = {"1kg_v3": unpack_members(EXP / "1kg.v3.tgz", REF_MEMBERS, REF_DIR),
           "gtex_v8": unpack_members(EXP / "GTEx_Analysis_v8_eQTL.tar", GTEX_MEMBERS, GTEX_DIR)}
    print(json.dumps(out, indent=1))


# ---------- reference panel ----------

class RefPanel:
    """PLINK binary reference (1000 Genomes phase 3 European, GRCh37)."""

    def __init__(self, prefix=REF_DIR / "EUR"):
        prefix = Path(prefix)
        self.n = sum(1 for _ in open(str(prefix) + ".fam", encoding="ascii"))
        self.chrom, self.rsid, self.a1, self.a2 = [], [], [], []
        pos = []
        self.by_rsid, self.by_pos = {}, {}
        with open(str(prefix) + ".bim", encoding="ascii") as handle:
            for i, line in enumerate(handle):
                c, rs, _, p, a1, a2 = line.split()
                c = norm_chrom(c)
                p = int(p)
                self.chrom.append(c)
                self.rsid.append(rs)
                self.a1.append(a1.upper())
                self.a2.append(a2.upper())
                pos.append(p)
                self.by_rsid.setdefault(rs, i)
                self.by_pos.setdefault((c, p), []).append(i)
        self.pos = np.array(pos, dtype=np.int64)
        self.m = len(pos)
        self.bpr = (self.n + 3) // 4
        bed = Path(str(prefix) + ".bed")
        with bed.open("rb") as handle:
            magic = handle.read(3)
        if magic != b"\x6c\x1b\x01":
            raise ValueError("Not a SNP-major PLINK .bed file")
        if bed.stat().st_size != 3 + self.m * self.bpr:
            raise ValueError(".bed size does not match .bim and .fam")
        self.bed = np.memmap(bed, dtype=np.uint8, mode="r", offset=3, shape=(self.m, self.bpr))

    def close(self):
        mm = getattr(self.bed, "_mmap", None)
        self.bed = None
        if mm is not None:
            mm.close()

    def dosages(self, idxs):
        """A1 allele dosage (0, 1, 2; NaN missing) for reference rows idxs, shape (k, n)."""
        raw = np.asarray(self.bed[np.asarray(idxs, dtype=np.int64)])
        codes = np.stack([(raw >> (2 * k)) & 3 for k in range(4)], axis=2).reshape(len(raw), -1)[:, : self.n]
        return np.array([2.0, np.nan, 1.0, 0.0])[codes]

    def a1_freq(self, idxs):
        return np.nanmean(self.dosages(idxs), axis=1) / 2

    def find(self, rsid=None, chrom=None, pos=None, alleles=None):
        """Reference row for a variant: by rsID, else by position with matching alleles (either strand)."""
        if rsid and rsid in self.by_rsid:
            return self.by_rsid[rsid]
        if chrom is None or pos is None:
            return None
        for i in self.by_pos.get((norm_chrom(chrom), int(pos)), []):
            ref = {self.a1[i], self.a2[i]}
            if alleles is None or set(alleles) == ref or {a.translate(COMP) for a in alleles} == ref:
                return i
        return None

    def freq_of(self, i, allele):
        """Frequency of allele at reference row i, or None if the alleles do not match."""
        f = float(self.a1_freq([i])[0])
        a1, a2 = self.a1[i], self.a2[i]
        if allele == a1 or (allele.translate(COMP) == a1 and not palindromic(a1, a2)):
            return f
        if allele == a2 or (allele.translate(COMP) == a2 and not palindromic(a1, a2)):
            return 1 - f
        return None


def standardize(d):
    mean = np.nanmean(d, axis=1, keepdims=True)
    d = np.where(np.isnan(d), mean, d) - mean
    norm = np.sqrt((d ** 2).sum(axis=1, keepdims=True))
    norm[norm == 0] = np.inf  # monomorphic in the panel: correlation 0
    return d / norm


def clump(cands, panel, r2_max, kb):
    """Greedy clumping by ascending p (ties by position); returns (kept, dropped_absent_from_reference)."""
    absent = [c for c in cands if c.get("ref") is None]
    present = sorted((c for c in cands if c.get("ref") is not None),
                     key=lambda c: (-c["neglog10p"], c["chrom"], c["pos"]))
    kept = []
    by_chrom = {}
    for c in present:
        by_chrom.setdefault(c["chrom"], []).append(c)
    for items in by_chrom.values():
        z = standardize(panel.dosages([c["ref"] for c in items]))
        pos = np.array([c["pos"] for c in items], dtype=np.int64)
        alive = np.ones(len(items), dtype=bool)
        for i in range(len(items)):
            if not alive[i]:
                continue
            kept.append(items[i])
            alive[i] = False
            near = np.nonzero(alive & (np.abs(pos - pos[i]) <= kb * 1000))[0]
            if near.size:
                r2 = (z[near] @ z[i]) ** 2
                alive[near[r2 > r2_max + 1e-9]] = False  # tolerance for floating-point noise
    kept.sort(key=lambda c: (-c["neglog10p"], c["chrom"], c["pos"]))
    return kept, absent


def r2_between(panel, i, others):
    z = standardize(panel.dosages([i] + list(others)))
    return (z[1:] @ z[0]) ** 2


# ---------- exposure files ----------

def parse_float(x):
    try:
        v = float(x)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def scan(key, p_max=None, rsids=frozenset(), chroms=None):
    """One pass over an exposure file: rows with p <= p_max (optionally on given chromosomes) or listed rsIDs."""
    layout, filename = FILES[key]
    spec = DATASETS[layout]
    path = EXP / filename
    rows = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\r\n").split(spec["sep"])
        try:
            col = {name: header.index(c) for name, c in spec["cols"].items()}
        except ValueError as exc:
            raise ValueError(f"{filename}: expected column missing ({exc})") from None
        width = max(col.values()) + 1
        for line in handle:
            parts = line.rstrip("\r\n").split(spec["sep"])
            if len(parts) < width:
                continue
            rsid = parts[col["rsid"]] if "rsid" in col else None
            wanted = rsid in rsids if rsid else False
            if "neglog10p" in col:
                nl = parse_float(parts[col["neglog10p"]])
            elif "log10p" in col:
                v = parse_float(parts[col["log10p"]])
                nl = -v if v is not None else None
            else:
                p = parse_float(parts[col["p"]])
                nl = -math.log10(p) if p and p > 0 else (math.inf if p == 0 else None)
            if nl is None:
                continue
            passes = p_max is not None and nl >= -math.log10(p_max)
            if not (passes or wanted):
                continue
            chrom = norm_chrom(parts[col["chrom"]])
            if passes and not wanted and chroms is not None and chrom not in chroms:
                continue
            beta, se = parse_float(parts[col["beta"]]), parse_float(parts[col["se"]])
            if beta is None or not se:
                continue
            rows.append({"rsid": rsid if rsid not in (None, "", "NA", ".") else None, "chrom": chrom,
                         "pos": int(float(parts[col["pos"]])), "ea": parts[col["ea"]].upper(),
                         "oa": parts[col["oa"]].upper(), "beta": beta, "se": se,
                         "eaf": parse_float(parts[col["eaf"]]) if "eaf" in col else None,
                         "p": parse_float(parts[col["p"]]), "neglog10p": nl})
    return rows


def attach_reference(rows, panel):
    """Fill rsID (from the reference when the file has positions only) and reference row."""
    for r in rows:
        i = panel.find(r["rsid"], r["chrom"], r["pos"], {r["ea"], r["oa"]})
        if i is not None and r["rsid"] is None:
            # position match: alleles were checked in find()
            r["rsid"] = panel.rsid[i]
        if i is not None and r["rsid"] != panel.rsid[i]:
            i = None
        r["ref"] = i
    return rows


def dedupe(rows):
    """Keep the most significant row per rsID (multi-allelic sites); rows without rsID stay."""
    best = {}
    for r in rows:
        k = r["rsid"] or (r["chrom"], r["pos"], r["ea"], r["oa"])
        if k not in best or r["neglog10p"] > best[k]["neglog10p"]:
            best[k] = r
    return list(best.values())


def variant_record(r, **extra):
    out = {k: r[k] for k in ("rsid", "chrom", "pos", "ea", "oa", "eaf", "beta", "se", "p")}
    out["neglog10p"] = round(r["neglog10p"], 4) if math.isfinite(r["neglog10p"]) else None
    out["F"] = round(f_stat(r["beta"], r["se"]), 2)
    out.update(extra)
    return out


# ---------- claim definitions ----------

def load_rules():
    rules = {}
    for name in ("rules_1.json", "rules_2.json"):
        for paper in json.loads((RD / name).read_text(encoding="utf-8")):
            for pc in paper.get("per_claim", []):
                rules[pc["claim_id"]] = (paper, pc)
    return rules


ZEN_TISSUE = {"Whole_Blood_GTEx_V8": "Whole_Blood", "Uterus_GTEx_V8": "Uterus"}
LI_READINGS = [("gene region +/- 100 Mb (as printed)", 100_000_000), ("gene region +/- 50 Mb (100 Mb total)", 50_000_000),
               ("gene region +/- 100 kb (kb/Mb misprint)", 100_000)]
ARU_READINGS = [("aptamer CTSH.8465.52.3 (GCST90240624)", "sun2018_CTSH_8465_52_3"),
                ("aptamer CTSH.3737.6.3 (GCST90240625)", "sun2018_CTSH_3737_6_3")]
WEI = ["WeiY2024ITT-%d" % i for i in range(1, 7)]
GUO = ["Guo2025Med-%d" % i for i in range(1, 6)]
GAO = ["Gao2025Med-%d" % i for i in range(1, 4)]
ZENG = ["Zeng2024Gene-%d" % i for i in (1, 2, 3, 9, 10, 11)]
LI = ["Li2026tg-%d" % i for i in range(1, 8)]
# Published variant list conflicts with the analysed set (protocol amendment 2026-10-30: sensitivity only).
LIST_CONFLICTS = {"Gao2025Med-3": "Table 1 lists 7 SLAMF1 variants; Supplementary Tables 2-3 (IVW) report 3"}
# Reported count used for the check: the figure or table the claim's estimate came from.
LI_REPORTED_SOURCE = "Fig 5 nsnp (the claim's estimate is from Fig 5)"


def gene_region_37(gene, pc):
    region = pc.get("gene_region")
    if region:
        chrom, span = region.removeprefix("chr").split(":")
        start, end = span.split("-")
        return norm_chrom(chrom), int(start), int(end), "paper S4 file (read as GRCh37)"
    g = json.loads(GENES.read_text(encoding="utf-8"))["genes"][gene]["GRCh37"]
    return norm_chrom(g["chrom"]), int(g["start"]), int(g["end"]), "Ensembl GRCh37"


def provenance(key):
    _, filename = FILES[key]
    prov = EXP / (filename + ".provenance.json")
    return json.loads(prov.read_text(encoding="utf-8"))["sha256"] if prov.exists() else None


def zhao_key(gene):
    key = f"zhao2023_{gene}"
    if key not in FILES:
        raise KeyError(f"No Zhao 2023 file for {gene}")
    return key


# ---------- route 2: stated rule ----------

def apply_rule(rows, panel, p_max, r2, kb, f_min=None, palin_eaf=None, region=None):
    """Window -> p threshold -> (MAF: none stated here) -> clumping -> F filter. Returns kept, log."""
    log = {}
    cands = [r for r in rows if r["neglog10p"] >= -math.log10(p_max)]
    if region:
        chrom, lo, hi = region
        cands = [r for r in cands if r["chrom"] == chrom and lo <= r["pos"] <= hi]
    cands = dedupe(cands)
    log["pass_p"] = len(cands)
    if palin_eaf:
        lo_f, hi_f = palin_eaf
        keep = []
        dropped = 0
        for r in cands:
            if palindromic(r["ea"], r["oa"]):
                eaf = r["eaf"]
                if eaf is None and r.get("ref") is not None:
                    eaf = panel.freq_of(r["ref"], r["ea"])
                if eaf is None or lo_f < eaf < hi_f:
                    dropped += 1
                    continue
            keep.append(r)
        log["palindromic_dropped"] = dropped
        cands = keep
    kept, absent = clump(cands, panel, r2, kb)
    log["absent_from_reference"] = len(absent)
    log["after_clumping"] = len(kept)
    if f_min is not None:
        before = len(kept)
        kept = [r for r in kept if f_stat(r["beta"], r["se"]) > f_min]
        log["dropped_F"] = before - len(kept)
    return kept, log


def route2_claims(rules, panel, cache):
    out = []
    # WeiY 2024: Zhao proteins, genome-wide
    for cid in WEI:
        paper, pc = rules[cid]
        key = zhao_key(pc["gene"])
        rows = cache[key]
        kept, log = apply_rule(rows, panel, paper["p_threshold"], paper["clump_r2"], paper["clump_kb"],
                               f_min=paper["f_stat_min"])
        out.append(result(cid, pc, "stated_rule", key, kept, log,
                          params={"p": paper["p_threshold"], "r2": paper["clump_r2"], "kb": paper["clump_kb"],
                                  "F_min": paper["f_stat_min"], "window": "genome-wide",
                                  "ld_reference": "1000 Genomes EUR (unstated in paper)"}))
    # Xu 2025: TRAIL
    cid = "Xu2025JCMM-1"
    paper, pc = rules[cid]
    key = "aholaolli2017_TRAIL"
    kept, log = apply_rule(cache[key], panel, paper["p_threshold"], paper["clump_r2"], paper["clump_kb"])
    out.append(result(cid, pc, "stated_rule", key, kept, log,
                      params={"p": paper["p_threshold"], "r2": paper["clump_r2"], "kb": paper["clump_kb"],
                              "window": "genome-wide", "ld_reference": "1000 Genomes EUR (unstated in paper)"}))
    # Aru 2025: CTSH, two aptamer readings
    cid = "Aru2025cath-1"
    paper, pc = rules[cid]
    readings = []
    for label, key in ARU_READINGS:
        kept, log = apply_rule(cache[key], panel, paper["p_threshold"], paper["clump_r2"], TWO_SAMPLE_MR_KB,
                               f_min=paper["f_stat_min"], palin_eaf=(0.42, 0.58))
        readings.append((label, key, kept, log))
    out.append(choose_reading(cid, pc, readings,
                              params={"p": paper["p_threshold"], "r2": paper["clump_r2"],
                                      "kb": f"{TWO_SAMPLE_MR_KB} (TwoSampleMR default; unstated)",
                                      "F_min": paper["f_stat_min"],
                                      "palindromic": "drop A/T, C/G with 0.42 < EAF < 0.58; EAF from 1000 Genomes EUR "
                                                     "because the exposure file has none",
                                      "window": "genome-wide", "ld_reference": "1000 Genomes EUR"}))
    # Li 2026: GLGC TG, drug-target regions, three window readings
    key = "glgc2021_logTG_EUR"
    for cid in LI:
        paper, pc = rules[cid]
        chrom, start, end, region_source = gene_region_37(pc["gene"], pc)
        readings = []
        for label, half in LI_READINGS:
            kept, log = apply_rule(cache[key], panel, paper["p_threshold"], paper["clump_r2"], paper["clump_kb"],
                                   region=(chrom, max(0, start - half), end + half))
            readings.append((label, key, kept, log))
        out.append(choose_reading(cid, pc, readings,
                                  params={"p": paper["p_threshold"], "r2": paper["clump_r2"], "kb": paper["clump_kb"],
                                          "gene_region_GRCh37": f"{chrom}:{start}-{end}",
                                          "gene_region_source": region_source,
                                          "n_reported_source": LI_REPORTED_SOURCE,
                                          "ld_reference": "1000 Genomes EUR (unstated in paper)"}))
    return out


def choose_reading(cid, pc, readings, params):
    reported = pc.get("n_instruments_reported")
    matches = [r for r in readings if len(r[2]) == reported]
    chosen = matches[0] if len(matches) == 1 else readings[0]
    label, key, kept, log = chosen
    res = result(cid, pc, "stated_rule", key, kept, log, params=params)
    res["readings"] = [{"reading": lab, "dataset": k, "n_rebuilt": len(kp), "log": lg} for lab, k, kp, lg in readings]
    res["reading_used"] = label
    if len(matches) != 1:
        res["reasons"].append(f"{len(matches)} readings reproduce the reported count; first reading used")
        res["analysis_role"] = "sensitivity_only"
    return res


def result(cid, pc, route, key, kept, log, params, extra_variant=None):
    reported = pc.get("n_instruments_reported")
    res = {"claim_id": cid, "gene": pc["gene"], "route": route, "dataset": key, "dataset_file": FILES[key][1],
           "dataset_sha256": provenance(key), "parameters": params, "n_reported": reported,
           "n_rebuilt": len(kept), "steps": log, "reasons": [],
           "variants": [variant_record(r, **(extra_variant or {})) for r in kept]}
    if route == "stated_rule":
        res["checks"] = {"count_matches": len(kept) == reported}
        if len(kept) != reported:
            res["reasons"].append(f"rebuilt {len(kept)} variants vs {reported} reported")
    res["analysis_role"] = "primary_rederived" if not res["reasons"] else "sensitivity_only"
    return res


# ---------- route 1: known variants ----------

def route1_zhao(rules, panel, cache, scanner):
    reported_f = json.loads((RD / "reported_f.json").read_text(encoding="utf-8"))["claims"]
    out = []
    for cid in GUO + GAO:
        paper, pc = rules[cid]
        listed = pc.get("instrument_rsids")
        key = zhao_key(pc["gene"])
        rows = {r["rsid"]: r for r in dedupe([r for r in cache[key] if r["rsid"] in set(listed)])}
        kept, missing, proxies = [], [], []
        for rs in listed:
            if rs in rows:
                kept.append((rows[rs], {}))
            else:
                missing.append(rs)
        for rs in missing:
            proxy = find_proxy(rs, key, panel, scanner)
            if proxy:
                kept.append(proxy)
                proxies.append({"missing": rs, "proxy": proxy[0]["rsid"], "r2": proxy[1]["proxy_r2"]})
        res = {"claim_id": cid, "gene": pc["gene"], "route": "known_variants", "dataset": key,
               "dataset_file": FILES[key][1], "dataset_sha256": provenance(key),
               "parameters": {"listed_variants": listed, "source": paper.get("source_location")},
               "n_reported": pc.get("n_instruments_reported"), "n_listed": len(listed),
               "n_rebuilt": len(kept), "reasons": [], "proxies": proxies,
               "variants": [variant_record(r, **x) for r, x in kept]}
        checks = {"all_found": not missing, "missing": missing}
        if cid in reported_f:
            fmap = {x["rsid"]: x["F_reported"] for x in reported_f[cid]}
            comp = []
            for v in res["variants"]:
                rs = v.get("proxy_for") or v["rsid"]
                if rs in fmap:
                    comp.append({"rsid": rs, "F_reported": fmap[rs], "F_rebuilt": v["F"],
                                 "within_10pct": abs(v["F"] - fmap[rs]) <= 0.10 * fmap[rs]})
            checks["F_comparison"] = comp
            checks["F_match"] = bool(comp) and all(c["within_10pct"] for c in comp) and len(comp) == len(fmap)
            if not checks["F_match"]:
                res["reasons"].append("reported F statistics not reproduced within 10%")
        if missing:
            res["reasons"].append(f"{len(missing)} listed variants not in the exposure file"
                                  + (f" ({len(proxies)} replaced by LD proxies)" if proxies else ""))
        if cid in LIST_CONFLICTS:
            res["reasons"].append("published list conflicts with the analysed set: " + LIST_CONFLICTS[cid])
        if pc.get("n_instruments_reported") is not None and pc["n_instruments_reported"] != len(listed):
            res["reasons"].append(f"paper reports {pc['n_instruments_reported']} variants but lists {len(listed)}")
        res["checks"] = checks
        res["analysis_role"] = "primary_rederived" if not res["reasons"] else "sensitivity_only"
        out.append(res)
    return out


def find_proxy(rsid, key, panel, scanner):
    """Best LD proxy (r2 >= 0.8, within 500 kb, 1000 Genomes EUR) present in the exposure file."""
    i = panel.by_rsid.get(rsid)
    if i is None:
        return None
    chrom, pos = panel.chrom[i], panel.pos[i]
    near = [j for j in range(max(0, i - 5000), min(panel.m, i + 5000))
            if j != i and panel.chrom[j] == chrom and abs(int(panel.pos[j]) - int(pos)) <= PROXY_KB * 1000]
    if not near:
        return None
    r2 = r2_between(panel, i, near)
    cands = {panel.rsid[j]: float(v) for j, v in zip(near, r2) if v >= PROXY_R2}
    if not cands:
        return None
    rows = dedupe([r for r in scanner(key, rsids=frozenset(cands)) if r["rsid"] in cands])
    if not rows:
        return None
    best = max(rows, key=lambda r: (cands[r["rsid"]], r["neglog10p"]))
    return best, {"proxy_for": rsid, "proxy_r2": round(cands[best["rsid"]], 4)}


ENSEMBL_CACHE = RD / "ensembl_grch38_top_snps.json"


def ensembl_grch38(rsids):
    """GRCh38 position and alleles from Ensembl; the raw answer is kept so reruns use the same positions."""
    cached = json.loads(ENSEMBL_CACHE.read_text(encoding="utf-8")) if ENSEMBL_CACHE.exists() else None
    if cached and set(cached["requested"]) == set(rsids):
        data = cached["response"]
    else:
        for attempt in range(6):
            req = urllib.request.Request(f"{ENSEMBL}/variation/homo_sapiens", method="POST",
                                         data=json.dumps({"ids": sorted(rsids)}).encode(),
                                         headers={**UA, "Content-Type": "application/json", "Accept": "application/json"})
            try:
                with urllib.request.urlopen(req, timeout=120) as handle:
                    data = json.load(handle)
                break
            except (OSError, ValueError) as exc:
                if attempt == 5:
                    raise
                print(f"Ensembl lookup: {exc}; retrying", file=sys.stderr)
                time.sleep(15 * (attempt + 1))
        ENSEMBL_CACHE.write_text(json.dumps({"retrieved_at": time.strftime("%Y-%m-%d"), "endpoint": f"{ENSEMBL}/variation/homo_sapiens",
                                             "requested": sorted(rsids), "response": data}, indent=1) + "\n", encoding="utf-8")
    out = {}
    for rs, rec in data.items():
        for m in rec.get("mappings", []):
            if m.get("assembly_name") == "GRCh38" and m.get("seq_region_name") in {str(c) for c in range(1, 23)} | {"X"}:
                out[rs] = {"chrom": m["seq_region_name"], "pos": int(m["start"]), "alleles": m["allele_string"].split("/")}
                break
    return out


def route1_gtex(rules):
    claims = [(cid, *rules[cid]) for cid in ZENG]
    loc = ensembl_grch38({pc["top_snp"] for _, _, pc in claims})
    wanted = {}
    for _, _, pc in claims:
        wanted.setdefault(ZEN_TISSUE[pc["exposure_id"]], set()).add(pc["probe"])
    pairs = {t: gtex_rows(t, probes) for t, probes in wanted.items()}
    out = []
    for cid, paper, pc in claims:
        tissue = ZEN_TISSUE[pc["exposure_id"]]
        rs = pc["top_snp"]
        probe = pc["probe"]
        res = {"claim_id": cid, "gene": pc["gene"], "route": "known_variants", "dataset": "gtex_v8_eqtl",
               "dataset_file": f"GTEx_Analysis_v8_eQTL.tar > {tissue}.v8.signif_variant_gene_pairs.txt.gz",
               "dataset_sha256": json.loads((EXP / "GTEx_Analysis_v8_eQTL.tar.provenance.json").read_text())["sha256"],
               "parameters": {"tissue": tissue, "probe": probe, "listed_variants": [rs],
                              "source": "Zeng 2024 Table S1 (SMR top SNP)"},
               "n_reported": pc.get("n_instruments_reported"), "n_listed": 1, "reasons": [], "variants": []}
        where = loc.get(rs)
        hit = None
        if where:
            prefix = f"chr{where['chrom']}_{where['pos']}_"
            for row in pairs[tissue][probe]:
                vid = row["variant_id"]
                if vid.startswith(prefix):
                    _, _, ref, alt, _ = vid.split("_")
                    if {ref, alt} <= set(where["alleles"]):
                        hit = row
                        break
        top = gtex_top(tissue, probe)
        if hit:
            _, pos, ref, alt, _ = hit["variant_id"].split("_")
            beta, se, p = float(hit["slope"]), float(hit["slope_se"]), float(hit["pval_nominal"])
            res["variants"] = [{"rsid": rs, "chrom": where["chrom"], "pos_GRCh38": int(pos), "ea": alt, "oa": ref,
                                "maf": float(hit["maf"]), "beta": beta, "se": se, "p": p,
                                "F": round(f_stat(beta, se), 2), "gtex_variant_id": hit["variant_id"]}]
            if p >= SMR_P:
                res["reasons"].append(f"top SNP exposure p = {p:.2g}, not below the SMR default 5e-8")
        else:
            res["reasons"].append("listed top SNP not among GTEx v8 significant pairs for this gene and tissue"
                                  if where else "rsID not found in Ensembl GRCh38")
        res["n_rebuilt"] = len(res["variants"])
        res["checks"] = {"all_found": bool(hit), "gtex_top_variant": top,
                         "listed_is_gtex_top": bool(top and hit and top["variant_id"] == hit["variant_id"])}
        res["analysis_role"] = "primary_rederived" if not res["reasons"] else "sensitivity_only"
        out.append(res)
    return out


def gtex_rows(tissue, probes):
    """Significant variant-gene pairs for the given genes (Ensembl IDs without version) in one tissue."""
    by_gene = {p: [] for p in probes}
    with gzip.open(GTEX_DIR / f"{tissue}.v8.signif_variant_gene_pairs.txt.gz", "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        gi = header.index("gene_id")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            gene = parts[gi].split(".")[0]
            if gene in by_gene:
                by_gene[gene].append(dict(zip(header, parts)))
    return by_gene


def gtex_top(tissue, probe):
    with gzip.open(GTEX_DIR / f"{tissue}.v8.egenes.txt.gz", "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            if row["gene_id"].split(".")[0] == probe:
                return {"variant_id": row["variant_id"], "rsid": row.get("rs_id_dbSNP151_GRCh38p7"),
                        "pval_nominal": float(row["pval_nominal"]), "qval": float(row["qval"])}
    return None


# ---------- 1000 Genomes check ----------

def spotcheck():
    panel = RefPanel()
    rng = random.Random(SEED)
    pool = [i for i in range(panel.m) if panel.rsid[i].startswith("rs") and panel.chrom[i] not in ("X", "Y")
            and len(panel.a1[i]) == 1 and len(panel.a2[i]) == 1]
    sample = sorted(rng.sample(pool, 200))
    freqs = panel.a1_freq(sample)
    ids = [panel.rsid[i] for i in sample]
    ens = {}
    for start in range(0, len(ids), 25):
        batch = ids[start:start + 25]
        for attempt in range(5):
            req = urllib.request.Request(f"{ENSEMBL}/variation/homo_sapiens?pops=1", method="POST",
                                         data=json.dumps({"ids": batch}).encode(),
                                         headers={**UA, "Content-Type": "application/json", "Accept": "application/json"})
            try:
                with urllib.request.urlopen(req, timeout=300) as handle:
                    ens.update(json.load(handle))
                break
            except (OSError, ValueError) as exc:
                if attempt == 4:
                    raise
                print(f"Ensembl batch {start // 25 + 1}: {exc}; retrying", file=sys.stderr)
                time.sleep(10 * (attempt + 1))
        time.sleep(1)
    rows, found, agree, allele_mismatch = [], 0, 0, 0
    for i, f, rs in zip(sample, freqs, ids):
        rec = ens.get(rs)
        eur = {p["allele"]: p["frequency"] for p in (rec or {}).get("populations", [])
               if p.get("population") == "1000GENOMES:phase_3:EUR"}
        row = {"rsid": rs, "chrom": panel.chrom[i], "pos_GRCh37": int(panel.pos[i]), "a1": panel.a1[i],
               "a2": panel.a2[i], "a1_freq_panel": round(float(f), 4)}
        if eur:
            found += 1
            if panel.a1[i] in eur:
                row["a1_freq_ensembl"] = eur[panel.a1[i]]
            elif panel.a2[i] in eur:
                row["a1_freq_ensembl"] = round(1 - eur[panel.a2[i]], 4)
            if "a1_freq_ensembl" in row and {panel.a1[i], panel.a2[i]} >= set(eur):
                row["diff"] = round(abs(row["a1_freq_ensembl"] - row["a1_freq_panel"]), 4)
                agree += row["diff"] <= 0.01
            else:
                allele_mismatch += 1
                row["allele_mismatch"] = sorted(eur)
        rows.append(row)
    passed = found > 0 and agree / found >= 0.98
    out = {"checked_at": time.strftime("%Y-%m-%d"), "rule": "PROTOCOL.md section 3, fixed details: 1000 Genomes file check",
           "structure": {"people": panel.n, "variants": panel.m, "bed_bytes_ok": True},
           "sample_seed": SEED, "sampled": len(sample), "found_in_ensembl": found, "within_0.01": agree,
           "allele_mismatch": allele_mismatch, "share_within_0.01": round(agree / found, 4) if found else None,
           "passed": passed and panel.n == 503, "variants": rows}
    SPOT.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in out.items() if k != "variants"}, indent=1))
    return out


# ---------- build ----------

def build():
    check = json.loads(SPOT.read_text(encoding="utf-8")) if SPOT.exists() else None
    if not check or not check.get("passed"):
        raise SystemExit("1000 Genomes reference has not passed the spot-check; run spotcheck first")
    rules = load_rules()
    gtex = route1_gtex(rules)  # network lookup first, so an Ensembl outage fails before the long file reads
    print("GTEx claims:", [(r["claim_id"], r["n_rebuilt"]) for r in gtex], flush=True)
    t0 = time.time()
    panel = RefPanel()
    print(f"reference loaded: {panel.m} variants, {panel.n} people ({time.time() - t0:.0f} s)", flush=True)

    # One pass per exposure file, collecting what every claim needs from it.
    needs = {}
    for cid in WEI:
        needs.setdefault(zhao_key(rules[cid][1]["gene"]), {"p": None, "rsids": set(), "chroms": None})
        needs[zhao_key(rules[cid][1]["gene"])]["p"] = rules[cid][0]["p_threshold"]
    for cid in GUO + GAO:
        n = needs.setdefault(zhao_key(rules[cid][1]["gene"]), {"p": None, "rsids": set(), "chroms": None})
        n["rsids"] |= set(rules[cid][1]["instrument_rsids"])
    needs["aholaolli2017_TRAIL"] = {"p": rules["Xu2025JCMM-1"][0]["p_threshold"], "rsids": set(), "chroms": None}
    for _, key in ARU_READINGS:
        needs[key] = {"p": rules["Aru2025cath-1"][0]["p_threshold"], "rsids": set(), "chroms": None}
    li_chroms = {gene_region_37(rules[c][1]["gene"], rules[c][1])[0] for c in LI}
    needs["glgc2021_logTG_EUR"] = {"p": rules[LI[0]][0]["p_threshold"], "rsids": set(), "chroms": li_chroms}

    cache = {}
    for key, n in needs.items():
        t0 = time.time()
        rows = scan(key, p_max=n["p"], rsids=frozenset(n["rsids"]), chroms=n["chroms"])
        cache[key] = attach_reference(rows, panel)
        print(f"{key}: {len(rows)} rows kept ({time.time() - t0:.0f} s)", flush=True)

    def scanner(key, rsids):
        return attach_reference(scan(key, rsids=rsids), panel)

    results = gtex + route1_zhao(rules, panel, cache, scanner) + route2_claims(rules, panel, cache)
    results.sort(key=lambda r: r["claim_id"])
    out = {"built_at": time.strftime("%Y-%m-%d"), "protocol": "PROTOCOL.md section 3, re-derivation (fixed details 2026-10-31)",
           "ld_reference": {"file": "1kg.v3.tgz > EUR.bed/bim/fam", "check": "results/ld_reference_check.json"},
           "claims": results}
    OUT.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    summary = [{k: r.get(k) for k in ("claim_id", "gene", "route", "n_reported", "n_rebuilt", "analysis_role")}
               | {"reading_used": r.get("reading_used"), "reasons": r["reasons"]} for r in results]
    counts = {}
    for r in results:
        counts[(r["route"], r["analysis_role"])] = counts.get((r["route"], r["analysis_role"]), 0) + 1
    REPORT.write_text(json.dumps({"built_at": out["built_at"],
                                  "counts": {f"{a} / {b}": n for (a, b), n in sorted(counts.items())},
                                  "claims": summary}, indent=1) + "\n", encoding="utf-8")
    for s in summary:
        print(s["claim_id"], s["route"], s["n_reported"], s["n_rebuilt"], s["analysis_role"], "; ".join(s["reasons"]))
    print(json.dumps(counts and {f"{a} / {b}": n for (a, b), n in sorted(counts.items())}, indent=1))


if __name__ == "__main__":
    commands = {"unpack": unpack, "spotcheck": spotcheck, "build": build}
    if len(sys.argv) != 2 or sys.argv[1] not in commands:
        raise SystemExit(__doc__)
    commands[sys.argv[1]]()
