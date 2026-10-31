"""Blinded access to GWAS summary statistics before the protocol freeze.

Every row within 1 Mb of any claimed gene is dropped before the caller sees it, so pre-freeze work
(file checks, positive controls) cannot look at the regions the replication will test. Every read is
appended to results/blinding_log.jsonl with the file, purpose, rows read and rows masked.

Exemptions (data/reference/blinding_exemptions.json) unmask one gene in one file only, each with a
written reason, for example the reproduction check in a claim's own discovery file.
"""
import bisect
import csv
import gzip
import io
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GENES = ROOT / "data" / "reference" / "claimed_genes.json"
EXEMPTIONS = ROOT / "data" / "reference" / "blinding_exemptions.json"
LOG = ROOT / "results" / "blinding_log.jsonl"
WINDOW = 1_000_000
PAR1_END = 2_781_479  # pseudoautosomal region 1 ends here (GRCh38; GRCh37 ends earlier, so this is conservative)


def norm_chrom(chrom):
    c = str(chrom).strip().upper()
    c = c[3:] if c.startswith("CHR") else c
    return {"23": "X", "24": "Y", "25": "X", "XY": "X", "PAR": "X", "M": "MT", "26": "MT"}.get(c, c)


class Mask:
    def __init__(self, build, genes_path=None, exemptions_path=None, file_key=None, window=WINDOW):
        # Resolve paths at call time (not import time) so the module-level settings are always the ones used.
        genes_path = genes_path or GENES
        exemptions_path = exemptions_path or EXEMPTIONS
        data = json.loads(Path(genes_path).read_text(encoding="utf-8"))
        if data.get("unresolved"):
            raise RuntimeError(f"Refusing to read: genes without positions {data['unresolved']}")
        exempt = set()
        if file_key and Path(exemptions_path).exists():
            for e in json.loads(Path(exemptions_path).read_text(encoding="utf-8")):
                if e["file_key"] == file_key:
                    if not e.get("reason"):
                        raise RuntimeError("Exemption without a reason")
                    exempt.add(e["gene"])
        self.exempt = sorted(exempt)
        intervals = {}
        for sym, builds in data["genes"].items():
            if sym in exempt:
                continue
            g = builds[build]
            chroms = [norm_chrom(g["chrom"])]
            if chroms[0] in {"X", "Y"} and g["start"] <= PAR1_END:
                chroms = ["X", "Y"]
            for c in chroms:
                intervals.setdefault(c, []).append((max(0, g["start"] - window), g["end"] + window))
        self.merged = {}
        for c, ivs in intervals.items():
            ivs.sort()
            out = [list(ivs[0])]
            for s, e in ivs[1:]:
                if s <= out[-1][1]:
                    out[-1][1] = max(out[-1][1], e)
                else:
                    out.append([s, e])
            self.merged[c] = ([s for s, _ in out], [e for _, e in out])

    def masked(self, chrom, pos):
        c = norm_chrom(chrom)
        if c not in self.merged:
            return False
        starts, ends = self.merged[c]
        i = bisect.bisect_right(starts, pos) - 1
        return i >= 0 and pos <= ends[i]


def _log(entry):
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(entry) + "\n")


def read_blinded(path, build, chrom_col, pos_col, purpose, file_key, delimiter="\t", keep=None, row_filter=None):
    """Yield rows (dicts) outside every masked window. Logs the read when the generator is exhausted.

    keep: optional list of columns to return (saves memory).
    row_filter: optional predicate row_filter(parts, idx) on the split line of an unmasked row, where idx maps
    column name to position; it runs after masking, so it can never see a masked row.
    """
    mask = Mask(build, genes_path=GENES, exemptions_path=EXEMPTIONS, file_key=file_key)
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    read = masked = returned = 0
    started = datetime.now(timezone.utc).isoformat()
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        header_line = handle.readline()
        header = next(csv.reader(io.StringIO(header_line), delimiter=delimiter))
        header = [h.lstrip("#").strip() for h in header]
        idx = {h: i for i, h in enumerate(header)}
        ci, pi = idx[chrom_col], idx[pos_col]
        for line in handle:
            read += 1
            parts = line.rstrip("\r\n").split(delimiter)
            try:
                if mask.masked(parts[ci], int(float(parts[pi]))):
                    masked += 1
                    continue
            except (ValueError, IndexError):
                masked += 1  # unparseable position: treat as masked, never shown
                continue
            if row_filter and not row_filter(parts, idx):
                continue
            returned += 1
            row = dict(zip(header, parts))
            yield {k: row.get(k) for k in keep} if keep else row
    _log({"started": started, "finished": datetime.now(timezone.utc).isoformat(), "file": Path(path).name,
          "file_key": file_key, "build": build, "purpose": purpose, "rows_read": read, "rows_masked": masked,
          "rows_returned": returned, "exemptions": mask.exempt, "window_bp": WINDOW})
