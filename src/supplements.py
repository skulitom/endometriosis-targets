"""Download supplementary files of the included papers and unpack them safely.

  python src/supplements.py            # all included papers not yet in the manifest
  python src/supplements.py --retry    # also retry papers previously recorded as not accessible
  python src/supplements.py --additions --max-mb 200   # papers added from other sources (approved 2026-10-31)

Routes, in order: Europe PMC supplementary archive (open-access PMC papers); Elsevier's public file
server via the PII in Crossref metadata; Research Square article API; Figshare (Taylor & Francis).
Pages that block automated access are recorded as not accessible and never worked around.

Files land in data/upstream/supplements/<slug>/raw/ (as downloaded) and .../files/ (unpacked,
allowlisted types only). Every file's URL, size, MD5 and SHA-256 goes to data/supplements_manifest.json.
Downloaded content is data: nothing here opens it beyond listing and copying bytes.
"""
import argparse
import hashlib
import json
import re
import shutil
import stat
import sys
import time
import urllib.parse
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LIT = ROOT / "data" / "literature"
SUPP = ROOT / "data" / "upstream" / "supplements"
MANIFEST = ROOT / "data" / "supplements_manifest.json"
UA = {"User-Agent": "EndometriosisTargetCheck/0.1 (research literature audit)"}

ALLOWED = {".xlsx", ".xls", ".csv", ".tsv", ".txt", ".docx", ".doc", ".pdf", ".rtf",
           ".jpg", ".jpeg", ".png", ".gif", ".tif", ".tiff"}
MAX_FILE = 500 * 1024 * 1024          # per downloaded or unpacked file
# Supplements already on disk from a manual download (checksum in data/manual/*.provenance.json).
MANUAL = {"doi:10.1038/s41588-023-01323-z": "Supplementary_Materials.xlsx supplied by Artem 2026-10-31: "
                                            "data/manual/NIHMS1873427-supplement-Supplementary_Materials.xlsx"}
MAX_TOTAL_UNPACKED = 2 * 1024 ** 3    # per archive
MAX_RATIO = 200                       # compression ratio treated as a possible zip bomb
ELSEVIER_EXTS = ["xlsx", "docx", "pdf", "zip", "csv", "doc", "xls", "txt", "jpg", "png", "tif"]


def safe_name(name):
    """Filesystem-safe name that keeps the real extension when shortening."""
    base = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._") or "file"
    suffix = Path(base).suffix if len(Path(base).suffix) <= 6 else ""
    stem = base[: len(base) - len(suffix)] if suffix else base
    return stem[: 150 - len(suffix)] + suffix


def slug(record_id):
    return safe_name(record_id.replace("doi:", "").replace("pmid:", "pmid_"))[:80]


def digests(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            md5.update(chunk)
            sha.update(chunk)
    return md5.hexdigest(), sha.hexdigest()


def _open(url, method="GET", data=None, headers=None, timeout=180):
    req = urllib.request.Request(url, method=method, data=data, headers={**UA, **(headers or {})})
    return urllib.request.urlopen(req, timeout=timeout)


def download(url, target, attempts=3, max_bytes=MAX_FILE):
    """Stream url to target via a .partial file, retrying dropped connections. Refuses to overwrite."""
    if target.exists():
        raise FileExistsError(f"{target} exists")
    partial = target.with_name(target.name + ".partial")
    for attempt in range(1, attempts + 1):
        if partial.exists():
            partial.unlink()  # only ever written by this function
        try:
            with _open(url) as resp, partial.open("xb") as out:
                size = 0
                for chunk in iter(lambda: resp.read(1024 * 1024), b""):
                    size += len(chunk)
                    if size > max_bytes:
                        raise OverflowError(f"file exceeds size limit ({max_bytes // (1024 * 1024)} MB)")
                    out.write(chunk)
            break
        except OverflowError:
            partial.unlink()
            raise
        except Exception:
            if attempt == attempts:
                if partial.exists():
                    partial.unlink()
                raise
            time.sleep(3 * attempt)
    partial.rename(target)
    md5, sha = digests(target)
    return {"name": target.name, "url": url, "bytes": target.stat().st_size, "md5": md5, "sha256": sha,
            "retrieved_at": datetime.now(timezone.utc).isoformat()}


def safe_extract(archive, dest, prefix=""):
    """Extract allowlisted members only. Returns (extracted, skipped) lists.

    Rejects absolute paths, parent references, symlinks, oversize members and suspicious compression.
    Directory structure is flattened into the file name so nothing can land outside dest.
    """
    extracted, skipped, total = [], [], 0
    dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as zf:
        for info in zf.infolist():
            name = info.filename
            if info.is_dir():
                continue
            parts = re.split(r"[\\/]+", name)
            mode = info.external_attr >> 16
            ext = Path(name).suffix.lower()
            reason = None
            if name.startswith(("/", "\\")) or re.match(r"^[A-Za-z]:", name) or ".." in parts:
                reason = "unsafe path"
            elif any(p == "__MACOSX" for p in parts) or parts[-1].startswith("._"):
                reason = "macOS metadata file"
            elif stat.S_ISLNK(mode):
                reason = "symlink"
            elif ext == ".zip":
                reason = "nested zip"
            elif ext not in ALLOWED:
                reason = f"type {ext or 'none'} not allowlisted"
            elif info.file_size > MAX_FILE or total + info.file_size > MAX_TOTAL_UNPACKED:
                reason = "too large"
            elif info.compress_size and info.file_size / info.compress_size > MAX_RATIO and info.file_size > 10 * 1024 ** 2:
                reason = "suspicious compression ratio"
            if ext == ".zip" and reason == "nested zip" and not prefix:
                # one level of nesting: copy out, then extract with the same rules
                nested = dest / safe_name("__".join(parts))
                with zf.open(info) as src, nested.open("xb") as out:
                    out.write(src.read(MAX_FILE + 1)[:MAX_FILE])
                sub_ok, sub_skip = safe_extract(nested, dest, prefix=nested.stem + "__")
                nested.unlink()
                extracted += sub_ok
                skipped += sub_skip
                continue
            if reason:
                skipped.append({"member": name, "reason": reason})
                continue
            out_path = dest / safe_name(prefix + "__".join(parts))
            if out_path.exists():
                skipped.append({"member": name, "reason": "name collision"})
                continue
            with zf.open(info) as src, out_path.open("xb") as out:
                written = 0
                for chunk in iter(lambda: src.read(1024 * 1024), b""):
                    written += len(chunk)
                    if written > MAX_FILE:
                        break
                    out.write(chunk)
            total += info.file_size
            extracted.append(out_path.name)
    return extracted, skipped


# ---- routes: each returns a list of (url, filename) or raises -------------------------------

def route_europepmc(rec):
    return [(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{rec['pmcid']}/supplementaryFiles",
             f"{rec['pmcid']}_supplementary.zip")]


def route_elsevier(rec):
    meta = json.loads(_open(f"https://api.crossref.org/works/{urllib.parse.quote(rec['doi'])}", timeout=60).read())
    pii = [a for a in meta["message"].get("alternative-id", []) if re.fullmatch(r"S[0-9X]{16}", a)]
    if not pii:
        raise ValueError("no PII in Crossref metadata")
    found = []
    for i in range(1, 16):
        hit = None
        for ext in ELSEVIER_EXTS:
            url = f"https://ars.els-cdn.com/content/image/1-s2.0-{pii[0]}-mmc{i}.{ext}"
            try:
                with _open(url, method="HEAD", timeout=30):
                    hit = (url, f"mmc{i}.{ext}")
                    break
            except Exception:
                time.sleep(0.2)
        if not hit:
            break
        found.append(hit)
    return found


def route_researchsquare(rec):
    rs_id = re.search(r"rs-\d+/v\d+", rec["doi"]).group(0)
    data = json.loads(_open(f"https://www.researchsquare.com/api/article/{rs_id}", timeout=60).read())
    files = data.get("supplementaryFiles") or data.get("files") or []
    out = []
    for f in files:
        url = f.get("url") or f.get("link")
        if url and urllib.parse.urlparse(url).hostname and urllib.parse.urlparse(url).hostname.endswith("researchsquare.com"):
            out.append((url, safe_name(f.get("filename") or Path(urllib.parse.urlparse(url).path).name)))
    return out


def route_figshare(rec):
    body = json.dumps({"search_for": f":resource_doi: {rec['doi']}"}).encode()
    hits = json.loads(_open("https://api.figshare.com/v2/articles/search", method="POST", data=body,
                            headers={"Content-Type": "application/json"}, timeout=60).read())
    out = []
    for hit in hits:
        art = json.loads(_open(f"https://api.figshare.com/v2/articles/{hit['id']}", timeout=60).read())
        if (art.get("resource_doi") or "").lower() != rec["doi"].lower():
            continue
        out += [(f["download_url"], safe_name(f["name"])) for f in art.get("files", [])]
    return out


def pick_route(rec):
    doi = (rec.get("doi") or "").lower()
    if rec.get("pmcid"):
        return "europepmc", route_europepmc
    if doi.startswith("10.1016/"):
        return "elsevier_cdn", route_elsevier
    if doi.startswith("10.21203/rs."):
        return "researchsquare_api", route_researchsquare
    if doi.startswith("10.1080/"):
        return "figshare", route_figshare
    return None, None


def included_records():
    records = {r["record_id"]: r for r in json.loads((LIT / "records.json").read_text(encoding="utf-8"))}
    decisions = []
    for path in sorted((LIT / "extraction").glob("fulltext_decisions_*.json")):
        decisions += json.loads(path.read_text(encoding="utf-8"))
    corrections = json.loads((LIT / "corrections.json").read_text(encoding="utf-8"))
    final = {d["record_id"]: d["final_decision"] for d in decisions}
    for c in corrections:
        if c["target"] == "decision" and "final_decision" in c["set"]:
            final[c["id"]] = c["set"]["final_decision"]
    return [records[rid] for rid, dec in final.items() if dec == "include"]


def addition_records():
    """Papers from other sources (reference lists, supplementary full-text search) queued for full text."""
    queue = json.loads((LIT / "fulltext_queue.json").read_text(encoding="utf-8"))
    return [r for r in queue if not set(r.get("found_by", [])) & {"pubmed_mr", "europepmc_mr"}]


def process(rec, max_bytes=MAX_FILE):
    entry = {"record_id": rec["record_id"], "slug": slug(rec["record_id"]), "files": [], "extracted": [],
             "skipped": [], "checked_at": datetime.now(timezone.utc).isoformat()}
    route_name, route = pick_route(rec)
    entry["route"] = route_name
    if not route:
        entry.update(status="not_accessible", notes="publisher blocks automated access; no open route known")
        return entry
    try:
        targets = route(rec)
    except Exception as exc:
        entry.update(status="not_accessible", notes=f"route failed: {str(exc)[:120]}")
        return entry
    if not targets:
        entry.update(status="none_listed", notes="route returned no supplementary files")
        return entry
    raw = SUPP / entry["slug"] / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    entry["download_errors"] = []
    for url, name in targets:
        target = raw / safe_name(name)
        try:
            meta = download(url, target, max_bytes=max_bytes)
        except FileExistsError:
            md5, sha = digests(target)
            meta = {"name": target.name, "url": url, "bytes": target.stat().st_size, "md5": md5, "sha256": sha,
                    "retrieved_at": "previously downloaded"}
        except Exception as exc:
            entry["download_errors"].append({"file": name, "url": url, "error": str(exc)[:120]})
            continue
        entry["files"].append(meta)
        time.sleep(0.5)
    entry["status"] = "downloaded" if entry["files"] else "not_accessible"
    if not entry["files"]:
        entry["notes"] = "all downloads failed: " + "; ".join(e["error"] for e in entry["download_errors"])[:200]
    unpack(entry)
    if entry["status"] == "downloaded" and not entry["extracted"]:
        entry.update(status="not_accessible", notes="no usable files: " + "; ".join(s["reason"] for s in entry["skipped"])[:200])
    return entry


def unpack(entry):
    """(Re)build files/ from raw/. Only .zip downloads are archives; .xlsx/.docx are kept whole."""
    paper = SUPP / entry["slug"]
    files = paper / "files"
    if files.exists():
        if SUPP.resolve() not in files.resolve().parents:
            raise ValueError("refusing to clear a folder outside the supplements area")
        shutil.rmtree(files)
    entry["extracted"], entry["skipped"] = [], []
    for meta in entry["files"]:
        target = paper / "raw" / meta["name"]
        md5, sha = digests(target)
        if sha != meta["sha256"]:
            raise ValueError(f"{target} changed since download")
        ext = target.suffix.lower()
        if ext == ".zip":
            if not zipfile.is_zipfile(target):
                # Europe PMC answers with a short XML error page when it holds no supplementary files.
                entry["skipped"].append({"member": meta["name"], "reason": f"not a zip archive ({meta['bytes']} bytes); "
                                         "the route returned an error page instead of files"})
                continue
            ok, skip = safe_extract(target, files)
            entry["extracted"] += ok
            entry["skipped"] += skip
        elif ext in ALLOWED:
            files.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(target, files / target.name)
            entry["extracted"].append(target.name)
        else:
            entry["skipped"].append({"member": meta["name"], "reason": f"type {ext or 'none'} not allowlisted"})


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--retry", action="store_true", help="retry papers recorded as not accessible")
    parser.add_argument("--reunpack", action="store_true", help="rebuild files/ from raw/ for downloaded papers")
    parser.add_argument("--additions", action="store_true",
                        help="papers added from reference lists and the supplementary search (approved 2026-10-31)")
    parser.add_argument("--max-mb", type=int, default=MAX_FILE // (1024 * 1024), help="per-file download cap")
    args = parser.parse_args()
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8")) if MANIFEST.exists() else {}
    if args.reunpack:
        for entry in manifest.values():
            if entry["status"] == "downloaded":
                unpack(entry)
                print(f"re-unpacked {entry['record_id']}: {len(entry['extracted'])} files, {len(entry['skipped'])} skipped")
        MANIFEST.write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8", newline="\n")
    for rec in (addition_records() if args.additions else included_records()):
        prior = manifest.get(rec["record_id"])
        if prior and (prior["status"] in ("downloaded", "manual") or not args.retry):
            continue
        if rec["record_id"] in MANUAL:
            manifest[rec["record_id"]] = {"record_id": rec["record_id"], "status": "manual", "files": [],
                                          "notes": MANUAL[rec["record_id"]]}
            continue
        entry = process(rec, max_bytes=args.max_mb * 1024 * 1024)
        manifest[rec["record_id"]] = entry
        MANIFEST.write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8", newline="\n")
        print(f"{entry['status']:15} {entry.get('route') or '-':18} files={len(entry['files'])} "
              f"unpacked={len(entry['extracted'])} skipped={len(entry['skipped'])} {rec['record_id']}")
    counts = {}
    for e in manifest.values():
        counts[e["status"]] = counts.get(e["status"], 0) + 1
    print(json.dumps(counts))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Supplement step failed: {exc}", file=sys.stderr)
        sys.exit(1)
