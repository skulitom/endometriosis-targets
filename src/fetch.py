"""Download one allowlisted public source file and verify it.

Nothing is downloaded by default. Existing files are re-verified, never overwritten.

Checks:
- byte count must match when the manifest gives one;
- MD5 must match the source's published MD5;
- files whose source publishes no checksum (md5 null) need "checksum_policy": "none_published"
  and an expected byte count; their SHA-256 is pinned in the provenance file at first download
  and checked on every later verification;
- https only, unless the manifest records an approval for plain http ("plain_http_approved")
  and the caller passes allow_http (command line: --allow-http).
"""
import argparse
import gzip
import hashlib
import http.client
import json
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "downloads.json"


def checksums(path):
    hashes = {"md5": hashlib.md5(), "sha256": hashlib.sha256()}
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            for digest in hashes.values():
                digest.update(chunk)
    return {name: digest.hexdigest() for name, digest in hashes.items()}


def gzip_intact(path):
    """Decompress the whole gzip stream: the trailer CRC-32 and length catch truncation or corruption."""
    with gzip.open(path, "rb") as handle:
        while handle.read(1 << 24):
            pass


def verify(path, record, pinned_sha256=None):
    if "expected_bytes" in record and path.stat().st_size != record["expected_bytes"]:
        raise ValueError("Unexpected byte count; file must not enter analysis")
    actual = checksums(path)
    if record.get("md5") is None:
        if record.get("checksum_policy") != "none_published" or "expected_bytes" not in record:
            raise ValueError("No source MD5: manifest must say checksum_policy none_published and give expected_bytes")
        if pinned_sha256 is not None and actual["sha256"] != pinned_sha256:
            raise ValueError("SHA-256 differs from the one pinned at first download; file must not enter analysis")
        if pinned_sha256 is None and str(record.get("filename", path.name)).endswith((".gz", ".tgz")):
            try:
                gzip_intact(path)
            except (OSError, EOFError) as exc:
                raise ValueError(f"gzip integrity check failed ({exc}); file must not enter analysis") from None
            actual["gzip_crc"] = "whole stream decompressed; CRC-32 and length matched"
    elif actual["md5"] != record["md5"]:
        raise ValueError("Source MD5 mismatch; file must not enter analysis")
    return actual


def download(record, partial, resume, attempts=8):
    """Write the URL to partial; with resume, continue an existing partial file using HTTP ranges."""
    for attempt in range(attempts):
        have = partial.stat().st_size if partial.exists() else 0
        headers = {"User-Agent": "EndometriosisTargetCheck/0.1"}
        if have:
            headers["Range"] = f"bytes={have}-"
        request = urllib.request.Request(record["url"], headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                if have:
                    expected = f"bytes {have}-"
                    if response.status != 206 or not response.headers.get("Content-Range", "").startswith(expected):
                        raise ValueError("Server did not honour the range request; cannot resume safely")
                with partial.open("ab" if have else "xb") as handle:
                    for chunk in iter(lambda: response.read(1024 * 1024), b""):
                        handle.write(chunk)
        except (OSError, http.client.HTTPException) as exc:
            if not resume or attempt == attempts - 1:
                raise
            print(f"attempt {attempt + 1} interrupted at {partial.stat().st_size if partial.exists() else 0} bytes: {exc}",
                  file=sys.stderr)
            time.sleep(min(60, 5 * (attempt + 1)))
            continue
        size = partial.stat().st_size
        if not resume or "expected_bytes" not in record or size >= record["expected_bytes"]:
            return attempt + 1
        print(f"attempt {attempt + 1} ended early at {size} bytes; resuming", file=sys.stderr)
    raise ValueError("Download still incomplete after retries")


def check_url(record, allow_http):
    scheme = urllib.parse.urlparse(record["url"]).scheme
    if scheme == "https":
        return
    if scheme == "http" and record.get("plain_http_approved") and allow_http:
        return
    raise ValueError(f"Refusing {scheme}:// download (https only unless approved in the manifest and --allow-http given)")


def fetch(key, manifest, target_dir, allow_http=False, resume=False):
    record = manifest[key]
    if Path(record["filename"]).name != record["filename"]:
        raise ValueError("Invalid manifest basename")
    subdir = record.get("subdir", "")
    if subdir and (Path(subdir).is_absolute() or ".." in Path(subdir).parts):
        raise ValueError("Invalid manifest subdir")
    target_dir = target_dir / subdir if subdir else target_dir
    target = target_dir / record["filename"]
    partial = target.with_name(target.name + ".partial")
    provenance_path = target.with_name(target.name + ".provenance.json")
    if target.exists():
        pinned = None
        if provenance_path.exists():
            pinned = json.loads(provenance_path.read_text(encoding="utf-8")).get("sha256")
        return {"status": "existing_verified", "path": str(target), **verify(target, record, pinned)}
    if partial.exists() and not resume:
        raise ValueError(f"Partial file exists; inspect or delete it, or pass --resume: {partial}")
    check_url(record, allow_http)
    target_dir.mkdir(parents=True, exist_ok=True)
    resumed_from = partial.stat().st_size if partial.exists() else 0
    attempts = download(record, partial, resume)
    digests = verify(partial, record)
    partial.rename(target)
    provenance = {"key": key, "url": record["url"], "source": record["source"],
                  "retrieved_at": datetime.now(timezone.utc).isoformat(),
                  "bytes": target.stat().st_size, **digests,
                  "integrity": ("source MD5 matched" if record.get("md5")
                                else "no source checksum: byte count matched; SHA-256 pinned here")}
    if urllib.parse.urlparse(record["url"]).scheme == "http":
        provenance["transport"] = f"plain http ({record['plain_http_approved']})"
    if resumed_from or attempts > 1:
        provenance["resumed"] = {"started_from_bytes": resumed_from, "connections": attempts}
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8", newline="\n")
    return {"status": "downloaded_verified", "path": str(target), **provenance}


def main():
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--key", choices=sorted(manifest))
    group.add_argument("--list", action="store_true", help="show allowlisted files and sizes")
    parser.add_argument("--allow-http", action="store_true",
                        help="permit plain http for manifest entries that record an approval")
    parser.add_argument("--resume", action="store_true",
                        help="continue an interrupted download with HTTP range requests (checks still apply)")
    args = parser.parse_args()
    if args.list:
        for key, record in sorted(manifest.items()):
            size = record.get("expected_bytes")
            print(f"{key:32} {size / 1e6:9.1f} MB  {record['source']}" if size else f"{key:32}  size unknown")
        return
    print(json.dumps(fetch(args.key, manifest, ROOT / "data" / "upstream",
                           allow_http=args.allow_http, resume=args.resume)))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Download not completed: {exc}", file=sys.stderr)
        sys.exit(1)
