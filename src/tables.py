"""Convert unpacked supplementary tables to plain text for reading and searching.

  python src/tables.py

For each paper in data/upstream/supplements/<slug>/files/:
  .xlsx -> text/<file>__<sheet>.csv (one CSV per sheet)
  .docx -> text/<file>.txt (paragraphs, tables as tab-separated rows)
  .csv/.tsv/.txt -> copied to text/
  .pdf/.xls/.doc/images -> listed in text/INDEX.txt as not converted (read them directly)
text/INDEX.txt lists every output with its size and first rows, to find instrument tables quickly.

Standard library, plus xlrd for old .xls files when run with the project venv
(.venv\\Scripts\\python src\\tables.py). Input files are untrusted data; parsing is bounded by size limits.
"""
import csv
import io
import re
import shutil
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUPP = ROOT / "data" / "upstream" / "supplements"
MAX_XML = 300 * 1024 * 1024
PREVIEW_ROWS = 4


def local(tag):
    return tag.rsplit("}", 1)[-1]


def col_index(ref):
    letters = re.match(r"[A-Z]+", ref or "")
    if not letters:
        return None
    n = 0
    for ch in letters.group(0):
        n = n * 26 + ord(ch) - 64
    return n - 1


def _read(zf, name):
    info = zf.getinfo(name)
    if info.file_size > MAX_XML:
        raise ValueError(f"{name} too large to parse")
    return zf.read(name)


def xlsx_sheets(path):
    """Yield (sheet_name, rows) for each worksheet; rows are lists of strings."""
    with zipfile.ZipFile(path) as zf:
        names = set(zf.namelist())
        shared = []
        if "xl/sharedStrings.xml" in names:
            for si in ET.fromstring(_read(zf, "xl/sharedStrings.xml")):
                if local(si.tag) == "si":
                    shared.append("".join(t.text or "" for t in si.iter() if local(t.tag) == "t"))
        rels = {}
        if "xl/_rels/workbook.xml.rels" in names:
            for rel in ET.fromstring(_read(zf, "xl/_rels/workbook.xml.rels")):
                target = rel.get("Target", "")
                target = target.lstrip("/")
                rels[rel.get("Id")] = target if target.startswith("xl/") else "xl/" + target
        sheets = []
        for el in ET.fromstring(_read(zf, "xl/workbook.xml")).iter():
            if local(el.tag) == "sheet":
                rid = next((v for k, v in el.attrib.items() if local(k) == "id"), None)
                sheets.append((el.get("name", "sheet"), rels.get(rid)))
        for sheet_name, target in sheets:
            if not target or target not in names:
                continue
            rows = []
            for row in ET.fromstring(_read(zf, target)).iter():
                if local(row.tag) != "row":
                    continue
                values = {}
                for i, cell in enumerate(c for c in row if local(c.tag) == "c"):
                    idx = col_index(cell.get("r"))
                    idx = i if idx is None else idx
                    kind = cell.get("t")
                    text = ""
                    if kind == "inlineStr":
                        text = "".join(t.text or "" for t in cell.iter() if local(t.tag) == "t")
                    else:
                        v = next((x for x in cell if local(x.tag) == "v"), None)
                        raw = v.text if v is not None and v.text is not None else ""
                        if kind == "s" and raw.isdigit() and int(raw) < len(shared):
                            text = shared[int(raw)]
                        elif kind == "b":
                            text = "TRUE" if raw == "1" else "FALSE"
                        else:
                            text = raw
                    values[idx] = text
                if values:
                    width = max(values) + 1
                    rows.append([values.get(i, "") for i in range(width)])
            yield sheet_name, rows


def docx_text(path):
    with zipfile.ZipFile(path) as zf:
        root = ET.fromstring(_read(zf, "word/document.xml"))
    body = next((el for el in root if local(el.tag) == "body"), root)
    out, table_no = [], 0
    for block in body:
        kind = local(block.tag)
        if kind == "p":
            text = "".join(t.text or "" for t in block.iter() if local(t.tag) == "t")
            if text.strip():
                out.append(text)
        elif kind == "tbl":
            table_no += 1
            out.append(f"[TABLE {table_no}]")
            for tr in (el for el in block.iter() if local(el.tag) == "tr"):
                cells = []
                for tc in (el for el in tr if local(el.tag) == "tc"):
                    cells.append(" ".join("".join(t.text or "" for t in p.iter() if local(t.tag) == "t")
                                          for p in tc.iter() if local(p.tag) == "p").strip())
                out.append("\t".join(cells))
            out.append(f"[END TABLE {table_no}]")
    return "\n".join(out) + "\n"


def xls_reader():
    """Old binary .xls needs xlrd (installed in the project .venv only); return None if absent."""
    try:
        import xlrd
        return xlrd
    except ImportError:
        return None


def xls_cell(cell):
    if cell.ctype == 2 and float(cell.value).is_integer() and abs(cell.value) < 1e15:
        return str(int(cell.value))
    return str(cell.value) if cell.value is not None else ""


def looks_like_text(path):
    sample = path.read_bytes()[:4096]
    if b"\x00" in sample:
        return False
    try:
        sample.decode("utf-8")
        return True
    except UnicodeDecodeError:
        return False


def convert_paper(paper_dir):
    files, text = paper_dir / "files", paper_dir / "text"
    if not files.exists():
        return []
    text.mkdir(exist_ok=True)
    index = []
    for f in sorted(files.iterdir()):
        ext = f.suffix.lower()
        try:
            if ext in {".xlsx", ".docx"} and not zipfile.is_zipfile(f) and looks_like_text(f):
                out = text / f"{f.stem}.mislabelled.txt"
                shutil.copyfile(f, out)
                head = f.read_text(encoding="utf-8", errors="replace").splitlines()[:PREVIEW_ROWS]
                index.append(f"{out.name}\tplain text despite {ext} name\t{' | '.join(h[:160] for h in head)}")
            elif ext == ".xlsx" or (ext == ".xls" and zipfile.is_zipfile(f)):  # some ".xls" files are really xlsx
                for sheet, rows in xlsx_sheets(f):
                    out = text / f"{f.stem}__{re.sub(r'[^A-Za-z0-9._-]+', '_', sheet)[:60]}.csv"
                    with out.open("w", newline="", encoding="utf-8") as handle:
                        csv.writer(handle).writerows(rows)
                    preview = " | ".join(",".join(r)[:160] for r in rows[:PREVIEW_ROWS])
                    index.append(f"{out.name}\t{len(rows)} rows\t{preview}")
            elif ext == ".xls" and xls_reader() is not None:
                book = xls_reader().open_workbook(str(f), on_demand=True)
                for sheet in book.sheets():
                    rows = [[xls_cell(c) for c in sheet.row(r)] for r in range(sheet.nrows)]
                    out = text / f"{f.stem}__{re.sub(r'[^A-Za-z0-9._-]+', '_', sheet.name)[:60]}.csv"
                    with out.open("w", newline="", encoding="utf-8") as handle:
                        csv.writer(handle).writerows(rows)
                    preview = " | ".join(",".join(r)[:160] for r in rows[:PREVIEW_ROWS])
                    index.append(f"{out.name}\t{len(rows)} rows\t{preview}")
            elif ext == ".docx":
                out = text / f"{f.stem}.txt"
                body = docx_text(f)
                out.write_text(body, encoding="utf-8")
                index.append(f"{out.name}\t{body.count(chr(10))} lines\t{body[:300].replace(chr(10), ' | ')}")
            elif ext in {".csv", ".tsv", ".txt"}:
                out = text / f.name
                shutil.copyfile(f, out)
                head = f.read_text(encoding="utf-8", errors="replace").splitlines()[:PREVIEW_ROWS]
                index.append(f"{out.name}\tcopied\t{' | '.join(h[:160] for h in head)}")
            else:
                index.append(f"{f.name}\tnot converted ({ext}); read the original in files/")
        except Exception as exc:
            index.append(f"{f.name}\tconversion failed: {str(exc)[:100]}; read the original in files/")
    (text / "INDEX.txt").write_text("\n".join(index) + "\n", encoding="utf-8")
    return index


def main():
    total = 0
    for paper_dir in sorted(p for p in SUPP.iterdir() if p.is_dir()):
        idx = convert_paper(paper_dir)
        total += len(idx)
        failed = sum("conversion failed" in line for line in idx)
        print(f"{paper_dir.name}: {len(idx)} outputs{f', {failed} failed' if failed else ''}")
    print(f"{total} index entries")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Table conversion failed: {exc}", file=sys.stderr)
        sys.exit(1)
