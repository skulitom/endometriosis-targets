# Endometriosis drug-target check

> **Work in progress, not peer reviewed.** No replication analysis has been run yet, so this repository does not say which targets hold up. The extracted claims are AI-assisted and not yet verified by a person. Nothing here is medical advice or a recommendation to use any drug.

Many studies use human genetics (Mendelian randomization) to propose drug targets for endometriosis. This project collects every such published claim and checks, under a protocol fixed in advance:
- whether the claim holds up in independent patient data;
- whether the genetic signals for the gene and for endometriosis share a cause;
- whether the gene is active in the cell types involved in the disease.

Null and failed claims are reported alongside the ones that hold up.

Start with [GOAL.md](GOAL.md), [PROTOCOL.md](PROTOCOL.md) and [STATUS.md](STATUS.md).

Current state: Gates A and B (claim inventory and data checks) are nearly finished.
- **Claims:** 488 claims from 60 papers. 133 of them come from 15 papers added from reference lists and a supplementary full-text search.
- **Checks:** all 189 headline claims were double-extracted and adjudicated.
- **Instruments:** extracted and verified for 128 claims, and re-derived from public exposure data for 29 more.

No replication analysis has been run and no target is endorsed.

## Files

- `GOAL.md`: outcome-based goal, completion criteria, decision gates, working boundaries.
- `PROTOCOL.md`: analysis plan, draft v0.1. Becomes binding only when frozen (`results/freeze.json`).
- `EXTRACTION.md`: field-by-field rules for extracting claims from papers.
- `EVIDENCE.md`: sources and dataset decisions so far.
- `STATUS.md`: current checkpoint, open questions, next action.
- `data/datasets.json`: every planned data source with role, build, sample sizes, access and independence notes.
- `data/downloads.json`: allowlisted download URLs with source checksums.
- `data/literature/`: screening and extraction decisions, extracted claims (`claims.json`), logged corrections, instrument specs.
- `data/literature/public/`: the search log and every record screened, with its decision and reason. Identifiers and titles only; abstracts can be re-fetched from PubMed or Europe PMC.
- `src/public_export.py`: writes `data/literature/public/` from the working files.
- `src/fetch.py`: downloads one allowlisted file, verifies size and the source MD5 (or, where the source publishes none, size plus a full gzip integrity test, with SHA-256 pinned at first download); https only unless plain http is approved in the manifest and `--allow-http` is given; `--resume` continues an interrupted download; never overwrites.
- `src/literature.py`: runs the stored PubMed and Europe PMC queries and de-duplicates results.
- `src/claims.py`: claim schema, validation, independence check, summary and CSV export.
- `src/screening.py`: merges screening batches (database records and, separately, records from other sources with the lead's re-checks), builds the full-text queue and the human-review samples, keeping any answers already given.
- `src/supp_batches.py`: builds screening batches for the supplementary full-text search from titles, abstracts and full-text excerpts.
- `src/fulltext.py`: saves the open PMC full text of papers added from other sources as plain text (tables as tab-separated rows), with the raw XML and its checksum, so extractors read local files only.
- `src/refscan.py`: reference-list scan for missed studies; collects included papers' reference lists (Europe PMC, Crossref, PMC full text), matches them to the search records, lists unmatched endometriosis genetic-causal studies for screening, and logs additions with the paper they came from (`data/literature/refscan/`).
- `src/inventory.py`: merges full-text decisions and claims, applies logged corrections, writes the claim table and flow counts.
- `src/supplements.py`: downloads included papers' supplementary files through open routes, unpacks them safely, logs checksums in `data/supplements_manifest.json`.
- `src/tables.py`: converts supplementary spreadsheets and Word files to text with an index per paper.
- `src/second.py`: compares the blind second extraction with the first under fixed agreement rules.
- `src/reconcile.py`: turns agreements, adjudications and not-found results into `corrections_second.json`.
- `src/instruments.py`: extracts per-variant instruments from supplement tables using the mapping specs in `data/literature/instruments/specs/`, and recomputes each estimate as a check (`--check <spec>` for one spec).
- `src/verify_tables.py`: re-reads original .xlsx files with openpyxl and compares every extracted value (project venv).
- `src/genes.py`: looks up GRCh37 and GRCh38 positions of every claimed gene from Ensembl (`data/reference/claimed_genes.json`).
- `src/blinding.py`: the only way GWAS files are read before the freeze; drops rows within 1 Mb of any claimed gene and logs each read.
- `src/gateb.py`: positive controls at known endometriosis loci in each replication source.
- `src/reproduce.py`: recomputes a published estimate from FinnGen R10 and the paper's own instrument.
- `src/tier_dryrun.py`: applies the replication-source and tier rules to the claim table without outcome data (eligible sources, primary source, best possible tier).
- `src/rederive.py`: rebuilds instruments for claims that publish no per-variant effects, from the stated exposure files (known variants, or the paper's own selection rule with LD clumping against 1000 Genomes European); also checks the 1000 Genomes file against Ensembl (project venv, needs numpy).
- `.venv/`: project-only Python environment with xlrd 2.0.1 (old .xls), openpyxl 3.1.5 (cross-check) and numpy 2.3.4 (LD for re-derivation). Nothing is installed system-wide.
- `results/`: inventory flow counts, claim summary and a review CSV.
- `tests/`: synthetic software tests (no biological data); `test_rederive.py` needs the venv.

## Run

Python 3.13 standard library only so far.

```powershell
python -m unittest discover -s tests -v
python src/fetch.py --list
python src/literature.py search   # re-runs the dated searches
python src/literature.py merge
python src/screening.py           # merge screening batches
python -I src/refscan.py fetch    # reference-list scan: fetch, match, then report after screening
python -I src/refscan.py match
python -I src/refscan.py report
python -I src/literature.py supplementary   # one-off full-text search (already run 2026-10-31)
python -I src/literature.py merge           # records.json (databases) + records_supplementary.json (other sources)
python -I src/refscan.py passages           # full-text excerpts for screening
python -I src/supp_batches.py               # screening input batches
python src/inventory.py           # rebuild claims.json from raw extraction plus corrections.json
python src/claims.py validate data/literature/claims.json
.venv\Scripts\python -I -m unittest discover -s tests -p "test_rederive.py"
.venv\Scripts\python -I src\rederive.py unpack      # after the exposure downloads
.venv\Scripts\python -I src\rederive.py spotcheck
.venv\Scripts\python -I src\rederive.py build
```

Large downloads (0.8 to 1.6 GB each; the newest single-cell atlas is 39 GB) are not started automatically:

```powershell
python src/fetch.py --key rahmioglu2023_meta
```

## What is not in this repository, and how to get it

Other people's data and text are never redistributed here. Each file is listed with its source and checksum, so anyone can fetch the same bytes:

| Not included | Where it is listed | How to get it |
|---|---|---|
| GWAS summary statistics, exposure QTL files, 1000 Genomes reference (about 16 GB) | `data/downloads.json` | `python src/fetch.py --key <key>` |
| Supplementary files of included papers | `data/supplements_manifest.json` | `python src/supplements.py` (and `--additions`) |
| Full texts of added papers | `data/upstream/fulltext/*/provenance.json` (not included; regenerated) | `python -I src/fulltext.py` |
| Rahmioglu 2023 supplement (publisher site blocks scripts) | `data/manual/*.provenance.json` (MD5, SHA-256) | download by hand from the paper's page |
| Search results with abstracts, screening inputs with full-text excerpts | `data/literature/public/` | `python src/literature.py search`, `python -I src/refscan.py passages` |

Small per-variant extracts that the analysis needs are included, as in the supplement of any MR paper, with their source named next to them:
- instrument tables copied from papers' supplements (`data/literature/instruments/`);
- re-derived instruments (`data/literature/instruments/rederive/`);
- positive-control effects at 14 variants (`results/gateb_positive_controls.json`).

## Licence and reuse

- Code (`src/`, `tests/`): MIT licence, see `LICENSE`.
- Documents, decisions and the claim table (`*.md`, `data/`, `results/`): Creative Commons Attribution 4.0, see `LICENSE-docs.md`.
- Original authors own their data and findings; cite the original papers for any claim you use.

## How this was made

The project lead directed the work and makes every protocol decision. The decisions are recorded with dates in `PROTOCOL.md` and `STATUS.md`. Literature search, screening, claim extraction, adjudication and code were done with Claude (Anthropic), an AI model, under that direction:
- every AI decision is logged with its reason;
- raw outputs are kept unedited;
- every change goes through a logged correction file.

Human verification is recorded per claim (`human_verified`), and so far no claim has been verified by a person.
