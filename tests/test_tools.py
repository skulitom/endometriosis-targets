"""Synthetic SOFTWARE tests. These records contain no biological observations."""
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claims import independent, validate_claim, validate_table  # noqa: E402
from fetch import fetch, verify  # noqa: E402
from literature import merge  # noqa: E402
from screening import exclusion_sample, merge_decisions  # noqa: E402


from inventory import apply_corrections, build  # noqa: E402
from supplements import safe_extract, safe_name  # noqa: E402
import zipfile  # noqa: E402


class SupplementTests(unittest.TestCase):
    def test_safe_extract_blocks_traversal_symlink_and_scripts(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            archive = d / "a.zip"
            with zipfile.ZipFile(archive, "w") as zf:
                zf.writestr("../evil.txt", "x")
                zf.writestr("C:/abs.txt", "x")
                zf.writestr("run.bat", "echo hi")
                link = zipfile.ZipInfo("link.txt")
                link.external_attr = (0o120777 << 16)
                zf.writestr(link, "target")
                zf.writestr("tables/S1 Table.xlsx", "data")
            ok, skipped = safe_extract(archive, d / "out")
            self.assertEqual(ok, ["tables__S1_Table.xlsx"])
            reasons = sorted(s["reason"] for s in skipped)
            self.assertIn("unsafe path", reasons)
            self.assertIn("symlink", reasons)
            self.assertTrue(any("not allowlisted" in r for r in reasons))
            self.assertFalse((d / "evil.txt").exists())

    def test_nested_zip_extracted_one_level_with_same_rules(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            inner = d / "inner.zip"
            with zipfile.ZipFile(inner, "w") as zf:
                zf.writestr("S2.csv", "a,b")
                zf.writestr("../up.csv", "a,b")
            outer = d / "outer.zip"
            with zipfile.ZipFile(outer, "w") as zf:
                zf.write(inner, "inner.zip")
            ok, skipped = safe_extract(outer, d / "out")
            self.assertEqual(ok, ["inner__S2.csv"])
            self.assertEqual([s["reason"] for s in skipped], ["unsafe path"])

    def test_xlsx_and_docx_conversion(self):
        from tables import docx_text, xlsx_sheets
        ns = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
        rns = 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
        with tempfile.TemporaryDirectory() as directory:
            x = Path(directory) / "t.xlsx"
            with zipfile.ZipFile(x, "w") as zf:
                zf.writestr("xl/workbook.xml", f'<workbook {ns} {rns}><sheets><sheet name="S1" sheetId="1" r:id="rId1"/></sheets></workbook>')
                zf.writestr("xl/_rels/workbook.xml.rels", '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                            '<Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
                zf.writestr("xl/sharedStrings.xml", f'<sst {ns}><si><t>SNP</t></si><si><t>rs123</t></si></sst>')
                zf.writestr("xl/worksheets/sheet1.xml", f'<worksheet {ns}><sheetData>'
                            '<row r="1"><c r="A1" t="s"><v>0</v></c><c r="C1" t="inlineStr"><is><t>beta</t></is></c></row>'
                            '<row r="2"><c r="A2" t="s"><v>1</v></c><c r="C2"><v>0.12</v></c></row></sheetData></worksheet>')
            sheets = list(xlsx_sheets(x))
            self.assertEqual(sheets, [("S1", [["SNP", "", "beta"], ["rs123", "", "0.12"]])])
            d = Path(directory) / "t.docx"
            w = 'xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"'
            with zipfile.ZipFile(d, "w") as zf:
                zf.writestr("word/document.xml", f'<w:document {w}><w:body><w:p><w:r><w:t>Table S1</w:t></w:r></w:p>'
                            '<w:tbl><w:tr><w:tc><w:p><w:r><w:t>Gene</w:t></w:r></w:p></w:tc><w:tc><w:p><w:r><w:t>OR</w:t></w:r></w:p></w:tc></w:tr></w:tbl>'
                            '</w:body></w:document>')
            self.assertEqual(docx_text(d), "Table S1\n[TABLE 1]\nGene\tOR\n[END TABLE 1]\n")

    def test_safe_name_strips_separators(self):
        self.assertNotIn("/", safe_name("../a/b c.xlsx"))
        self.assertTrue(safe_name("../a/b c.xlsx").endswith(".xlsx"))

    def test_safe_name_keeps_extension_when_shortening(self):
        long_name = "x" * 300 + ".xlsx"
        self.assertTrue(safe_name(long_name).endswith(".xlsx"))
        self.assertLessEqual(len(safe_name(long_name)), 150)

    def test_macos_metadata_skipped(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            archive = d / "m.zip"
            with zipfile.ZipFile(archive, "w") as zf:
                zf.writestr("__MACOSX/._S1.xlsx", "junk")
                zf.writestr("._S2.csv", "junk")
                zf.writestr("S3.csv", "a,b")
            ok, skipped = safe_extract(archive, d / "out")
            self.assertEqual(ok, ["S3.csv"])
            self.assertEqual({s["reason"] for s in skipped}, {"macOS metadata file"})


class SecondExtractionTests(unittest.TestCase):
    def test_rounding_is_agreement(self):
        from second import compare
        diffs, _ = compare(make_claim(estimate=1.357, ci_low=1.10, ci_high=1.30 * 1.3, p_value=5.57e-9),
                           {"estimate": 1.36, "ci_low": 1.1, "ci_high": 1.69, "p_value": 5.6e-9})
        self.assertEqual(diffs, {})

    def test_real_differences_flagged(self):
        from second import compare
        diffs, _ = compare(make_claim(), {"estimate": 0.8, "ci_low": 0.7, "ci_high": 0.9, "p_value": 0.001,
                                          "direction": "protective", "outcome_cohorts": ["UK Biobank"]})
        self.assertEqual(set(diffs), {"estimate", "ci_low", "ci_high", "direction", "outcome_cohorts"})

    def test_missing_versus_present_is_a_difference(self):
        from second import same_p
        self.assertFalse(same_p(None, 0.01))
        self.assertTrue(same_p(None, None))
        self.assertTrue(same_p(1.0e-5, 1.1e-5))
        self.assertFalse(same_p(1e-5, 1e-4))


class BlindingTests(unittest.TestCase):
    def write_genes(self, d, genes, unresolved=()):
        p = Path(d) / "genes.json"
        p.write_text(json.dumps({"genes": genes, "unresolved": list(unresolved)}), encoding="utf-8")
        return p

    def test_mask_window_par_and_chrom_names(self):
        from blinding import Mask
        with tempfile.TemporaryDirectory() as d:
            genes = {"G1": {"GRCh37": {"chrom": "6", "start": 5_000_000, "end": 5_100_000}},
                     "PARG": {"GRCh37": {"chrom": "X", "start": 1_400_000, "end": 1_500_000}}}
            m = Mask("GRCh37", genes_path=self.write_genes(d, genes), exemptions_path=Path(d) / "none.json")
            self.assertTrue(m.masked("chr6", 4_000_000))
            self.assertTrue(m.masked("6", 6_100_000))
            self.assertFalse(m.masked("6", 6_100_001))
            self.assertTrue(m.masked("23", 1_450_000))
            self.assertTrue(m.masked("Y", 1_450_000))   # pseudoautosomal gene masked on both X and Y
            self.assertFalse(m.masked("7", 5_000_000))

    def test_reader_never_returns_masked_rows_and_logs(self):
        import gzip
        import blinding
        with tempfile.TemporaryDirectory() as d:
            genes = {"G1": {"GRCh37": {"chrom": "1", "start": 5_000_000, "end": 5_001_000}}}
            gp = self.write_genes(d, genes)
            f = Path(d) / "gwas.tsv.gz"
            with gzip.open(f, "wt") as h:
                h.write("chromosome\tbase_pair_location\tp_value\n1\t5000500\t1e-20\n1\t9000000\t1e-9\n2\t5000500\t0.5\n")
            old = (blinding.GENES, blinding.EXEMPTIONS, blinding.LOG)
            blinding.GENES, blinding.EXEMPTIONS, blinding.LOG = gp, Path(d) / "none.json", Path(d) / "log.jsonl"
            try:
                rows = list(blinding.read_blinded(f, "GRCh37", "chromosome", "base_pair_location", "test", "k",
                                                  row_filter=lambda parts, idx: float(parts[idx["p_value"]]) < 1e-8))
                log = [json.loads(x) for x in (Path(d) / "log.jsonl").read_text().splitlines()]
            finally:
                blinding.GENES, blinding.EXEMPTIONS, blinding.LOG = old
            self.assertEqual([r["base_pair_location"] for r in rows], ["9000000"])
            self.assertEqual((log[0]["rows_read"], log[0]["rows_masked"], log[0]["rows_returned"]), (3, 1, 1))

    def test_refuses_unresolved_genes(self):
        from blinding import Mask
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(RuntimeError, "Refusing"):
                Mask("GRCh37", genes_path=self.write_genes(d, {}, ["X1 (GRCh37)"]))

    def test_exemption_applies_to_named_file_only(self):
        from blinding import Mask
        with tempfile.TemporaryDirectory() as d:
            genes = {"G1": {"GRCh38": {"chrom": "1", "start": 10_000_000, "end": 10_010_000}}}
            ex = Path(d) / "ex.json"
            ex.write_text(json.dumps([{"file_key": "disc", "gene": "G1", "reason": "reproduction"}]), encoding="utf-8")
            gp = self.write_genes(d, genes)
            self.assertFalse(Mask("GRCh38", genes_path=gp, exemptions_path=ex, file_key="disc").masked("1", 10_000_500))
            self.assertTrue(Mask("GRCh38", genes_path=gp, exemptions_path=ex, file_key="other").masked("1", 10_000_500))


class GateBTests(unittest.TestCase):
    def test_alignment_and_binomial(self):
        from gateb import align, binom_upper
        lead = {"ea": "A", "oa": "G", "eaf": 0.3}
        self.assertEqual(align(lead, "A", "G", 0.1, None)[0], 0.1)
        self.assertEqual(align(lead, "G", "A", 0.1, None)[0], -0.1)
        self.assertEqual(align(lead, "T", "C", 0.1, None)[0], 0.1)       # strand flip
        self.assertIsNone(align(lead, "A", "C", 0.1, None)[0])
        self.assertIsNone(align({"ea": "A", "oa": "T", "eaf": 0.5}, "A", "T", 0.1, None)[0])
        self.assertAlmostEqual(binom_upper(10, 10), 1 / 1024)
        self.assertAlmostEqual(binom_upper(0, 5), 1.0)


class InstrumentTests(unittest.TestCase):
    def test_number_parsing(self):
        from instruments import parse_number
        self.assertAlmostEqual(parse_number("1.2E-05"), 1.2e-5)
        self.assertAlmostEqual(parse_number("−0.5"), -0.5)
        self.assertAlmostEqual(parse_number("1.2 × 10^-5"), 1.2e-5)
        self.assertIsNone(parse_number("NA"))
        with self.assertRaises(ValueError):
            parse_number("abc")

    def test_wald_ivw_and_allele_flip(self):
        import math as m
        from instruments import recompute
        one = [{"beta_exposure": 0.5, "beta_outcome": 0.1, "se_outcome": 0.02}]
        self.assertAlmostEqual(recompute(one)[0], 0.2)
        flipped = [{"beta_exposure": 0.5, "beta_outcome": -0.1, "se_outcome": 0.02,
                    "effect_allele": "A", "other_allele": "G", "effect_allele_outcome": "G", "other_allele_outcome": "A"}]
        self.assertAlmostEqual(recompute(flipped)[0], 0.2)
        two = [{"beta_exposure": 0.5, "beta_outcome": 0.1, "se_outcome": 0.02},
               {"beta_exposure": 0.25, "beta_outcome": 0.05, "se_outcome": 0.02}]
        self.assertAlmostEqual(recompute(two)[0], 0.2)
        self.assertAlmostEqual(recompute(two)[1], 1 / m.sqrt(0.5 ** 2 / 0.02 ** 2 + 0.25 ** 2 / 0.02 ** 2))
        unalignable = [{"beta_exposure": 0.5, "beta_outcome": 0.1, "se_outcome": 0.02,
                        "effect_allele": "A", "other_allele": "G", "effect_allele_outcome": "C", "other_allele_outcome": "T"}]
        self.assertIsNone(recompute(unalignable))

    def test_apply_spec_extracts_checks_and_recomputes(self):
        import math as m
        from instruments import apply_spec
        with tempfile.TemporaryDirectory() as directory:
            supp = Path(directory) / "paper"
            (supp / "text").mkdir(parents=True)
            (supp / "text" / "S1.csv").write_text(
                "Table S1 title\nSNP,exposure,ea,oa,bx,sx,by,sy\n"
                "rs1,GENE1,A,G,0.5,0.05,0.1,0.02\nrs2,GENE1,C,T,0.25,0.05,0.05,0.02\n"
                "rs3,OTHER,A,G,0.3,0.05,0.0,0.02\n", encoding="utf-8")
            spec = {"record_id": "pmid:1", "tables": {"T1": {"file": "text/S1.csv", "header_row": 2, "columns": {
                "rsid": "SNP", "effect_allele": "ea", "other_allele": "oa", "beta_exposure": "bx",
                "se_exposure": "sx", "beta_outcome": "by", "se_outcome": "sy"}}},
                "claims": {"T1": {"status": "mapped", "table": "T1",
                                  "filters": [{"column": "exposure", "op": "equals", "value": "GENE1"}], "expected_n": 2},
                           "T2": {"status": "not_available", "reason": "no table"}}}
            claims = {"T1": make_claim(estimate=round(m.exp(0.2), 4), ci_low=1.0, ci_high=1.5, n_instruments=3),
                      "T2": make_claim(claim_id="T2")}
            out = apply_spec(spec, supp, claims)
            self.assertEqual(out["T1"]["status"], "extracted")
            self.assertEqual([v["rsid"] for v in out["T1"]["variants"]], ["rs1", "rs2"])
            self.assertEqual(out["T1"]["variants"][0]["source"], {"file": "text/S1.csv", "line": 3})
            self.assertTrue(out["T1"]["recompute"]["agree"])
            self.assertTrue(any("reports 3" in f for f in out["T1"]["flags"]))
            self.assertEqual(out["T2"]["status"], "not_available")

    def test_spec_cannot_escape_supplement_folder(self):
        from instruments import apply_spec
        with tempfile.TemporaryDirectory() as directory:
            supp = Path(directory) / "paper"
            supp.mkdir()
            (Path(directory) / "outside.csv").write_text("SNP\nrs1\n", encoding="utf-8")
            spec = {"record_id": "x", "tables": {"T1": {"file": "../outside.csv", "columns": {
                "rsid": "SNP", "effect_allele": "SNP", "beta_exposure": "SNP"}}}, "claims": {}}
            with self.assertRaisesRegex(ValueError, "escapes"):
                apply_spec(spec, supp, {})


class ReconcileTests(unittest.TestCase):
    def test_agree_disagree_and_not_found(self):
        from reconcile import build as reconcile
        claims = {"A": make_claim(claim_id="A"), "B": make_claim(claim_id="B", estimate=None, ci_low=None, ci_high=None),
                  "C": make_claim(claim_id="C")}
        comparison = {"claims": [
            {"claim_id": "A", "status": "agree"},
            {"claim_id": "B", "status": "disagree", "differences": {"estimate": [None, 1.3]}},
            {"claim_id": "C", "status": "not_found", "notes": "paywalled"}]}
        adjudicated = {"B": {"claim_id": "B", "fields": [{"field": "estimate", "decision": "second", "value": 1.3,
                                                          "source_inconsistent": True, "note": "abstract says 1.2"}],
                             "overall_note": None}}
        out = {c["id"]: c for c in reconcile(claims, comparison, {"A": {"instruments_location": "S2"}}, adjudicated)}
        self.assertEqual(out["A"]["set"]["instruments_location"], "S2")
        self.assertEqual(out["B"]["expect"], {"estimate": None})
        self.assertEqual(out["B"]["set"]["estimate"], 1.3)
        self.assertIn("abstract says 1.2", out["B"]["set"]["source_issues"])
        self.assertTrue(out["C"]["set"]["second_extraction"].startswith("not found"))

    def test_missing_adjudication_fails(self):
        from reconcile import build as reconcile
        with self.assertRaisesRegex(ValueError, "No adjudication"):
            reconcile({"B": make_claim(claim_id="B")}, {"claims": [{"claim_id": "B", "status": "disagree",
                                                                   "differences": {"estimate": [1.2, 1.3]}}]}, {}, {})


class CorrectionTests(unittest.TestCase):
    def test_correction_applies_when_expectation_matches(self):
        _, claims = apply_corrections([], [make_claim()], [
            {"target": "claim", "id": "T1", "expect": {"direction": "risk_increasing"},
             "set": {"direction": "null_or_unclear"}, "why": "test"}])
        self.assertEqual(claims[0]["direction"], "null_or_unclear")

    def test_stale_correction_fails(self):
        with self.assertRaisesRegex(ValueError, "expected"):
            apply_corrections([], [make_claim()], [
                {"target": "claim", "id": "T1", "expect": {"direction": "protective"},
                 "set": {"direction": "null_or_unclear"}, "why": "test"}])

    def test_correction_needs_reason(self):
        with self.assertRaisesRegex(ValueError, "reason"):
            apply_corrections([], [make_claim()], [{"target": "claim", "id": "T1", "set": {}, "why": ""}])


class InventoryTests(unittest.TestCase):
    screened = [{"record_id": "pmid:1", "screening": {"decision": "include", "reason_code": "I"}},
                {"record_id": "pmid:2", "screening": {"decision": "unsure", "reason_code": "U"}},
                {"record_id": "pmid:3", "screening": {"decision": "exclude", "reason_code": "E2"}}]

    def decisions(self, n1=1):
        return [{"record_id": "pmid:1", "final_decision": "include", "reason_code": "I", "full_text_accessed": True, "n_claims": n1},
                {"record_id": "pmid:2", "final_decision": "exclude", "reason_code": "E6", "full_text_accessed": True,
                 "n_claims": 0, "duplicate_of": "pmid:1"}]

    def test_flow_counts(self):
        flow, summ = build(self.screened, self.decisions(), [make_claim(paper_id="pmid:1")])
        self.assertEqual(flow["papers_included"], 1)
        self.assertEqual(flow["full_text_excluded"], {"E6": 1})
        self.assertEqual(summ["claims"], 1)

    def test_claim_count_must_match_decision(self):
        with self.assertRaisesRegex(ValueError, "n_claims"):
            build(self.screened, self.decisions(n1=2), [make_claim(paper_id="pmid:1")])

    def test_claims_from_excluded_paper_rejected(self):
        with self.assertRaisesRegex(ValueError, "excluded"):
            build(self.screened, self.decisions(), [make_claim(paper_id="pmid:2", pmid="2")])

    def test_missing_decision_rejected(self):
        with self.assertRaisesRegex(ValueError, "no full-text decision"):
            build(self.screened, self.decisions()[:1], [make_claim(paper_id="pmid:1")])


class ScreeningTests(unittest.TestCase):
    records = [{"record_id": "a"}, {"record_id": "b"}]

    def test_every_record_needs_one_decision(self):
        with self.assertRaisesRegex(ValueError, "missing"):
            merge_decisions(self.records, [{"record_id": "a", "decision": "include", "reason_code": "I"}])

    def test_code_must_match_decision(self):
        rows = [{"record_id": "a", "decision": "exclude", "reason_code": "I"},
                {"record_id": "b", "decision": "include", "reason_code": "I"}]
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            merge_decisions(self.records, rows)

    def test_exclusion_sample_has_minimum_and_is_reproducible(self):
        merged = [{"record_id": str(i), "screening": {"decision": "exclude"}} for i in range(30)]
        first = [m["record_id"] for m in exclusion_sample(merged)]
        self.assertEqual(len(first), 10)
        self.assertEqual(first, [m["record_id"] for m in exclusion_sample(merged)])


def make_claim(**overrides):
    claim = {
        "claim_id": "T1", "paper_id": "pmid:1", "doi": None, "pmid": "1", "first_author": "Example",
        "year": 2025, "journal": None, "preprint": False, "gene_symbol": "GENE1", "exposure_type": "pQTL",
        "exposure_source": "synthetic", "exposure_tissue": None, "outcome_source": "synthetic",
        "outcome_release": None, "outcome_cases": 10, "outcome_controls": 20, "outcome_cohorts": ["FinnGen"],
        "instrument_scope": "cis", "n_instruments": 1, "instruments": None, "instruments_location": None, "method": "Wald ratio",
        "estimate_scale": "OR", "estimate": 1.2, "ci_low": 1.1, "ci_high": 1.3, "p_value": 0.001,
        "direction": "risk_increasing", "coloc_reported": False, "coloc_pph4": None,
        "replication_reported": False, "replication_details": None, "experimental_validation": None,
        "headline": True, "source_location": "Table 2", "full_text_accessed": True,
        "extraction_by": "test", "second_extraction": None, "human_verified": False, "notes": None,
    }
    claim.update(overrides)
    return claim


class ClaimTests(unittest.TestCase):
    def test_valid_claim_passes(self):
        self.assertEqual(validate_claim(make_claim()), [])

    def test_common_direction_flips_per_lower_and_drops_unclear(self):
        from claims import common_direction
        self.assertEqual(common_direction(make_claim()), 1)
        self.assertEqual(common_direction(make_claim(effect_scale="per_lower")), -1)
        self.assertEqual(common_direction(make_claim(estimate=0.8, ci_low=0.7, ci_high=0.9, direction="protective",
                                                     effect_scale="per_lower")), 1)
        self.assertIsNone(common_direction(make_claim(effect_scale="unclear")))
        self.assertIsNone(common_direction(make_claim(direction="null_or_unclear")))
        self.assertTrue(any("effect_scale" in p for p in validate_claim(make_claim(effect_scale="per_SD"))))

    def test_discovery_conflicts_group_smr_with_qtl_type(self):
        from tier_dryrun import discovery_direction_conflicts
        claims = [make_claim(claim_id="A", exposure_type="pQTL"),
                  make_claim(claim_id="B", exposure_type="SMR_pQTL", estimate=0.8, ci_low=0.7, ci_high=0.9,
                             direction="protective"),
                  make_claim(claim_id="C", exposure_type="eQTL", estimate=0.8, ci_low=0.7, ci_high=0.9,
                             direction="protective", effect_scale="per_lower")]
        out = discovery_direction_conflicts(claims)
        self.assertEqual([(x["exposure_type"], x["higher_level_raises_risk"], x["higher_level_lowers_risk"])
                          for x in out], [("pQTL", ["A"], ["B"])])

    def test_rederived_instruments_merge_and_validate(self):
        from inventory import merge_rederived
        rebuilt = {"claims": [{"claim_id": "T1", "route": "stated_rule", "dataset": "d", "dataset_file": "f.gz",
                               "dataset_sha256": "0" * 64, "n_reported": 1, "n_rebuilt": 1,
                               "analysis_role": "primary_rederived", "reasons": [], "reading_used": None,
                               "variants": [{"rsid": "rs1", "chrom": "1", "pos": 5, "ea": "A", "oa": "G",
                                             "eaf": 0.2, "beta": 0.3, "se": 0.05, "p": 1e-9}]}]}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "rederived.json"
            path.write_text(json.dumps(rebuilt), encoding="utf-8")
            claims = merge_rederived([make_claim(instruments_status="not_available")], path)
            self.assertEqual(claims[0]["instruments_status"], "rederived")
            self.assertEqual(claims[0]["instruments"][0]["beta_exposure"], 0.3)
            self.assertEqual(validate_claim(claims[0]), [])
            claims[0]["instruments_check"] = {}
            self.assertTrue(any("re-derivation record" in p for p in validate_claim(claims[0])))
            with self.assertRaisesRegex(ValueError, "already has extracted"):
                merge_rederived([make_claim(instruments_status="extracted")], path)

    def test_direction_must_match_estimate(self):
        problems = validate_claim(make_claim(estimate=0.8, ci_low=0.7, ci_high=0.9))
        self.assertTrue(any("risk_increasing" in p for p in problems))

    def test_estimate_outside_interval_rejected(self):
        problems = validate_claim(make_claim(estimate=1.5))
        self.assertTrue(any("confidence interval" in p for p in problems))

    def test_missing_source_location_rejected(self):
        problems = validate_claim(make_claim(source_location="  "))
        self.assertTrue(any("source_location" in p for p in problems))

    def test_boolean_not_accepted_as_number(self):
        problems = validate_claim(make_claim(outcome_cases=True))
        self.assertTrue(any("outcome_cases" in p for p in problems))

    def test_unknown_cohort_name_rejected(self):
        problems = validate_claim(make_claim(outcome_cohorts=["Finngen R9"]))
        self.assertTrue(any("unrecognised cohorts" in p for p in problems))

    def test_duplicate_ids_reported(self):
        report = validate_table([make_claim(), make_claim()])
        self.assertIn("duplicate claim_id", report["T1"])


class IndependenceTests(unittest.TestCase):
    def test_finngen_vs_ukb_independent(self):
        self.assertTrue(independent(["FinnGen"], ["UK Biobank"]))

    def test_same_cohort_not_independent(self):
        self.assertFalse(independent(["FinnGen"], ["FinnGen"]))

    def test_rahmioglu_vs_finngen_near_independent_only_through_nfbc(self):
        from claims import independence_note
        self.assertTrue(independent(["FinnGen"], ["Rahmioglu2023_meta"]))
        self.assertTrue(independent(["Rahmioglu2023_meta"], ["FinnGen"]))
        self.assertIn("NFBC", independence_note(["FinnGen"], ["Rahmioglu2023_meta"]))
        self.assertEqual(independence_note(["FinnGen"], ["UK Biobank"]), "independent")
        self.assertIsNone(independence_note(["FinnGen", "UK Biobank"], ["Rahmioglu2023_meta"]))
        self.assertFalse(independent(["Estonian Biobank", "FinnGen"], ["Rahmioglu2023_meta"]))
        self.assertFalse(independent(["Rahmioglu2023_meta"], ["UK Biobank"]))

    def test_primary_source_most_cases_and_combines_only_disjoint_sources(self):
        from tier_dryrun import eligible_sources, primary_sources
        finngen_claim = make_claim(outcome_cohorts=["FinnGen"])
        names = [e["source"] for e in eligible_sources(finngen_claim)]
        self.assertEqual(sorted(names), ["Pan-UKBB females", "Rahmioglu 2023 (public, no 23andMe)"])
        self.assertEqual([p["source"] for p in primary_sources(eligible_sources(finngen_claim))],
                         ["Rahmioglu 2023 (public, no 23andMe)"])  # shares UK Biobank with Pan-UKBB
        ukb_claim = make_claim(outcome_cohorts=["UK Biobank"])
        self.assertEqual([p["source"] for p in primary_sources(eligible_sources(ukb_claim))], ["FinnGen R12"])
        mvp_claim = make_claim(outcome_cohorts=["Million Veteran Program"])
        self.assertEqual([p["source"] for p in primary_sources(eligible_sources(mvp_claim))],
                         ["Rahmioglu 2023 (public, no 23andMe)", "FinnGen R12"])  # disjoint apart from NFBC

    def test_rahmioglu_overlaps_its_listed_cohorts(self):
        for cohort in ("Estonian Biobank", "deCODE", "other", "BioBank Japan", "UK Biobank"):
            self.assertFalse(independent([cohort], ["Rahmioglu2023_meta"]))
        self.assertTrue(independent(["Million Veteran Program"], ["Rahmioglu2023_meta"]))

    def test_unknown_cohort_never_independent(self):
        self.assertFalse(independent(["unknown"], ["UK Biobank"]))


class LiteratureTests(unittest.TestCase):
    def test_merge_by_doi_and_pmid(self):
        a = [{"pmid": "1", "doi": "10.1/x", "title": "Endometriosis pQTL", "abstract": "short"}]
        b = [{"pmid": "1", "doi": None, "title": "Endometriosis pQTL", "abstract": "a longer abstract"},
             {"pmid": None, "doi": "10.1/X", "title": "Same paper, upper-case DOI", "abstract": ""}]
        merged = merge([("pubmed", a), ("europepmc", b)])
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0]["found_by"], ["europepmc", "pubmed"])
        self.assertEqual(merged[0]["abstract"], "a longer abstract")
        self.assertTrue(merged[0]["gene_level_hint"])


class FetchTests(unittest.TestCase):
    def test_corruption_and_size_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "x.bin"
            path.write_bytes(b"original")
            record = {"md5": hashlib.md5(b"original").hexdigest(), "expected_bytes": 8}
            self.assertEqual(len(verify(path, record)["sha256"]), 64)
            path.write_bytes(b"corrupt!")
            with self.assertRaisesRegex(ValueError, "MD5"):
                verify(path, record)
            path.write_bytes(b"short")
            with self.assertRaisesRegex(ValueError, "byte count"):
                verify(path, record)

    def test_path_traversal_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = {"bad": {"filename": "../evil", "url": "https://example.invalid", "md5": "0", "source": "x"}}
            with self.assertRaisesRegex(ValueError, "basename"):
                fetch("bad", manifest, Path(directory))

    def test_leftover_partial_blocks_retry(self):
        with tempfile.TemporaryDirectory() as directory:
            (Path(directory) / "f.gz.partial").write_bytes(b"half")
            manifest = {"k": {"filename": "f.gz", "url": "https://example.invalid", "md5": "0", "source": "x"}}
            with self.assertRaisesRegex(ValueError, "Partial"):
                fetch("k", manifest, Path(directory))

    def test_no_source_md5_needs_policy_size_and_pinned_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "x.bin"
            path.write_bytes(b"original")
            with self.assertRaisesRegex(ValueError, "none_published"):
                verify(path, {"md5": None, "expected_bytes": 8})
            with self.assertRaisesRegex(ValueError, "none_published"):
                verify(path, {"md5": None, "checksum_policy": "none_published"})
            record = {"md5": None, "checksum_policy": "none_published", "expected_bytes": 8}
            pinned = verify(path, record)["sha256"]
            self.assertEqual(verify(path, record, pinned)["sha256"], pinned)
            path.write_bytes(b"corrupt!")
            with self.assertRaisesRegex(ValueError, "pinned"):
                verify(path, record, pinned)

    def test_truncated_gzip_without_source_md5_rejected(self):
        import gzip as gz
        with tempfile.TemporaryDirectory() as directory:
            whole = gz.compress(b"line\n" * 10000)
            path = Path(directory) / "x.gz"
            record = {"filename": "x.gz", "md5": None, "checksum_policy": "none_published"}
            path.write_bytes(whole)
            self.assertIn("gzip_crc", verify(path, dict(record, expected_bytes=len(whole))))
            path.write_bytes(whole[:-6])
            with self.assertRaisesRegex(ValueError, "gzip integrity"):
                verify(path, dict(record, expected_bytes=len(whole) - 6))

    def test_plain_http_needs_manifest_approval_and_flag(self):
        with tempfile.TemporaryDirectory() as directory:
            plain = {"filename": "f.gz", "url": "http://example.invalid/f.gz", "md5": "0", "source": "x"}
            with self.assertRaisesRegex(ValueError, "https only"):
                fetch("k", {"k": plain}, Path(directory), allow_http=True)
            approved = dict(plain, plain_http_approved="test")
            with self.assertRaisesRegex(ValueError, "https only"):
                fetch("k", {"k": approved}, Path(directory))
            ftp = dict(plain, url="ftp://example.invalid/f.gz", plain_http_approved="test")
            with self.assertRaisesRegex(ValueError, "https only"):
                fetch("k", {"k": ftp}, Path(directory), allow_http=True)

    def test_subdir_traversal_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            record = {"filename": "f.gz", "url": "https://example.invalid", "md5": "0", "source": "x",
                      "subdir": "../outside"}
            with self.assertRaisesRegex(ValueError, "subdir"):
                fetch("k", {"k": record}, Path(directory))


class GateBLocusTests(unittest.TestCase):
    def test_claims_near_positive_control_leads_are_capped_and_unknown_positions_count_as_near(self):
        from tier_dryrun import GATEB_LABEL, ceiling, gate_b_lead_near
        genes = {"NEAR": {"GRCh37": {"chrom": "2", "start": 11_000_000, "end": 11_100_000}},
                 "FAR": {"GRCh37": {"chrom": "2", "start": 40_000_000, "end": 40_100_000}}}
        leads = [{"rsid": "rs1", "chrom": "2", "pos": 11_721_535}]
        claim = {"gene_symbol": "NEAR", "instruments_status": "extracted", "instruments_check": {}, "instrument_scope": "cis"}
        self.assertEqual(gate_b_lead_near(claim, genes, leads), ["rs1"])
        self.assertEqual(gate_b_lead_near({**claim, "gene_symbol": "FAR"}, genes, leads), [])
        self.assertEqual(gate_b_lead_near({**claim, "gene_symbol": "MISSING"}, genes, leads), ["position unknown"])
        self.assertEqual(ceiling(claim, [{"source": "x"}], genes, leads), ("Tier 2", [GATEB_LABEL]))
        self.assertEqual(ceiling({**claim, "gene_symbol": "FAR"}, [{"source": "x"}], genes, leads), ("Tier 1", []))


class SupplementarySourceTests(unittest.TestCase):
    def test_supplementary_records_drop_those_already_found(self):
        from literature import supplementary_only
        main = [{"doi": "10.1/a", "pmid": "1", "title": "Alpha study"}]
        supp = [{"doi": "10.1/A", "pmid": None, "title": "x"}, {"doi": None, "pmid": "1", "title": "y"},
                {"doi": None, "pmid": None, "title": "Alpha <i>study</i>"}, {"doi": "10.1/b", "pmid": "2", "title": "Beta"}]
        self.assertEqual([r["doi"] for r in supplementary_only(main, supp)], ["10.1/b"])

    def test_review_sample_keeps_human_answers(self):
        from screening import write_sample
        rec = {"record_id": "r1", "title": "t", "screening": {"reason_code": "E1", "reason": "x"}}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "s.json"
            write_sample(path, [rec])
            rows = json.loads(path.read_text(encoding="utf-8"))
            rows[0].update(human_agrees=False, note="disagree")
            path.write_text(json.dumps(rows), encoding="utf-8")
            write_sample(path, [rec])
            row = json.loads(path.read_text(encoding="utf-8"))[0]
        self.assertEqual((row["human_agrees"], row["note"]), (False, "disagree"))

    def test_lead_checks_override_only_the_expected_screener_row(self):
        from screening import apply_lead_checks
        rows = [{"record_id": "r1", "decision": "exclude", "reason_code": "E3", "reason": "a"}]
        check = {"record_id": "r1", "screener_decision": "exclude", "screener_reason_code": "E3",
                 "decision": "include", "reason_code": "I", "reason": "b"}
        out = apply_lead_checks([dict(r) for r in rows], [check])[0]
        self.assertEqual((out["decision"], out["screener_reason_code"], out["lead_checked"]), ("include", "E3", True))
        with self.assertRaises(ValueError):
            apply_lead_checks([dict(r) for r in rows], [{**check, "screener_reason_code": "E1"}])

    def test_own_analysis_flag(self):
        from refscan import OWN_ANALYSIS
        self.assertTrue(OWN_ANALYSIS.search("In this study, we performed a two-sample Mendelian randomization analysis."))
        self.assertTrue(OWN_ANALYSIS.search("SMR analysis was performed with eQTLGen."))
        self.assertIsNone(OWN_ANALYSIS.search("Previous Mendelian randomization studies suggested a link."))


class RefScanTests(unittest.TestCase):
    def _raw(self, d, name, citing, routes):
        path = d / f"{name}.json"
        path.write_text(json.dumps({"paper": {"record_id": citing}, "routes": {
            k: {"status": "ok", "references": v} for k, v in routes.items()}}), encoding="utf-8")
        return path

    def test_routes_and_papers_merge_into_one_reference(self):
        from refscan import norm_doi, unique_references
        self.assertEqual(norm_doi("https://doi.org/10.1/ABC"), "10.1/abc")
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            a = self._raw(d, "a", "doi:p1", {
                "europepmc": [{"title": "A Synthetic Study.", "pmid": "111", "doi": None}],
                "crossref": [{"title": "A synthetic study", "doi": "10.1/x", "pmid": None}]})
            b = self._raw(d, "b", "doi:p2", {"crossref": [{"title": None, "doi": "10.1/x", "pmid": None}]})
            refs = unique_references([a, b])
        self.assertEqual(len(refs), 1)
        r = next(iter(refs.values()))
        self.assertEqual((r["pmid"], r["doi"], r["cited_by"]), ("111", "10.1/x", {"doi:p1", "doi:p2"}))

    def test_match_record_by_doi_pmid_or_title(self):
        from refscan import match_record, record_index
        idx = record_index([{"record_id": "doi:10.1/y", "doi": "10.1/y", "pmid": "5", "title": "Some Title: here"}])
        self.assertEqual(match_record({"doi": "10.1/y"}, idx), "doi:10.1/y")
        self.assertEqual(match_record({"pmid": "5"}, idx), "doi:10.1/y")
        self.assertEqual(match_record({"title": "some title here"}, idx), "doi:10.1/y")
        self.assertIsNone(match_record({"title": "other"}, idx))

    def test_candidate_decisions_must_cover_candidates_with_valid_codes(self):
        from refscan import merge_decisions
        cands = [{"ref_id": "doi:1"}, {"ref_id": "doi:2"}]
        ok = [{"ref_id": "doi:1", "decision": "exclude", "reason_code": "E1", "reason": "x"},
              {"ref_id": "doi:2", "decision": "include", "reason_code": "I", "reason": "y"}]
        self.assertEqual(len(merge_decisions(cands, ok)), 2)
        with self.assertRaises(ValueError):
            merge_decisions(cands, ok[:1])
        with self.assertRaises(ValueError):
            merge_decisions(cands, [ok[0], {**ok[1], "reason_code": "E2"}])


if __name__ == "__main__":
    unittest.main()
