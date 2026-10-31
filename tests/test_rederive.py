"""Synthetic tests for src/rederive.py (needs numpy: run with the project venv).

  .venv\\Scripts\\python -I -m unittest discover -s tests -p "test_rederive.py" -v
"""
import io
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

try:
    import numpy as np
except ImportError:  # system Python without numpy: skip
    np = None

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))


def write_plink(prefix, genotypes, variants):
    """genotypes: list per variant of A1 dosages (0/1/2, None missing); variants: (chrom, rsid, pos, a1, a2)."""
    n = len(genotypes[0])
    Path(str(prefix) + ".fam").write_text("".join(f"f{i} i{i} 0 0 2 -9\n" for i in range(n)))
    Path(str(prefix) + ".bim").write_text("".join(f"{c}\t{rs}\t0\t{p}\t{a1}\t{a2}\n" for c, rs, p, a1, a2 in variants))
    code = {2: 0b00, None: 0b01, 1: 0b10, 0: 0b11}
    out = bytearray(b"\x6c\x1b\x01")
    for row in genotypes:
        for start in range(0, n, 4):
            byte = 0
            for k, g in enumerate(row[start:start + 4]):
                byte |= code[g] << (2 * k)
            out.append(byte)
    Path(str(prefix) + ".bed").write_bytes(bytes(out))


@unittest.skipIf(np is None, "numpy not installed (use the project venv)")
class RederiveTests(unittest.TestCase):
    def setUp(self):
        import rederive
        self.rd = rederive
        self.tmp = tempfile.TemporaryDirectory()
        base = [0, 1, 2, 1, 0, 2, 1, 0, 1, 2]
        other = [2, 1, 2, 1, 1, 1, 1, 1, 1, 1]  # uncorrelated with base (r = 0)
        self.prefix = Path(self.tmp.name) / "EUR"
        write_plink(self.prefix,
                    [base, base, other, [None] + base[1:], [2 - g for g in base]],
                    [("1", "rs1", 1000, "A", "G"), ("1", "rs2", 2000, "C", "T"), ("1", "rs3", 3000, "A", "C"),
                     ("1", "rs4", 50_000_000, "G", "T"), ("2", "rs5", 1000, "A", "T")])
        self.panel = rederive.RefPanel(self.prefix)

    def tearDown(self):
        self.panel.close()
        self.tmp.cleanup()

    def test_bed_decoding_and_frequency(self):
        d = self.panel.dosages([0, 3])
        self.assertEqual(d.shape, (2, 10))
        self.assertEqual(list(d[0]), [0, 1, 2, 1, 0, 2, 1, 0, 1, 2])
        self.assertTrue(np.isnan(d[1, 0]))
        self.assertAlmostEqual(float(self.panel.a1_freq([0])[0]), 10 / 20)

    def test_bed_size_mismatch_rejected(self):
        self.panel.close()
        bed = Path(str(self.prefix) + ".bed")
        bed.write_bytes(bed.read_bytes() + b"\x00")
        with self.assertRaisesRegex(ValueError, "size"):
            self.rd.RefPanel(self.prefix)
        self.panel = self.rd.RefPanel.__new__(self.rd.RefPanel)
        self.panel.bed = None

    def test_find_by_position_checks_alleles_on_either_strand(self):
        self.assertEqual(self.panel.find(None, "1", 1000, {"A", "G"}), 0)
        self.assertEqual(self.panel.find(None, "chr1", 1000, {"T", "C"}), 0)
        self.assertIsNone(self.panel.find(None, "1", 1000, {"A", "C"}))
        self.assertEqual(self.panel.find("rs3"), 2)

    def test_clumping_removes_correlated_within_distance_only(self):
        def cand(rs, chrom, pos, nl):
            return {"rsid": rs, "chrom": chrom, "pos": pos, "neglog10p": nl, "ref": self.panel.find(rs)}
        cands = [cand("rs2", "1", 2000, 9), cand("rs1", "1", 1000, 10), cand("rs3", "1", 3000, 8),
                 cand("rs4", "1", 50_000_000, 7), cand("rsX", "1", 4000, 20)]
        kept, absent = self.rd.clump(cands, self.panel, 0.001, 10_000)
        self.assertEqual([c["rsid"] for c in kept], ["rs1", "rs3", "rs4"])  # rs2 = rs1; rs4 too far (49 Mb)
        self.assertEqual([c["rsid"] for c in absent], ["rsX"])
        kept, _ = self.rd.clump(cands[:4], self.panel, 0.001, 50_000)
        self.assertEqual([c["rsid"] for c in kept], ["rs1", "rs3"])  # rs4 now within distance and correlated
        kept, _ = self.rd.clump(cands[:2], self.panel, 1.0, 10_000)
        self.assertEqual(len(kept), 2)  # r2 must exceed the threshold

    def test_freq_of_handles_flips_and_palindromes(self):
        self.assertAlmostEqual(self.panel.freq_of(0, "A"), 0.5)
        self.assertAlmostEqual(self.panel.freq_of(0, "T"), 0.5)  # complement of A1 on a non-palindromic SNP
        self.assertAlmostEqual(self.panel.freq_of(0, "C"), 0.5)  # complement of A2
        self.assertIsNone(self.panel.freq_of(0, "AT"))
        self.assertAlmostEqual(self.panel.freq_of(4, "T"), 1 - float(self.panel.a1_freq([4])[0]))

    def test_unpack_takes_wanted_regular_files_by_basename(self):
        tar_path = Path(self.tmp.name) / "a.tgz"
        with tarfile.open(tar_path, "w:gz") as tar:
            for name, data in [("x/../../EUR.bim", b"bim"), ("x/AFR.bim", b"no"), ("x/EUR.fam", b"fam")]:
                info = tarfile.TarInfo(name)
                info.size = len(data)
                tar.addfile(info, io.BytesIO(data))
        dest = Path(self.tmp.name) / "out"
        self.rd.unpack_members(tar_path, {"EUR.bim", "EUR.fam"}, dest)
        self.assertEqual(sorted(p.name for p in dest.iterdir()), ["EUR.bim", "EUR.fam"])
        self.assertEqual((dest / "EUR.bim").read_bytes(), b"bim")
        with tarfile.open(tar_path, "w:gz") as tar:
            info = tarfile.TarInfo("EUR.bed")
            info.type = tarfile.SYMTYPE
            info.linkname = "/etc/passwd"
            tar.addfile(info)
        with self.assertRaisesRegex(ValueError, "regular file"):
            self.rd.unpack_members(tar_path, {"EUR.bed"}, dest)

    def test_rule_order_window_threshold_clump_and_F(self):
        rows = [{"rsid": "rs1", "chrom": "1", "pos": 1000, "ea": "A", "oa": "G", "beta": 0.5, "se": 0.05,
                 "eaf": 0.3, "p": 1e-20, "neglog10p": 20, "ref": 0},
                {"rsid": "rs2", "chrom": "1", "pos": 2000, "ea": "C", "oa": "T", "beta": 0.4, "se": 0.05,
                 "eaf": 0.3, "p": 1e-15, "neglog10p": 15, "ref": 1},
                {"rsid": "rs3", "chrom": "1", "pos": 3000, "ea": "A", "oa": "C", "beta": 0.02, "se": 0.0045,
                 "eaf": 0.3, "p": 1e-5, "neglog10p": 5, "ref": 2},
                {"rsid": "rs5", "chrom": "2", "pos": 1000, "ea": "A", "oa": "T", "beta": 0.5, "se": 0.05,
                 "eaf": 0.5, "p": 1e-20, "neglog10p": 20, "ref": 4}]
        kept, log = self.rd.apply_rule(rows, self.panel, 1e-4, 0.001, 10_000, f_min=20, palin_eaf=(0.42, 0.58))
        self.assertEqual([r["rsid"] for r in kept], ["rs1"])
        self.assertEqual(log, {"pass_p": 4, "palindromic_dropped": 1, "absent_from_reference": 0,
                               "after_clumping": 2, "dropped_F": 1})
        kept, _ = self.rd.apply_rule(rows, self.panel, 1e-4, 0.001, 10_000, region=("1", 2500, 3500))
        self.assertEqual([r["rsid"] for r in kept], ["rs3"])


if __name__ == "__main__":
    unittest.main()
