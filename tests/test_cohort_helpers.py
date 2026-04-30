import os
import tempfile
import unittest

from pathlib import Path

from dviont.scripts.dviont import parse_reads_list, resolve_sample_vcf, write_cohort_helper_scripts


class CohortHelperTests(unittest.TestCase):
    def test_parse_reads_list_skips_comments_and_blanks(self):
        with tempfile.NamedTemporaryFile("w", delete=False) as h:
            h.write("# comment\n\nS1\t/a.fastq.gz\nS2\t/b.fastq.gz\n")
            path = h.name
        try:
            rows = parse_reads_list(path)
            self.assertEqual(rows, [("S1", "/a.fastq.gz"), ("S2", "/b.fastq.gz")])
        finally:
            os.unlink(path)

    def test_resolve_vcf(self):
        self.assertTrue(resolve_sample_vcf("calls/S1", "S1", "clean").endswith("S1.dviont.clean.vcf.gz"))
        self.assertTrue(resolve_sample_vcf("calls/S1", "S1", "clair3_raw").endswith("S1.clair3.raw.vcf.gz"))
        self.assertTrue(resolve_sample_vcf("calls/S1", "S1", "legacy_merge").endswith("S1.dviont.legacy.vcf.gz"))

    def test_write_helper_scripts(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            align = out / "alignments" / "cohort.snp_alignment.fasta"
            align.parent.mkdir(parents=True)
            align.write_text(">S1\nACGT\n")
            write_cohort_helper_scripts(out, align)
            g = out / "run_gubbins.sh"
            m = out / "run_masked_snp_dists.sh"
            self.assertTrue(g.exists())
            self.assertTrue(m.exists())
            self.assertIn("GUBBINS_ENV=\"/path/to/gubbins/env\"", g.read_text())
            self.assertIn("cohort.masked_snp_distance_matrix.tsv", m.read_text())


if __name__ == "__main__":
    unittest.main()
