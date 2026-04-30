import os
import tempfile
import unittest

from dviont.scripts.dviont import parse_reads_list, resolve_sample_vcf


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


if __name__ == "__main__":
    unittest.main()
