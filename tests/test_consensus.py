import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from dviont.scripts.dviont import write_consensus


@unittest.skipUnless(shutil.which("bcftools"), "bcftools is required")
class ConsensusTest(unittest.TestCase):
    def test_alt_ref_and_missing_genotypes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference.fa"
            reference.write_text(">chr1\nACGTA\n")
            vcf = root / "merged.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "##contig=<ID=chr1,length=5>\n"
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\talt\tref\tmissing\n"
                "chr1\t3\t.\tG\tT\t.\tPASS\t.\tGT\t1/1\t0/0\t./.\n"
            )
            compressed = root / "merged.vcf.gz"
            subprocess.run(
                ["bcftools", "view", "-Oz", "-o", str(compressed), str(vcf)],
                check=True, capture_output=True,
            )
            subprocess.run(
                ["bcftools", "index", str(compressed)],
                check=True, capture_output=True,
            )
            for sample, expected in (
                ("alt", "ACTTA"), ("ref", "ACGTA"), ("missing", "ACNTA")
            ):
                with self.subTest(sample=sample):
                    output = root / f"{sample}.fa"
                    length = write_consensus(compressed, reference, sample, output)
                    self.assertEqual(length, 5)
                    self.assertEqual(output.read_text(), f">{sample}\n{expected}\n")


if __name__ == "__main__":
    unittest.main()
