import csv
import os
import subprocess
from collections import defaultdict, deque
from dataclasses import dataclass

import pysam
from Bio import SeqIO


@dataclass
class CleanConfig:
    min_af_review: float = 0.8
    min_dp_review: int = 10
    min_qual_review: float = 15.0
    dense_window_bp: int = 100
    dense_min_variants: int = 5


class VCFProcessor:
    def __init__(self, vcf_file, ref_fmt, output_dir, sample, reference_fasta=None, genbank_file=None, clean_config=None):
        self.vcf_file = vcf_file
        self.ref_fmt = ref_fmt
        self.output_dir = output_dir
        self.sample = sample
        self.reference_fasta = reference_fasta
        self.genbank_file = genbank_file
        self.clean_config = clean_config or CleanConfig()
        self.genbank_dict = self.load_genbank(self.genbank_file) if self.genbank_file else {}

    def load_genbank(self, genbank_file):
        genbank_dict = {}
        for record in SeqIO.parse(genbank_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS" and 'locus_tag' in feature.qualifiers:
                    locus_tag = feature.qualifiers['locus_tag'][0]
                    genbank_dict[locus_tag] = {
                        'product_id': feature.qualifiers.get('protein_id', [''])[0],
                        'product': feature.qualifiers.get('product', [''])[0],
                    }
        return genbank_dict

    def _run(self, cmd):
        subprocess.run(cmd, check=True)

    def normalize_sort_index(self):
        os.makedirs(self.output_dir, exist_ok=True)
        norm = os.path.join(self.output_dir, f"{self.sample}.tmp.norm.vcf.gz")
        sorted_vcf = os.path.join(self.output_dir, f"{self.sample}.tmp.norm.sorted.vcf.gz")
        self._run(["bcftools", "norm", "-f", os.path.abspath(self.reference_fasta), "-m", "-any", "-Oz", "-o", norm, self.vcf_file])
        self._run(["bcftools", "sort", "-Oz", "-o", sorted_vcf, norm])
        self._run(["bcftools", "index", "-f", sorted_vcf])
        return sorted_vcf

    def _is_dense(self, records):
        by_chrom = defaultdict(list)
        for i, rec in enumerate(records):
            by_chrom[rec.chrom].append((rec.pos, i))
        dense_idxs = set()
        w = self.clean_config.dense_window_bp
        k = self.clean_config.dense_min_variants
        for _, arr in by_chrom.items():
            dq = deque()
            for pos, idx in arr:
                dq.append((pos, idx))
                while dq and pos - dq[0][0] > w:
                    dq.popleft()
                if len(dq) >= k:
                    dense_idxs.update(i for _, i in dq)
        return dense_idxs

    def annotate_clean(self, input_vcf, output_vcf):
        invcf = pysam.VariantFile(input_vcf)
        header = invcf.header.copy()
        header.info.add("DVIONT_STATUS", 1, "String", "DviONT status: PASS/REVIEW/EXCLUDE")
        header.info.add("DVIONT_FLAGS", ".", "String", "DviONT QC flags")
        records = [r for r in invcf]
        dense_idxs = self._is_dense(records)

        outvcf = pysam.VariantFile(output_vcf, "wz", header=header)
        report_rows = []
        for i, rec in enumerate(records):
            flags = []
            filt = list(rec.filter.keys())
            qual = rec.qual if rec.qual is not None else 0.0
            gt = rec.samples[self.sample].get("GT") if self.sample in rec.samples else None
            dp = rec.samples[self.sample].get("DP") if self.sample in rec.samples else None
            ad = rec.samples[self.sample].get("AD") if self.sample in rec.samples else None
            af = rec.samples[self.sample].get("AF") if self.sample in rec.samples else rec.info.get("AF")
            af_val = af[0] if isinstance(af, (tuple, list)) and af else (float(af) if af is not None else None)

            if af_val is not None and af_val < self.clean_config.min_af_review:
                flags.append("LOW_AF")
            if dp is not None and dp < self.clean_config.min_dp_review:
                flags.append("LOW_DP")
            if qual < self.clean_config.min_qual_review:
                flags.append("LOW_QUAL")
            if gt and len(gt) >= 2 and gt[0] != gt[1]:
                flags.append("HET_LIKE")
            if len(rec.alts or []) > 1:
                flags.append("MULTIALLELIC")
            if any(str(a).startswith("<") for a in (rec.alts or [])):
                flags.append("SYMBOLIC_ALT")
            if filt and filt != ["PASS"]:
                flags.append("NONPASS_FILTER")
            if i in dense_idxs:
                flags.append("DENSE_REGION")

            status = "PASS"
            if "SYMBOLIC_ALT" in flags:
                status = "EXCLUDE"
            elif flags:
                status = "REVIEW"

            rec.info["DVIONT_STATUS"] = status
            rec.info["DVIONT_FLAGS"] = ",".join(sorted(set(flags))) if flags else ""
            outvcf.write(rec)

            variant_type = "SNP" if len(rec.ref) == 1 and all(len(a) == 1 for a in (rec.alts or [])) else "INDEL"
            report_rows.append([self.sample, rec.chrom, rec.pos, variant_type, rec.ref, ",".join(rec.alts or []), qual, ";".join(filt) if filt else "PASS", str(gt), dp, str(ad), af_val, status, rec.info.get("DVIONT_FLAGS", "")])
        outvcf.close()
        pysam.tabix_index(output_vcf, preset="vcf", force=True)
        return report_rows

    def write_report(self, rows):
        output_path = os.path.join(self.output_dir, f"{self.sample}.dviont.report.tsv")
        header = ["SAMPLE", "CHROM", "POS", "TYPE", "REF", "ALT", "QUAL", "FILTER", "GT", "DP", "AD", "AF", "DVIONT_STATUS", "DVIONT_FLAGS"]
        with open(output_path, "w", newline="") as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t')
            writer.writerow(header)
            writer.writerows(rows)
        return output_path
