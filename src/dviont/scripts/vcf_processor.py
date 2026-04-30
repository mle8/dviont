import csv
import os
import subprocess
from collections import defaultdict, deque
from dataclasses import dataclass
from typing import List, Tuple

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
    def __init__(self, vcf_file, output_dir, sample, reference_fasta, genbank_file=None, clean_config=None):
        self.vcf_file = vcf_file
        self.output_dir = output_dir
        self.sample = sample
        self.reference_fasta = reference_fasta
        self.genbank_file = genbank_file
        self.clean_config = clean_config or CleanConfig()
        self.genbank_dict = self._load_genbank(genbank_file) if genbank_file else {}

    def _run(self, cmd):
        subprocess.run(cmd, check=True)

    def _load_genbank(self, genbank_file):
        out = {}
        for record in SeqIO.parse(genbank_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS" and "locus_tag" in feature.qualifiers:
                    locus_tag = feature.qualifiers["locus_tag"][0]
                    out[locus_tag] = {
                        "product_id": feature.qualifiers.get("protein_id", [""])[0],
                        "product": feature.qualifiers.get("product", [""])[0],
                    }
        return out

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
        for arr in by_chrom.values():
            dq = deque()
            for pos, idx in arr:
                dq.append((pos, idx))
                while dq and pos - dq[0][0] > self.clean_config.dense_window_bp:
                    dq.popleft()
                if len(dq) >= self.clean_config.dense_min_variants:
                    dense_idxs.update(j for _, j in dq)
        return dense_idxs

    def _ann_extract(self, rec):
        ann = rec.info.get("ANN")
        if not ann:
            return [""] * 8
        first = ann[0].split("|") if isinstance(ann, (list, tuple)) else str(ann).split("|")
        annot = first[1] if len(first) > 1 else ""
        impact = first[2] if len(first) > 2 else ""
        gene = first[3] if len(first) > 3 else ""
        locus_tag = first[4] if len(first) > 4 else ""
        hgvs_c = first[9] if len(first) > 9 else ""
        hgvs_p = first[10] if len(first) > 10 else ""
        product_id = self.genbank_dict.get(locus_tag, {}).get("product_id", "")
        product = self.genbank_dict.get(locus_tag, {}).get("product", "")
        return [annot, impact, gene, locus_tag, hgvs_c, hgvs_p, product_id, product]

    def annotate_clean(self, input_vcf, output_vcf) -> List[List]:
        invcf = pysam.VariantFile(input_vcf)
        header = invcf.header.copy()
        if "DVIONT_STATUS" not in header.info:
            header.info.add("DVIONT_STATUS", 1, "String", "DviONT status: PASS/REVIEW/EXCLUDE")
        if "DVIONT_FLAGS" not in header.info:
            header.info.add("DVIONT_FLAGS", 1, "String", "DviONT QC flags")

        records = [r for r in invcf]
        dense_idxs = self._is_dense(records)
        outvcf = pysam.VariantFile(output_vcf, "wz", header=header)
        rows = []

        for i, rec in enumerate(records):
            flags = []
            filt = list(rec.filter.keys())
            qual = rec.qual if rec.qual is not None else 0.0
            sample_data = rec.samples[self.sample] if self.sample in rec.samples else next(iter(rec.samples.values()))
            gt = sample_data.get("GT")
            dp = sample_data.get("DP")
            ad = sample_data.get("AD")
            af = sample_data.get("AF")
            if af is None:
                af = rec.info.get("AF")
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

            status = "PASS" if not flags else "REVIEW"
            if "SYMBOLIC_ALT" in flags:
                status = "EXCLUDE"

            rec.info["DVIONT_STATUS"] = status
            rec.info["DVIONT_FLAGS"] = ",".join(sorted(set(flags))) if flags else ""
            outvcf.write(rec)

            vartype = "SNP" if len(rec.ref) == 1 and all(len(a) == 1 for a in (rec.alts or [])) else "INDEL"
            rows.append([
                self.sample, rec.chrom, rec.pos, vartype, rec.ref, ",".join(rec.alts or []), qual,
                ";".join(filt) if filt else "PASS", str(gt), dp, str(ad), af_val, status,
                rec.info.get("DVIONT_FLAGS", ""), *self._ann_extract(rec)
            ])

        outvcf.close()
        pysam.tabix_index(output_vcf, preset="vcf", force=True)
        return rows

    def write_report(self, rows):
        output = os.path.join(self.output_dir, f"{self.sample}.dviont.report.tsv")
        header = ["SAMPLE", "CHROM", "POS", "TYPE", "REF", "ALT", "QUAL", "FILTER", "GT", "DP", "AD", "AF", "DVIONT_STATUS", "DVIONT_FLAGS", "ANNOT", "IMPACT", "GENE", "LOCUS_TAG", "HGVS.c", "HGVS.p", "PRODUCT_ID", "PRODUCT"]
        with open(output, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(header)
            w.writerows(rows)
        return output

    def process_clean(self, output_vcf):
        sorted_vcf = self.normalize_sort_index()
        rows = self.annotate_clean(sorted_vcf, output_vcf)
        report = self.write_report(rows)
        return output_vcf, report
