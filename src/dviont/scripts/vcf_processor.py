import csv
import os
import gzip

class VCFProcessor:
    def __init__(self, vcf_file, ref_fmt, output_dir, sample, reference_fasta=None, genbank_file=None, clean_config=None):
        self.vcf_file = vcf_file
        self.ref_fmt = ref_fmt.lower()
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

    def parse_vcf(self):
        """Parse the VCF file and generate the summary report."""
        header = ["CHROM", "POS", "TYPE", "REF", "ALT", "EVIDENCE"]
        if self.ref_fmt == "genbank":
            header += ["ANNOT", "IMPACT", "GENE", "LOCUS_TAG", "HGVS.c", "HGVS.p", "PRODUCT_ID", "PRODUCT"]

        rows = []

        open_func = gzip.open if self.vcf_file.endswith(".gz") else open
        with open_func(self.vcf_file, 'rt') as vcf:
            for line in vcf:
                if line.startswith("#"):
                    continue

                fields = line.strip().split("\t")
                # Minimum 8 columns for a valid VCF
                if len(fields) < 8:
                    continue

                chrom, pos, _, ref, alt, _, _, info, *rest = fields

                # Determine TYPE
                if len(ref) == 1 and len(alt) == 1:
                    variant_type = "SNP"
                elif len(ref) == len(alt):
                    variant_type = "MNP"
                else:
                    variant_type = "INDEL"

                # Extract EVIDENCE
                evidence = "NA"
                if rest:
                    sample_fields = rest[-1].split(":")  # FORMAT: GT:GQ:DP:AD:AF
                    if len(sample_fields) >= 4:
                        dp = sample_fields[2]
                        ad_values = sample_fields[3].split(",")  # e.g. "15,19"
                        if len(ad_values) >= 2:
                            ref_count, alt_count = ad_values[0], ad_values[1]
                            evidence = f"ALT:{dp}/{alt_count};REF:{dp}/{ref_count}"

                row = [chrom, pos, variant_type, ref, alt, evidence]

                # Only add GenBank annotations if requested
                if self.ref_fmt == "genbank" and "ANN=" in info:
                    info_fields = info.split("|")
                    annot = info_fields[1] if len(info_fields) > 1 else ""
                    impact = info_fields[2] if len(info_fields) > 2 else ""
                    gene = info_fields[3] if len(info_fields) > 3 else ""
                    locus_tag = info_fields[4] if len(info_fields) > 4 else ""

                    if annot == "intergenic_region":
                        gene = ""
                    elif gene == locus_tag:
                        gene = "hyp"

                    hgvs_c = info_fields[9] if len(info_fields) > 9 else ""
                    hgvs_p = info_fields[10] if len(info_fields) > 10 else ""

                    product_id, product = "", ""
                    if locus_tag in self.genbank_dict:
                        product_id = self.genbank_dict[locus_tag]["product_id"]
                        product = self.genbank_dict[locus_tag]["product"]

                    row += [annot, impact, gene, locus_tag, hgvs_c, hgvs_p, product_id, product]

                rows.append(row)

        # Write to TSV file
        self.write_tsv(header, rows)

    def write_tsv(self, header, rows):
        """Write the summary data to a TSV file."""
        os.makedirs(self.output_dir, exist_ok=True)
        output_path = os.path.join(self.output_dir, f"{self.sample}_dviont_report.tsv")

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
