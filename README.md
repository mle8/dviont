# dviONT

dviONT is a bacterial ONT variant workflow where **Clair3 performs variant calling** and dviONT performs conservative VCF standardization, lightweight QC annotation, reporting, cohort SNP alignment generation, and SNP distance analysis.

## VCF modes
- `clair3_raw`: preserves Clair3 `merge_output.vcf.gz` as `sample.clair3.raw.vcf.gz`.
- `clean` (default): standardizes and annotates variants as `sample.dviont.clean.vcf.gz` + `sample.dviont.report.tsv`.
- `legacy_merge`: historical mode output `sample.dviont.legacy.vcf.gz` for benchmarking.

## Important interpretation notes
- Dense region flags in single-sample mode are **not recombination calls**.
- Dense regions can reflect multiple processes (mapping instability, structural effects, divergence, etc.).
- Recombination-aware masking is only done in **cohort mode with Gubbins**.

## Single-sample call
```bash
dviont call \
  -o sample_out \
  -r ref.fasta \
  -i sample.fastq.gz \
  -t 16 \
  -m r1041_e82_400bps_sup_v430_bacteria_finetuned \
  -p /path/to/model \
  --preset ont-q20 \
  -s SAMPLE \
  --vcf-mode clean
```

## Cohort mode
`samples.tsv` format:
```text
sample_id<TAB>reads_path
```
Allows `#` comment lines and blank lines.

Example:
```bash
dviont cohort \
  --ref ref.fasta \
  --reads-list samples.tsv \
  --out cohort_out \
  --threads 16 \
  --model-name r1041_e82_400bps_sup_v430_bacteria_finetuned \
  --model-path /path/to/model \
  --preset ont-q20 \
  --vcf-mode clean \
  --cohort-vcf-source clean \
  --recombination none
```

`--cohort-vcf-source` options:
- `clean` (default): `calls/{sample}/{sample}.dviont.clean.vcf.gz`
- `clair3_raw`: `calls/{sample}/{sample}.clair3.raw.vcf.gz`
- `legacy_merge`: `calls/{sample}/{sample}.dviont.legacy.vcf.gz`

`--recombination` options:
- `none` (default): produce unmasked SNP distances only.
- `gubbins`: run Gubbins and also generate masked SNP distances.

## Cohort outputs
```text
cohort_out/
  calls/
  cohort_vcfs/
    cohort_merged.vcf.gz
    cohort_merged.norm.vcf.gz
    cohort_merged.snps.vcf.gz
  alignments/
    cohort.snp_alignment.fasta
    consensus_snps/
    alignment_lengths.txt
  distances/
    cohort.unmasked_snp_distance_matrix.tsv
    cohort.masked_snp_distance_matrix.tsv   # if gubbins enabled
  gubbins/
    cohort.filtered_polymorphic_sites.fasta
    cohort.final_tree.tre
    cohort.recombination_predictions.gff
```

`cohort.snp_alignment.fasta` is a **SNP-only alignment**, not a full whole-genome alignment.

## Environments
- `environment.yml`: dviONT + core tools (includes `snp-dists`).
- `environment-gubbins.yml`: optional separate environment for Gubbins.

## Examples
See `examples/` for call/cohort and LSF workflow templates.
