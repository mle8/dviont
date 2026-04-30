<img src="src/dviont/docs/DVIONT.png" width="50%">

# DviONT

DviONT (DNA Variant Identification using ONT) is a bacterial long-read variant pipeline by William Shropshire. Reads are aligned with minimap2 (or winnowmap), variants are called with Clair3, and DviONT standardizes VCF output for reporting and downstream comparison.

The default Clair3 model is `r1041_e82_400bps_sup_v430_bacteria_finetuned`.

## Features

- FASTA or GenBank reference input.
- ONT long-read alignment with minimap2 or winnowmap.
- Clair3 variant calling with selectable `--vcf-mode`:
  - `clean` (default): normalize/sort/index + conservative QC annotations.
  - `clair3_raw`: preserve Clair3 merged output for benchmarking.
  - `legacy_merge`: retain historical DviONT legacy merge output.
- Optional SnpEff/GenBank report fields when annotation is available.
- Cohort SNP alignment from per-sample VCFs.
- SNP distance matrix with `snp-dists`.
- Optional Gubbins masking in cohort mode.

## Installation

> Tested primarily on Linux/HPC environments (RHEL-like systems).

### 1) Create environment

```bash
git clone https://github.com/wshropshire/dviont
cd dviont

conda create -n dviont_env python=3.10
conda activate dviont_env

pip install build
python -m build
pip install ./dist/dviont-0.3.1.tar.gz

# Install pipeline dependencies
conda env update --name dviont_env --file environment.yml
```

### 2) Optional: separate Gubbins environment

```bash
conda env create -n dviont_gubbins -f environment-gubbins.yml
```

### 3) Clair3 models

```bash
./src/dviont/bin/download_clair3_models --output-dir models
```

If no model list is provided, the downloader uses default models including `r1041_e82_400bps_sup_v430_bacteria_finetuned`.

### Required executables

DviONT expects these tools on `PATH`:
- `minimap2` (or `winnowmap` when selected)
- `samtools`
- `bcftools`
- `run_clair3.sh`
- `snp-dists` (cohort distance matrix)
- `run_gubbins.py` (only when `--recombination gubbins`)

## Single-sample usage

### Preferred subcommand

```bash
dviont call \
  -o out/SAMPLE1 \
  -r ref.fasta \
  -i SAMPLE1.fastq.gz \
  -s SAMPLE1 \
  -t 16 \
  -m r1041_e82_400bps_sup_v430_bacteria_finetuned \
  -p /path/to/clair3/models \
  --preset ont-q20 \
  --aligner minimap2 \
  --vcf-mode clean
```

### Backward-compatible single command style

```bash
dviont \
  -o out/SAMPLE1 \
  -r ref.fasta \
  -i SAMPLE1.fastq.gz \
  -s SAMPLE1 \
  --vcf-mode clean
```

## Cohort usage

`samples.tsv` format:

```text
sample_id<TAB>reads_path
```

Example:

```bash
dviont cohort \
  --ref ref.fasta \
  --reads-list examples/samples.tsv \
  --out cohort_out \
  --threads 16 \
  --model-name r1041_e82_400bps_sup_v430_bacteria_finetuned \
  --model-path /path/to/clair3/models \
  --vcf-mode clean \
  --cohort-vcf-source clean \
  --recombination none
```

Use `--recombination gubbins` to run Gubbins masking after SNP alignment generation.

## Outputs

### Single sample

- `sample.clair3.raw.vcf.gz` (+ index)
- `sample.dviont.clean.vcf.gz` (+ index) in clean mode
- `sample.dviont.legacy.vcf.gz` (+ index) in legacy mode
- `sample.dviont.report.tsv`
- `sample.consensus.fasta`

### Cohort

- `cohort_vcfs/cohort_merged.vcf.gz`
- `cohort_vcfs/cohort_merged.norm.vcf.gz`
- `cohort_vcfs/cohort_merged.snps.vcf.gz`
- `alignments/cohort.snp_alignment.fasta`
- `distances/cohort.unmasked_snp_distance_matrix.tsv`
- `distances/cohort.masked_snp_distance_matrix.tsv` (when Gubbins is run)

## Notes

- Clean mode is conservative: it standardizes and annotates VCFs without aggressive default filtering.
- Dense variant regions are flagged as `DENSE_REGION`; this is not a direct recombination call.
- Recombination-aware masking is optional and only performed in cohort mode via Gubbins.

## Citation / Authors

DviONT was developed by **William Shropshire**.

## Examples
See `examples/` for call/cohort and LSF workflow templates.
dviONT v0.3.1
