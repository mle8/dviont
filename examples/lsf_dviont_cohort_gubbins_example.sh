#!/bin/bash
#BSUB -q medium
#BSUB -n 16
#BSUB -M 64000
#BSUB -R "rusage[mem=64000]"

source ~/.bashrc
conda activate dviont_env

PROJECT_DIR=/path/to/project
REF=${PROJECT_DIR}/ref.fasta
READS_LIST=${PROJECT_DIR}/samples.tsv
OUT=${PROJECT_DIR}/cohort_out
CLAIR3_MODEL_PATH=/path/to/clair3/model

dviont cohort \
  --ref "$REF" \
  --reads-list "$READS_LIST" \
  --out "$OUT" \
  --threads 16 \
  --model-name r1041_e82_400bps_sup_v430_bacteria_finetuned \
  --model-path "$CLAIR3_MODEL_PATH" \
  --preset ont-q20 \
  --vcf-mode clean \
  --cohort-vcf-source clean \
  --recombination gubbins
