#!/bin/bash
#BSUB -q <queue>
#BSUB -n 16
#BSUB -M 64000
source ~/.bashrc
conda activate dviont_env
PROJECT_DIR=/path/to/project
REF_GENOME=$PROJECT_DIR/ref.fasta
LONG_READS=$PROJECT_DIR/reads
OUTPUT_BASE=$PROJECT_DIR/calls
MODEL_NAME=r1041_e82_400bps_sup_v430_bacteria_finetuned
MODEL_PATH=/path/to/model
MAX_PARALLEL=4
