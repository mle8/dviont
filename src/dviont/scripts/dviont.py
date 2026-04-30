#!/usr/bin/env python3
import argparse
import logging
import os
import shutil
import subprocess
from types import SimpleNamespace

from Bio import SeqIO

from .clair3_module import run_clair3
from .directory_management import PipelineManager
from .extract_fasta_and_gbk import extract_fasta_and_gbk
from .minimap2 import run_minimap2_alignment
from .ref_format import determine_ref_format
from .vcf_processor import CleanConfig, VCFProcessor
from .version import __version__

VCF_SOURCE_SUFFIX = {
    "clean": "{sample}.dviont.clean.vcf.gz",
    "clair3_raw": "{sample}.clair3.raw.vcf.gz",
    "legacy_merge": "{sample}.dviont.legacy.vcf.gz",
}


def require_tool(name):
    if shutil.which(name) is None:
        raise RuntimeError(f"Required executable not found in PATH: {name}")


def ensure_vcf_index(vcf_path):
    if not (os.path.exists(vcf_path + ".tbi") or os.path.exists(vcf_path + ".csi")):
        subprocess.run(["bcftools", "index", "-f", vcf_path], check=True)


def parse_reads_list(path):
    samples = []
    with open(path) as handle:
        for lineno, line in enumerate(handle, start=1):
            striped = line.strip()
            if not striped or striped.startswith("#"):
                continue
            parts = striped.split("\t")
            if len(parts) != 2:
                raise ValueError(f"Invalid samples.tsv format at line {lineno}: expected <sample_id><TAB><reads_path>")
            samples.append((parts[0], parts[1]))
    if not samples:
        raise ValueError("No samples were found in reads list after removing comments/blank lines")
    return samples


def resolve_sample_vcf(call_dir, sample, source):
    suffix = VCF_SOURCE_SUFFIX[source].format(sample=sample)
    return os.path.join(call_dir, suffix)


def add_common_call_args(parser):
    parser.add_argument("-o", "--output_dir", required=True)
    parser.add_argument("-r", "--ref", required=True)
    parser.add_argument("-i", "--reads", required=True)
    parser.add_argument("-t", "--threads", type=int, default=2)
    parser.add_argument("-m", "--model_name", default="r1041_e82_400bps_sup_v430_bacteria_finetuned")
    parser.add_argument("-s", "--sample", default="SAMPLE")
    parser.add_argument("-p", "--model_path", default=None)
    parser.add_argument("--preset", default="ont-q20")
    parser.add_argument("--vcf-mode", choices=["clair3_raw", "clean", "legacy_merge"], default="clean")
    parser.add_argument("--min-af-review", type=float, default=0.8)
    parser.add_argument("--min-dp-review", type=int, default=10)
    parser.add_argument("--min-qual-review", type=float, default=15.0)
    parser.add_argument("--dense-window-bp", type=int, default=100)
    parser.add_argument("--dense-min-variants", type=int, default=5)


def parse_args():
    parser = argparse.ArgumentParser(description="dviONT")
    sub = parser.add_subparsers(dest="command")
    call_parser = sub.add_parser("call")
    add_common_call_args(call_parser)

    cohort = sub.add_parser("cohort")
    cohort.add_argument("--ref", required=True)
    cohort.add_argument("--reads-list", required=True)
    cohort.add_argument("--out", required=True)
    cohort.add_argument("--threads", type=int, default=2)
    cohort.add_argument("--model-name", default="r1041_e82_400bps_sup_v430_bacteria_finetuned")
    cohort.add_argument("--model-path", default=None)
    cohort.add_argument("--preset", default="ont-q20")
    cohort.add_argument("--vcf-mode", choices=["clair3_raw", "clean", "legacy_merge"], default="clean")
    cohort.add_argument("--cohort-vcf-source", choices=["clean", "clair3_raw", "legacy_merge"], default="clean")
    cohort.add_argument("--recombination", choices=["none", "gubbins"], default="none")

    add_common_call_args(parser)
    parser.add_argument("-v", "--version", action="version", version=f"dviONT v{__version__}")
    return parser.parse_args()


def run_call(args):
    for exe in ["bcftools", "samtools", "minimap2", "run_clair3.sh"]:
        require_tool(exe)
    pm = PipelineManager(args.output_dir, args.sample)
    ref_dir = pm.create_output_directory()
    ref_fmt = determine_ref_format(args.ref)
    fasta_out = extract_fasta_and_gbk(args.ref, ref_dir, ref_fmt, args.output_dir)
    bam_output = run_minimap2_alignment(fasta_out, args.reads, args.threads, args.output_dir, args.sample)
    merge_vcf, clair3_dir = run_clair3(args.output_dir, fasta_out, bam_output, args.threads, args.model_name, args.sample, args.model_path)

    raw_vcf = os.path.join(args.output_dir, f"{args.sample}.clair3.raw.vcf.gz")
    shutil.copy2(merge_vcf, raw_vcf)
    for ext in (".tbi", ".csi"):
        if os.path.exists(merge_vcf + ext):
            shutil.copy2(merge_vcf + ext, raw_vcf + ext)
    ensure_vcf_index(raw_vcf)

    if args.vcf_mode == "clair3_raw":
        return raw_vcf
    if args.vcf_mode == "legacy_merge":
        legacy = os.path.join(args.output_dir, f"{args.sample}.dviont.legacy.vcf.gz")
        subprocess.run(["bcftools", "norm", "-m", "-any", raw_vcf, "-Oz", "-o", legacy], check=True)
        ensure_vcf_index(legacy)
        return legacy

    processor = VCFProcessor(vcf_file=raw_vcf, ref_fmt=ref_fmt, output_dir=args.output_dir, sample=args.sample, reference_fasta=fasta_out,
                             clean_config=CleanConfig(args.min_af_review, args.min_dp_review, args.min_qual_review, args.dense_window_bp, args.dense_min_variants))
    normed = processor.normalize_sort_index()
    clean_vcf = os.path.join(args.output_dir, f"{args.sample}.dviont.clean.vcf.gz")
    rows = processor.annotate_clean(normed, clean_vcf)
    processor.write_report(rows)
    return clean_vcf


def build_cohort(args):
    require_tool("bcftools")
    require_tool("snp-dists")
    for exe in ["minimap2", "samtools", "run_clair3.sh"]:
        require_tool(exe)
    if args.recombination == "gubbins":
        require_tool("run_gubbins.py")

    samples = parse_reads_list(args.reads_list)
    calls_dir = os.path.join(args.out, "calls")
    cohort_vcfs = os.path.join(args.out, "cohort_vcfs")
    align_dir = os.path.join(args.out, "alignments")
    consensus_dir = os.path.join(align_dir, "consensus_snps")
    dist_dir = os.path.join(args.out, "distances")
    gubbins_dir = os.path.join(args.out, "gubbins")
    for d in [calls_dir, cohort_vcfs, align_dir, consensus_dir, dist_dir]:
        os.makedirs(d, exist_ok=True)

    for sample, reads in samples:
        sample_out = os.path.join(calls_dir, sample)
        os.makedirs(sample_out, exist_ok=True)
        call_args = SimpleNamespace(output_dir=sample_out, ref=args.ref, reads=reads, threads=args.threads, model_name=args.model_name,
                                    sample=sample, model_path=args.model_path, preset=args.preset, vcf_mode=args.vcf_mode,
                                    min_af_review=0.8, min_dp_review=10, min_qual_review=15.0, dense_window_bp=100, dense_min_variants=5)
        run_call(call_args)

    selected_vcfs = []
    for sample, _ in samples:
        call_dir = os.path.join(calls_dir, sample)
        vcf = resolve_sample_vcf(call_dir, sample, args.cohort_vcf_source)
        if not os.path.exists(vcf):
            raise FileNotFoundError(f"Expected cohort VCF source missing for sample {sample}: {vcf}")
        ensure_vcf_index(vcf)
        selected_vcfs.append(vcf)

    merged = os.path.join(cohort_vcfs, "cohort_merged.vcf.gz")
    merged_norm = os.path.join(cohort_vcfs, "cohort_merged.norm.vcf.gz")
    merged_snps = os.path.join(cohort_vcfs, "cohort_merged.snps.vcf.gz")
    subprocess.run(["bcftools", "merge", "-m", "none", "-Oz", "-o", merged, *selected_vcfs], check=True)
    ensure_vcf_index(merged)
    subprocess.run(["bcftools", "norm", "-m", "-any", merged, "-Oz", "-o", merged_norm], check=True)
    ensure_vcf_index(merged_norm)
    subprocess.run(["bcftools", "view", "-i", 'strlen(REF)=1 && strlen(ALT)=1 && ALT!="*"', merged_norm, "-Oz", "-o", merged_snps], check=True)
    ensure_vcf_index(merged_snps)

    sample_list = subprocess.check_output(["bcftools", "query", "-l", merged_snps], text=True).strip().splitlines()
    alignment_file = os.path.join(align_dir, "cohort.snp_alignment.fasta")
    with open(alignment_file, "w") as outfa:
        for sample in sample_list:
            seq = subprocess.check_output(["bcftools", "consensus", "-f", args.ref, "-s", sample, merged_snps], text=True)
            lines = [l.strip() for l in seq.splitlines() if l.strip()]
            sequence = "".join([l for l in lines if not l.startswith(">")])
            sample_fa = os.path.join(consensus_dir, f"{sample}.fa")
            with open(sample_fa, "w") as sf:
                sf.write(f">{sample}\n{sequence}\n")
            outfa.write(f">{sample}\n{sequence}\n")

    records = list(SeqIO.parse(alignment_file, "fasta"))
    if len(records) != len(sample_list):
        raise RuntimeError(f"Alignment record count mismatch: expected {len(sample_list)} found {len(records)}")
    lengths = sorted({len(r.seq) for r in records})
    with open(os.path.join(align_dir, "alignment_lengths.txt"), "w") as h:
        for l in lengths:
            h.write(f"{l}\n")
    if len(lengths) != 1:
        raise RuntimeError(f"SNP alignment sequences have non-identical lengths: {lengths}")

    with open(os.path.join(dist_dir, "cohort.unmasked_snp_distance_matrix.tsv"), "w") as h:
        subprocess.run(["snp-dists", alignment_file], check=True, stdout=h)

    if args.recombination == "gubbins":
        os.makedirs(gubbins_dir, exist_ok=True)
        prefix = os.path.join(gubbins_dir, "cohort")
        subprocess.run(["run_gubbins.py", "--prefix", prefix, "--threads", str(args.threads), alignment_file], check=True)
        masked = f"{prefix}.filtered_polymorphic_sites.fasta"
        if not os.path.exists(masked):
            raise FileNotFoundError(f"Gubbins completed but expected output missing: {masked}")
        with open(os.path.join(dist_dir, "cohort.masked_snp_distance_matrix.tsv"), "w") as h:
            subprocess.run(["snp-dists", masked], check=True, stdout=h)


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO)
    if args.command == "cohort":
        build_cohort(args)
    else:
        run_call(args)


if __name__ == "__main__":
    main()
