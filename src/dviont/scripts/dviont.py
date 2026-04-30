#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess
from pathlib import Path



def check_executable(exe, hint=None):
    if shutil.which(exe) is None:
        msg = f"Required executable not found on PATH: {exe}"
        if hint:
            msg += f" ({hint})"
        raise RuntimeError(msg)


def parse_reads_list(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            sample, reads = line.split("\t", 1)
            rows.append((sample, reads))
    return rows


def resolve_sample_vcf(sample_outdir, sample, source):
    if source == "clean":
        return os.path.join(sample_outdir, f"{sample}.dviont.clean.vcf.gz")
    if source == "clair3_raw":
        return os.path.join(sample_outdir, f"{sample}.clair3.raw.vcf.gz")
    return os.path.join(sample_outdir, f"{sample}.dviont.legacy.vcf.gz")


def run_cmd(cmd):
    subprocess.run(cmd, check=True)


def run_call(args):
    from .clair3_module import run_clair3
    from .extract_fasta_and_gbk import extract_fasta_and_gbk
    from .minimap2 import run_minimap2_alignment
    from .ref_format import determine_ref_format
    from .vcf_processor import VCFProcessor, CleanConfig
    from .winnowmap import run_winnowmap_alignment
    for exe in ["bcftools", "samtools", "minimap2", "run_clair3.sh"]:
        check_executable(exe)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    ref_fmt = determine_ref_format(args.ref)
    fasta = extract_fasta_and_gbk(args.ref, str(out / "reference"), ref_fmt, str(out))
    if args.aligner == "winnowmap":
        bam = run_winnowmap_alignment(fasta, args.reads, args.threads, args.output_dir, args.sample, preset=args.preset)
    else:
        bam = run_minimap2_alignment(fasta, args.reads, args.threads, args.output_dir, args.sample, preset=args.preset)

    clair3_vcf, _consensus = run_clair3(args.output_dir, fasta, bam, args.threads, args.model_name, args.sample, args.model_path)
    raw_vcf = out / f"{args.sample}.clair3.raw.vcf.gz"
    shutil.copy2(clair3_vcf, raw_vcf)
    run_cmd(["bcftools", "index", "-f", str(raw_vcf)])

    if args.vcf_mode == "clair3_raw":
        return
    if args.vcf_mode == "legacy_merge":
        # keep compatibility: if merge_vcfs output exists, preserve naming contract
        legacy = out / f"{args.sample}.dviont.legacy.vcf.gz"
        shutil.copy2(clair3_vcf, legacy)
        run_cmd(["bcftools", "index", "-f", str(legacy)])
        return

    clean_out = out / f"{args.sample}.dviont.clean.vcf.gz"
    processor = VCFProcessor(
        vcf_file=str(raw_vcf),
        output_dir=args.output_dir,
        sample=args.sample,
        reference_fasta=fasta,
        genbank_file=args.ref if ref_fmt == "genbank" else None,
        clean_config=CleanConfig(
            min_af_review=args.min_af_review,
            min_dp_review=args.min_dp_review,
            min_qual_review=args.min_qual_review,
            dense_window_bp=args.dense_window_bp,
            dense_min_variants=args.dense_min_variants,
        ),
    )
    processor.process_clean(str(clean_out))


def run_cohort(args):
    for exe in ["bcftools", "snp-dists"]:
        check_executable(exe)
    if args.recombination == "gubbins":
        check_executable("run_gubbins.py")

    out = Path(args.out)
    calls_dir = out / "calls"
    cohort_vcfs = out / "cohort_vcfs"
    alns = out / "alignments"
    dist = out / "distances"
    for p in [calls_dir, cohort_vcfs, alns / "consensus_snps", dist]:
        p.mkdir(parents=True, exist_ok=True)

    samples = parse_reads_list(args.reads_list)
    vcfs = []
    for sample, reads in samples:
        sample_out = calls_dir / sample
        call_args = argparse.Namespace(**vars(args), output_dir=str(sample_out), reads=reads, sample=sample)
        run_call(call_args)
        vcfs.append(resolve_sample_vcf(str(sample_out), sample, args.cohort_vcf_source))

    for vcf in vcfs:
        run_cmd(["bcftools", "index", "-f", vcf])
    merged = cohort_vcfs / "cohort_merged.vcf.gz"
    run_cmd(["bcftools", "merge", "-m", "none", "-Oz", "-o", str(merged), *vcfs])
    run_cmd(["bcftools", "index", "-f", str(merged)])
    norm = cohort_vcfs / "cohort_merged.norm.vcf.gz"
    run_cmd(["bcftools", "norm", "-m", "-any", str(merged), "-Oz", "-o", str(norm)])
    run_cmd(["bcftools", "index", "-f", str(norm)])
    snps = cohort_vcfs / "cohort_merged.snps.vcf.gz"
    run_cmd(["bcftools", "view", "-i", 'strlen(REF)=1 && strlen(ALT)=1 && ALT!="*"', str(norm), "-Oz", "-o", str(snps)])
    run_cmd(["bcftools", "index", "-f", str(snps)])

    seq_paths = []
    for sample, _ in samples:
        fa = alns / "consensus_snps" / f"{sample}.fasta"
        with open(fa, "w") as h:
            subprocess.run(["bcftools", "consensus", "-f", args.ref, "-s", sample, str(snps)], check=True, stdout=h)
        txt = fa.read_text().splitlines()
        txt[0] = f">{sample}"
        fa.write_text("\n".join(txt) + "\n")
        seq_paths.append(fa)

    alignment = alns / "cohort.snp_alignment.fasta"
    with open(alignment, "w") as out_h:
        lengths = set()
        for fa in seq_paths:
            content = fa.read_text()
            out_h.write(content)
            seq = "".join([x for x in content.splitlines() if not x.startswith(">")])
            lengths.add(len(seq))
        if len(lengths) != 1:
            raise RuntimeError("Consensus SNP FASTA sequence lengths differ across samples.")

    with open(dist / "cohort.unmasked_snp_distance_matrix.tsv", "w") as h:
        subprocess.run(["snp-dists", str(alignment)], check=True, stdout=h)


def build_parser():
    parser = argparse.ArgumentParser(description="DviONT")
    sub = parser.add_subparsers(dest="command")
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-r", "--ref", required=True)
    common.add_argument("-t", "--threads", type=int, default=2)
    common.add_argument("-m", "--model-name", default="r1041_e82_400bps_sup_v430_bacteria_finetuned")
    common.add_argument("-p", "--model-path", default=None)
    common.add_argument("--aligner", choices=["minimap2", "winnowmap"], default="minimap2")
    common.add_argument("--preset", default="ont-q20")
    common.add_argument("--vcf-mode", choices=["clair3_raw", "clean", "legacy_merge"], default="clean")
    common.add_argument("--min-af-review", type=float, default=0.8)
    common.add_argument("--min-dp-review", type=int, default=10)
    common.add_argument("--min-qual-review", type=float, default=15)
    common.add_argument("--dense-window-bp", type=int, default=100)
    common.add_argument("--dense-min-variants", type=int, default=5)

    call = sub.add_parser("call", parents=[common])
    call.add_argument("-o", "--output-dir", required=True)
    call.add_argument("-i", "--reads", required=True)
    call.add_argument("-s", "--sample", required=True)

    cohort = sub.add_parser("cohort", parents=[common])
    cohort.add_argument("--reads-list", required=True)
    cohort.add_argument("--out", required=True)
    cohort.add_argument("--recombination", choices=["none", "gubbins"], default="none")
    cohort.add_argument("--cohort-vcf-source", choices=["clean", "clair3_raw", "legacy_merge"], default="clean")

    # backward-compatible flat mode
    parser.add_argument("-o", "--output_dir")
    parser.add_argument("-i", "--reads")
    parser.add_argument("-s", "--sample", default="SAMPLE")
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "call":
        run_call(args)
    elif args.command == "cohort":
        run_cohort(args)
    else:
        if args.output_dir and args.reads and args.ref:
            run_call(args)
        else:
            parser.print_help()


if __name__ == "__main__":
    main()
