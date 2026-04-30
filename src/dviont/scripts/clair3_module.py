import logging
import os
import shutil
import subprocess


def run_clair3(output_dir, ref, bam_output, threads=2, model_name="r1041_e82_400bps_sup_v430_bacteria_finetuned", sample="", model_path=None):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)
    default_model_path = os.path.join(base_dir, "models", model_name)
    parallel_path = os.path.join(base_dir, "binaries", "parallel")
    model_path = model_path or default_model_path

    output_dir = os.path.abspath(output_dir)
    clair3_output_dir = os.path.join(output_dir, "clair3")
    os.makedirs(clair3_output_dir, exist_ok=True)

    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Clair3 model path does not exist: {model_path}")

    clair3_cmd = [
        "run_clair3.sh",
        f"--bam_fn={os.path.abspath(bam_output)}",
        f"--ref_fn={os.path.abspath(ref)}",
        f"--threads={threads}",
        "--platform=ont",
        f"--parallel={parallel_path}",
        f"--model_path={os.path.abspath(model_path)}",
        f"--output={clair3_output_dir}",
        f"--sample_name={sample}",
        "--include_all_ctgs",
        "--haploid_precise",
        "--no_phasing_for_fa",
        "--enable_long_indel",
    ]
    logging.info("Running Clair3: %s", " ".join(clair3_cmd))
    subprocess.run(clair3_cmd, check=True)

    merge_vcf = os.path.join(clair3_output_dir, "merge_output.vcf.gz")
    if not os.path.exists(merge_vcf):
        raise FileNotFoundError(f"Expected Clair3 output missing: {merge_vcf}")
    if not (os.path.exists(f"{merge_vcf}.tbi") or os.path.exists(f"{merge_vcf}.csi")):
        subprocess.run(["bcftools", "index", "-f", merge_vcf], check=True)

    tmp_dir = os.path.join(clair3_output_dir, "tmp")
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    return merge_vcf, clair3_output_dir
