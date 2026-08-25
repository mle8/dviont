import logging
import os
import shutil
import subprocess


def _ensure_repeats(fasta_file, repeats_path):
    """Create Winnowmap's weighted repetitive k-mer list with meryl."""
    if os.path.exists(repeats_path):
        return repeats_path
    if shutil.which("meryl") is None:
        raise RuntimeError("meryl is required when --aligner winnowmap is selected")

    database = os.path.join(os.path.dirname(repeats_path), "ref_k15.meryl")
    logging.info("Generating Winnowmap weighted k-mers with meryl")
    subprocess.run(
        ["meryl", "count", "k=15", "output", database, fasta_file], check=True
    )
    with open(repeats_path, "w") as repeats:
        subprocess.run(
            ["meryl", "print", "greater-than", "distinct=0.9998", database],
            check=True,
            stdout=repeats,
        )
    return repeats_path


def run_winnowmap_alignment(
    fasta_file, reads, threads, output_dir, sample, preset="ont-q20"
):
    """Align reads with Winnowmap and return a sorted, indexed BAM."""
    if shutil.which("winnowmap") is None:
        raise RuntimeError("winnowmap is required when --aligner winnowmap is selected")

    bam_output = os.path.join(output_dir, f"{sample}_aln_sort.bam")
    repeats = _ensure_repeats(
        fasta_file, os.path.join(output_dir, "repeat_k15.txt")
    )
    x_value = "map-ont" if preset.startswith("ont-") else "map-pb"
    winnowmap_cmd = [
        "winnowmap", "-W", repeats, "-t", str(threads), "-ax", x_value,
        fasta_file, reads,
    ]
    sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", bam_output]
    logging.info("Running Winnowmap with command: %s", " ".join(winnowmap_cmd))

    with subprocess.Popen(winnowmap_cmd, stdout=subprocess.PIPE) as aligner:
        subprocess.run(sort_cmd, stdin=aligner.stdout, check=True)
        aligner.stdout.close()
        if aligner.wait() != 0:
            raise subprocess.CalledProcessError(aligner.returncode, winnowmap_cmd)

    subprocess.run(["samtools", "index", bam_output], check=True)
    return bam_output
