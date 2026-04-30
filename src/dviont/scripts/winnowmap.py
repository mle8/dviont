import subprocess
import os
import logging
import shutil

def _ensure_repeats(fasta_file, repeats_path, threads):
    """
    Generate a weighted k-mer list (repeat_k15.txt) with meryl if it doesn't exist.
    """
    if os.path.exists(repeats_path):
        return repeats_path
    if not shutil.which("meryl"):
        raise RuntimeError("❌ meryl not found. Install it or run via container.")
    db = os.path.join(os.path.dirname(repeats_path), "ref_k15.meryl")

    logging.info(f"🔹 Generating weighted minimizer list at {repeats_path} using meryl...")
    subprocess.run(["meryl", "count", "k=15", "output", db, fasta_file], check=True)
    with open(repeats_path, "w") as fh:
        subprocess.run(["meryl", "print", "greater-than", "distinct=0.9998", db],
                       check=True, stdout=fh)
    return repeats_path


def run_winnowmap_alignment(fasta_file, reads, threads, output_dir, sample,
                             preset="ont-q20", repeats="repeat_k15.txt"):
    """
    Run Winnowmap2 alignment and sort the output BAM file with sample name.

    Args:
        fasta_file (str): Path to reference FASTA.
        reads (str): Path to FASTQ reads.
        threads (int): Number of threads.
        output_dir (str): Directory where BAM will be saved.
        sample (str): Sample name for BAM prefix.
        preset (str): One of ['ont-legacy','ont-q20','pb-clr','pb-hifi','asm'].
        repeats (str): Path (or filename) for repeat k-mer list (default: repeat_k15.txt).

    Returns:
        str: Path to the sorted BAM file if successful, else None.
    """
    bam_output = os.path.join(output_dir, f"{sample}_aln_sort.bam")

    # Map dviont presets to Winnowmap -x
    if preset in ("ont-legacy", "ont-q20"):
        x_value = "map-ont"
    elif preset in ("pb-clr", "pb-hifi", "asm"):
        x_value = "map-pb"
    else:
        x_value = "map-ont"

    # Prepare repeats file
    repeats_path = repeats if os.path.isabs(repeats) else os.path.join(output_dir, repeats)
    repeats_path = _ensure_repeats(fasta_file, repeats_path, threads)

    winnowmap_cmd = [
        "winnowmap", "-W", repeats_path, "-t", str(threads), "-ax", x_value, fasta_file, reads
    ]
    samtools_sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", bam_output]

    try:
        logging.info(f"🔹 Running Winnowmap with command: {' '.join(winnowmap_cmd)}")

        with subprocess.Popen(winnowmap_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as wmap_proc, \
             subprocess.Popen(samtools_sort_cmd, stdin=wmap_proc.stdout, stderr=subprocess.PIPE, text=True) as sort_proc:

            for line in wmap_proc.stderr:
                logging.info(f"Winnowmap2: {line.strip()}")

            wmap_proc.stdout.close()
            logging.info(f"🔹 Running samtools sort with command: {' '.join(samtools_sort_cmd)}")

            wmap_return = wmap_proc.wait()
            sort_stderr = sort_proc.communicate()[1]

            if wmap_return != 0:
                logging.error(f"❌ Winnowmap failed with error code {wmap_return}.")
                return None
            if sort_proc.returncode != 0:
                logging.error(f"❌ Samtools sorting failed: {sort_stderr}")
                return None

        subprocess.run(["samtools", "index", bam_output], check=True)
        logging.info(f"✅ Winnowmap alignment and Samtools sort complete: {bam_output}")
        return bam_output

    except subprocess.CalledProcessError as e:
        logging.error(f"❌ Error running Winnowmap2 or Samtools: {getattr(e, 'stderr', e)}")
        return None
    except Exception as e:
        logging.error(f"❌ Unexpected error: {e}")
        return None
