import os
import shutil
import logging
from Bio import SeqIO
import subprocess
import gzip


def _sanitize_fasta_in_place(fasta_path: str) -> None:
    """
    Rewrite FASTA to remove blank/whitespace-only lines and CRLF artifacts.
    Preserves headers (trimmed), removes all whitespace inside sequence lines,
    and uppercases sequence.

    This fixes cases where empty rows between records can break tools like
    bcftools consensus (variants not applied).
    """
    tmp_path = fasta_path + ".tmp"

    with open(fasta_path, "r", newline=None) as fin, open(tmp_path, "w", newline="\n") as fout:
        for raw in fin:
            line = raw.rstrip("\n").rstrip("\r")

            # Skip empty/whitespace-only lines (key fix)
            if not line.strip():
                continue

            if line.startswith(">"):
                # Preserve full header content, just trim ends
                fout.write(line.strip() + "\n")
            else:
                # Remove ANY whitespace inside sequence lines
                seq = "".join(line.split()).upper()
                if seq:
                    fout.write(seq + "\n")

    os.replace(tmp_path, fasta_path)


def extract_fasta_and_gbk(reference, ref_dir, ref_fmt, output_dir):
    """Extract FASTA and GenBank files, ensuring correct formatting (no blank lines in ref.fa)."""
    try:
        if not os.path.exists(reference):
            logging.error(f"Reference file not found: {reference}")
            return None

        os.makedirs(os.path.join(ref_dir, "genomes"), exist_ok=True)
        os.makedirs(os.path.join(ref_dir, "ref"), exist_ok=True)

        fasta_out = os.path.join(ref_dir, "ref.fa")  # same output as original

        if ref_fmt == "genbank":
            # Define paths
            genes_gbk_path = os.path.join(ref_dir, "ref", "genes.gbk")
            genes_gbk_gz_path = genes_gbk_path + ".gz"

            # Copy and gzip the GenBank file
            shutil.copy(reference, genes_gbk_path)
            with open(genes_gbk_path, "rb") as f_in, gzip.open(genes_gbk_gz_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

            # Remove uncompressed file after gzipping
            os.remove(genes_gbk_path)
            logging.info(f"Copied and gzipped GenBank file: {genes_gbk_gz_path}")

            # Convert GenBank to FASTA
            with open(fasta_out, "w") as fasta_file:
                records = SeqIO.parse(reference, "genbank")
                for record in records:
                    SeqIO.write(record, fasta_file, "fasta")

            logging.info(f"Converted GenBank to FASTA: {fasta_out}")

        elif ref_fmt == "fasta":
            # Copy FASTA file
            shutil.copy(reference, fasta_out)
            logging.info(f"Copied FASTA file to {fasta_out}")

        else:
            logging.error(f"Unknown ref_fmt: {ref_fmt}. Expected 'genbank' or 'fasta'.")
            return None

        # sanitize ref.fa in place (removes blank lines between records, etc.)
        try:
            _sanitize_fasta_in_place(fasta_out)
        except Exception as e:
            logging.error(f"❌ Error sanitizing FASTA {fasta_out}: {e}")
            return None

        # Ensure `fasta_out` is symlinked to `genomes/`
        genomes_fasta_out = os.path.join(ref_dir, "genomes", "ref.fa")
        if not os.path.exists(genomes_fasta_out):
            os.symlink(os.path.relpath(fasta_out, start=os.path.join(ref_dir, "genomes")), genomes_fasta_out)
            logging.info(f"Created relative symlink: {genomes_fasta_out} -> {fasta_out}")
        else:
            logging.info(f"Symlink already exists: {genomes_fasta_out}")

        # Run samtools faidx (remove stale indices first)
        for ext in (".fai", ".gzi"):
            p = fasta_out + ext
            if os.path.exists(p):
                try:
                    os.remove(p)
                except OSError:
                    pass

        try:
            subprocess.run(["samtools", "faidx", fasta_out], capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running samtools faidx: {e.stderr}")
            return None

        return fasta_out

    except Exception as e:
        logging.error(f"Error in extract_fasta_and_gbk: {e}")
        return None

