#!/usr/bin/env python3
"""Run MAFFT alignments on cluster FASTA files."""

import argparse
import os
import subprocess
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

MAFFT_CMD = {
    "default": ["mafft", "--localpair", "--maxiterate", "1000", "--anysymbol"],
    "fast": ["mafft", "--anysymbol"],
}


def run_mafft(input_fasta, output_fasta, seq_count, threads):
    if seq_count > 500:
        cmd = MAFFT_CMD["fast"] + ["--thread", str(threads), input_fasta]
    else:
        cmd = MAFFT_CMD["default"] + ["--thread", str(threads), input_fasta]
    with open(output_fasta, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, check=True)


def main():
    parser = argparse.ArgumentParser(description="Run MAFFT alignments on cluster FASTAs.")
    parser.add_argument("--res", required=True)
    parser.add_argument("--seed", required=True)
    parser.add_argument("--cluster_base", default="results/cluster_fastas")
    parser.add_argument("--out_base", default="results/mafft_alignments")
    parser.add_argument("--num_parallel_jobs", type=int, default=4)
    parser.add_argument("--cores_per_alignment", type=int, default=4)
    args = parser.parse_args()

    clusters_dir = os.path.join(args.cluster_base, f"res{args.res}_seed{args.seed}")
    alignments_dir = os.path.join(args.out_base, f"res{args.res}_seed{args.seed}")
    os.makedirs(alignments_dir, exist_ok=True)

    if not os.path.isdir(clusters_dir):
        raise FileNotFoundError(f"Missing cluster FASTA directory: {clusters_dir}")

    cluster_files = []
    for fname in os.listdir(clusters_dir):
        if fname.endswith(".fa"):
            cluster_path = os.path.join(clusters_dir, fname)
            seq_count = sum(1 for _ in SeqIO.parse(cluster_path, "fasta"))
            if seq_count < 2:
                continue
            base = os.path.splitext(fname)[0]
            out_file = os.path.join(alignments_dir, f"{base}.mafft.aln.fa")
            if not os.path.exists(out_file):
                cluster_files.append((cluster_path, out_file, seq_count))

    with ThreadPoolExecutor(max_workers=max(args.num_parallel_jobs, 1)) as executor:
        futures = [executor.submit(run_mafft, c[0], c[1], c[2], args.cores_per_alignment) for c in cluster_files]
        for f in tqdm(as_completed(futures), total=len(futures), desc="MAFFT Alignments"):
            f.result()

    print(f"Done. Alignments written to {alignments_dir}")


if __name__ == "__main__":
    main()
