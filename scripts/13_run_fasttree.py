#!/usr/bin/env python3
"""Run FastTree on MAFFT or MUSCLE alignments."""

import argparse
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm


def run_fasttree(input_fasta, output_tree, fasttree_exe="FastTree"):
    cmd = [fasttree_exe, "-nt", "-gtr", input_fasta]
    with open(output_tree, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, check=True)


def main():
    parser = argparse.ArgumentParser(description="Run FastTree GTR on nucleotide alignments.")
    parser.add_argument("--res", required=True)
    parser.add_argument("--seed", required=True)
    parser.add_argument("--alignment", choices=["mafft", "muscle"], default="mafft")
    parser.add_argument("--align_base", default=None, help="Base alignment directory. Default depends on --alignment.")
    parser.add_argument("--out_base", default="results/fasttree")
    parser.add_argument("--cores", type=int, default=8, help="Number of parallel FastTree jobs")
    parser.add_argument("--fasttree", default="FastTree")
    args = parser.parse_args()

    if args.align_base is None:
        args.align_base = "results/mafft_alignments" if args.alignment == "mafft" else "results/muscle_alignments"

    suffix = "mafft.aln.fa" if args.alignment == "mafft" else "muscle.aln.fa"
    alignments_dir = os.path.join(args.align_base, f"res{args.res}_seed{args.seed}")
    trees_dir = os.path.join(args.out_base, f"{args.alignment}_res{args.res}_seed{args.seed}")
    os.makedirs(trees_dir, exist_ok=True)

    if not os.path.isdir(alignments_dir):
        raise FileNotFoundError(f"Missing alignment directory: {alignments_dir}")

    jobs = []
    for fname in os.listdir(alignments_dir):
        if fname.endswith(suffix):
            aln = os.path.join(alignments_dir, fname)
            base = fname.replace(f".{suffix}", "")
            tree = os.path.join(trees_dir, f"{base}.tree")
            if not os.path.exists(tree):
                jobs.append((aln, tree))

    with ThreadPoolExecutor(max_workers=args.cores) as executor:
        futures = [executor.submit(run_fasttree, aln, tree, args.fasttree) for aln, tree in jobs]
        for f in tqdm(as_completed(futures), total=len(futures), desc="FastTree GTR"):
            f.result()

    print(f"Done. Trees written to {trees_dir}")


if __name__ == "__main__":
    main()
