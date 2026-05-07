#!/usr/bin/env python3
"""Generate one FASTA file per Louvain cluster."""

import argparse
import os
from collections import defaultdict
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description="Generate cluster FASTA files from Louvain assignments.")
    parser.add_argument("--res", required=True, help="Louvain resolution")
    parser.add_argument("--seed", required=True, help="Louvain seed")
    parser.add_argument("--species_fastas_dir", default="results/region_level_orthology/species_1kb_fastas")
    parser.add_argument("--cluster_dir", default="results/region_level_orthology/louvain_clusters")
    parser.add_argument("--out_base", default="results/cluster_fastas")
    args = parser.parse_args()

    cluster_file = os.path.join(args.cluster_dir, f"louvain_clusters_res{args.res}_seed{args.seed}.txt")
    out_dir = os.path.join(args.out_base, f"res{args.res}_seed{args.seed}")
    os.makedirs(out_dir, exist_ok=True)

    if not os.path.exists(cluster_file):
        raise FileNotFoundError(f"Missing cluster file: {cluster_file}")

    region_to_cluster = {}
    with open(cluster_file) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            region_to_cluster[parts[0]] = parts[1]

    clusters = defaultdict(list)
    for fasta_file in os.listdir(args.species_fastas_dir):
        if not fasta_file.endswith((".fa", ".fasta")):
            continue
        fasta_path = os.path.join(args.species_fastas_dir, fasta_file)
        for rec in SeqIO.parse(fasta_path, "fasta"):
            if rec.id in region_to_cluster:
                clusters[region_to_cluster[rec.id]].append(rec)

    for clust, records in clusters.items():
        out_file = os.path.join(out_dir, f"cluster_{clust}.fa")
        SeqIO.write(records, out_file, "fasta")

    print(f"Wrote {len(clusters)} cluster FASTA files to {out_dir}")


if __name__ == "__main__":
    main()
