#!/usr/bin/env python3
"""Root and resolve unrooted region trees using OrthoFinder code."""

import os
import glob
import argparse
import multiprocessing as mp
from tqdm import tqdm
import sys
import csv
from argparse import Namespace
import traceback


def load_orthofinder(orthofinder_src):
    sys.path.append(orthofinder_src)
    from orthofinder.tools import tree as tree_lib
    from orthofinder.gene_tree_inference import trees2ologs_of as om1
    return tree_lib, om1


def root_and_resolve_tree(args):
    tree_fn, species_tree_fn, out_dir, orthofinder_src = args
    try:
        tree_lib, om1 = load_orthofinder(orthofinder_src)
        species_tree_rooted = tree_lib.Tree(species_tree_fn)
        tree = tree_lib.Tree(tree_fn)

        def gene_to_species(name):
            return name.split("|")[1] if "|" in name else name

        args_ns = Namespace(separator="|", gene_to_species_func=gene_to_species)
        GeneToSpecies = om1.GetGeneToSpeciesMap(args_ns)

        root = om1.GetRoot(tree, species_tree_rooted, GeneToSpecies)
        if root is None:
            return (os.path.basename(tree_fn), "no_root")

        if root != tree:
            tree.set_outgroup(root)

        resolved_tree = om1.Resolve(tree, GeneToSpecies)
        out_fn = os.path.join(out_dir, os.path.basename(tree_fn))
        resolved_tree.write(outfile=out_fn, format=1)
        return (os.path.basename(tree_fn), "ok")

    except Exception as e:
        print(f"\n[ERROR] in {tree_fn}:\n{traceback.format_exc()}")
        return (os.path.basename(tree_fn), f"error: {str(e)}")


def main():
    parser = argparse.ArgumentParser(description="Root and resolve unrooted region trees using a rooted species tree.")
    parser.add_argument("-t", "--trees_dir", required=True, help="Directory with unrooted tree files")
    parser.add_argument("-s", "--species_tree", required=True, help="Rooted species tree in Newick format")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory for rooted/resolved trees")
    parser.add_argument("-n", "--threads", type=int, default=4)
    parser.add_argument("--tree_suffix", default=".tree", help="Tree file suffix, e.g. .tree or .treefile")
    parser.add_argument("--orthofinder_src", default="scripts/OrthoFinder/src", help="Path to OrthoFinder/src")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    if not os.path.isdir(args.orthofinder_src):
        raise FileNotFoundError(f"Missing OrthoFinder source directory: {args.orthofinder_src}")

    tree_files = glob.glob(os.path.join(args.trees_dir, f"*{args.tree_suffix}"))
    if not tree_files:
        raise FileNotFoundError(f"No *{args.tree_suffix} files found in {args.trees_dir}")

    print(f"[+] Found {len(tree_files)} trees")
    arg_list = [(fn, args.species_tree, args.out_dir, args.orthofinder_src) for fn in tree_files]

    with mp.Pool(processes=args.threads) as pool:
        results = list(tqdm(pool.imap_unordered(root_and_resolve_tree, arg_list), total=len(tree_files), desc="Rooting & Resolving Trees"))

    summary_fn = os.path.join(args.out_dir, "tree_status.csv")
    with open(summary_fn, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["tree_file", "status"])
        writer.writerows(results)

    ok = sum(1 for _, status in results if status == "ok")
    no_root = sum(1 for _, status in results if status == "no_root")
    errors = [(fn, status) for fn, status in results if status.startswith("error")]
    print(f"Completed: {ok} trees")
    print(f"No root: {no_root}")
    print(f"Errors: {len(errors)}")
    print(f"Summary: {summary_fn}")


if __name__ == "__main__":
    main()
