#!/usr/bin/env python3

"""
17_validate_trees_with_liftover.py

Validate rooted region trees using LiftOver-supported region pairs.

Purpose
-------
For each rooted region tree, this script compares tree-derived cross-species
region pairs against LiftOver-supported graph edges.

This is useful for checking whether regions grouped together in the inferred
tree are also supported by LiftOver-based mappings.

Input
-----
1. One rooted/resolved tree file OR a directory of tree files.
2. A LiftOver graph edge file.

Expected tree format
--------------------
Newick tree with leaf names formatted like:

    region_id|species

Example:

    1_1000_2000|homo_sapiens

Expected LiftOver graph format
------------------------------
Tab-separated file with two columns:

    regionA|speciesA    regionB|speciesB

Example:

    4_106025048_106025953|bos_taurus    16_6825707_6826526|canis_lupus_familiaris

Output
------
A summary TSV file with:

    tree_file
    num_tree_pairs
    num_liftover_pairs_in_tree
    TP
    FP
    FN
    precision
    recall
    f1
    liftover_pairs_present_together
    fraction_liftover_pairs_present_together

Example
-------
python scripts/17_validate_trees_with_liftover.py \
    --tree_file results/rooted_trees/fasttree_mafft_res1000_seed10/cluster_92502.tree \
    --liftover_graph results/liftover/graph_outputs/Graph_AllSpecies_t0.5.txt \
    --out_file results/validation/cluster_92502_validation.tsv

Batch mode:

python scripts/17_validate_trees_with_liftover.py \
    --tree_dir results/rooted_trees/fasttree_mafft_res1000_seed10 \
    --tree_suffix ".tree" \
    --liftover_graph results/liftover/graph_outputs/Graph_AllSpecies_t0.5.txt \
    --out_file results/validation/all_tree_validation.tsv
"""

import argparse
import os
import glob
from itertools import combinations

import pandas as pd
import newick


def get_all_leaf_names(node):
    """
    Recursively collect leaf names under a tree node.

    Internal nodes without names are ignored.
    """
    names = set()

    if not node.descendants:
        if node.name and node.name != "NoName":
            names.add(node.name)
        return names

    for child in node.descendants:
        names.update(get_all_leaf_names(child))

    return names


def get_internal_leaf_sets(tree):
    """
    For every internal node, collect the set of leaf names below it.

    These sets are used to generate region pairs implied by the tree.
    """
    leaf_sets = []

    def helper(node):
        if not node.descendants:
            if node.name and node.name != "NoName":
                return {node.name}
            return set()

        leaves = set()
        for child in node.descendants:
            leaves.update(helper(child))

        if len(leaves) > 1:
            leaf_sets.append(leaves)

        return leaves

    helper(tree)
    return leaf_sets


def species_from_name(name):
    """
    Extract species from a region name formatted as region|species.
    """
    if "|" not in name:
        return None
    return name.split("|")[-1]


def get_cross_species_pairs_from_tree(tree):
    """
    Generate all cross-species region pairs implied by internal nodes
    in the tree.

    Only pairs where the two leaves come from different species are kept.
    """
    node_sets = get_internal_leaf_sets(tree)

    all_pairs = set()

    for ns in node_sets:
        nodes_species = []

        for node in ns:
            sp = species_from_name(node)
            if sp is not None:
                nodes_species.append((node, sp))

        for (n1, s1), (n2, s2) in combinations(nodes_species, 2):
            if s1 != s2:
                all_pairs.add(tuple(sorted([n1, n2])))

    all_names = set()
    for ns in node_sets:
        all_names.update(ns)

    return all_pairs, node_sets, all_names


def read_tree(tree_file):
    """
    Read a Newick tree file using the newick Python package.
    """
    with open(tree_file) as f:
        tree_str = f.read().strip()

    trees = newick.loads(tree_str)
    if len(trees) == 0:
        raise ValueError(f"No trees found in {tree_file}")

    return trees[0]


def read_liftover_pairs(liftover_graph):
    """
    Read LiftOver graph file and return a dataframe with columns A and B.
    """
    df = pd.read_csv(
        liftover_graph,
        sep="\t",
        header=None,
        usecols=[0, 1],
        names=["A", "B"]
    )

    df = df.dropna()
    return df


def validate_one_tree(tree_file, lift_df):
    """
    Validate one tree against LiftOver graph edges.
    """
    tree = read_tree(tree_file)

    tree_pairs, node_sets, all_names = get_cross_species_pairs_from_tree(tree)

    # Restrict LiftOver pairs to only nodes present in this tree.
    lift_sub = lift_df[
        lift_df["A"].isin(all_names) & lift_df["B"].isin(all_names)
    ].copy()

    liftover_pairs = set(
        tuple(sorted(pair))
        for pair in lift_sub[["A", "B"]].itertuples(index=False, name=None)
    )

    tp = tree_pairs & liftover_pairs
    fp = tree_pairs - liftover_pairs
    fn = liftover_pairs - tree_pairs

    precision = len(tp) / (len(tp) + len(fp)) if (len(tp) + len(fp)) > 0 else 0.0
    recall = len(tp) / (len(tp) + len(fn)) if (len(tp) + len(fn)) > 0 else 0.0
    f1 = (
        2 * precision * recall / (precision + recall)
        if (precision + recall) > 0
        else 0.0
    )

    # Additional check:
    # How many LiftOver pairs are present together under the same internal node?
    liftover_pairs_present_together = set()

    for pair in liftover_pairs:
        for leaf_set in node_sets:
            if pair[0] in leaf_set and pair[1] in leaf_set:
                liftover_pairs_present_together.add(pair)
                break

    fraction_present_together = (
        len(liftover_pairs_present_together) / len(liftover_pairs)
        if len(liftover_pairs) > 0
        else 0.0
    )

    return {
        "tree_file": os.path.basename(tree_file),
        "tree_path": tree_file,
        "num_tree_pairs": len(tree_pairs),
        "num_liftover_pairs_in_tree": len(liftover_pairs),
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "liftover_pairs_present_together": len(liftover_pairs_present_together),
        "fraction_liftover_pairs_present_together": fraction_present_together,
    }


def collect_tree_files(args):
    """
    Collect one tree file or all tree files from a directory.
    """
    if args.tree_file:
        return [args.tree_file]

    if args.tree_dir:
        pattern = os.path.join(args.tree_dir, f"*{args.tree_suffix}")
        tree_files = sorted(glob.glob(pattern))
        if len(tree_files) == 0:
            raise FileNotFoundError(f"No tree files found with pattern: {pattern}")
        return tree_files

    raise ValueError("Provide either --tree_file or --tree_dir")


def main():
    parser = argparse.ArgumentParser(
        description="Validate rooted region trees using LiftOver-supported region pairs."
    )

    parser.add_argument(
        "--tree_file",
        default=None,
        help="Single rooted tree file in Newick format."
    )

    parser.add_argument(
        "--tree_dir",
        default=None,
        help="Directory containing rooted tree files."
    )

    parser.add_argument(
        "--tree_suffix",
        default=".tree",
        help="Tree file suffix for batch mode. Default: .tree"
    )

    parser.add_argument(
        "--liftover_graph",
        required=True,
        help="Two-column LiftOver graph edge file."
    )

    parser.add_argument(
        "--out_file",
        required=True,
        help="Output summary TSV file."
    )

    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out_file), exist_ok=True)

    print(f"Reading LiftOver graph: {args.liftover_graph}")
    lift_df = read_liftover_pairs(args.liftover_graph)
    print(f"Loaded {len(lift_df)} LiftOver graph edges")

    tree_files = collect_tree_files(args)
    print(f"Validating {len(tree_files)} tree(s)")

    rows = []
    for tree_file in tree_files:
        print(f"Validating: {tree_file}")
        try:
            rows.append(validate_one_tree(tree_file, lift_df))
        except Exception as e:
            rows.append({
                "tree_file": os.path.basename(tree_file),
                "tree_path": tree_file,
                "num_tree_pairs": 0,
                "num_liftover_pairs_in_tree": 0,
                "TP": 0,
                "FP": 0,
                "FN": 0,
                "precision": 0.0,
                "recall": 0.0,
                "f1": 0.0,
                "liftover_pairs_present_together": 0,
                "fraction_liftover_pairs_present_together": 0.0,
                "error": str(e),
            })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.out_file, sep="\t", index=False)

    print(f"Saved validation results to: {args.out_file}")


if __name__ == "__main__":
    main()
