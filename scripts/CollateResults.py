#!/usr/bin/env python3
"""Collate reciprocal LiftOver-overlap mappings for one species pair."""

import argparse
import os


def read_forward(path, s1, s2):
    pairs = set()
    if not os.path.exists(path):
        return pairs
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            n1 = parts[3] + "|" + s1
            n2 = parts[7] + "|" + s2
            pairs.add((n1, n2))
    return pairs


def read_reverse(path, s1, s2):
    pairs = set()
    if not os.path.exists(path):
        return pairs
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            n1 = parts[7] + "|" + s1
            n2 = parts[3] + "|" + s2
            pairs.add((n1, n2))
    return pairs


def main():
    parser = argparse.ArgumentParser(description="Collate reciprocal edges from two overlap files.")
    parser.add_argument("mapFile1")
    parser.add_argument("mapFile2")
    parser.add_argument("S1")
    parser.add_argument("S2")
    parser.add_argument("--allow_one_direction", action="store_true", help="If one file is missing, output available one-direction edges.")
    args = parser.parse_args()

    set1 = read_forward(args.mapFile1, args.S1, args.S2)
    set2 = read_reverse(args.mapFile2, args.S1, args.S2)

    if args.allow_one_direction and (not set1 or not set2):
        out = set1 if set1 else set2
    else:
        out = set1.intersection(set2)

    for n1, n2 in sorted(out):
        print(f"{n1}\t{n2}")


if __name__ == "__main__":
    main()
