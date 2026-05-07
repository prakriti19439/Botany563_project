#!/usr/bin/env python3
"""
07_graph_stats.py

Create LiftOver/overlap/reciprocal count summary tables from the region-level
LiftOver workflow outputs.

The original version had server-specific BED paths. This version uses relative
paths and command-line arguments.
"""

import argparse
import os
import subprocess
from pathlib import Path
import pandas as pd


def run_shell_command(command: str) -> str:
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        return "0"
    return result.stdout.strip()


def count_lines(path: Path) -> str:
    if not path.exists():
        return "0"
    return run_shell_command(f"wc -l < {path}")


def count_unique_col(path: Path, col: int) -> str:
    if not path.exists():
        return "0"
    return run_shell_command(f"awk -F'\\t' '{{print ${col}}}' {path} | sort -u | wc -l")


def main():
    parser = argparse.ArgumentParser(description="Summarize LiftOver, overlap, and reciprocal graph counts.")
    parser.add_argument("--base_dir", default=".", help="Project base directory")
    parser.add_argument("--bin_dir", default="data/bins", help="Directory containing *_1000bp.txt bins")
    parser.add_argument("--liftover_dir", default="results/liftover", help="LiftOver result directory")
    parser.add_argument("--out_file", default="results/liftover/results_stats/combined_all_stats.csv")
    parser.add_argument("--res", default="1kb")
    parser.add_argument("--thresholds", nargs="+", default=["t0.10", "t0.50"])
    args = parser.parse_args()

    base = Path(args.base_dir)
    bin_dir = base / args.bin_dir
    liftover_dir = base / args.liftover_dir
    stitch_dir = liftover_dir / "liftOver_stitching_outputs"
    overlap_dir = liftover_dir / "overlap_outputs"
    graph_dir = liftover_dir / "graph_outputs"

    rows = ["Human", "pig", "acomys_dimidiatus", "acomys_russatus", "macaca_fas", "rattus_norvegicus", "mouse"]
    columns = rows[:]

    species_map = {
        "Human": "homo_sapiens",
        "pig": "sus_scrofa",
        "mouse": "mus_musculus",
        "rattus_norvegicus": "rattus_norvegicus",
        "acomys_dimidiatus": "acomys_dimidiatus",
        "acomys_russatus": "acomys_russatus",
        "macaca_fas": "macaca_fas",
    }

    bin_name = {
        "Human": "homo_sapiens_1000bp.txt",
        "pig": "sus_scrofa_1000bp.txt",
        "mouse": "mus_musculus_1000bp.txt",
        "rattus_norvegicus": "rattus_norvegicus_1000bp.txt",
        "acomys_dimidiatus": "acomys_dimidiatus_1000bp.txt",
        "acomys_russatus": "acomys_russatus_1000bp.txt",
        "macaca_fas": "macaca_fas_1000bp.txt",
    }

    all_lines = []

    for thresh in args.thresholds:
        # LiftOver counts
        all_lines.append([f"Liftover counts ({args.res}, {thresh})"])
        all_lines.append([""] + columns)
        for row in rows:
            vals = []
            for col in columns:
                if row == col:
                    vals.append(count_lines(bin_dir / bin_name[row]))
                else:
                    fname = stitch_dir / f"{species_map[row]}_to_{species_map[col]}_{thresh}.txt"
                    vals.append(count_unique_col(fname, 4))
            all_lines.append([row] + vals)

        # Overlap counts
        all_lines.append([])
        all_lines.append([f"Overlap counts ({args.res}, {thresh})"])
        all_lines.append([""] + columns)
        for row in rows:
            vals = []
            for col in columns:
                if row == col:
                    vals.append(count_lines(bin_dir / bin_name[row]))
                else:
                    fname = overlap_dir / f"{row}_to_{col}_{thresh}_overlap.txt"
                    vals.append(count_unique_col(fname, 4))
            all_lines.append([row] + vals)

        # Reciprocal graph counts
        all_lines.append([])
        all_lines.append([f"Reciprocal counts ({args.res}, {thresh})"])
        all_lines.append([""] + columns)
        for row in rows:
            vals = []
            for col in columns:
                if row == col:
                    vals.append(count_lines(bin_dir / bin_name[row]))
                else:
                    fname = graph_dir / f"Graph_{row}_to_{col}_{thresh}.txt"
                    vals.append(count_unique_col(fname, 1))
            all_lines.append([row] + vals)
        all_lines.append([])

    out_file = base / args.out_file
    out_file.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(all_lines).to_csv(out_file, index=False, header=False)
    print(f"Wrote {out_file}")


if __name__ == "__main__":
    main()
