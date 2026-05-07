#!/bin/bash

# ----------------------------------------------------------------------
# 05_run_liftover.sh
#
# Purpose:
#   Run the LiftOver-based region mapping workflow for 1000 bp genomic bins.
#
# Method:
#   For every ordered species pair S1 -> S2:
#     Step 1: Map 1000 bp regions from S1 to S2 using liftOver.
#     Step 2: Intersect mapped S2 regions with the universe of 1000 bp bins
#             in S2 using bedtools intersect.
#     Step 3: Collate reciprocal mappings between species pairs to create
#             graph edges for the region orthology graph.
#
# Main outputs:
#   results/liftover/liftOver_stitching_outputs/
#   results/liftover/overlap_outputs/
#   results/liftover/graph_outputs/
#   results/liftover/results_stats/
#
# Notes:
#   - The main analysis uses 1000 bp bins.
#   - This wrapper assumes the three project scripts below exist:
#       scripts/run_liftover/RunStitchingLiftOver.sh
#       scripts/run_liftover/OverlapLiftOverResults.sh
#       scripts/run_liftover/CompileStitchedGraphs.sh
#   - These scripts were adapted from the original gene-region LiftOver
#     workflow and applied here to uniform 1000 bp genomic bins.
# ----------------------------------------------------------------------

set -e

# -----------------------
# Paths
# -----------------------

BASE_DIR=$(pwd)

BIN_DIR=${BASE_DIR}/data/bins
CHAIN_DIR=${BASE_DIR}/results/chain_files_res/chain_files

LIFTOVER_DIR=${BASE_DIR}/results/liftover
STITCH_OUT=${LIFTOVER_DIR}/liftOver_stitching_outputs
OVERLAP_OUT=${LIFTOVER_DIR}/overlap_outputs
GRAPH_OUT=${LIFTOVER_DIR}/graph_outputs
STATS_OUT=${LIFTOVER_DIR}/results_stats

SCRIPT_DIR=${BASE_DIR}/scripts/run_liftover

mkdir -p ${STITCH_OUT}
mkdir -p ${OVERLAP_OUT}
mkdir -p ${GRAPH_OUT}
mkdir -p ${STATS_OUT}

# -----------------------
# Parameters
# -----------------------

# Main bin size used in the project
BIN_SIZE=1000

# LiftOver minMatch thresholds tested
THRESHOLDS="0.10 0.50"

# Species names should match bin names and chain-file naming convention.
# Expected bin files:
#   data/bins/${species}_1000bp.txt
SPECIES_LIST="acomys_dimidiatus acomys_russatus homo_sapiens mus_musculus sus_scrofa rattus_norvegicus macaca_fas"

# -----------------------
# Input checks
# -----------------------

for script in RunStitchingLiftOver.sh OverlapLiftOverResults.sh CompileStitchedGraphs.sh
do
    if [ ! -e ${SCRIPT_DIR}/${script} ]; then
        echo "ERROR: Missing script ${SCRIPT_DIR}/${script}"
        exit 1
    fi
done

for SPC in ${SPECIES_LIST}
do
    if [ ! -e ${BIN_DIR}/${SPC}_${BIN_SIZE}bp.txt ]; then
        echo "ERROR: Missing bin file: ${BIN_DIR}/${SPC}_${BIN_SIZE}bp.txt"
        exit 1
    fi
done

# ----------------------------------------------------------------------
# Step 1: Run LiftOver mapping
# ----------------------------------------------------------------------
# This step maps 1000 bp bins from one species to another.
#
# Example output:
#   results/liftover/liftOver_stitching_outputs/
#       bos_taurus_to_canis_lupus_familiaris_t0.10.txt
#
# Output format:
#   column 1: chromosome number in S2 genome
#   column 2: start of mapped region in S2 genome
#   column 3: end of mapped region in S2 genome
#   column 4: input region from S1 genome
#   column 5: mapped region number based on multiplicity
#
# Example:
#   16  6825195 6826107 4_106025048_106025953 1
# ----------------------------------------------------------------------

echo "Step 1: Running LiftOver mapping for ${BIN_SIZE} bp bins"

bash ${SCRIPT_DIR}/RunStitchingLiftOver.sh \
    ${BIN_DIR} \
    ${CHAIN_DIR} \
    ${STITCH_OUT} \
    ${BIN_SIZE} \
    "${THRESHOLDS}" \
    "${SPECIES_LIST}"

# ----------------------------------------------------------------------
# Step 2: Overlap mapped regions with target-species uniform bins
# ----------------------------------------------------------------------
# This step uses bedtools intersect to overlap mapped regions in S2 with
# the 1000 bp bin universe for S2.
#
# Example output:
#   results/liftover/overlap_outputs/
#       cow_to_dog_t0.10_overlap.txt
#
# Output format:
#   column 1: chromosome number of mapped region in S2 genome
#   column 2: start of mapped region in S2 genome
#   column 3: end of mapped region in S2 genome
#   column 4: input region id from S1 genome
#   column 5: chromosome of overlapped 1000 bp region in S2 genome
#   column 6: start of overlapped region in S2 genome
#   column 7: end of overlapped region in S2 genome
#   column 8: overlapped region id in S2 genome
#   column 9: number of base pairs overlapped
#
# Example:
#   16 6825195 6826107 4_106025048_106025953 16 6825707 6826526 16_6825707_6826526 400
# ----------------------------------------------------------------------

echo "Step 2: Overlapping lifted regions with target species ${BIN_SIZE} bp bins"

bash ${SCRIPT_DIR}/OverlapLiftOverResults.sh \
    ${STITCH_OUT} \
    ${BIN_DIR} \
    ${OVERLAP_OUT} \
    ${BIN_SIZE} \
    "${THRESHOLDS}" \
    "${SPECIES_LIST}"

# ----------------------------------------------------------------------
# Step 3: Collate reciprocal mappings into graph edges
# ----------------------------------------------------------------------
# This step extracts reciprocal region mappings from the overlap outputs.
#
# For each species pair:
#   - Take column 4 and column 8 from S1 -> S2 overlap output.
#   - Take column 4 and column 8 from S2 -> S1 overlap output, reversing
#     the order for comparison.
#   - Keep region pairs supported in both directions.
#
# Output graph format:
#   region_from_S1|species1    region_from_S2|species2
#
# Example:
#   4_106025048_106025953|bos_taurus    16_6825707_6826526|canis_lupus_familiaris
#
# Note:
#   For each species pair, there may be two graph files containing the same
#   reciprocal edges but in different ordering.
# ----------------------------------------------------------------------

echo "Step 3: Compiling reciprocal LiftOver graph edges"

bash ${SCRIPT_DIR}/CompileStitchedGraphs.sh \
    ${OVERLAP_OUT} \
    ${GRAPH_OUT} \
    ${STATS_OUT} \
    ${BIN_SIZE} \
    "${THRESHOLDS}" \
    "${SPECIES_LIST}"

echo "Done running LiftOver workflow."
echo
echo "Main outputs:"
echo "  LiftOver mapped regions: ${STITCH_OUT}"
echo "  Overlap outputs:         ${OVERLAP_OUT}"
echo "  Graph outputs:           ${GRAPH_OUT}"
echo "  Summary statistics:      ${STATS_OUT}"
