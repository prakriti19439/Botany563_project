#!/bin/bash

# ----------------------------------------------------------------------
# 12_run_muscle.sh
#
# Purpose:
#   Run MUSCLE v5 multiple sequence alignment for each Louvain cluster FASTA.
#
# Input:
#   results/cluster_fastas/res${RES}_seed${SEED}/cluster_*.fa
#
# Output:
#   results/muscle_alignments/res${RES}_seed${SEED}/cluster_*.muscle.aln.fa
#
# Notes:
#   - MUSCLE is used as an alternative/complement to MAFFT.
#   - MUSCLE can be slower than MAFFT but may improve alignment quality
#     for difficult clusters.
# ----------------------------------------------------------------------

set -e

# -----------------------
# Parameters
# -----------------------

BASE_DIR=$(pwd)

CLUSTER_FASTA_DIR=${BASE_DIR}/results/cluster_fastas
OUT_BASE=${BASE_DIR}/results/muscle_alignments

# Louvain settings tested in the project
RESOLUTIONS="100 1000"
SEEDS="10 100"

# MUSCLE executable
MUSCLE=${MUSCLE:-muscle}

mkdir -p ${OUT_BASE}

# -----------------------
# Run MUSCLE
# -----------------------

for RES in ${RESOLUTIONS}
do
    for SEED in ${SEEDS}
    do
        INDIR=${CLUSTER_FASTA_DIR}/res${RES}_seed${SEED}
        OUTDIR=${OUT_BASE}/res${RES}_seed${SEED}

        mkdir -p ${OUTDIR}

        if [ ! -d ${INDIR} ]; then
            echo "WARNING: input directory not found: ${INDIR}"
            continue
        fi

        echo "Running MUSCLE for res=${RES}, seed=${SEED}"

        for FA in ${INDIR}/cluster_*.fa
        do
            if [ ! -e ${FA} ]; then
                echo "No cluster FASTA files found in ${INDIR}"
                continue
            fi

            BASE=$(basename ${FA} .fa)
            OUT=${OUTDIR}/${BASE}.muscle.aln.fa

            if [ -s ${OUT} ]; then
                echo "Skipping existing alignment: ${OUT}"
                continue
            fi

            echo "Aligning ${BASE}"

            ${MUSCLE} \
                -align ${FA} \
                -output ${OUT}
        done
    done
done

echo "Done running MUSCLE alignments."
