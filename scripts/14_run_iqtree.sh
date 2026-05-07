#!/bin/bash

# ----------------------------------------------------------------------
# 14_run_iqtree.sh
#
# Purpose:
#   Run IQ-TREE on cluster multiple sequence alignments to infer
#   maximum-likelihood phylogenetic trees and estimate bootstrap support.
#
# Input:
#   results/mafft_alignments/res${RES}_seed${SEED}/cluster_*.mafft.aln.fa
#   or
#   results/muscle_alignments/res${RES}_seed${SEED}/cluster_*.muscle.aln.fa
#
# Output:
#   results/iqtree/${ALIGNMENT}_res${RES}_seed${SEED}/cluster_*.treefile
#   results/iqtree/${ALIGNMENT}_res${RES}_seed${SEED}/cluster_*.iqtree
#   results/iqtree/${ALIGNMENT}_res${RES}_seed${SEED}/cluster_*.log
#
# Notes:
#   - IQ-TREE provides maximum-likelihood trees.
#   - ModelFinder is used with -m MFP.
#   - Ultrafast bootstrap support is estimated with -B.
# ----------------------------------------------------------------------

set -e

# -----------------------
# Parameters
# -----------------------

BASE_DIR=$(pwd)

# Choose alignment type: mafft or muscle
ALIGNMENT=${ALIGNMENT:-mafft}

MAFFT_DIR=${BASE_DIR}/results/mafft_alignments
MUSCLE_DIR=${BASE_DIR}/results/muscle_alignments
OUT_BASE=${BASE_DIR}/results/iqtree

RESOLUTIONS="100 1000"
SEEDS="10 100"

# IQ-TREE executable
IQTREE=${IQTREE:-iqtree2}

# Number of ultrafast bootstrap replicates
BOOTSTRAPS=${BOOTSTRAPS:-1000}

# Number of threads
THREADS=${THREADS:-AUTO}

mkdir -p ${OUT_BASE}

# -----------------------
# Select input directory and suffix
# -----------------------

if [ "${ALIGNMENT}" == "mafft" ]; then
    IN_BASE=${MAFFT_DIR}
    SUFFIX="mafft.aln.fa"
elif [ "${ALIGNMENT}" == "muscle" ]; then
    IN_BASE=${MUSCLE_DIR}
    SUFFIX="muscle.aln.fa"
else
    echo "ERROR: ALIGNMENT must be either 'mafft' or 'muscle'"
    exit 1
fi

# -----------------------
# Run IQ-TREE
# -----------------------

for RES in ${RESOLUTIONS}
do
    for SEED in ${SEEDS}
    do
        INDIR=${IN_BASE}/res${RES}_seed${SEED}
        OUTDIR=${OUT_BASE}/${ALIGNMENT}_res${RES}_seed${SEED}

        mkdir -p ${OUTDIR}

        if [ ! -d ${INDIR} ]; then
            echo "WARNING: input directory not found: ${INDIR}"
            continue
        fi

        echo "Running IQ-TREE for ${ALIGNMENT}, res=${RES}, seed=${SEED}"

        for ALN in ${INDIR}/cluster_*.${SUFFIX}
        do
            if [ ! -e ${ALN} ]; then
                echo "No alignment files found in ${INDIR}"
                continue
            fi

            BASE=$(basename ${ALN} .${SUFFIX})
            PREFIX=${OUTDIR}/${BASE}

            if [ -s ${PREFIX}.treefile ]; then
                echo "Skipping existing IQ-TREE output: ${PREFIX}.treefile"
                continue
            fi

            echo "Inferring IQ-TREE for ${BASE}"

            ${IQTREE} \
                -s ${ALN} \
                -m MFP \
                -B ${BOOTSTRAPS} \
                -T ${THREADS} \
                --prefix ${PREFIX}
        done
    done
done

echo "Done running IQ-TREE."
