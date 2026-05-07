#!/bin/bash

# ----------------------------------------------------------------------
# 15_run_mrbayes.sh
#
# Purpose:
#   Run MrBayes on selected cluster alignments to estimate Bayesian
#   phylogenetic trees and posterior probabilities.
#
# Input:
#   results/mafft_alignments/res${RES}_seed${SEED}/cluster_*.mafft.aln.fa
#   or
#   results/muscle_alignments/res${RES}_seed${SEED}/cluster_*.muscle.aln.fa
#
# Output:
#   results/mrbayes/${ALIGNMENT}_res${RES}_seed${SEED}/cluster_*.nex
#   results/mrbayes/${ALIGNMENT}_res${RES}_seed${SEED}/cluster_*.con.tre
#
# Notes:
#   - MrBayes is computationally expensive.
#   - It is recommended to run this script on selected representative
#     clusters rather than all clusters.
#   - This script uses seqmagick to convert FASTA alignment to NEXUS.
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
OUT_BASE=${BASE_DIR}/results/mrbayes

RESOLUTIONS=${RESOLUTIONS:-"1000"}
SEEDS=${SEEDS:-"10"}

# Optional: provide a file with cluster IDs to run, one per line.
# Example:
#   cluster_92502
#   cluster_12345
#
# If CLUSTER_LIST is empty, all clusters in the input directory are run.
CLUSTER_LIST=${CLUSTER_LIST:-}

# MrBayes parameters
NGEN=${NGEN:-1000000}
SAMPLEFREQ=${SAMPLEFREQ:-1000}
NCHAINS=${NCHAINS:-4}
PRINTFREQ=${PRINTFREQ:-1000}
DIAGNFREQ=${DIAGNFREQ:-10000}
BURNIN=${BURNIN:-250}

# Executables
MB=${MB:-mb}
SEQMAGICK=${SEQMAGICK:-seqmagick}

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
# Run MrBayes
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

        echo "Running MrBayes for ${ALIGNMENT}, res=${RES}, seed=${SEED}"

        # Decide which alignments to run.
        if [ -n "${CLUSTER_LIST}" ]; then
            FILES=""
            while read CLUSTER
            do
                FILES="${FILES} ${INDIR}/${CLUSTER}.${SUFFIX}"
            done < ${CLUSTER_LIST}
        else
            FILES=${INDIR}/cluster_*.${SUFFIX}
        fi

        for ALN in ${FILES}
        do
            if [ ! -e ${ALN} ]; then
                echo "Missing alignment: ${ALN}"
                continue
            fi

            BASE=$(basename ${ALN} .${SUFFIX})
            NEX=${OUTDIR}/${BASE}.nex

            if [ -s ${OUTDIR}/${BASE}.con.tre ]; then
                echo "Skipping existing MrBayes consensus tree for ${BASE}"
                continue
            fi

            echo "Preparing MrBayes input for ${BASE}"

            # Convert FASTA alignment to NEXUS.
            ${SEQMAGICK} convert ${ALN} ${NEX}

            # Append MrBayes block.
            cat >> ${NEX} << EOF

begin mrbayes;
    set autoclose=yes nowarn=yes;
    lset nst=6 rates=gamma;
    mcmc ngen=${NGEN} samplefreq=${SAMPLEFREQ} nchains=${NCHAINS} printfreq=${PRINTFREQ} diagnfreq=${DIAGNFREQ};
    sump burnin=${BURNIN};
    sumt burnin=${BURNIN};
end;
EOF

            echo "Running MrBayes for ${BASE}"

            cd ${OUTDIR}
            ${MB} ${BASE}.nex
            cd ${BASE_DIR}
        done
    done
done

echo "Done running MrBayes."
