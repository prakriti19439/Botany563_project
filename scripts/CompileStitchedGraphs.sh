#!/bin/bash
set -e

# ----------------------------------------------------------------------
# CompileStitchedGraphs.sh
# Collate reciprocal LiftOver overlap files into pairwise graph edge files.
# ----------------------------------------------------------------------

BASE_DIR=${BASE_DIR:-$(pwd)}
ODIR=${1:-${BASE_DIR}/results/liftover/overlap_outputs}
GDIR=${2:-${BASE_DIR}/results/liftover/graph_outputs}
STATS_DIR=${3:-${BASE_DIR}/results/liftover/results_stats}
BIN_SIZE=${4:-1000}
THRESHOLDS=${5:-"0.10 0.50"}
SPECIES_LIST=${6:-"homo_sapiens sus_scrofa acomys_dimidiatus acomys_russatus macaca_fas rattus_norvegicus mus_musculus"}
COLLATE=${COLLATE:-${BASE_DIR}/scripts/run_liftover/CollateResults.py}
JOBS=${JOBS:-10}
ALLOW_ONE_DIRECTION=${ALLOW_ONE_DIRECTION:-0}

mkdir -p "${GDIR}" "${STATS_DIR}"

if [ ! -s "${COLLATE}" ]; then
    echo "ERROR: Missing CollateResults.py: ${COLLATE}"
    exit 1
fi

declare -A cell
cell[homo_sapiens]=Human
cell[sus_scrofa]=pig
cell[mus_musculus]=mouse
cell[rattus_norvegicus]=rat
cell[acomys_dimidiatus]=acomys_dimidiatus
cell[acomys_russatus]=acomys_russatus
cell[macaca_fas]=macaca_fas

COMMANDS=${GDIR}/collateCommands.txt
rm -f "${COMMANDS}"

for S1 in ${SPECIES_LIST}; do
    for S2 in ${SPECIES_LIST}; do
        if [ "${S1}" == "${S2}" ]; then
            continue
        fi

        for T in ${THRESHOLDS}; do
            THRESH="t${T}"
            OUT1=${ODIR}/${cell[$S1]}_to_${cell[$S2]}_${THRESH}_overlap.txt
            OUT2=${ODIR}/${cell[$S2]}_to_${cell[$S1]}_${THRESH}_overlap.txt
            GOUT=${GDIR}/Graph_${cell[$S1]}_to_${cell[$S2]}_${THRESH}.txt

            if [ "${ALLOW_ONE_DIRECTION}" == "1" ]; then
                echo "python ${COLLATE} ${OUT1} ${OUT2} ${S1} ${S2} --allow_one_direction > ${GOUT}" >> "${COMMANDS}"
            else
                echo "python ${COLLATE} ${OUT1} ${OUT2} ${S1} ${S2} > ${GOUT}" >> "${COMMANDS}"
            fi
        done
    done
done

cat "${COMMANDS}" | xargs -d'\n' --max-procs=${JOBS} -I CMD bash -c CMD

echo "Done compiling graph outputs in ${GDIR}"
