#!/bin/bash
set -e

# ----------------------------------------------------------------------
# 08_run_louvain.sh
# Run Louvain clustering on the combined region orthology graph.
# ----------------------------------------------------------------------

BASE_DIR=${BASE_DIR:-$(pwd)}
INPUT_FILE=${INPUT_FILE:-${BASE_DIR}/results/liftover/graph_outputs/Graph_AllSpecies_t0.5_wtcol1.txt}
OUTDIR=${OUTDIR:-${BASE_DIR}/results/region_level_orthology/louvain_clusters}
LOUVAIN_SCRIPT=${LOUVAIN_SCRIPT:-${BASE_DIR}/scripts/louvain_with_adj_list.py}

seed_values=${SEEDS:-"10 100"}
res_values=${RESOLUTIONS:-"100 1000"}
JOBS=${JOBS:-4}

mkdir -p "${OUTDIR}"

if [ ! -s "${INPUT_FILE}" ]; then
    echo "ERROR: missing graph input: ${INPUT_FILE}"
    exit 1
fi

if [ ! -s "${LOUVAIN_SCRIPT}" ]; then
    echo "ERROR: missing Louvain script: ${LOUVAIN_SCRIPT}"
    exit 1
fi

running=0
for seed in ${seed_values}; do
    for res in ${res_values}; do
        output_file=${OUTDIR}/louvain_clusters_res${res}_seed${seed}.txt
        echo "Running Louvain: res=${res}, seed=${seed}"
        python "${LOUVAIN_SCRIPT}" "${INPUT_FILE}" "${res}" "${output_file}" "${seed}" &

        running=$((running + 1))
        if [ "${running}" -ge "${JOBS}" ]; then
            wait
            running=0
        fi
    done
done
wait

echo "Done running Louvain clustering."
