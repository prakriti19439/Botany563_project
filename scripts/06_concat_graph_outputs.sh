#!/bin/bash
set -e

# ----------------------------------------------------------------------
# 06_concat_graph_outputs.sh
# Concatenate pairwise reciprocal LiftOver graph files into one master graph.
# ----------------------------------------------------------------------

BASE_DIR=${BASE_DIR:-$(pwd)}
GDIR=${GDIR:-${BASE_DIR}/results/liftover/graph_outputs}
OUT=${OUT:-${GDIR}/Graph_AllSpecies_t0.5.txt}
WEIGHTED_OUT=${WEIGHTED_OUT:-${GDIR}/Graph_AllSpecies_t0.5_wtcol1.txt}
THRESH=${THRESH:-t0.50}

# Common names are used in graph file names produced by CompileStitchedGraphs.sh.
declare -A cell
cell[homo_sapiens]=Human
cell[sus_scrofa]=pig
cell[mus_musculus]=mouse
cell[rattus_norvegicus]=rat
cell[acomys_dimidiatus]=acomys_dimidiatus
cell[acomys_russatus]=acomys_russatus
cell[macaca_fas]=macaca_fas

species_list=(homo_sapiens sus_scrofa acomys_dimidiatus acomys_russatus macaca_fas rattus_norvegicus mus_musculus)

mkdir -p "${GDIR}"
: > "${OUT}"

for ((i=0; i<${#species_list[@]}; i++)); do
    for ((j=i+1; j<${#species_list[@]}; j++)); do
        S1=${species_list[i]}
        S2=${species_list[j]}
        FILE=${GDIR}/Graph_${cell[$S1]}_to_${cell[$S2]}_${THRESH}.txt

        if [ -s "${FILE}" ]; then
            echo "Adding ${FILE}"
            cat "${FILE}" >> "${OUT}"
        else
            echo "WARNING: missing or empty graph file: ${FILE}"
        fi
    done
done

# Add weight column for Louvain code that expects three columns.
awk -v OFS='\t' '{print $1, $2, 1}' "${OUT}" > "${WEIGHTED_OUT}"

echo "Wrote master graph: ${OUT}"
echo "Wrote weighted graph: ${WEIGHTED_OUT}"
