#!/bin/bash
set -e

# ------------------------------------------------------------
# 04_genBins.sh
# Generate uniform genomic bins for Acomys, macaque, and rat rn7
# ------------------------------------------------------------
# Main analysis uses 1000 bp bins.
# 5000 bp and 10000 bp bins are also generated for optional analyses.
# ------------------------------------------------------------

BASE_DIR=${BASE_DIR:-$(pwd)}
CHROM_DIR=${CHROM_DIR:-${BASE_DIR}/data/chrom_sizes}
BIN_DIR=${BIN_DIR:-${BASE_DIR}/data/bins}
GENBINS=${GENBINS:-${BASE_DIR}/scripts/genBins.R}

mkdir -p "${CHROM_DIR}" "${BIN_DIR}"
cd "${CHROM_DIR}"

# -----------------------
# Download chromosome-size files for Acomys species
# -----------------------

wget -c https://ftp.ensembl.org/pub/rapid-release/species/Acomys_russatus/GCA_903995435.1/ensembl/genome/Acomys_russatus-GCA_903995435.1-chromosomes.tsv.gz
wget -c https://ftp.ensembl.org/pub/rapid-release/species/Acomys_dimidiatus/GCA_907164435.1/ensembl/genome/Acomys_dimidiatus-GCA_907164435.1-chromosomes.tsv.gz

gunzip -f Acomys_russatus-GCA_903995435.1-chromosomes.tsv.gz
gunzip -f Acomys_dimidiatus-GCA_907164435.1-chromosomes.tsv.gz

sort -k1,1n Acomys_russatus-GCA_903995435.1-chromosomes.tsv \
| awk '{print "chr"$1"\t"$2}' \
> acomys_russatus.chrom.sizes

sort -k1,1n Acomys_dimidiatus-GCA_907164435.1-chromosomes.tsv \
| awk '{print "chr"$1"\t"$2}' \
> acomys_dimidiatus.chrom.sizes

# ------------------------------------------------------------
# Note for macaque and rat rn7
# ------------------------------------------------------------
# For macaca_fas and rattus_norvegicus, ready-to-use chrom.sizes
# files were not found directly on the Ensembl website. These files
# should be prepared manually from NCBI/GenBank assembly information.
#
# Expected files:
#   data/chrom_sizes/macaca_fas.chrom.sizes
#   data/chrom_sizes/rattus_norvegicus.chrom.sizes
# ------------------------------------------------------------

cd "${BASE_DIR}"

# -----------------------
# Generate uniform bins
# -----------------------

declare -A CHRSIZES
CHRSIZES[acomys_russatus]=acomys_russatus.chrom.sizes
CHRSIZES[acomys_dimidiatus]=acomys_dimidiatus.chrom.sizes
CHRSIZES[macaca_fas]=macaca_fas.chrom.sizes
CHRSIZES[rattus_norvegicus]=rattus_norvegicus.chrom.sizes

for SPC in "${!CHRSIZES[@]}"; do
    echo "Processing species: ${SPC}"

    CHROM_FILE=${CHROM_DIR}/${CHRSIZES[$SPC]}
    if [ ! -s "${CHROM_FILE}" ]; then
        echo "ERROR: Missing chromosome-size file: ${CHROM_FILE}"
        exit 1
    fi

    for SIZE in 1000 5000 10000; do
        echo "Generating ${SIZE} bp bins for ${SPC}"

        Rscript "${GENBINS}" \
            "${CHROM_FILE}" \
            "${SIZE}" \
            "${BIN_DIR}/${SPC}_${SIZE}bp.txt"
    done
done

echo "Done generating uniform bins."
