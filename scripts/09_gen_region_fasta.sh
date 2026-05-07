#!/bin/bash
set -e

# ----------------------------------------------------------------------
# 09_gen_region_fasta.sh
# Extract 1000 bp region FASTA sequences for each species.
# ----------------------------------------------------------------------

BASE_DIR=${BASE_DIR:-$(pwd)}
BED_DIR=${BED_DIR:-${BASE_DIR}/data/bins}
GENOME_DIR=${GENOME_DIR:-${BASE_DIR}/data/genome_assemblies}
OUTDIR=${OUTDIR:-${BASE_DIR}/results/region_level_orthology/species_1kb_fastas}
BEDTOOLS=${BEDTOOLS:-bedtools}

mkdir -p "${OUTDIR}"

declare -A fasta
fasta[acomys_dimidiatus]="${GENOME_DIR}/Acomys_dimidiatus-GCA_907164435.1-softmasked.fa"
fasta[acomys_russatus]="${GENOME_DIR}/Acomys_russatus-GCA_903995435.1-softmasked.fa"
fasta[homo_sapiens]="${GENOME_DIR}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
fasta[mus_musculus]="${GENOME_DIR}/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa"
fasta[rattus_norvegicus]="${GENOME_DIR}/Rattus_norvegicus.mRatBN7.2.dna_sm.toplevel.fa"
fasta[sus_scrofa]="${GENOME_DIR}/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
fasta[macaca_fas]="${GENOME_DIR}/Macaca_fascicularis.Macaca_fascicularis_5.0.dna_sm.toplevel.fa"

species_list=(homo_sapiens sus_scrofa mus_musculus acomys_dimidiatus acomys_russatus macaca_fas rattus_norvegicus)

for sp in "${species_list[@]}"; do
    echo "Processing ${sp}"

    bedfile=${BED_DIR}/${sp}_1000bp.txt
    genome=${fasta[$sp]}

    if [ ! -s "${bedfile}" ]; then
        echo "ERROR: missing BED/bin file: ${bedfile}"
        exit 1
    fi

    if [ ! -s "${genome}" ]; then
        echo "ERROR: missing genome FASTA: ${genome}"
        exit 1
    fi

    # Make temporary BED with region ID + species name in the 4th column.
    # The region ID becomes: chrom_start_end|species
    awk -v sp="${sp}" 'BEGIN{OFS="\t"} {
        gsub(/^chr/, "", $1);
        region=$1"_"$2"_"$3;
        print $1, $2, $3, region"|"sp;
    }' "${bedfile}" \
    | ${BEDTOOLS} getfasta -fi "${genome}" -bed - -fo "${OUTDIR}/${sp}_1kb_regions.fa" -nameOnly

done

echo "Done generating species 1kb FASTA files."
