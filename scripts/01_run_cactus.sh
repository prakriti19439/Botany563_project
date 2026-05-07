#!/bin/bash

# Run Progressive Cactus for kidney species WGA.
# Input: data/kidneySpecies_rn7.txt
# Output: results/cactus_alignment/kidneySpecies_rn7.hal

# module load apptainer/cactus-2

mkdir -p results/cactus_alignment
mkdir -p results/cactus_alignment/workingDir_kidney_temp
mkdir -p scratch/cactus_kidney/cactus_alignment/coord_dir

export TMPDIR=$(pwd)/results/cactus_alignment/workingDir_kidney_temp

cd results/cactus_alignment

cactus \
    --binariesMode local \
    --writeLogs \
    --logDebug \
    --maxCores 50 \
    --maxMemory 500G \
    ./jobstore_kidney \
    data/kidneySpecies_rn7.txt \
    ../results/cactus_alignment/kidneySpecies_rn7.hal \
    --coordinationDir ../scratch/cactus_kidney/cactus_alignment/coord_dir
