#!/bin/bash

# Download genome assemblies for region-level phylogenomics.
# Output: data/genome_assemblies/

mkdir -p data/genome_assemblies
cd data/genome_assemblies

wget -c http://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget -c http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz
wget -c http://ftp.ensembl.org/pub/release-102/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.gz
wget -c http://ftp.ensembl.org/pub/release-102/fasta/macaca_fascicularis/dna/Macaca_fascicularis.Macaca_fascicularis_5.0.dna_sm.toplevel.fa.gz
wget -c http://ftp.ensembl.org/pub/release-112/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna_sm.toplevel.fa.gz
wget -c https://ftp.ensembl.org/pub/rapid-release/species/Acomys_russatus/GCA_903995435.1/ensembl/genome/Acomys_russatus-GCA_903995435.1-softmasked.fa.gz
wget -c https://ftp.ensembl.org/pub/rapid-release/species/Acomys_dimidiatus/GCA_907164435.1/ensembl/genome/Acomys_dimidiatus-GCA_907164435.1-softmasked.fa.gz
