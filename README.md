# Region-level Phylogenomics Pipeline


This repository contains a reproducible workflow for inferring rooted evolutionary trees from uniform 1 kb coding and non-coding mammalian genomic regions. The workflow starts from raw genome assemblies, performs whole-genome alignment and LiftOver-based coordinate mapping, constructs a region orthology graph, clusters genomic regions into putative orthogroups, aligns sequences, infers trees, roots/resolves the trees, and validates selected trees using LiftOver-supported region pairs.

---

## 1. Project Goal

Traditional phylogenomics pipelines usually infer trees from genes or proteins. This project instead treats uniform 1 kb genomic regions as evolutionary units so that both coding and non-coding regions can be included.

In this project, **region orthology** refers to genomic regions inferred to descend from a shared ancestral locus through whole-genome alignment and reciprocal coordinate mapping.

The main biological question is:

> Can coding and non-coding genomic regions, especially regulatory regions, retain sufficient evolutionary signal to infer meaningful rooted evolutionary trees across mammals?

---

## 2. High-Level Workflow

```text
Soft-masked genome FASTA files
↓
Whole-genome alignment with Progressive Cactus
↓
HAL alignment file
↓
Pairwise chain files
↓
LiftOver coordinate mappings
↓
Reciprocal region orthology graph
↓
Louvain clustering
↓
Cluster-specific FASTA files
↓
Multiple sequence alignment with MAFFT and/or MUSCLE
↓
Tree inference with FastTree, IQ-TREE, and MrBayes
↓
OrthoFinder-inspired rooting and resolution
↓
Rooted region-level evolutionary trees
↓
Validation using LiftOver-supported region pairs
```

---

## 3. Species and Assembly Naming

The scripts use the following species labels.

| Script label | Species | Assembly/source name |
|---|---|---|
| `homo_sapiens` | Human | GRCh38 / hg38 |
| `mus_musculus` | Mouse | GRCm38 / mm10 |
| `sus_scrofa` | Pig | Sscrofa11.1 / susScr11 |
| `rattus_norvegicus` | Rat | mRatBN7.2 / rn7 |
| `macaca_fas` | Macaca fascicularis | MacFas5 |
| `acomys_dimidiatus` | Acomys dimidiatus | GCA_907164435.1 / mAcoDim1 |
| `acomys_russatus` | Acomys russatus | GCA_903995435.1 / mAcoRus1.1 |

Important note: this workflow uses **Macaca fascicularis**, not Macaca mulatta.

---

## 4. Repository Organization

Current repository structure:

```text
Botany563_project/
│
├── README.md
├── notebook-log.md
├── .gitignore
│
├── data/
│
├── results/
│
└── scripts/
    ├── 00_download_genomes.sh
    ├── 01_run_cactus.sh
    ├── 02_MakeAllChains.sh
    ├── 03_RunCommands.sh
    ├── 04_genBins.R
    ├── 04_genBins.sh
    ├── 05_run_liftover.sh
    ├── 06_concat_graph_outputs.sh
    ├── 07_graph_stats.py
    ├── 08_run_louvain.sh
    ├── 09_gen_region_fasta.sh
    ├── 10_gen_cluster_fastas.py
    ├── 11_gen_mafft_alignment_args.py
    ├── 12_run_muscle.sh
    ├── 13_run_fasttree.py
    ├── 14_run_iqtree.sh
    ├── 15_run_mrbayes.sh
    ├── 16_gen_resolved_rooted_trees.py
    ├── 17_validate_trees_with_liftover.py
    ├── CollateResults.py
    ├── CompileStitchedGraphs.sh
    ├── OverlapLiftOverResults.sh
    └── RunStitchingLiftOver.sh
```

Recommended output folders created during the workflow:

```text
data/genome_assemblies/
data/chrom_sizes/
data/bins/

results/cactus_alignment/
results/chain_files_res/
results/liftover/
results/region_level_orthology/
results/cluster_fastas/
results/mafft_alignments/
results/muscle_alignments/
results/fasttree/
results/iqtree/
results/mrbayes/
results/rooted_trees/
results/validation/
```

---

## 5. Software Requirements

The workflow requires the following software.

| Software | Purpose |
|---|---|
| `bash` | Shell scripts |
| `wget` | Download genome/chromosome-size files |
| Progressive Cactus | Whole-genome alignment |
| HAL tools: `hal2fasta`, `halStats`, `halLiftover` | Extract genomes and generate pairwise mappings from HAL |
| UCSC tools: `faToTwoBit`, `pslPosTarget`, `axtChain`, `liftOver` | 2bit conversion, PSL processing, chain generation, coordinate mapping |
| `bedtools` | Intersect mapped regions with uniform genomic bins |
| R | Run `04_genBins.R` |
| R package `plyr` | Used by `04_genBins.R` |
| Python 3 | Python workflow scripts |
| Python packages: `pandas`, `BioPython`, `newick`, `tqdm` | Data processing, FASTA/tree parsing, progress bars |
| `scikit-network` or compatible Louvain script dependency | Louvain clustering if using `louvain_with_adj_list.py` |
| MAFFT | Multiple sequence alignment |
| MUSCLE v5 | Alternative multiple sequence alignment |
| FastTree | Approximate maximum-likelihood tree inference |
| IQ-TREE 2 | Maximum-likelihood tree inference and bootstrap support |
| MrBayes | Bayesian phylogenetic inference |
| Seqmagick | FASTA-to-NEXUS conversion for MrBayes |
| OrthoFinder source code | Rooting/resolving tree logic |

Example environment setup:

```bash
mamba create -n region_trees \
    python=3.10 pandas biopython tqdm newick numpy scipy networkx matplotlib \
    -c conda-forge -c bioconda

mamba activate region_trees

mamba install -c bioconda bedtools samtools mafft muscle fasttree iqtree mrbayes seqmagick ucsc-liftover
```

Cactus may need a separate install or cluster module:

```bash
module load apptainer/cactus-2
```

---

# 7. Reproducible Workflow

---

## Step 1. Download soft-masked genome assemblies

Script:

```text
scripts/00_download_genomes.sh
```

Purpose:

Downloads soft-masked genome FASTA files into:

```text
data/genome_assemblies/
```

Run:

```bash
bash scripts/00_download_genomes.sh
```

No command-line arguments are required.

Main output:

```text
data/genome_assemblies/
```

The script downloads genome assemblies for:

- human
- mouse
- pig
- Macaca fascicularis
- rat
- Acomys russatus
- Acomys dimidiatus

---

## Step 2. Run Progressive Cactus whole-genome alignment

Script:

```text
scripts/01_run_cactus.sh
```

Purpose:

Runs Progressive Cactus and creates a HAL whole-genome alignment.

Run:

```bash
bash scripts/01_run_cactus.sh
```

No command-line arguments are required, but the script expects the Cactus input file:

```text
data/kidneySpecies_rn7.txt
```

Expected output:

```text
results/cactus_alignment/kidneySpecies_rn7.hal
```

Important settings used in the script:

```text
--binariesMode local
--writeLogs
--logDebug
--maxCores 50
--maxMemory 500G
```

Comment:

This is the most computationally expensive step. It may require a high-memory cluster node.

---

## Step 3. Generate pairwise chain-file commands from the HAL file

Script:

```text
scripts/02_MakeAllChains.sh
```

Purpose:

Generates pairwise chain-file commands from the Progressive Cactus HAL file.

Run:

```bash
bash scripts/02_MakeAllChains.sh
```

No command-line arguments are required.

The script expects:

```text
results/cactus_alignment/kidneySpecies_rn7.hal
thirdparty/hal/bin/
thirdparty/ucsc_tools/
```

The script creates:

```text
commands.txt
commands2.txt
```

The script also creates intermediate files in:

```text
results/chain_files_res/species_files/
results/chain_files_res/results_halLiftover/
results/chain_files_res/chain_files/
```

---

## Step 4. Run generated chain-file commands

Script:

```text
scripts/03_RunCommands.sh
```

Purpose:

Runs the commands generated by `02_MakeAllChains.sh`.

Run:

```bash
bash scripts/03_RunCommands.sh
```

No command-line arguments are required.

The script runs:

```bash
commands.txt
commands2.txt
```

with:

```text
--max-procs=3
```

Expected output:

```text
results/chain_files_res/results_halLiftover/*.psl
results/chain_files_res/chain_files/*.chain
```

---

## Step 5. Generate uniform genomic bins

Scripts:

```text
scripts/04_genBins.sh
scripts/04_genBins.R
```

Purpose:

Generates 1000 bp bins for selected species.

Run:

```bash
bash scripts/04_genBins.sh
```

The shell script does not require command-line arguments, but it supports optional environment variables:

```bash
BASE_DIR=$(pwd)
CHROM_DIR=${BASE_DIR}/data/chrom_sizes
BIN_DIR=${BASE_DIR}/data/bins
GENBINS=${BASE_DIR}/scripts/genBins.R
```

To override defaults:

```bash
BASE_DIR=$(pwd) \
CHROM_DIR=data/chrom_sizes \
BIN_DIR=data/bins \
GENBINS=scripts/04_genBins.R \
bash scripts/04_genBins.sh
```

```

Expected output:

```text
data/bins/acomys_russatus_1000bp.txt
data/bins/acomys_dimidiatus_1000bp.txt
data/bins/macaca_fas_1000bp.txt
data/bins/rattus_norvegicus_1000bp.txt
```

The R script itself takes three positional arguments:

```bash
Rscript scripts/04_genBins.R <chrom_sizes_file> <bin_size> <output_file>
```

Example:

```bash
Rscript scripts/04_genBins.R \
    data/chrom_sizes/rattus_norvegicus.chrom.sizes \
    1000 \
    data/bins/rattus_norvegicus_1000bp.txt
```

---

## Step 6. Run LiftOver mapping workflow

Wrapper script:

```text
scripts/05_run_liftover.sh
```

Helper scripts:

```text
scripts/RunStitchingLiftOver.sh
scripts/OverlapLiftOverResults.sh
scripts/CompileStitchedGraphs.sh
scripts/CollateResults.py
```

Purpose:

Maps 1 kb genomic bins between species, overlaps mapped regions with target-species bins, and compiles reciprocal region mappings into graph edges.

Run:

```bash
bash scripts/05_run_liftover.sh
```

No command-line arguments are required by the wrapper.


Important parameters:

```bash
BIN_SIZE=1000
THRESHOLDS="0.10 0.50"
SPECIES_LIST="acomys_dimidiatus acomys_russatus homo_sapiens mus_musculus sus_scrofa rattus_norvegicus macaca_fas"
```

Expected output:

```text
results/liftover/liftOver_stitching_outputs/
results/liftover/overlap_outputs/
results/liftover/graph_outputs/
results/liftover/results_stats/
```

---

## Step 7. Concatenate pairwise reciprocal graph files

Script:

```text
scripts/06_concat_graph_outputs.sh
```

Purpose:

Concatenates pairwise reciprocal LiftOver graph files into one master graph file and adds a third weight column for Louvain.

Run:

```bash
bash scripts/06_concat_graph_outputs.sh
```

No command-line arguments are required, but optional environment variables can override defaults:

```bash
BASE_DIR=$(pwd)
GDIR=${BASE_DIR}/results/liftover/graph_outputs
OUT=${GDIR}/Graph_AllSpecies_t0.5.txt
WEIGHTED_OUT=${GDIR}/Graph_AllSpecies_t0.5_wtcol1.txt
THRESH=t0.50
```

Expected output:

```text
results/liftover/graph_outputs/Graph_AllSpecies_t0.5.txt
results/liftover/graph_outputs/Graph_AllSpecies_t0.5_wtcol1.txt
```

---

## Step 8. Summarize LiftOver, overlap, and reciprocal graph statistics

Script:

```text
scripts/07_graph_stats.py
```

Purpose:

Creates count summary tables for LiftOver, overlap, and reciprocal graph results.

Run with defaults:

```bash
python scripts/07_graph_stats.py
```

Available arguments:

```text
--base_dir       Default: .
--bin_dir        Default: data/bins
--liftover_dir   Default: results/liftover
--out_file       Default: results/liftover/results_stats/combined_all_stats.csv
--res            Default: 1kb
--thresholds     Default: t0.10 t0.50
```

Recommended command:

```bash
python scripts/07_graph_stats.py \
    --base_dir . \
    --bin_dir data/bins \
    --liftover_dir results/liftover \
    --out_file results/liftover/results_stats/combined_all_stats.csv \
    --res 1kb \
    --thresholds t0.10 t0.50
```

Expected output:

```text
results/liftover/results_stats/combined_all_stats.csv
```

---

## Step 9. Run Louvain clustering

Script:

```text
scripts/08_run_louvain.sh
```

Purpose:

Runs Louvain clustering on the combined weighted region orthology graph.

Run:

```bash
bash scripts/08_run_louvain.sh
```

No command-line arguments are required, but environment variables can override defaults:

```bash
INPUT_FILE=results/liftover/graph_outputs/Graph_AllSpecies_t0.5_wtcol1.txt
OUTDIR=results/region_level_orthology/louvain_clusters
LOUVAIN_SCRIPT=scripts/louvain_with_adj_list.py
SEEDS="10 100"
RESOLUTIONS="100 1000"
JOBS=4
```

Example with explicit values:

```bash
INPUT_FILE=results/liftover/graph_outputs/Graph_AllSpecies_t0.5_wtcol1.txt \
OUTDIR=results/region_level_orthology/louvain_clusters \
LOUVAIN_SCRIPT=scripts/louvain_with_adj_list.py \
SEEDS="10 100" \
RESOLUTIONS="100 1000" \
JOBS=4 \
bash scripts/08_run_louvain.sh
```

Expected output:

```text
results/region_level_orthology/louvain_clusters/louvain_clusters_res100_seed10.txt
results/region_level_orthology/louvain_clusters/louvain_clusters_res100_seed100.txt
results/region_level_orthology/louvain_clusters/louvain_clusters_res1000_seed10.txt
results/region_level_orthology/louvain_clusters/louvain_clusters_res1000_seed100.txt
```

---

## Step 10. Generate species-level 1 kb FASTA files

Script:

```text
scripts/09_gen_region_fasta.sh
```

Purpose:

Extracts DNA sequences for each 1 kb genomic bin using `bedtools getfasta`.

Run:

```bash
bash scripts/09_gen_region_fasta.sh
```

No command-line arguments are required, but environment variables can override defaults:

```bash
BASE_DIR=$(pwd)
BED_DIR=${BASE_DIR}/data/bins
GENOME_DIR=${BASE_DIR}/data/genome_assemblies
OUTDIR=${BASE_DIR}/results/region_level_orthology/species_1kb_fastas
BEDTOOLS=bedtools
```

Example with explicit values:

```bash
BED_DIR=data/bins \
GENOME_DIR=data/genome_assemblies \
OUTDIR=results/region_level_orthology/species_1kb_fastas \
BEDTOOLS=bedtools \
bash scripts/09_gen_region_fasta.sh
```

Expected output:

```text
results/region_level_orthology/species_1kb_fastas/*_1kb_regions.fa
```

FASTA header format:

```text
region_id|species
```

Example:

```text
1_1000_2000|homo_sapiens
```

This species suffix is required for the rooting step.

---

## Step 11. Generate one FASTA file per Louvain cluster

Script:

```text
scripts/10_gen_cluster_fastas.py
```

Purpose:

Uses Louvain cluster assignments and species-level FASTA files to create one FASTA file per cluster.

Required arguments:

```text
--res
--seed
```

Optional arguments:

```text
--species_fastas_dir    Default: results/region_level_orthology/species_1kb_fastas
--cluster_dir           Default: results/region_level_orthology/louvain_clusters
--out_base              Default: results/cluster_fastas
```

Example:

```bash
python scripts/10_gen_cluster_fastas.py \
    --res 1000 \
    --seed 10 \
    --species_fastas_dir results/region_level_orthology/species_1kb_fastas \
    --cluster_dir results/region_level_orthology/louvain_clusters \
    --out_base results/cluster_fastas
```

Expected output:

```text
results/cluster_fastas/res1000_seed10/cluster_*.fa
```

---

## Step 12. Run MAFFT alignments

Script:

```text
scripts/11_gen_mafft_alignment_args.py
```

Purpose:

Runs MAFFT on cluster-specific FASTA files.

Required arguments:

```text
--res
--seed
```

Optional arguments:

```text
--cluster_base           Default: results/cluster_fastas
--out_base               Default: results/mafft_alignments
--num_parallel_jobs      Default: 4
--cores_per_alignment    Default: 4
```

Example:

```bash
python scripts/11_gen_mafft_alignment_args.py \
    --res 1000 \
    --seed 10 \
    --cluster_base results/cluster_fastas \
    --out_base results/mafft_alignments \
    --num_parallel_jobs 8 \
    --cores_per_alignment 4
```

Expected output:

```text
results/mafft_alignments/res1000_seed10/cluster_*.mafft.aln.fa
```

Method details:

- If a cluster has more than 500 sequences, the script uses faster MAFFT settings.
- Smaller clusters use `--localpair --maxiterate 1000 --anysymbol`.

---

## Step 13. Run MUSCLE alignments

Script:

```text
scripts/12_run_muscle.sh
```

Purpose:

Runs MUSCLE v5 as an alternative alignment method.

Run:

```bash
bash scripts/12_run_muscle.sh
```

No command-line arguments are required.

To use a custom MUSCLE executable:

```bash
MUSCLE=/path/to/muscle bash scripts/12_run_muscle.sh
```

Expected output:

```text
results/muscle_alignments/res100_seed10/
results/muscle_alignments/res100_seed100/
results/muscle_alignments/res1000_seed10/
results/muscle_alignments/res1000_seed100/
```

---

## Step 14. Run FastTree

Script:

```text
scripts/13_run_fasttree.py
```

Purpose:

Runs FastTree on MAFFT or MUSCLE alignments.

Required arguments:

```text
--res
--seed
```

Optional arguments:

```text
--alignment    Choices: mafft, muscle. Default: mafft
--align_base   Default depends on --alignment
--out_base     Default: results/fasttree
--cores        Default: 8
--fasttree     Default: FastTree
```

Example for MAFFT alignments:

```bash
python scripts/13_run_fasttree.py \
    --res 1000 \
    --seed 10 \
    --alignment mafft \
    --align_base results/mafft_alignments \
    --out_base results/fasttree \
    --cores 8 \
    --fasttree FastTree
```

Example for MUSCLE alignments:

```bash
python scripts/13_run_fasttree.py \
    --res 1000 \
    --seed 10 \
    --alignment muscle \
    --align_base results/muscle_alignments \
    --out_base results/fasttree \
    --cores 8 \
    --fasttree FastTree
```

Expected output:

```text
results/fasttree/mafft_res1000_seed10/cluster_*.tree
results/fasttree/muscle_res1000_seed10/cluster_*.tree
```

---

## Step 15. Run IQ-TREE

Script:

```text
scripts/14_run_iqtree.sh
```

Purpose:

Runs IQ-TREE to infer maximum-likelihood trees with model selection and ultrafast bootstrap support.

Run with defaults:

```bash
bash scripts/14_run_iqtree.sh
```

Optional environment variables:

```bash
ALIGNMENT=mafft
IQTREE=iqtree2
BOOTSTRAPS=1000
THREADS=AUTO
`

Run on MAFFT alignments:

```bash
ALIGNMENT=mafft \
BOOTSTRAPS=1000 \
THREADS=AUTO \
IQTREE=iqtree2 \
bash scripts/14_run_iqtree.sh
```

Run on MUSCLE alignments:

```bash
ALIGNMENT=muscle \
BOOTSTRAPS=1000 \
THREADS=AUTO \
IQTREE=iqtree2 \
bash scripts/14_run_iqtree.sh
```

Expected output:

```text
results/iqtree/mafft_res1000_seed10/cluster_*.treefile
results/iqtree/mafft_res1000_seed10/cluster_*.iqtree
results/iqtree/mafft_res1000_seed10/cluster_*.log
```

---

## Step 16. Run MrBayes

Script:

```text
scripts/15_run_mrbayes.sh
```

Purpose:

Runs MrBayes on selected cluster alignments to estimate Bayesian trees and posterior probabilities.

Run with defaults:

```bash
bash scripts/15_run_mrbayes.sh
```

Optional environment variables:

```bash
ALIGNMENT=mafft
RESOLUTIONS="1000"
SEEDS="10"
CLUSTER_LIST=
NGEN=1000000
SAMPLEFREQ=1000
NCHAINS=4
PRINTFREQ=1000
DIAGNFREQ=10000
BURNIN=250
MB=mb
SEQMAGICK=seqmagick
```

Recommended: run only selected representative clusters.

Create a selected cluster list:

```bash
echo "cluster_92502" > selected_clusters.txt
```

Run MrBayes only on selected clusters:

```bash
ALIGNMENT=mafft \
RESOLUTIONS="1000" \
SEEDS="10" \
CLUSTER_LIST=selected_clusters.txt \
NGEN=1000000 \
SAMPLEFREQ=1000 \
NCHAINS=4 \
BURNIN=250 \
MB=mb \
SEQMAGICK=seqmagick \
bash scripts/15_run_mrbayes.sh
```

Expected output:

```text
results/mrbayes/mafft_res1000_seed10/cluster_*.nex
results/mrbayes/mafft_res1000_seed10/cluster_*.con.tre
```

---

## Step 17. Root and resolve unrooted trees

Script:

```text
scripts/16_gen_resolved_rooted_trees.py
```

Purpose:

Roots and resolves unrooted region trees using OrthoFinder code and a rooted species tree.

Required arguments:

```text
-t / --trees_dir
-s / --species_tree
-o / --out_dir
```

Optional arguments:

```text
-n / --threads        Default: 4
--tree_suffix         Default: .tree
--orthofinder_src     Default: scripts/OrthoFinder/src
```

Example for FastTree output:

```bash
python scripts/16_gen_resolved_rooted_trees.py \
    --trees_dir results/fasttree/mafft_res1000_seed10 \
    --species_tree data/species_tree_rooted.tre \
    --out_dir results/rooted_trees/fasttree_mafft_res1000_seed10 \
    --threads 10 \
    --tree_suffix .tree \
    --orthofinder_src scripts/OrthoFinder/src
```

Example for IQ-TREE output:

```bash
python scripts/16_gen_resolved_rooted_trees.py \
    --trees_dir results/iqtree/mafft_res1000_seed10 \
    --species_tree data/species_tree_rooted.tre \
    --out_dir results/rooted_trees/iqtree_mafft_res1000_seed10 \
    --threads 10 \
    --tree_suffix .treefile \
    --orthofinder_src scripts/OrthoFinder/src
```

Expected output:

```text
results/rooted_trees/*/cluster_*.tree
results/rooted_trees/*/tree_status.csv
```

Important:

Tree leaf names must include the species name after `|`, for example:

```text
1_1000_2000|homo_sapiens
```

---

## Step 18. Validate rooted trees using LiftOver-supported pairs

Script:

```text
scripts/17_validate_trees_with_liftover.py
```

Purpose:

Compares tree-derived cross-species region pairs against LiftOver-supported graph edges.

Required arguments:

```text
--liftover_graph
--out_file
```

Provide either:

```text
--tree_file
```

or:

```text
--tree_dir
```

Optional argument:

```text
--tree_suffix    Default: .tree
```

Validate one tree:

```bash
python scripts/17_validate_trees_with_liftover.py \
    --tree_file results/rooted_trees/fasttree_mafft_res1000_seed10/cluster_92502.tree \
    --liftover_graph results/liftover/graph_outputs/Graph_AllSpecies_t0.5.txt \
    --out_file results/validation/cluster_92502_validation.tsv
```

Validate all trees in a directory:

```bash
python scripts/17_validate_trees_with_liftover.py \
    --tree_dir results/rooted_trees/fasttree_mafft_res1000_seed10 \
    --tree_suffix .tree \
    --liftover_graph results/liftover/graph_outputs/Graph_AllSpecies_t0.5.txt \
    --out_file results/validation/all_tree_validation.tsv
```

Expected output columns:

```text
tree_file
tree_path
num_tree_pairs
num_liftover_pairs_in_tree
TP
FP
FN
precision
recall
f1
liftover_pairs_present_together
fraction_liftover_pairs_present_together
```

Interpretation:

- High recall means the inferred tree recovered most LiftOver-supported pairs present in that tree.
- Lower precision is expected because trees imply many indirect evolutionary relationships, while LiftOver is a stricter direct coordinate-mapping reference.

---

# 8. Full Command Order

Run from the repository root:

```bash
cd Botany563_project
```

Recommended full order:

```bash
# 1. Download genome assemblies
bash scripts/00_download_genomes.sh

# 2. Run Progressive Cactus WGA
bash scripts/01_run_cactus.sh

# 3. Generate chain-file commands
bash scripts/02_MakeAllChains.sh

# 4. Run generated chain-file commands
bash scripts/03_RunCommands.sh

# 5. Generate chromosome sizes and uniform bins
GENBINS=scripts/04_genBins.R bash scripts/04_genBins.sh

# 6. Run LiftOver workflow
bash scripts/05_run_liftover.sh

# 7. Concatenate graph outputs
bash scripts/06_concat_graph_outputs.sh

# 8. Summarize graph statistics
python scripts/07_graph_stats.py \
    --base_dir . \
    --bin_dir data/bins \
    --liftover_dir results/liftover \
    --out_file results/liftover/results_stats/combined_all_stats.csv \
    --res 1kb \
    --thresholds t0.10 t0.50

# 9. Run Louvain clustering
INPUT_FILE=results/liftover/graph_outputs/Graph_AllSpecies_t0.5_wtcol1.txt \
OUTDIR=results/region_level_orthology/louvain_clusters \
LOUVAIN_SCRIPT=scripts/louvain_with_adj_list.py \
SEEDS="10 100" \
RESOLUTIONS="100 1000" \
JOBS=4 \
bash scripts/08_run_louvain.sh

# 10. Generate species-level FASTA files
BED_DIR=data/bins \
GENOME_DIR=data/genome_assemblies \
OUTDIR=results/region_level_orthology/species_1kb_fastas \
BEDTOOLS=bedtools \
bash scripts/09_gen_region_fasta.sh

# 11. Generate cluster FASTA files
python scripts/10_gen_cluster_fastas.py \
    --res 1000 \
    --seed 10 \
    --species_fastas_dir results/region_level_orthology/species_1kb_fastas \
    --cluster_dir results/region_level_orthology/louvain_clusters \
    --out_base results/cluster_fastas

# 12. Run MAFFT alignments
python scripts/11_gen_mafft_alignment_args.py \
    --res 1000 \
    --seed 10 \
    --cluster_base results/cluster_fastas \
    --out_base results/mafft_alignments \
    --num_parallel_jobs 8 \
    --cores_per_alignment 4

# 13. Optional: run MUSCLE alignments
bash scripts/12_run_muscle.sh

# 14. Run FastTree
python scripts/13_run_fasttree.py \
    --res 1000 \
    --seed 10 \
    --alignment mafft \
    --align_base results/mafft_alignments \
    --out_base results/fasttree \
    --cores 8 \
    --fasttree FastTree

# 15. run IQ-TREE
ALIGNMENT=mafft \
BOOTSTRAPS=1000 \
THREADS=AUTO \
IQTREE=iqtree2 \
bash scripts/14_run_iqtree.sh

# 16. run MrBayes on selected clusters
echo "cluster_92502" > selected_clusters.txt

ALIGNMENT=mafft \
RESOLUTIONS="1000" \
SEEDS="10" \
CLUSTER_LIST=selected_clusters.txt \
NGEN=1000000 \
SAMPLEFREQ=1000 \
NCHAINS=4 \
BURNIN=250 \
MB=mb \
SEQMAGICK=seqmagick \
bash scripts/15_run_mrbayes.sh

# 17. Root and resolve FastTree trees
python scripts/16_gen_resolved_rooted_trees.py \
    --trees_dir results/fasttree/mafft_res1000_seed10 \
    --species_tree data/species_tree_rooted.tre \
    --out_dir results/rooted_trees/fasttree_mafft_res1000_seed10 \
    --threads 10 \
    --tree_suffix .tree \
    --orthofinder_src scripts/OrthoFinder/src

# 18. Validate rooted trees
python scripts/17_validate_trees_with_liftover.py \
    --tree_dir results/rooted_trees/fasttree_mafft_res1000_seed10 \
    --tree_suffix .tree \
    --liftover_graph results/liftover/graph_outputs/Graph_AllSpecies_t0.5.txt \
    --out_file results/validation/all_tree_validation.tsv
```

---

# 9. Biological Interpretation

This workflow tests whether non-coding and coding genomic regions retain enough evolutionary information to reconstruct local region-level evolutionary relationships.

Overall,
- Human and Macaca fascicularis should often cluster together.
- Mouse and rat should often cluster together.
- Acomys dimidiatus and Acomys russatus should often cluster together.
- Accessible chromatin regions show stronger consistency with LiftOver-supported mappings.
- Region trees may differ from the global species tree because local genomic regions can have different evolutionary histories.

---

# 10. Limitations

- LiftOver is conservative and may miss deeply diverged orthologous regions.
- Repetitive regions can produce ambiguous mappings.
- Region-level trees may reflect local genomic evolution rather than the full species tree.
- Some clusters may include paralogous or ambiguously mapped regions.
- MrBayes is computationally expensive and is best used on selected representative clusters.

---

# 11. Future Directions

Possible extensions:

- Compare coding vs non-coding region trees.
- Compare accessible vs inaccessible chromatin regions.
- Test genomic language model embeddings for region orthology inference.
- Use graph neural networks for region-level orthology prediction.
