declare -A cell
cell[homo_sapiens]=Human
cell[sus_scrofa]=pig
cell[mus_musculus]=mouse
cell[rattus_norvegicus]=rat
cell[acomys_dimidiatus]=acomys_dimidiatus
cell[acomys_russatus]=acomys_russatus
cell[macaca_fas]=macaca_fas

if [ -e overlapCommands.txt ]; then
    rm overlapCommands.txt
fi

beddir=/mnt/dv/wid/projects5/Roy-LongRangeEvolution/Wilson_Data/signal_processing/bins

# Bedtools executable
bexe=/mnt/dv/wid/projects2/Roy-common/programs/thirdparty/bedtools2/bin/bedtools

mkdir -p overlap_outputs
ODIR=overlap_outputs

SDIR=liftOver_stitching_outputs

species_list=(homo_sapiens sus_scrofa acomys_dimidiatus acomys_russatus macaca_fas rattus_norvegicus)

for S1 in "${species_list[@]}"; do
    for S2 in "${species_list[@]}"; do
        if [[ $S1 == $S2 ]]; then
            continue
        fi

        # Choose bed file path depending on species
        case $S2 in
            acomys_dimidiatus|acomys_russatus|macaca_fas|rattus_norvegicus)
                F2=/mnt/dv/wid/projects7/Roy-singlecell2/multispecies_singlecell/pg_work/data/bins/${S2}_1000bp.txt
                ;;
            *)
                F2=${beddir}/${cell[$S2]}_1000bp.txt
                ;;
        esac

        # LiftOver outputs
        OUT2=${SDIR}/${S1}_to_${S2}_t0.10.txt
        OUT3=${SDIR}/${S1}_to_${S2}_t0.50.txt

        # Add commands
        printf "${bexe} intersect -wo -a <(cut -f1-4 -d$'\t' ${OUT2}) -b <(sed 's/chr//g' ${F2}) > ${ODIR}/${cell[$S1]}_to_${cell[$S2]}_t0.10_overlap.txt\n" >> overlapCommands.txt
        printf "${bexe} intersect -wo -a <(cut -f1-4 -d$'\t' ${OUT3}) -b <(sed 's/chr//g' ${F2}) > ${ODIR}/${cell[$S1]}_to_${cell[$S2]}_t0.50_overlap.txt\n" >> overlapCommands.txt
    done
done

# Run commands in parallel
cat overlapCommands.txt | xargs -d'\n' --max-procs=100 -I CMD bash -c CMD
