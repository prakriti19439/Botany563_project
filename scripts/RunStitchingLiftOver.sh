## adapted from /mnt/dv/wid/projects5/Roy-LongRangeEvolution/ENSEMBL_Release_102_multiple_alignments/full_uniform_bin_halLiftover/RunStitchingLiftOver.sh

exe=/mnt/dv/wid/projects5/Roy-LongRangeEvolution/ENSEMBL_Release_102_multiple_alignments/full_uniform_bin_halLiftover/liftOver/liftOver
beddir=/mnt/dv/wid/projects5/Roy-LongRangeEvolution/Wilson_Data/signal_processing/bins

declare -A cell
cell[homo_sapiens]=Human
cell[sus_scrofa]=pig
cell[rattus_norvegicus]=rat
cell[acomys_dimidiatus]=acomys_dimidiatus
cell[acomys_russatus]=acomys_russatus
cell[macaca_fas]=macaca_fas
cell[rattus_norvegicus]=rattus_norvegicus

if [ -e stitchingMapCommands.txt ]; then
    rm stitchingMapCommands.txt
fi

if [[ ! -e liftOver_stitching_outputs ]]; then
    mkdir liftOver_stitching_outputs
fi

# Loop over all species for mapping
species_list=(homo_sapiens sus_scrofa acomys_dimidiatus acomys_russatus macaca_fas rattus_norvegicus)

for S1 in "${species_list[@]}"; do
    for S2 in "${species_list[@]}"; do
        if [ "$S1" == "$S2" ]; then
            continue
        fi

        printf "Mapping ${S1} to ${S2}\n"

        # Choose appropriate bed file
        case $S1 in
            acomys_dimidiatus|acomys_russatus|macaca_fas|rattus_norvegicus)
                # For these species, use bins path
                F1=/mnt/dv/wid/projects7/Roy-singlecell2/multispecies_singlecell/pg_work/data/bins/${S1}_1000bp.txt
                ;;
            *)
                # Default bed file
                F1=${beddir}/${cell[$S1]}_1000bp.txt
                ;;
        esac

        # Chain file
        CFILE=/mnt/dv/wid/projects7/Roy-singlecell2/multispecies_singlecell/pg_work/results/chain_files_res/chain_files/${S1}_to_${S2}.linearGap_medium.chain

        # Output files for two thresholds
        OUT1=liftOver_stitching_outputs/${S1}_to_${S2}_t0.10.txt
        OUT2=liftOver_stitching_outputs/${S1}_to_${S2}_t0.50.txt

        # Add commands to stitchingMapCommands.txt
        printf "eval ${exe} -multiple -minMatch=0.10 <(sed 's/chr//g' ${F1}) ${CFILE} ${OUT1} tmp.txt\n" >> stitchingMapCommands.txt
        printf "eval ${exe} -multiple -minMatch=0.50 <(sed 's/chr//g' ${F1}) ${CFILE} ${OUT2} tmp.txt\n" >> stitchingMapCommands.txt

    done
done

# Run the commands in parallel
cat stitchingMapCommands.txt | xargs -d'\n' --max-procs=100 -I CMD bash -c CMD
