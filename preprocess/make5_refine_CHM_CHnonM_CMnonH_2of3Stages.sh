#!/bin/bash
dirPATH=${HOME}/CHMsInOtherContexts/CellStateTransition
# shPATH=${HOME}/bin/utilities
shPATH_Y=/mnt/Storage/home/yanghui/scripts/shell
rPATH=${HOME}/bin/utilities
anPATH=${HOME}/../yanghui/annotations

CHM_organize(){
    # D: Feb-17-2022 13:47 Thu
    mkdir -p ${dirPATH}/CHMOrganization;cd ${dirPATH}/CHMOrganization
    CHM_group(){
        name=${1} # Spermatogenesis
        CHM_PATH=${2} # ${dirPATH}/Spermatogenesis_GSE137744/CHMs
        CHM_BIN=${3} # 200bp
        CHM_STAGE=${4} # "US DS PS RS"
        C_class=${5} # CHM CHnonM CMnonH

        mkdir -p ${dirPATH}/CHMOrganization/${name};cd ${dirPATH}/CHMOrganization/${name}

        refine(){
            echo ${C_class}
            touch merged.${C_class}.bed
            for stage in ${CHM_STAGE}
            do
                cut -f 1-3 ${CHM_PATH}/CHMs/${stage}_${CHM_BIN}.${C_class}.bed >> merged.${C_class}.bed
            done # for stage end

            sort -k1,1 -k2,2n merged.${C_class}.bed | mergeBed -i - -d 2000 -c 2,3 -o median,median | cut -f 1,4-5 | \
                grep -v -w "0" > tmp && mv tmp merged.${C_class}.bed #

            echo -ne "#Chrom\tStart\tEnd" > title.txt
            cat merged.${C_class}.bed > merged.${C_class}.existence.txt
            for stage in ${CHM_STAGE}
            do
                echo -ne "\t${stage}" >> title.txt
                intersectBed -c -a merged.${C_class}.bed -b ${CHM_PATH}/CHMs/${stage}_${CHM_BIN}.${C_class}.bed | \
                    cut -f 4 | awk 'BEGIN{FS=OFS="\t"}{if($1==0){print 0} else{print 1}}' | \
                    paste merged.${C_class}.existence.txt - > tmp && mv tmp merged.${C_class}.existence.txt
            done # for stage end
            echo -ne "\n" >> title.txt
            cat title.txt merged.${C_class}.existence.txt > tmp && mv tmp merged.${C_class}.existence.txt
        }
        # refine

        category(){
            gen_stable_bed(){
                # for PGC :
                # tail -n +2 merged.${C_class}.existence.txt | awk 'BEGIN{FS=OFS="\t"}{if($4==1 && ($5==1 || $6==1)){print $1, $2, $3;}}' > stable.${C_class}s.bed # 6,261
                # for others :
                # tail -n +2 merged.${C_class}.existence.txt | awk 'BEGIN{FS=OFS="\t"} {t=0; for(i=4;i<=NF;i++){t+=$i}; if(t>=3/4*(NF-3)){print $0}}' > stable.${C_class}s.bed #
                tail -n +2 merged.${C_class}.existence.txt | awk 'BEGIN{FS=OFS="\t"} {t=0; for(i=4;i<=NF;i++){t+=$i}; if(t>=2/3*(NF-3)){print $0}}' > stable.${C_class}s.bed #
                tail -n +2 merged.${C_class}.existence.txt | intersectBed -v -a - -b stable.${C_class}s.bed > stage-specific.${C_class}s.bed #
            }
            # gen_stable_bed

            filterAccordingToSignal(){
                methyl(){
                    for type in stable stage-specific
                    do
                        for stage in ${CHM_STAGE}
                        do
                            bash ${shPATH_Y}/averageMethylInRegionMultipleThreads.sh \
                                ${type}.${C_class}s.bed \
                                ${CHM_PATH}/PreparedBeforeCallCHM/${stage}.methyl.sam.G.bed \
                                ${type}.methyl_${stage}.${C_class}.txt \
                                1
                            sort -k1,1 -k2,2n ${type}.methyl_${stage}.${C_class}.txt > tmp && mv tmp ${type}.methyl_${stage}.${C_class}.txt
                        done # for stage end
                    done # for type end
                }
                methyl

                K9(){
                    for type in stable stage-specific
                    do
                        for stage in ${CHM_STAGE}
                        do
                            bash ${shPATH_Y}/averageSignalInRegion.sh \
                                ${type}.${C_class}s.bed \
                                ${CHM_PATH}/PreparedBeforeCallCHM/${stage}.H3K9me3.rmDup.bw \
                                1 \
                                ${type}.H3K9me3_${stage}.${C_class}.txt
                            sort -k1,1 -k2,2n ${type}.H3K9me3_${stage}.${C_class}.txt > tmp && mv tmp ${type}.H3K9me3_${stage}.${C_class}.txt
                        done # for stage end
                    done # for type end
                }
                K9

                merged(){
                    for type in stable stage-specific
                    do
                        paste ${type}.methyl_*.txt ${type}.H3K9me3_*.txt | \
                            grep -viw NA | grep -w "^chr[0-9]\{1,2\}" | cut -f 1-3 |sort -k1,1 -k2,2n > ${type}.filteredBySignal.bed #
                    done # for type end
                }
                # merged
            }
            filterAccordingToSignal

            nonCHM_CGIs(){
                cat *.filteredBySignal.bed | intersectBed -v -a ${anPATH}/mm10/mm10.cpgIslandExtUnmasked.bed -b - | \
                    grep -w "^chr[0-9]\{1,2\}" > nonCHM_CGIs.bed #
                STAGE1=$(echo ${CHM_STAGE} | cut -d " " -f 1)

                # DNA methylation in STAGE1
                bash ${shPATH_Y}/averageMethylInRegionMultipleThreads.sh \
                    nonCHM_CGIs.bed \
                    ${CHM_PATH}/PreparedBeforeCallCHM/${STAGE1}.methyl.sam.G.bed \
                    nonCHM_CGIs.methyl_${STAGE1}.txt \
                    1
                sort -k1,1 -k2,2n nonCHM_CGIs.methyl_${STAGE1}.txt > tmp && mv tmp nonCHM_CGIs.methyl_${STAGE1}.txt

                # H3K9me3 in STAGE2
                bash ${shPATH_Y}/averageSignalInRegion.sh \
                    nonCHM_CGIs.bed \
                    ${CHM_PATH}/PreparedBeforeCallCHM/${STAGE1}.H3K9me3.rmDup.bw \
                    1 \
                    nonCHM_CGIs.H3K9me3_${STAGE1}.txt
                sort -k1,1 -k2,2n nonCHM_CGIs.H3K9me3_${STAGE1}.txt > tmp && mv tmp nonCHM_CGIs.H3K9me3_${STAGE1}.txt

                paste nonCHM_CGIs.methyl_${STAGE1}.txt nonCHM_CGIs.H3K9me3_${STAGE1}.txt | \
                    grep -viw "NA" | awk 'BEGIN{FS=OFS="\t"}{if($4>=0.5 && $8<0.3){print $1, $2, $3}}' > nonCHM_CGIs_HMR_euch.bed

            }
            # nonCHM_CGIs
        }
        category
    }

    # CHM_group HeartDevelopment ${dirPATH}/HeartDevelopment 200bp "E10.5 E11.5 E12.5 E13.5 E14.5 E15.5 E16.5 P0" CHM
    CHM_group HeartDevelopment ${dirPATH}/HeartDevelopment 200bp "E10.5 E11.5 E12.5 E13.5 E14.5 E15.5 E16.5 P0" CHnonM
    CHM_group HeartDevelopment ${dirPATH}/HeartDevelopment 200bp "E10.5 E11.5 E12.5 E13.5 E14.5 E15.5 E16.5 P0" CMnonH
    # CHM_group Spermatogenesis ${dirPATH}/Spermatogenesis 200bp "US DS RS PS" CHM
    CHM_group Spermatogenesis ${dirPATH}/Spermatogenesis 200bp "US DS RS PS" CHnonM
    CHM_group Spermatogenesis ${dirPATH}/Spermatogenesis 200bp "US DS RS PS" CMnonH
    # CHM_group RetinalDevelopment ${dirPATH}/RetinalDevelopment 200bp "E14.5 E17.5 P0 P3 P7 P10 P14 P21" CHM
    CHM_group RetinalDevelopment ${dirPATH}/RetinalDevelopment 200bp "E14.5 E17.5 P0 P3 P7 P10 P14 P21" CHnonM
    CHM_group RetinalDevelopment ${dirPATH}/RetinalDevelopment 200bp "E14.5 E17.5 P0 P3 P7 P10 P14 P21" CMnonH
    # CHM_group EarlyEmbryogenesis ${dirPATH}/EarlyEmbryogenesis 200bp "2cell 8cell Morula ICM" CHM
    CHM_group EarlyEmbryogenesis ${dirPATH}/EarlyEmbryogenesis 200bp "2cell 8cell Morula ICM" CHnonM
    CHM_group EarlyEmbryogenesis ${dirPATH}/EarlyEmbryogenesis 200bp "2cell 8cell Morula ICM" CMnonH
    # CHM_group PGCsDevelopment ${dirPATH}/PGCsDevelopment 200bp "E10.5 E13.5_female E13.5_male" CHM
    CHM_group PGCsDevelopment ${dirPATH}/PGCsDevelopment 200bp "E10.5 E13.5_female E13.5_male" CHnonM
    CHM_group PGCsDevelopment ${dirPATH}/PGCsDevelopment 200bp "E10.5 E13.5_female E13.5_male" CMnonH
    # CHM_group LiverDevelopment ${dirPATH}/LiverDevelopment 200bp "E11.5 E12.5 E13.5 E14.5 E15.5 E16.5 P0" CHM
    CHM_group LiverDevelopment ${dirPATH}/LiverDevelopment 200bp "E11.5 E12.5 E13.5 E14.5 E15.5 E16.5 P0" CHnonM
    CHM_group LiverDevelopment ${dirPATH}/LiverDevelopment 200bp "E11.5 E12.5 E13.5 E14.5 E15.5 E16.5 P0" CMnonH
}
CHM_organize
