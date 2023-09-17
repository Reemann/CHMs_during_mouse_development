#!/bin/bash

dirPATH=${HOME}/CHMsInOtherContexts/CellStateTransition
dirPATH_Y=/mnt/Storage/home/yanghui/imprinting/result.2021/Exploring/CHMsInOtherContexts/CellStateTransition
shPATH=/mnt/Storage/home/yanghui/scripts/shell
rPATH=/mnt/Storage/home/yanghui/scripts/R
dataPATH_CHM_earlyEmbryo=/mnt/Storage/home/yanghui/imprinting/result.2021/Filtering/CpGrichedHighMethylH3K9me3/Hybrid/ChromHMM/Manage
anPATH=/mnt/Storage/home/yanghui/annotations

prepare_ln(){
    for process in PGCsDevelopment 'RetinalDevelopment/GSE87064_Mouse' Spermatogenesis_GSE137744;do
        name1=${process%%_*}
        name=${name1%%/*}
        cd ${dirPATH}/${name}/CHMs
        ln -sf ${dirPATH_Y}/${process}/CHMs/Binarized_H3K9me3_* .
        ln -sf ${dirPATH_Y}/${process}/CHMs/Binarized_Methyl_* .
        ln -sf ${dirPATH_Y}/${process}/CHMs/Binarized_CpGNumber_mm10_200bp .
    done
    cd ${dirPATH}/EarlyEmbryogenesis/CHMs
    for stage in 2cell 8cell Morula ICM;do
        cd ${dirPATH}/EarlyEmbryogenesis/CHMs
        for chr in $(seq 19);do
            mkdir -p Binarized_H3K9me3_${stage}_200bp; cd Binarized_H3K9me3_${stage}_200bp
            ln -sf ${dataPATH_CHM_earlyEmbryo}/../Separate/${stage}/Binarized_H3K9me3_TreatOnlyPValue0.0001/${stage}_chr${chr}_binary.txt ${stage}_200bp_chr${chr}_binary.txt
            cd ..
            mkdir -p Binarized_Methyl_${stage}_200bp; cd Binarized_Methyl_${stage}_200bp
            ln -sf ${dataPATH_CHM_earlyEmbryo}/../Separate/${stage}/Binarized_Methyl0.5/${stage}_chr${chr}_binary.txt                     ${stage}_200bp_chr${chr}_binary.txt
            cd ..
            ln -sf ${dirPATH}/PGCsDevelopment/CHMs/Binarized_CpGNumber_mm10_200bp .
        done
        cd ${dirPATH}/EarlyEmbryogenesis/PreparedBeforeCallCHM
        ln -sf /mnt/Storage/home/yanghui/data.public/mm10/WGBS/EarlyEmbryos/sam.G.bed/${stage}.sam.G.bed ${stage}.WGBS.sam.G.bed
    done
}
# prepare_ln

callchm(){
	Name=$1
    process=$2
	OutDir=~/CHMsInOtherContexts/CellStateTransition/${process}/CHMs
	GenomeSiz=${anPATH}/mm10/mm10_euch.chrom.sizes
    cooocupancy(){
        cd ${OutDir}
        pwd
        # step4. binarize according to co-occupancy of CpG number, H3K9me3 and DNA methylation
        mkdir -p ${OutDir}/Binarized_Cooccupancy_${Name}.CHM
        mkdir -p ${OutDir}/Binarized_Cooccupancy_${Name}.CHnonM
        mkdir -p ${OutDir}/Binarized_Cooccupancy_${Name}.CMnonH
        for chrom in $(cut -f 1 ${GenomeSiz} | sort -u)
        do
            # CHM
            echo -e "${Name}\t${chrom}\nMerged" > ${OutDir}/Binarized_Cooccupancy_${Name}.CHM/${Name}_${chrom}_binary.txt
            paste Binarized_Methyl_${Name}/${Name}_${chrom}_binary.txt Binarized_H3K9me3_${Name}/${Name}_${chrom}_binary.txt Binarized_CpGNumber_mm10_200bp/mm10_${chrom}_binary.txt | \
            tail -n +3 | awk 'BEGIN{FS=OFS="\t"}{if($1==2){print 2} else{if($1==1 && $2==1 && $3==1){print 1} else{print 0}}}' >> ${OutDir}/Binarized_Cooccupancy_${Name}.CHM/${Name}_${chrom}_binary.txt #####
            # CH-nonM
            echo -e "${Name}\t${chrom}\nMerged" > ${OutDir}/Binarized_Cooccupancy_${Name}.CHnonM/${Name}_${chrom}_binary.txt
            paste Binarized_Methyl_${Name}/${Name}_${chrom}_binary.txt Binarized_H3K9me3_${Name}/${Name}_${chrom}_binary.txt Binarized_CpGNumber_mm10_200bp/mm10_${chrom}_binary.txt | \
            tail -n +3 | awk 'BEGIN{FS=OFS="\t"}{if($1==2){print 2} else{if($1==0 && $2==1 && $3==1){print 1} else{print 0}}}' >> ${OutDir}/Binarized_Cooccupancy_${Name}.CHnonM/${Name}_${chrom}_binary.txt #####
            # CM-nonH
            echo -e "${Name}\t${chrom}\nMerged" > ${OutDir}/Binarized_Cooccupancy_${Name}.CMnonH/${Name}_${chrom}_binary.txt
            paste Binarized_Methyl_${Name}/${Name}_${chrom}_binary.txt Binarized_H3K9me3_${Name}/${Name}_${chrom}_binary.txt Binarized_CpGNumber_mm10_200bp/mm10_${chrom}_binary.txt | \
            tail -n +3 | awk 'BEGIN{FS=OFS="\t"}{if($1==2){print 2} else{if($1==1 && $2==0 && $3==1){print 1} else{print 0}}}' >> ${OutDir}/Binarized_Cooccupancy_${Name}.CMnonH/${Name}_${chrom}_binary.txt #####

        done # for chrom end
    }
    cooocupancy
    LearnModel(){
        # step5. learn model
        for lm_Name in ${Name}.CHM ${Name}.CHnonM ${Name}.CMnonH;do
            pwd
            cd ${OutDir}
            java -Xmx50g -jar /mnt/Storage/home/yanghui/scripts/bin/ChromHMM/ChromHMM.jar LearnModel -b 200 -init random -p 6 Binarized_Cooccupancy_${lm_Name} Output_${lm_Name} 2 mm10
            cd ${OutDir}/Output_${lm_Name}
            State=$(tail -n +2 emissions_2_999.txt | awk 'BEGIN{FS=OFS="\t";S="E1";MAX=-1;}{if($2>MAX){MAX=$2;S=$1}}END{print "E"S}')
            tail -n +2 ${Name}_2_999_segments.bed | grep -w ${State} | cut -f 1-3 | sort -k1,1 -k2,2n | mergeBed -i - -d 2000 | awk 'BEGIN{FS=OFS="\t"}{if($3-$2>=600 && $2!=0){print }}' > ../${lm_Name}.bed
        done
    }
    LearnModel
}


unset DISPLAY

# for stage in E10.5;do
# for stage in E13.5_female E13.5_male;do
#     cd ~/CHMsInOtherContexts/CellStateTransition/PGCsDevelopment/CHMs
#     callchm ${stage}_200bp PGCsDevelopment > ~/CHMsInOtherContexts/CellStateTransition/PGCsDevelopment/CHMs/${stage}_pcar_CHM_CHnonM_CMnonH.log
# done

# for stage in E14.5 E17.5 P0 P3 P7 P10 P14 P21;do
#     cd ~/CHMsInOtherContexts/CellStateTransition/RetinalDevelopment/CHMs
#     callchm ${stage}_200bp RetinalDevelopment > ~/CHMsInOtherContexts/CellStateTransition/RetinalDevelopment/CHMs/${stage}_pcar_CHM_CHnonM_CMnonH.log
# done

# for stage in US DS PS RS;do
#     cd ~/CHMsInOtherContexts/CellStateTransition/Spermatogenesis/CHMs
#     callchm ${stage}_200bp Spermatogenesis > ~/CHMsInOtherContexts/CellStateTransition/Spermatogenesis/CHMs/${stage}_pcar_CHM_CHnonM_CMnonH.log 2>&1
# done

# export DISPLAY=:0
# for stage in 2cell 8cell Morula ICM;do
#     cd ~/CHMsInOtherContexts/CellStateTransition/EarlyEmbryogenesis/CHMs
#     callchm ${stage}_200bp EarlyEmbryogenesis > ~/CHMsInOtherContexts/CellStateTransition/EarlyEmbryogenesis/CHMs/${stage}_pcar_CHM_CHnonM_CMnonH.log 2>&1
# done

cd ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment/CHMs
ln -s ${dirPATH}/PGCsDevelopment/CHMs/Binarized_CpGNumber_mm10_200bp Binarized_CpGNumber_mm10_200bp
for stage in E11.5 E12.5 E15.5 E16.5;do
    cd ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment/CHMs
    callchm ${stage}_200bp LiverDevelopment > ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment/CHMs/${stage}_pcar_CHM_CHnonM_CMnonH.log 2>&1
done
