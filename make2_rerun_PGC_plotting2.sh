#!/bin/bash
cd ~/CHMsInOtherContexts/CellStateTransition/PGCsDevelopment/CHMs
ln -sf /mnt/Storage/home/yanghui/imprinting/result.2021/Exploring/CHMsInOtherContexts/CellStateTransition/PGCsDevelopment/CHMs/E10.5_200bp.CHM.bed .
ln -sf /mnt/Storage/home/yanghui/imprinting/result.2021/Exploring/CHMsInOtherContexts/CellStateTransition/PGCsDevelopment/CHMs/E13.5_female_200bp.CHM.bed .
ln -sf /mnt/Storage/home/yanghui/imprinting/result.2021/Exploring/CHMsInOtherContexts/CellStateTransition/PGCsDevelopment/CHMs/E13.5_male_200bp.CHM.bed .

cd ~/CHMsInOtherContexts/CellStateTransition/PGCsDevelopment/Merged
ln -sf /mnt/Storage/home/yanghui/imprinting/result.2021/Exploring/CHMsInOtherContexts/CellStateTransition/PGCsDevelopment/Merged/E10.5.H3K9me3.rmDup.bw .
ln -sf /mnt/Storage/home/yanghui/imprinting/result.2021/Exploring/CHMsInOtherContexts/CellStateTransition/PGCsDevelopment/Merged/E13.5_female.H3K9me3.rmDup.bw .
ln -sf /mnt/Storage/home/yanghui/imprinting/result.2021/Exploring/CHMsInOtherContexts/CellStateTransition/PGCsDevelopment/Merged/E13.5_male.H3K9me3.rmDup.bw .
ln -sf /mnt/Storage/home/yanghui/imprinting/result.2021/Exploring/CHMsInOtherContexts/CellStateTransition/PGCsDevelopment/Merged/E10.5.methyl.all.bw .
ln -sf /mnt/Storage/home/yanghui/imprinting/result.2021/Exploring/CHMsInOtherContexts/CellStateTransition/PGCsDevelopment/Merged/E13.5_female.methyl.all.bw .
ln -sf /mnt/Storage/home/yanghui/imprinting/result.2021/Exploring/CHMsInOtherContexts/CellStateTransition/PGCsDevelopment/Merged/E13.5_male.methyl.all.bw .
ln -sf /mnt/Storage/home/yanghui/imprinting/result.2021/Exploring/CHMsInOtherContexts/CellStateTransition/PGCsDevelopment/Merged/E10.5.methyl.sam.G.bed .
ln -sf /mnt/Storage/home/yanghui/imprinting/result.2021/Exploring/CHMsInOtherContexts/CellStateTransition/PGCsDevelopment/Merged/E13.5_male.methyl.sam.G.bed .
ln -sf /mnt/Storage/home/yanghui/imprinting/result.2021/Exploring/CHMsInOtherContexts/CellStateTransition/PGCsDevelopment/Merged/E13.5_female.methyl.sam.G.bed .

mkdir -p ~/CHMsInOtherContexts/CellStateTransition/Prepare;cd ~/CHMsInOtherContexts/CellStateTransition/Prepare
awk 'BEGIN{FS=OFS="\t"}{if($4<10){print $1, $2, $3;}}' ~/../yanghui/annotations/mm10/Bins/mm10.b1kb.euchr.CpGnumber > CpGs_low.bed # 1,722,250
awk 'BEGIN{FS=OFS="\t"}{if($4>=10 && $4<30){print $1, $2, $3;}}' ~/../yanghui/annotations/mm10/Bins/mm10.b1kb.euchr.CpGnumber > CpGs_intermediate.bed # 699,059
awk 'BEGIN{FS=OFS="\t"}{if($4>=30){print $1, $2, $3;}}' ~/../yanghui/annotations/mm10/Bins/mm10.b1kb.euchr.CpGnumber > CpGs_high.bed # 41,427

process_cor(){
    # D: Apr-14-2022 16:10 Thu
    name=${1} # Spermatogenesis_GSE137744
    CHM_STAGE=${2} # "US DS PS RS"
    CHM_STAGE_plot=${3} # "US;DS;PS;RS"

    shPATH=/mnt/Storage/home/yanghui/scripts/shell
    rPATH=/mnt/Storage/home/yanghui/scripts/R
    anPATH=/mnt/Storage/home/yanghui/annotations

    dirPATH=~/CHMsInOtherContexts/CellStateTransition
    
    cd ${dirPATH}/${name}/Merged
    echo -e "xValue\tyValue\tgroup" > cor.ggplot_basic.txt
    for group in low intermediate high
    do
        for stage in ${CHM_STAGE}
        do
            # DNA methylation level
            bash ${shPATH}/averageMethylInRegionMultipleThreads.sh \
                "${dirPATH}/Prepare/CpGs_${group}.bed" \
                "${stage}.methyl.sam.G.bed" \
                "CpGs_${group}.methyl_${stage}.txt" \
                "2"
            
            # H3K9me3 signal
            bash ${shPATH}/averageSignalInRegionPrimaryOrder.sh \
                "${dirPATH}/Prepare/CpGs_${group}.bed" \
                "${stage}.H3K9me3.rmDup.bw" \
                "1" \
                "CpGs_${group}.H3K9me3_${stage}.txt"
                
            Rscript ${rPATH}/calculateCor_for_ggplot.r \
                "CpGs_${group}.methyl_${stage}.txt" \
                "CpGs_${group}.H3K9me3_${stage}.txt" \
                "4;1" \
                "F" \
                "F" \
                "${stage}" \
                "${group}" \
                "cor.ggplot_basic.txt"
        done # for stage end
    done # for group end
    
    Rscript ${rPATH}/ggplot_barplot.r \
        "cor.ggplot_basic.txt" \
        "cor.ggplot_barplot.pdf" \
        "${name}" \
        "" \
        "Pearson's correlation coefficient__n(DNA methylation and H3K9me3)" \
        "-0.2;1;0.2" \
        "${CHM_STAGE_plot}" \
        "low;intermediate;high" \
        "dodge" \
        "4" \
        "3.5" \
        "-007" \
        "F"
}


process_cor PGCsDevelopment "E9.5 E10.5 E13.5_female E13.5_male" "E9.5;E10.5;E13.5_female;E13.5_male"
