#!/bin/bash
dirPATH=${HOME}/CHMsInOtherContexts/CellStateTransition
dirPATH_Y=/mnt/Storage/home/yanghui/imprinting/result.2021/Exploring/CHMsInOtherContexts/CellStateTransition
pyPATH=/mnt/Storage/home/yanghui/scripts/python
shPATH=/mnt/Storage/home/yanghui/scripts/shell
binPATH=/mnt/Storage/home/yanghui/scripts/bin
rPATH=/mnt/Storage/home/yanghui/scripts/R
anPATH=/mnt/Storage/home/yanghui/annotations
dataPATH_methyl_hybrid=/mnt/Storage/home/yanghui/data.public/mm10/WGBS/EarlyEmbryos/sam.G.bed
dataPATH_H3K9me3_hybrid=/mnt/Storage/home/yanghui/data.public/mm10/ChIP-seq/H3K9me3_GSE97778_SNPTraceablePreImplantation_Gao_2018_NCB/Processed
dataPATH_expr_earlyEmbryo=/mnt/Storage/home/yuzhaowei/projects/imprinting/data_mm10/public/RNA-seq/GSE70605_WCF_RNA-seq


###### ---------- Methyl & H3K9me3 signal (1kb bin) ----------

methyl(){
    # mkdir -p ${dirPATH}/CHMOrganization/1kbBins;
    cd ${dirPATH}/CHMOrganization/1kbBins

    echo -e "xValue\tyValue\tgroup" > methyl.ggplot_basic.txt
    echo -e "xValue\tyValue\tgroup" > H3K9me3.ggplot_basic.txt
    # PGC development
    for stage in PGC_E10.5 PGC_E13.5_female PGC_E13.5_male
    do
        cat ${dirPATH_Y}/PGCsDevelopment/Merged/CpGs_*.methyl_${stage}.txt | \
            awk -v XVALUE=${stage} 'BEGIN{FS=OFS="\t"}{print XVALUE, $4, "methyl";}' >> methyl.ggplot_basic.txt
        cat ${dirPATH_Y}/PGCsDevelopment/Merged/CpGs_*.H3K9me3_${stage}.txt | \
            awk -v XVALUE=${stage} 'BEGIN{FS=OFS="\t"}{print XVALUE, $1, "H3K9me3";}' >> H3K9me3.ggplot_basic.txt
    done # for stage end

    # Spermatogenesis_GSE137744
    for stage in US DS PS RS
    do
        cat ${dirPATH_Y}/Spermatogenesis_GSE137744/Merged/CpGs_*.methyl_${stage}.txt | \
            awk -v XVALUE=${stage} 'BEGIN{FS=OFS="\t"}{print XVALUE, $4, "methyl";}' >> methyl.ggplot_basic.txt
        cat ${dirPATH_Y}/Spermatogenesis_GSE137744/Merged/CpGs_*.H3K9me3_${stage}.txt | \
            awk -v XVALUE=${stage} 'BEGIN{FS=OFS="\t"}{print XVALUE, $1, "H3K9me3";}' >> H3K9me3.ggplot_basic.txt
    done # for stage end

    # RetinalDevelopment
    for stage in E14.5 E17.5 P0 P3 P7 P10 P14 P21
    do
        cat ${dirPATH_Y}/RetinalDevelopment/GSE87064_Mouse/Merged/CpGs_*.methyl_${stage}.txt | \
            awk -v XVALUE=${stage} 'BEGIN{FS=OFS="\t"}{print XVALUE, $4, "methyl";}' >> methyl.ggplot_basic.txt
        cat ${dirPATH_Y}/RetinalDevelopment/GSE87064_Mouse/Merged/CpGs_*.H3K9me3_${stage}.txt | \
            awk -v XVALUE=${stage} 'BEGIN{FS=OFS="\t"}{print XVALUE, $1, "H3K9me3";}' >> H3K9me3.ggplot_basic.txt
    done # for stage end

    # HeartDevelopment
    for stage in E10.5 E11.5 E12.5 E13.5 E14.5 E15.5 E16.5 P0
    do
        cat ${dirPATH}/HeartDevelopment/PreparedBeforeCallCHM/CpGs_*.methyl_${stage}.txt | \
            awk -v XVALUE=Heart_${stage} 'BEGIN{FS=OFS="\t"}{print XVALUE, $4, "methyl";}' >> methyl.ggplot_basic.txt
        cat ${dirPATH}/HeartDevelopment/PreparedBeforeCallCHM/CpGs_*.H3K9me3_${stage}.txt | \
            awk -v XVALUE=Heart_${stage} 'BEGIN{FS=OFS="\t"}{print XVALUE, $1, "H3K9me3";}' >> H3K9me3.ggplot_basic.txt
    done # for stage end

    # DNA methylation level
    Rscript ~/bin/utilities/ggplot_boxplot.r \
        methyl.ggplot_basic.txt \
        methyl.ggplot_boxplot.pdf \
        "DNA methylation level" \
        "" \
        "DNA methylation level" \
        "0" \
        "1" \
        "0.2" \
        "E10.5;E13.5_female;E13.5_male;US;DS;PS;RS;E14.5;E17.5;P0;P3;P7;P10;P14;P21;Heart_E10.5;Heart_E11.5;Heart_E12.5;Heart_E13.5;Heart_E14.5;Heart_E15.5;Heart_E16.5;Heart_P0" \
        "methyl" \
        "6" \
        "4"

    # H3K9me3
    Rscript ~/bin/utilities/ggplot_boxplot.r \
            H3K9me3.ggplot_basic.txt \
            H3K9me3.ggplot_boxplot.pdf \
            "H3K9me3" \
            "" \
            "H3K9me3" \
            "0" \
            "0.2" \
            "0.1" \
            "E10.5;E13.5_female;E13.5_male;US;DS;PS;RS;E14.5;E17.5;P0;P3;P7;P10;P14;P21;Heart_E10.5;Heart_E11.5;Heart_E12.5;Heart_E13.5;Heart_E14.5;Heart_E15.5;Heart_E16.5;Heart_P0" \
            "H3K9me3" \
            "6" \
            "4"

    cd ~/CHMsInOtherContexts/figures
    ln -sf ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/1kbBins/H3K9me3.ggplot_boxplot.pdf make7_boxplot_H3K9me3_global_signal.pdf
    ln -sf ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/1kbBins/methyl.ggplot_boxplot.pdf make7_boxplot_methyl_global_signal.pdf

}

# methyl

###### ---------- Methyl & H3K9me3 signal (1kb bin CpGrich) ----------
cd ${dirPATH}/CHMOrganization/1kbBins
echo -e "xValue\tyValue\tgroup" > CpGrich_methyl.ggplot_basic.txt
echo -e "xValue\tyValue\tgroup" > CpGrich_H3K9me3.ggplot_basic.txt
# PGC development
for stage in E10.5 E13.5_female E13.5_male
do
    cat ${dirPATH_Y}/PGCsDevelopment/Merged/CpGs_high.methyl_${stage}.txt | \
        awk -v XVALUE=${stage} 'BEGIN{FS=OFS="\t"}{print XVALUE, $4, "methyl";}' >> CpGrich_methyl.ggplot_basic.txt
    cat ${dirPATH_Y}/PGCsDevelopment/Merged/CpGs_high.H3K9me3_${stage}.txt | \
        awk -v XVALUE=${stage} 'BEGIN{FS=OFS="\t"}{print XVALUE, $1, "H3K9me3";}' >> CpGrich_H3K9me3.ggplot_basic.txt
done # for stage end

# Spermatogenesis_GSE137744
for stage in US DS PS RS
do
    cat ${dirPATH_Y}/Spermatogenesis_GSE137744/Merged/CpGs_high.methyl_${stage}.txt | \
        awk -v XVALUE=${stage} 'BEGIN{FS=OFS="\t"}{print XVALUE, $4, "methyl";}' >> CpGrich_methyl.ggplot_basic.txt
    cat ${dirPATH_Y}/Spermatogenesis_GSE137744/Merged/CpGs_high.H3K9me3_${stage}.txt | \
        awk -v XVALUE=${stage} 'BEGIN{FS=OFS="\t"}{print XVALUE, $1, "H3K9me3";}' >> CpGrich_H3K9me3.ggplot_basic.txt
done # for stage end

# RetinalDevelopment
for stage in E14.5 E17.5 P0 P3 P7 P10 P14 P21
do
    cat ${dirPATH_Y}/RetinalDevelopment/GSE87064_Mouse/Merged/CpGs_high.methyl_${stage}.txt | \
        awk -v XVALUE=${stage} 'BEGIN{FS=OFS="\t"}{print XVALUE, $4, "methyl";}' >> CpGrich_methyl.ggplot_basic.txt
    cat ${dirPATH_Y}/RetinalDevelopment/GSE87064_Mouse/Merged/CpGs_high.H3K9me3_${stage}.txt | \
        awk -v XVALUE=${stage} 'BEGIN{FS=OFS="\t"}{print XVALUE, $1, "H3K9me3";}' >> CpGrich_H3K9me3.ggplot_basic.txt
done # for stage end

# HeartDevelopment
for stage in E10.5 E11.5 E12.5 E13.5 E14.5 E15.5 E16.5 P0
do
    cat ${dirPATH}/HeartDevelopment/PreparedBeforeCallCHM/CpGs_high.methyl_${stage}.txt | \
        awk -v XVALUE=Heart_${stage} 'BEGIN{FS=OFS="\t"}{print XVALUE, $4, "methyl";}' >> CpGrich_methyl.ggplot_basic.txt
    cat ${dirPATH}/HeartDevelopment/PreparedBeforeCallCHM/CpGs_high.H3K9me3_${stage}.txt | \
        awk -v XVALUE=Heart_${stage} 'BEGIN{FS=OFS="\t"}{print XVALUE, $1, "H3K9me3";}' >> CpGrich_H3K9me3.ggplot_basic.txt
done # for stage end


# DNA methylation level
Rscript ~/bin/utilities/ggplot_boxplot.r \
      CpGrich_methyl.ggplot_basic.txt \
      CpGrich_methyl.ggplot_boxplot.pdf \
      "DNA methylation level" \
      "" \
      "DNA methylation level" \
      "0" \
      "1" \
      "0.2" \
      "E10.5;E13.5_female;E13.5_male;US;DS;PS;RS;E14.5;E17.5;P0;P3;P7;P10;P14;P21;Heart_E10.5;Heart_E11.5;Heart_E12.5;Heart_E13.5;Heart_E14.5;Heart_E15.5;Heart_E16.5;Heart_P0" \
      "methyl" \
      "6" \
      "4"

# H3K9me3
Rscript ~/bin/utilities/ggplot_boxplot.r \
        CpGrich_H3K9me3.ggplot_basic.txt \
        CpGrich_H3K9me3.ggplot_boxplot.pdf \
        "H3K9me3" \
        "" \
        "H3K9me3" \
        "0" \
        "0.75" \
        "0.25" \
        "E10.5;E13.5_female;E13.5_male;US;DS;PS;RS;E14.5;E17.5;P0;P3;P7;P10;P14;P21;Heart_E10.5;Heart_E11.5;Heart_E12.5;Heart_E13.5;Heart_E14.5;Heart_E15.5;Heart_E16.5;Heart_P0" \
        "H3K9me3" \
        "6" \
        "4"

cd ~/CHMsInOtherContexts/figures
ln -sf ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/1kbBins/CpGrich_H3K9me3.ggplot_boxplot.pdf make7_boxplot_H3K9me3_global_signal_CpGrich.pdf
ln -sf ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/1kbBins/CpGrich_methyl.ggplot_boxplot.pdf   make7_boxplot_methyl_global_signal_CpGrich.pdf





    H3K9me3_profile_DC(){
        # D: Apr-24-2022 18:15 Sun
        for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific DCSpecific
        do
            for stage in MPP CDP cDC pDC
            do
                bash ${shPATH}/profileSurroundingPeakCenter.sh \
                    ${dirPATH}/CHMOrganization/Universal_specific/${group}.bed \
                    10000 \
                    10000 \
                    100 \
                    ${stage}.H3K9me3.rmDup.bw \
                    ${stage} \
                    H3K9me3Surrounding${group} \
                    "-10kb;Center;10kb" \
                    "0;50.5;100" \
                    "H3K9me3 signal (RPM)" \
                    "0" \
                    "0.3" \
                    "Dark2"
            done # for stage end
        done # for group end
    }
    # H3K9me3_profile_DC

