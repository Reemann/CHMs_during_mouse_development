#!/bin/bash

dirPATH=${HOME}/CHMsInOtherContexts/CellStateTransition
pyPATH=/mnt/Storage/home/yanghui/scripts/python
shPATH=/mnt/Storage/home/yanghui/scripts/shell
binPATH=/mnt/Storage/home/yanghui/scripts/bin
rPATH=/mnt/Storage/home/yanghui/scripts/R
anPATH=/mnt/Storage/home/yanghui/annotations
dataPATH_methyl_hybrid=/mnt/Storage/home/yanghui/data.public/mm10/WGBS/EarlyEmbryos/sam.G.bed
dataPATH_H3K9me3_hybrid=/mnt/Storage/home/yanghui/data.public/mm10/ChIP-seq/H3K9me3_GSE97778_SNPTraceablePreImplantation_Gao_2018_NCB/Processed
dataPATH_expr_earlyEmbryo=/mnt/Storage/home/yuzhaowei/projects/imprinting/data_mm10/public/RNA-seq/GSE70605_WCF_RNA-seq

# D: Nov-06-2022 15:26 Sun
mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features/GenomeDistribution/Repeats;cd ${dirPATH}/CHMOrganization/Universal_specific/Features/GenomeDistribution/Repeats
simpleOverlap(){
    # D: Nov-06-2022 16:15 Sun
    # R: Almost CHMs located in repeats, especially in LTR
    echo -e "xValue\tyValue\tgroup" > repeats.ggplot_basic.txt
    echo -e "xValue\tyValue\tgroup" > LINE.ggplot_basic.txt
    echo -e "xValue\tyValue\tgroup" > SINE.ggplot_basic.txt
    echo -e "xValue\tyValue\tgroup" > LTR.ggplot_basic.txt
    for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific HeartSpecific LiverSpecific
    do
        cd ${dirPATH}/CHMOrganization/Universal_specific/Features/GenomeDistribution/Repeats
        N_TOTAL=$(cat ${dirPATH}/CHMOrganization/Universal_specific/${group}.CHM.bed | wc -l)
        for type in repeats LINE SINE LTR
        do
            N_TYPE=$(intersectBed -u -a ${dirPATH}/CHMOrganization/Universal_specific/${group}.CHM.bed -b ${anPATH}/mm10/Repeats/mm10.${type}.bed | wc -l)
            echo -e "${group}\t${N_TOTAL}\t${N_TYPE}\t${type}" | awk 'BEGIN{FS=OFS="\t"}{print $1, $3/$2*100, $4;}' >> ${type}.ggplot_basic.txt
        done # for type end
    done # for group end

    for type in repeats LINE SINE LTR
    do
        cd ${dirPATH}/CHMOrganization/Universal_specific/Features/GenomeDistribution/Repeats
        Rscript ${rPATH}/ggplot_barplot.r \
            "${type}.ggplot_basic.txt" \
            "${type}.ggplot_barplot.pdf" \
            "${type}" \
            "" \
            "Percentage of CHMs located in ${type}" \
            "0;100;20" \
            "Universal;EarlyEmbryoSpecific;PGCSpecific;SpermSpecific;RetinalSpecific;HeartSpecific;LiverSpecific" \
            "${type}" \
            "dodge" \
            "3" \
            "5" \
            "-007" \
            "F"
        cd ~/CHMsInOtherContexts/figures
        ln -s ${dirPATH}/CHMOrganization/Universal_specific/Features/GenomeDistribution/Repeats/${type}.ggplot_barplot.pdf make8_bar_CHM_repeats_overlapRatio_${type}.pdf
    done # for type end

}
# simpleOverlap

overlapDetails(){
    # D: Nov-06-2022 16:43 Sun
    mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features/GenomeDistribution/Repeats/OverlapDetails;cd ${dirPATH}/CHMOrganization/Universal_specific/Features/GenomeDistribution/Repeats/OverlapDetails
    # chr start end repeats LTR LINE SINE
    for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific HeartSpecific LiverSpecific
    do
        cut -f 1-3 ${dirPATH}/CHMOrganization/Universal_specific/${group}.CHM.bed > ${group}.overlapDetails.txt
        intersectBed -c -a ${group}.overlapDetails.txt -b ${anPATH}/mm10/Repeats/mm10.repeats.bed | \
            awk 'BEGIN{FS=OFS="\t"}{if($4>0){print $1, $2, $3, 1;} else{print $1, $2, $3, 0;}}' > tmp && mv tmp ${group}.overlapDetails.txt
        intersectBed -c -a ${group}.overlapDetails.txt -b ${anPATH}/mm10/Repeats/mm10.LTR.bed | \
            awk 'BEGIN{FS=OFS="\t"}{if($5>0){print $1, $2, $3, $4, 1;} else{print $1, $2, $3, $4, 0;}}' > tmp && mv tmp ${group}.overlapDetails.txt
        intersectBed -c -a ${group}.overlapDetails.txt -b ${anPATH}/mm10/Repeats/mm10.LINE.bed | \
            awk 'BEGIN{FS=OFS="\t"}{if($6>0){print $1, $2, $3, $4, $5, 1;} else{print $1, $2, $3, $4, $5, 0;}}' > tmp && mv tmp ${group}.overlapDetails.txt
        intersectBed -c -a ${group}.overlapDetails.txt -b ${anPATH}/mm10/Repeats/mm10.SINE.bed | \
            awk 'BEGIN{FS=OFS="\t"}{if($7>0){print $1, $2, $3, $4, $5, $6, 1;} else{print $1, $2, $3, $4, $5, $6, 0;}}' > tmp && mv tmp ${group}.overlapDetails.txt
        intersectBed -c -a ${group}.overlapDetails.txt -b ${anPATH}/mm10/Repeats/mm10.others.bed | \
            awk 'BEGIN{FS=OFS="\t"}{if($8>0){print $1, $2, $3, $4, $5, $6, $7, 1;} else{print $1, $2, $3, $4, $5, $6, $7, 0;}}' > tmp && mv tmp ${group}.overlapDetails.txt
    done # for group end
}
# overlapDetails

overlapGroups(){
    # D: Nov-07-2022 13:55 Mon
    mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features/GenomeDistribution/Repeats/OverlapDetails;cd ${dirPATH}/CHMOrganization/Universal_specific/Features/GenomeDistribution/Repeats/OverlapDetails
    piechart(){
        for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific HeartSpecific LiverSpecific
        do
            echo -e "value\tgroup" > ${group}.overlapDetails.ggplot_pieplot.txt
            cut -f 4- ${group}.overlapDetails.txt | sort | uniq -c | awk 'BEGIN{FS=" ";OFS="\t"}{print $1,"g_"$2$3$4$5$6}' >> ${group}.overlapDetails.ggplot_pieplot.txt
            Rscript ${rPATH}/ggplot_pie.r \
                ${group}.overlapDetails.ggplot_pieplot.txt \
                ${group}.overlapDetails \
                "g_00000;g_10001;g_10010;g_10011;g_10100;g_10101;g_10110;g_10111;g_11000;g_11001;g_11010;g_11011;g_11100;g_11101;g_11110;g_11111" \
                "#ffffff;#7fa718;#00007E;#0000FF;#007E00;#00FF00;#007E7E;#00FFFF;#7E0000;#FF0000;#7E007E;#FF00FF;#7E7E00;#FFFF00;#7E7E7E;#000000" \
                T \
                6 \
                6
        done # for group end
    }
    piechart
}
overlapGroups
