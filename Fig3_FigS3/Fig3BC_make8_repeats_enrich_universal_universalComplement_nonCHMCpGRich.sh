#!/bin/bash

# dirPATH=${HOME}/CHMsInOtherContexts/CellStateTransition
dirPATH=${HOME}/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Features/RegionEnrichment/Repeats
pyPATH=/mnt/Storage/home/yanghui/scripts/python
shPATH=/mnt/Storage/home/wangyiman/bin/utilities
# shPATH=/mnt/Storage/home/yanghui/scripts/shell
binPATH=/mnt/Storage/home/yanghui/scripts/bin
rPATH=/mnt/Storage/home/yanghui/scripts/R
rPATH_me=/mnt/Storage/home/wangyiman/bin/utilities
anPATH=/mnt/Storage/home/yanghui/annotations

CHMbedPATH=/mnt/Storage/home/wangyiman/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific

##### ---------- repeats subfamily enrichment ----------

##### 1. obtain enrichment score
cd ${dirPATH}
mkdir -p LTR_Subfamily LINE_Subfamily SINE_Subfamily
obtain_enrichment_score(){
    for group in Universal Universal_complementarySet NonCHMsCpGrich
    do
        case $group in
            "NonCHMsCpGrich")
                group_fi=Sequence/NonCHMsCpGrich.bed
                ;;
            "Universal")
                group_fi=Universal.CHM.bed
                ;;
            "Universal_complementarySet")
                group_fi=Universal_complementarySet.CHM.bed
                ;;
            *)
                group_fi=../Overlap/${group}_complementarySet.CHM.bed
                ;;
        esac

        for repeatType in LTR LINE SINE
        do
            bash ${shPATH}/regionEnrichment_${repeatType}_mm10_euch.sh \
                "${dirPATH}/${repeatType}_Subfamily" \
                "${CHMbedPATH}/${group_fi}" \
                "${group}"
        done # for repeatType end
    done # for group end
}
# obtain_enrichment_score


##### 2. prepare for sina plot
prepare_sina_plot(){
    for repeatType in LTR LINE SINE
    do
        cd ${dirPATH}/${repeatType}_Subfamily
        echo -e "xValue\tyValue\tgroup\tcol\tlabel" > enrichmentScore.ggplot_sina.txt
        for group in Universal Universal_complementarySet NonCHMsCpGrich
        do
            for fmlEnrich in $(ls ${group}.EnrichScore_${repeatType}_*.txt)
            do
                fmlName=$(echo ${fmlEnrich#*EnrichScore_${repeatType}_*})
                fmlName=$(echo ${fmlName%.txt})
                tail -n +2 ${fmlEnrich} | grep -v "inf" | \
                    awk -v XVALUE=${group} -v GROUP=${repeatType} -v COL=${fmlName} 'BEGIN{FS=OFS="\t"}{\
                    print XVALUE, $2, GROUP, COL, $1;}' >> enrichmentScore.ggplot_sina.txt
            done # for fmlEnrich end
        done # for group end
    done # for repeatType end
}
# prepare_sina_plot


##### 3. sina plot

LTR_sina(){
    cd ${dirPATH}/LTR_Subfamily
    Rscript ${rPATH_me}/ggplot_sina_colorfull.r \
        enrichmentScore.ggplot_sina.txt \
        enrichmentScore_filtered.ggplot_sina.pdf \
        "Enrichment score" \
        "-12;8;4" \
        "Universal;Universal_complementarySet;NonCHMsCpGrich" \
        "LTR" \
        "ERV1;ERVK;ERVL;ERVL-MaLR;Gypsy" \
        "#CE0013;#C7A609;#64C0AB;#A14C94;#1187CD" \
        "0.2" \
        "0.75" \
        "1" \
        "4.8" \
        "6.4"
}
# LTR_sina

LINE_sina(){
    cd ${dirPATH}/LINE_Subfamily
    Rscript ${rPATH_me}/ggplot_sina_colorfull.r \
        enrichmentScore.ggplot_sina.txt \
        enrichmentScore_filtered.ggplot_sina.pdf \
        "Enrichment score" \
        "-12;8;4" \
        "Universal;Universal_complementarySet;NonCHMsCpGrich" \
        "LINE" \
        "CR1;Dong-R4;Jockey;L1;L1-Tx1;L2;Penelope;RTE-BovB;RTE-X" \
        "#CE0013;#16557A;#C7A609;#87C232;#64C0AB;#A14C94;#8B7E75;#1187CD;#000000" \
        "0.2" \
        "0.75" \
        "1" \
        "4.8" \
        "5.4"
}
# LINE_sina

SINE_sina(){
    cd ${dirPATH}/SINE_Subfamily
    Rscript ${rPATH_me}/ggplot_sina_colorfull.r \
        enrichmentScore.ggplot_sina.txt \
        enrichmentScore_filtered.ggplot_sina.pdf \
        "Enrichment score" \
        "-12;8;4" \
        "Universal;Universal_complementarySet;NonCHMsCpGrich" \
        "SINE" \
        "5S-Deu-L2;Alu;B2;B4;ID;MIR;tRNA;tRNA-Deu;tRNA-RTE" \
        "#CE0013;#16557A;#C7A609;#87C232;#64C0AB;#A14C94;#8B7E75;#1187CD;#000000" \
        "0.2" \
        "0.75" \
        "1" \
        "4.8" \
        "5.4"
}
# SINE_sina

soft_link_sina_plot(){
    cd ${HOME}/CHMsInOtherContexts/figures
    for repeatType in LTR LINE SINE;do
        ln -sf ${dirPATH}/${repeatType}_Subfamily/enrichmentScore_filtered.ggplot_sina.pdf ./make8_sina_${repeatType}_overlapRatio.pdf
    done
}
# soft_link_sina_plot


##### 4. scatter plot

LTR_F(){
    cd ${dirPATH}/LTR_Subfamily
    # for group in Universal Universal_complementarySet NonCHMsCpGrich
    for group in Universal_complementarySet
    do
        echo -e "xValue\tyValue\tgroup\tshape\tsubfamily" > ${group}_enrichDetails.ggplot_scatter.txt # group: ERV1...; shape: UniversalCHMs
        # for fml in ERV1 ERVK ERVL ERVL-MaLR Gypsy
        for fml in ERV1 ERVK ERVL
        do
            paste ${group}.overlaped_LTR_${fml}.txt ${group}.EnrichScore_LTR_${fml}.txt | \
                tail -n +2 | awk -v GROUP=${fml} -v SHAPE=${group} 'BEGIN{FS=OFS="\t"}{if($1==$4 && $5!="-inf"){print $3*100, $5, GROUP, SHAPE, $1;}}' >> ${group}_enrichDetails.ggplot_scatter.txt
        done # for fml end
    
        Rscript ${rPATH}/ggplot_scatterplot_basic.r \
            ${group}_enrichDetails.ggplot_scatter.txt \
            enrichDetails \
            "Overlapped ratio (%)" \
            "Enrichment score" \
            "0;20;20" \
            "-12;8;4" \
            "ERV1;ERVK;ERVL" \
            "T;1;0" \
            "0.75" \
            "5" \
            "5" \
            "Set1"    
    
    done # for group end


}
LTR_F


LINE_F(){
    cd ${dirPATH}/LINE_Subfamily
    for group in Universal_complementarySet
    do
        echo -e "xValue\tyValue\tgroup\tshape\tsubfamily" > ${group}_enrichDetails.ggplot_scatter.txt # group: ERV1...; shape: UniversalCHMs
        for fml in CR1 L1
        do
            paste ${group}.overlaped_LINE_${fml}.txt ${group}.EnrichScore_LINE_${fml}.txt | \
                tail -n +2 | awk -v GROUP=${fml} -v SHAPE=${group} 'BEGIN{FS=OFS="\t"}{if($1==$4 && $5!="-inf"){print $3*100, $5, GROUP, SHAPE, $1;}}' >> ${group}_enrichDetails.ggplot_scatter.txt
        done # for fml end
    done # for group end
    
    Rscript ${rPATH}/ggplot_scatterplot_basic.r \
        ${group}_enrichDetails.ggplot_scatter.txt \
        enrichDetails \
        "Overlapped ratio (%)" \
        "Enrichment score" \
        "0;20;20" \
        "-12;8;4" \
        "CR1;L1" \
        "T;1;0" \
        "0.75" \
        "5" \
        "5" \
        "Set1"
}
# LINE_F


# SINE_F(){
#     cd ${dirPATH}/SINE_Subfamily
#     echo -e "xValue\tyValue\tgroup\tshape\tsubfamily" > enrichDetails.ggplot_scatter.txt # group: ERV1...; shape: UniversalCHMs
#     for group in Universal EarlyEmbryogenesis PGCsDevelopment Spermatogenesis RetinalDevelopment HeartDevelopment LiverDevelopment Universal_complementarySet NonCHMsCpGrich
#     do
#         for fml in 5S-Deu-L2 Alu B2 B4 ID MIR tRNA tRNA-Deu tRNA-RTE
#         do
#             paste ${group}.overlaped_SINE_${fml}.txt ${group}.EnrichScore_SINE_${fml}.txt | \
#                 tail -n +2 | awk -v GROUP=${fml} -v SHAPE=${group} 'BEGIN{FS=OFS="\t"}{if($1==$4 && $5!="-inf"){print $3*100, $5, GROUP, SHAPE, $1;}}' >> enrichDetails.ggplot_scatter.txt
#         done # for fml end
#     done # for group end

#     Rscript ${rPATH}/ggplot_scatterplot_basic.r \
#         enrichDetails.ggplot_scatter.txt \
#         enrichDetails \
#         "Overlapped ratio (%)" \
#         "Enrichment score" \
#         "0;100;20" \
#         "-12;8;4" \
#         "5S-Deu-L2;Alu;B2;B4;ID;MIR;tRNA;tRNA-Deu;tRNA-RTE" \
#         "T;1;0" \
#         "0.75" \
#         "5" \
#         "5" \
#         "Set3"
# }
# SINE_F

soft_link_scatter_plot(){
    cd ~/CHMsInOtherContexts/figures

    repeatType=LTR
    group=Universal_complementarySet
    ln -sf ${dirPATH}/${repeatType}_Subfamily/"Scatterplot of enrichDetails.pdf" ./make8_scatter_${group}_${repeatType}_overlapRatio.pdf

    repeatType=LINE
    group=Universal_complementarySet
    ln -sf ${dirPATH}/${repeatType}_Subfamily/"Scatterplot of enrichDetails.pdf" ./make8_scatter_${group}_${repeatType}_overlapRatio.pdf
}
soft_link_scatter_plot


##### 5. prepare for pieplot
prepare_pie_plot(){
    for repeatType in LTR LINE SINE
    do
        cd ${dirPATH}/${repeatType}_Subfamily
        for group in Universal Universal_complementarySet NonCHMsCpGrich
        do
            echo -e "value\tgroup" > ${group}.enriched.ggplot_pieplot.txt
            for fmlEnrich in $(ls ${group}.EnrichScore_${repeatType}_*.txt)
            do
                fmlName=$(echo ${fmlEnrich#*EnrichScore_${repeatType}_*})
                fmlName=$(echo ${fmlName%.txt})
                tail -n +2 ${fmlEnrich} | \
                    awk -v GROUP=${fmlName} 'BEGIN{FS=OFS="\t";TOTAL=0;}{\
                    if($2>=1){TOTAL+=1;}}END{print TOTAL, GROUP;}' >> ${group}.enriched.ggplot_pieplot.txt
            done # for fmlEnrich end
        done # for group end
    done # for repeatType end
}
# prepare_pie_plot


##### 6. pieplot

LTR_pie(){
    cd ${dirPATH}/LTR_Subfamily
    for group in Universal Universal_complementarySet NonCHMsCpGrich
    do
        Rscript ${rPATH}/ggplot_pie.r \
        ${group}.enriched.ggplot_pieplot.txt \
        ${group} \
        "ERV1;ERVK;ERVL;ERVL-MaLR;Gypsy" \
        "#CE0013;#C7A609;#64C0AB;#A14C94;#1187CD" \
        T \
        3.5 \
        3.5
    done # for group end
}
# LTR_pie

LINE_pie(){
    cd ${dirPATH}/LINE_Subfamily
    for group in Universal Universal_complementarySet NonCHMsCpGrich
    do
        Rscript ${rPATH}/ggplot_pie.r \
        ${group}.enriched.ggplot_pieplot.txt \
        ${group} \
        "CR1;Dong-R4;Jockey;L1;L1-Tx1;L2;Penelope;RTE-BovB;RTE-X" \
        "#CE0013;#16557A;#C7A609;#87C232;#64C0AB;#A14C94;#8B7E75;#1187CD;#000000" \
        T \
        3.5 \
        3.5
    done # for group end
}
# LINE_pie

SINE_pie(){
    cd ${dirPATH}/SINE_Subfamily
    for group in Universal Universal_complementarySet NonCHMsCpGrich
    do
        Rscript ${rPATH}/ggplot_pie.r \
        ${group}.enriched.ggplot_pieplot.txt \
        ${group} \
        "5S-Deu-L2;Alu;B2;B4;ID;MIR;tRNA;tRNA-Deu;tRNA-RTE" \
        "#CE0013;#16557A;#C7A609;#87C232;#64C0AB;#A14C94;#8B7E75;#1187CD;#000000" \
        T \
        3.5 \
        3.5
    done # for group end
}
# SINE_pie


soft_link_pie_plot(){
    cd ${HOME}/CHMsInOtherContexts/figures
    for repeatType in LTR LINE SINE;do
        for CHMType in Universal Universal_complementarySet NonCHMsCpGrich;do
            ln -sf "${dirPATH}/${repeatType}_Subfamily/Pie plot of ${CHMType}.pdf" ./make8_pie_${repeatType}_${CHMType}_enriched_subfamily.pdf
        done
    done
}
# soft_link_pie_plot
