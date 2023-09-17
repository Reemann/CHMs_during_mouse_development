#!/bin/bash
dirPATH=${HOME}/CHMsInOtherContexts/CellStateTransition
pyPATH=/mnt/Storage/home/yanghui/scripts/python
shPATH=/mnt/Storage/home/wangyiman/bin/utilities
binPATH=/mnt/Storage/home/yanghui/scripts/bin
rPATH=/mnt/Storage/home/yanghui/scripts/R
anPATH=/mnt/Storage/home/yanghui/annotations



###### ---------- region enrichment ----------
### universal CHM, universal-complementary CHM, non-CHM CpG-rich

# mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features;cd ${dirPATH}/CHMOrganization/Universal_specific/Features
# mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment;cd ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment

each_1(){
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

        bash ${shPATH}/regionEnrichment_userDefined.sh \
        "${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment" \
        "mm10_euch" \
        "${anPATH}/mm10/mm10.chrom.sizes" \
        "${dirPATH}/CHMOrganization/Universal_specific/${group_fi}" \
        "${anPATH}/mm10/RegionEnrichment/mm10_euch.CGI_Unmasked.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.promoter.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.geneBody.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.exon.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.intron.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.LTR.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.LINE.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.SINE.merged.bed;" \
        "${group}"
    done # for group end
}
# each_1



merge_and_plot_1(){
    cd ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment
    echo -e "xValue\tyValue\tgroup" > overlapRatio.ggplot_basic.txt
    echo -e "xValue\tyValue\tgroup" > regionEnrich.ggplot_basic.txt
    for group in Universal Universal_complementarySet NonCHMsCpGrich
    do
        tail -n +2 ${group}.overlaped.txt | awk -v GROUP=${group} 'BEGIN{FS=OFS="\t"}{print $1, $3*100, GROUP;}' >> overlapRatio.ggplot_basic.txt
        tail -n +2 ${group}.EnrichScore.txt | awk -v GROUP=${group} 'BEGIN{FS=OFS="\t"}{print $1, $2, GROUP;}' >> regionEnrich.ggplot_basic.txt
    done # for group end

    Rscript ${rPATH}/ggplot_barplot.r \
        "overlapRatio.ggplot_basic.txt" \
        "overlapRatio.ggplot_barplot.pdf" \
        "Overlapped ratio" \
        "" \
        "Percentage (%)" \
        "0;100;20" \
        "CGI_Unmasked;promoter;geneBody;exon;intron;LTR;LINE;SINE" \
        "Universal;Universal_complementarySet;NonCHMsCpGrich" \
        "dodge" \
        "10" \
        "3.5" \
        "-007" \
        "F"

    Rscript ${rPATH}/ggplot_barplot.r \
        "regionEnrich.ggplot_basic.txt" \
        "regionEnrich.ggplot_barplot.pdf" \
        "Enrichment score" \
        "" \
        "Enrichment score" \
        "-7;7;3" \
        "CGI_Unmasked;promoter;geneBody;exon;intron;LTR;LINE;SINE" \
        "Universal;Universal_complementarySet;NonCHMsCpGrich" \
        "dodge" \
        "10" \
        "3.5" \
        "-1;1" \
        "F"

    cd ${dirPATH}/../figures
    ln -sf ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/overlapRatio.ggplot_barplot.pdf make8_bar_overlapRatio.pdf
    ln -sf ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/regionEnrich.ggplot_barplot.pdf make8_bar_regionEnrich.pdf
}
# merge_and_plot_1


###### ---------- region enrichment (divided by A/B compartment) ----------
dataPATH_compartment=/mnt/Storage/home/wangyiman/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Features/ChromatinStates/Cavalli_Cell2017_ABcompartment_5tissue

each_cmpt(){
    cd ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/compartmentAB
    for group in UniversalCHM Universal_complementarySetCHM nonCHM
    do  
        case $group in
            "nonCHM")
                group_dir=nonCHM_intersect_table
                ;;
            *)
                group_dir=CHM_intersect_table
                ;;
        esac
        for cmpt in A B
        do
            bash ${shPATH}/regionEnrichment_userDefined.sh \
                "${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/compartmentAB" \
                "mm10_euch" \
                "${anPATH}/mm10/mm10_euch.chrom.sizes" \
                "${dataPATH_compartment}/${group_dir}/intersectOnly_${group}_5tissue_${cmpt}.bed" \
                "${anPATH}/mm10/RegionEnrichment/mm10_euch.CGI_Unmasked.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.promoter.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.geneBody.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.exon.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.intron.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.LTR.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.LINE.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.SINE.merged.bed" \
                "regionEnrich_${group}_${cmpt}"
        done # for cmpt end
    done # for group end
}
# each_cmpt


merge_and_plot_2(){
    cd ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/compartmentAB
    echo -e "xValue\tyValue\tgroup" > regionEnrich_cmpt.ggplot_basic.txt
    for group in UniversalCHM Universal_complementarySetCHM nonCHM
    do
        for cmpt in A B
        do
            tail -n +2 regionEnrich_${group}_${cmpt}.EnrichScore.txt | \
                awk -v GROUP="${group} ov ${cmpt}" 'BEGIN{FS=OFS="\t"}{\
                split($1, a, "_");\
                print a[1], $2, GROUP;}' >> regionEnrich_cmpt.ggplot_basic.txt
        done # for cmpt end
    done # for group end

    Rscript ${rPATH}/ggplot_barplot.r \
        "regionEnrich_cmpt.ggplot_basic.txt" \
        "regionEnrich_cmpt.ggplot_barplot.pdf" \
        "Region enrichment" \
        "" \
        "Enrichment score" \
        "-8;8.2;2" \
        "CGI;promoter;geneBody;exon;intron;LTR;LINE;SINE" \
        "UniversalCHM ov A;UniversalCHM ov B;Universal_complementarySetCHM ov A;Universal_complementarySetCHM ov B;nonCHM ov A;nonCHM ov B" \
        "dodge" \
        "7" \
        "3" \
        "1" \
        "F"

    cd ${dirPATH}/../figures
    ln -sf ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/compartmentAB/regionEnrich_cmpt.ggplot_barplot.pdf make8_bar_regionEnrich_compartmentAB.pdf
}

merge_and_plot_2



###### ---------- region enrichment ----------
### universal CHM, process-complementary CHM, non-CHM CpG-rich

# mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features;cd ${dirPATH}/CHMOrganization/Universal_specific/Features
# mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment;cd ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment

each(){
    for group in Universal EarlyEmbryogenesis PGCsDevelopment Spermatogenesis RetinalDevelopment HeartDevelopment LiverDevelopment Universal_complementarySet NonCHMsCpGrich
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

        bash ${shPATH}/regionEnrichment_userDefined.sh \
        "${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment" \
        "mm10_euch" \
        "${anPATH}/mm10/mm10.chrom.sizes" \
        "${dirPATH}/CHMOrganization/Universal_specific/${group_fi}" \
        "${anPATH}/mm10/RegionEnrichment/mm10_euch.CGI_Unmasked.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.promoter.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.geneBody.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.exon.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.intron.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.LTR.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.LINE.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.SINE.merged.bed;" \
        "${group}"
    done # for group end
}
# each



merge_and_plot(){
    cd ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment
    echo -e "xValue\tyValue\tgroup" > overlapRatio.ggplot_basic.txt
    echo -e "xValue\tyValue\tgroup" > regionEnrich.ggplot_basic.txt
    for group in Universal EarlyEmbryogenesis PGCsDevelopment Spermatogenesis RetinalDevelopment HeartDevelopment LiverDevelopment Universal_complementarySet NonCHMsCpGrich
    do
        tail -n +2 ${group}.overlaped.txt | awk -v GROUP=${group} 'BEGIN{FS=OFS="\t"}{print $1, $3*100, GROUP;}' >> overlapRatio.ggplot_basic.txt
        tail -n +2 ${group}.EnrichScore.txt | awk -v GROUP=${group} 'BEGIN{FS=OFS="\t"}{print $1, $2, GROUP;}' >> regionEnrich.ggplot_basic.txt
    done # for group end

    Rscript ${rPATH}/ggplot_barplot.r \
        "overlapRatio.ggplot_basic.txt" \
        "overlapRatio.ggplot_barplot.pdf" \
        "Overlapped ratio" \
        "" \
        "Percentage (%)" \
        "0;100;20" \
        "CGI_Unmasked;promoter;geneBody;exon;intron;LTR;LINE;SINE" \
        "Universal;EarlyEmbryogenesis;PGCsDevelopment;Spermatogenesis;RetinalDevelopment;HeartDevelopment;LiverDevelopment;Universal_complementarySet;NonCHMsCpGrich" \
        "dodge" \
        "10" \
        "3.5" \
        "-007" \
        "F"

    Rscript ${rPATH}/ggplot_barplot.r \
        "regionEnrich.ggplot_basic.txt" \
        "regionEnrich.ggplot_barplot.pdf" \
        "Enrichment score" \
        "" \
        "Enrichment score" \
        "-7;7;3" \
        "CGI_Unmasked;promoter;geneBody;exon;intron;LTR;LINE;SINE" \
        "Universal;EarlyEmbryogenesis;PGCsDevelopment;Spermatogenesis;RetinalDevelopment;HeartDevelopment;LiverDevelopment;Universal_complementarySet;NonCHMsCpGrich" \
        "dodge" \
        "10" \
        "3.5" \
        "-1;1" \
        "F"

    cd ${dirPATH}/../figures
    ln -sf ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/overlapRatio.ggplot_barplot.pdf make8_bar_overlapRatio.pdf
    ln -sf ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/regionEnrich.ggplot_barplot.pdf make8_bar_regionEnrich.pdf
}
# merge_and_plot
