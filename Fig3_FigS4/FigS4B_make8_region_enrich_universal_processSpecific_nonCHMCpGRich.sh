#!/bin/bash
dirPATH=${HOME}/CHMsInOtherContexts/CellStateTransition
pyPATH=/mnt/Storage/home/yanghui/scripts/python
shPATH=/mnt/Storage/home/wangyiman/bin/utilities
binPATH=/mnt/Storage/home/yanghui/scripts/bin
rPATH=/mnt/Storage/home/wangyiman/bin/utilities
anPATH=/mnt/Storage/home/yanghui/annotations

###### ---------- region enrichment ----------

mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features;cd ${dirPATH}/CHMOrganization/Universal_specific/Features
mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment;cd ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment

each(){
    for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific HeartSpecific LiverSpecific NonCHMsCpGrich
    do
        if [ $group == "NonCHMsCpGrich" ]
        then
            group_fi=Sequence/NonCHMsCpGrich.bed
        else
            group_fi=${group}.CHM.bed
        fi
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
    # D: Apr-17-2022 19:18 Sun
    cd ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment
    echo -e "xValue\tyValue\tgroup" > overlapRatio.ggplot_basic.txt
    echo -e "xValue\tyValue\tgroup" > regionEnrich.ggplot_basic.txt
    for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific HeartSpecific LiverSpecific NonCHMsCpGrich
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
        "Universal;EarlyEmbryoSpecific;PGCSpecific;SpermSpecific;RetinalSpecific;HeartSpecific;LiverSpecific;NonCHMsCpGrich" \
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
        "-6;6;3" \
        "CGI_Unmasked;promoter;geneBody;exon;intron;LTR;LINE;SINE" \
        "Universal;EarlyEmbryoSpecific;PGCSpecific;SpermSpecific;RetinalSpecific;HeartSpecific;LiverSpecific;NonCHMsCpGrich" \
        "dodge" \
        "10" \
        "3.5" \
        "-1;1" \
        "F"

    cd ${dirPATH}/../figures
    ln -s ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/overlapRatio.ggplot_barplot.pdf make8_bar_overlapRatio.pdf
    ln -s ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/regionEnrich.ggplot_barplot.pdf make8_bar_regionEnrich.pdf
}
# merge_and_plot

merge_and_plot_onlyProcessSpecific(){
    # D: Apr-17-2022 19:18 Sun
    cd ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment
    echo -e "xValue\tyValue\tgroup" > overlapRatio.ggplot_basic.processSpecific.txt
    echo -e "xValue\tyValue\tgroup" > regionEnrich.ggplot_basic.processSpecific.txt
    for group in EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific HeartSpecific LiverSpecific
    do
        tail -n +2 ${group}.overlaped.txt | awk -v GROUP=${group} 'BEGIN{FS=OFS="\t"}{print $1, $3*100, GROUP;}' >> overlapRatio.ggplot_basic.processSpecific.txt
        tail -n +2 ${group}.EnrichScore.txt | awk -v GROUP=${group} 'BEGIN{FS=OFS="\t"}{print $1, $2, GROUP;}' >> regionEnrich.ggplot_basic.processSpecific.txt
    done # for group end

    Rscript ${rPATH}/ggplot_barplot.r \
        "overlapRatio.ggplot_basic.processSpecific.txt" \
        "overlapRatio.ggplot_barplot.processSpecific.pdf" \
        "Overlapped ratio" \
        "" \
        "Percentage (%)" \
        "0;100;20" \
        "CGI_Unmasked;promoter;geneBody;exon;intron;LTR;LINE;SINE" \
        "EarlyEmbryoSpecific;PGCSpecific;SpermSpecific;RetinalSpecific;HeartSpecific;LiverSpecific" \
        "dodge" \
        "10" \
        "3.5" \
        "-007" \
        "F"

    Rscript ${rPATH}/ggplot_barplot.r \
        "regionEnrich.ggplot_basic.processSpecific.txt" \
        "regionEnrich.ggplot_barplot.processSpecific.pdf" \
        "Enrichment score" \
        "" \
        "Enrichment score" \
        "-6;6;3" \
        "CGI_Unmasked;promoter;geneBody;exon;intron;LTR;LINE;SINE" \
        "EarlyEmbryoSpecific;PGCSpecific;SpermSpecific;RetinalSpecific;HeartSpecific;LiverSpecific" \
        "dodge" \
        "10" \
        "3.5" \
        "-1;1" \
        "F"

    cd ${dirPATH}/../figures
    ln -sf ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/overlapRatio.ggplot_barplot.processSpecific.pdf make8_bar_overlapRatio_processSpecific.pdf
    ln -sf ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/regionEnrich.ggplot_barplot.processSpecific.pdf make8_bar_regionEnrich_processSpecific.pdf
}
merge_and_plot_onlyProcessSpecific
