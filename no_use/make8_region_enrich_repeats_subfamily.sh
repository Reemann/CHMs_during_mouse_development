#!/bin/bash
dirPATH=${HOME}/CHMsInOtherContexts/CellStateTransition
pyPATH=/mnt/Storage/home/yanghui/scripts/python
shPATH=/mnt/Storage/home/yanghui/scripts/shell
binPATH=/mnt/Storage/home/yanghui/scripts/bin
rPATH=/mnt/Storage/home/yanghui/scripts/R
rPATH_w=${HOME}/bin/utilities
anPATH=/mnt/Storage/home/yanghui/annotations
dataPATH_methyl_hybrid=/mnt/Storage/home/yanghui/data.public/mm10/WGBS/EarlyEmbryos/sam.G.bed
dataPATH_H3K9me3_hybrid=/mnt/Storage/home/yanghui/data.public/mm10/ChIP-seq/H3K9me3_GSE97778_SNPTraceablePreImplantation_Gao_2018_NCB/Processed
dataPATH_expr_earlyEmbryo=/mnt/Storage/home/yuzhaowei/projects/imprinting/data_mm10/public/RNA-seq/GSE70605_WCF_RNA-seq


LINE_F(){
    # D: Apr-17-2022 22:11 Sun
    mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LINE/SubFamily;cd ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LINE/SubFamily
    for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific HeartSpecific LiverSpecific
    do
        bash ${shPATH}/regionEnrichment_LINE_mm10.sh \
            ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LINE/SubFamily \
            ${dirPATH}/CHMOrganization/Universal_specific/${group}.CHM.bed \
            ${group}
    done # for group end
}
# LINE_F

LTR_F(){
    # D: Apr-17-2022 22:14 Sun
    mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LTR/SubFamily;cd ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LTR/SubFamily
    for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific HeartSpecific LiverSpecific
    do
        bash ${shPATH}/regionEnrichment_LTR_mm10.sh \
            ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LTR/SubFamily \
            ${dirPATH}/CHMOrganization/Universal_specific/${group}.CHM.bed \
            ${group}
    done # for group end
}
# LTR_F

LINE_F_sina(){
    # D: Apr-17-2022 22:58 Sun
    mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LINE/SubFamily;cd ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LINE/SubFamily
    echo -e "xValue\tyValue\tgroup\tcol\tlabel" > enrichmentScore.ggplot_sina.txt
    for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific HeartSpecific LiverSpecific
    do
        for fml in CR1 Dong-R4 Jockey L1 L1-Tx1 L2 Penelope RTE-BovB RTE-X
        do
            tail -n +2 ${group}.EnrichScore_LINE_${fml}.txt | \
                awk -v XVALUE=${group} -v COL=${fml} 'BEGIN{FS=OFS="\t"}{\
                if($2=="-inf"){{print XVALUE, "-8", "LINE", COL, "";}} else{\
                if($2>=3){print XVALUE, $2, "LINE", COL, $1;} else{print XVALUE, $2, "LINE", COL, "";}}}' >> enrichmentScore.ggplot_sina.txt
        done # for fml end
    done # for group end

    Rscript ${rPATH_w}/ggplot_sina_colorfull.r \
        enrichmentScore.ggplot_sina.txt \
        enrichmentScore.ggplot_sina.pdf \
        "Enrichment score" \
        "-8;8;2" \
        "Universal;EarlyEmbryoSpecific;PGCSpecific;SpermSpecific;RetinalSpecific;HeartSpecific;LiverSpecific" \
        "LINE" \
        "CR1;Dong-R4;Jockey;L1;L1-Tx1;L2;Penelope;RTE-BovB;RTE-X" \
        "#CE0013;#16557A;#C7A609;#87C232;#64C0AB;#A14C94;#8B7E75;#1187CD;#000000" \
        "0.2" \
        "0.75" \
        "1;3" \
        "5" \
        "5"

    echo -e "xValue\tyValue\tgroup\tcol\tlabel" > enrichmentScore.ggplot_sina_details.txt
    for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific HeartSpecific LiverSpecific
    do
        for fml in CR1 Dong-R4 Jockey L1 L1-Tx1 L2 Penelope RTE-BovB RTE-X
        do
            tail -n +2 ${group}.EnrichScore_LINE_${fml}.txt | \
                awk -v XVALUE=${group} -v COL=${fml} 'BEGIN{FS=OFS="\t"}{\
                if($2=="-inf"){{print XVALUE, "-8", "LINE", COL, $1;}} else{\
                print XVALUE, $2, "LINE", COL, $1;}}' >> enrichmentScore.ggplot_sina_details.txt
        done # for fml end
    done # for group end
    cd ~/CHMsInOtherContexts/figures
    ln -s ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LINE/SubFamily/enrichmentScore.ggplot_sina.pdf make8_sina_enrichmentScore_LINE_subfamily.pdf
}
# LINE_F_sina

LTR_F_sina(){
    # D: Apr-17-2022 23:34 Sun
    mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LTR/SubFamily;cd ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LTR/SubFamily
    echo -e "xValue\tyValue\tgroup\tcol\tlabel" > enrichmentScore.ggplot_sina.txt
    for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific HeartSpecific LiverSpecific
    do
        for fml in ERV1 "ERV1?" ERVK "ERVK?" ERVL ERVL-MaLR "ERVL?" Gypsy "Gypsy?"
        do
            tail -n +2 ${group}.EnrichScore_LTR_${fml}.txt | \
                awk -v XVALUE=${group} -v COL=${fml} 'BEGIN{FS=OFS="\t"}{\
                if($2=="-inf"){{print XVALUE, "-8", "LTR", COL, "";}} else{\
                if($2>=-10){print XVALUE, $2, "LTR", COL, $1;} else{print XVALUE, $2, "LINE", COL, "";}}}' >> enrichmentScore.ggplot_sina.txt
        done # for fml end
    done # for group end

    # Rscript ${rPATH_w}/ggplot_sina_colorfull.r \
    #     enrichmentScore.ggplot_sina.txt \
    #     enrichmentScore.ggplot_sina.pdf \
    #     "Enrichment score" \
    #     "-8;8;2" \
    #     "Universal;EarlyEmbryoSpecific;PGCSpecific;SpermSpecific;RetinalSpecific;HeartSpecific" \
    #     "LTR" \
    #     "ERV1;ERV1?;ERVK;ERVK?;ERVL;ERVL-MaLR;ERVL?;Gypsy;Gypsy?" \
    #     "#CE0013;#16557A;#C7A609;#87C232;#64C0AB;#A14C94;#8B7E75;#1187CD;#000000" \
    #     "0.2" \
    #     "0.75" \
    #     "1" \
    #     "5" \
    #     "5"

    Rscript ${rPATH_w}/ggplot_sina_colorfull.r \
        enrichmentScore.ggplot_sina.txt \
        enrichmentScore_filtered.ggplot_sina.pdf \
        "Enrichment score" \
        "-8;8;2" \
        "Universal;EarlyEmbryoSpecific;PGCSpecific;SpermSpecific;RetinalSpecific;HeartSpecific;LiverSpecific" \
        "LTR" \
        "ERV1;ERVK;ERVL;ERVL-MaLR;Gypsy" \
        "#CE0013;#C7A609;#64C0AB;#A14C94;#1187CD" \
        "0.2" \
        "0.75" \
        "1" \
        "5" \
        "5"

    cd ~/CHMsInOtherContexts/figures
    ln -s ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LINE/SubFamily/enrichmentScore_filtered.ggplot_sina.pdf make8_sina_enrichmentScore_LTR_subfamily.pdf
}
# LTR_F_sina

LINE_F_sina_enrich_pie(){
    # D: Apr-18-2022 10:59 Mon
    mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LINE/SubFamily;cd ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LINE/SubFamily
    for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific HeartSpecific LiverSpecific
    do
        echo -e "value\tgroup" > ${group}.enriched.basic_pieplot.txt
        for fml in CR1 Dong-R4 Jockey L1 L1-Tx1 L2 Penelope RTE-BovB RTE-X
        do
            tail -n +2 ${group}.EnrichScore_LINE_${fml}.txt | \
                awk -v GROUP=${fml} 'BEGIN{FS=OFS="\t";TOTAL=0;}{\
                if($2>=1){TOTAL+=1;}}END{print TOTAL, GROUP;}' >> ${group}.enriched.basic_pieplot.txt
        done # for fml end

        Rscript ${rPATH}/ggplot_pie.r \
            ${group}.enriched.basic_pieplot.txt \
            ${group}_LINE \
            "CR1;Dong-R4;Jockey;L1;L1-Tx1;L2;Penelope;RTE-BovB;RTE-X" \
            "#CE0013;#16557A;#C7A609;#87C232;#64C0AB;#A14C94;#8B7E75;#1187CD;#000000" \
            T \
            3.5 \
            3.5
    done # for group end
    mkdir -p ~/CHMsInOtherContexts/figures/make8_pie_LINE_subfamily_ESover1
    cd ~/CHMsInOtherContexts/figures/make8_pie_LINE_subfamily_ESover1
    ln -s ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LINE/SubFamily/Plot\ of*LINE* .
}
LINE_F_sina_enrich_pie

LTR_F_sina_enrich_pie(){
    echo 'no need'
    # D: Apr-18-2022 11:44 Mon
    # mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LTR/SubFamily;cd ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LTR/SubFamily
    # for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific HeartSpecific LiverSpecific
    # do
    #     echo -e "value\tgroup" > ${group}.enriched.basic_pieplot.txt
    #     for fml in ERV1 "ERV1?" ERVK "ERVK?" ERVL ERVL-MaLR "ERVL?" Gypsy "Gypsy?"
    #     do
    #         tail -n +2 ${group}.EnrichScore_LTR_${fml}.txt | \
    #             awk -v GROUP=${fml} 'BEGIN{FS=OFS="\t";TOTAL=0;}{\
    #             if($2>=1){TOTAL+=1;}}END{print TOTAL, GROUP;}' >> ${group}.enriched.basic_pieplot.txt
    #     done # for fml end

    #     Rscript ${rPATH}/ggplot_pie.r \
    #         ${group}.enriched.basic_pieplot.txt \
    #         ${group} \
    #         "ERV1;ERV1?;ERVK;ERVK?;ERVL;ERVL-MaLR;ERVL?;Gypsy;Gypsy?" \
    #         "#CE0013;#16557A;#C7A609;#87C232;#64C0AB;#A14C94;#8B7E75;#1187CD;#000000" \
    #         T \
    #         3.5 \
    #         3.5
    # done # for group end
}
# LTR_F_sina_enrich_pie



LTR_F_sina_enrich_pie_filtered(){
    # D: May-02-2022 19:37 Mon
    mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LTR/SubFamily;cd ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LTR/SubFamily
    for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific HeartSpecific LiverSpecific
    do
        echo -e "value\tgroup" > ${group}.enriched.basic_pieplot_filtered.txt
        for fml in ERV1 ERVK ERVL ERVL-MaLR Gypsy
        do
            tail -n +2 ${group}.EnrichScore_LTR_${fml}.txt | \
                awk -v GROUP=${fml} 'BEGIN{FS=OFS="\t";TOTAL=0;}{\
                if($2>=1){TOTAL+=1;}}END{print TOTAL, GROUP;}' >> ${group}.enriched.basic_pieplot_filtered.txt
        done # for fml end

        Rscript ${rPATH}/ggplot_pie.r \
            ${group}.enriched.basic_pieplot_filtered.txt \
            ${group}_LTR_filtered \
            "ERV1;ERVK;ERVL;ERVL-MaLR;Gypsy" \
            "#CE0013;#C7A609;#64C0AB;#A14C94;#1187CD" \
            T \
            3.5 \
            3.5
    done # for group end
    mkdir -p ~/CHMsInOtherContexts/figures/make8_pie_LTR_subfamily_ESover1
    cd ~/CHMsInOtherContexts/figures/make8_pie_LTR_subfamily_ESover1
    ln -s ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LTR/SubFamily/Plot\ of*LTR* .

}
LTR_F_sina_enrich_pie_filtered


