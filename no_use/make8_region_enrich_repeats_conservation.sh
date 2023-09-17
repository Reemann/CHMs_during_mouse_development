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

LINE_F_sina_enrich_group(){
    # D: Apr-18-2022 12:22 Mon
    # T: executed in 81ms, finished 12:36:47 2022-04-18
    mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LINE/SubFamily;cd ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LINE/SubFamily
    # cat *.EnrichScore_*.txt | grep -v "EnrichScore" | \
    #     awk 'BEGIN{FS=OFS="\t";}{if($2!="-inf" && $2>=1){print $1;}}' | sort -u > subfml_enriched.txt # 40
    #
    # cat *.EnrichScore_*.txt | grep -v "EnrichScore" | \
    #     awk 'BEGIN{FS=OFS="\t";}{if($2=="-inf" || $2<1){print $1;}}' | sort -u | grep -wvf subfml_enriched.txt - > subfml_nonenriched.txt # 201
    #
    # # T: executed in 334ms, finished 12:44:11 2022-04-18
    # for subfml in $(cut -f 1 subfml_enriched.txt)
    # do
    #     cat ${anPATH}/mm10/Repeats/LINE/*/mm10.${subfml}.bed | cut -f 1-3 >> subfml_enriched.bed #
    # done # for subfml

    # for subfml in $(cut -f 1 subfml_nonenriched.txt)
    # do
    #     cat ${anPATH}/mm10/Repeats/LINE/*/mm10.${subfml}.bed | cut -f 1-3 >> subfml_nonenriched.bed #
    # done # for subfml

    enriched_manage_for_table(){
        # D: Apr-19-2022 21:12 Tue
        # T: executed in 525ms, finished 21:18:13 2022-04-19
        echo -e "Process\tFamily\tSubfamily" > enriched.LINE.txt
        for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific HeartSpecific
        do
            for fml in CR1 Dong-R4 Jockey L1 L1-Tx1 L2 Penelope RTE-BovB RTE-X
            do
                tail -n +2 ${group}.EnrichScore_LINE_${fml}.txt | \
                    awk -v PROCESS=${group} -v FAMILY=${fml} 'BEGIN{FS=OFS="\t"}{\
                    if($2!="-inf" && $2>=1){print PROCESS, FAMILY, $1;}}' >> enriched.LINE.txt
            done # for fml end
        done # for group end
    }
    # enriched_manage_for_table

    enriched_overlap(){
        # D: Apr-19-2022 21:57 Tue
        Rscript ${HOME}/CHMsInOtherContexts/bin/make8_UpSetR.r
        cd ~/CHMsInOtherContexts/figures
        ln -s ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Features/RegionEnrichment/LINE/SubFamily/Overlap_UpSetR.pdf make8_UpSetR_LINE_overlap.pdf
    }
    # enriched_overlap

    phastCons60wayPlacental(){
        # D: Apr-18-2022 12:51 Mon
        echo -e "xValue\tyValue\tgroup" > phastCons60wayPlacental.ggplot_basic.txt
        for group in enriched nonenriched
        do
            bash ${shPATH}/averageSignalInRegion.sh \
                subfml_${group}.bed \
                ${anPATH}/mm10/Conservation/mm10.60way.phastCons60wayPlacental.bw \
                1 \
                ${group}.phastConsPlacental.txt
            awk -v GROUP=${group} 'BEGIN{FS=OFS="\t"}{print "LINE", $4, GROUP;}' \
                ${group}.phastConsPlacental.txt >> phastCons60wayPlacental.ggplot_basic.txt
        done # for group end

        # Violin plot
        # Rscript ${rPATH}/ggplot_violin.r \
        #     phastCons60wayPlacental.ggplot_basic.txt \
        #     phastCons60wayPlacental.ggplot_violin.pdf \
        #     LINE \
        #     "" \
        #     "Conservation score" \
        #     0 \
        #     0.6 \
        #     0.2 \
        #     "LINE" \
        #     "enriched;nonenriched" \
        #     4 \
        #     3

        # box plot
        Rscript ${rPATH}/ggplot_boxplot.r \
            phastCons60wayPlacental.ggplot_basic.txt \
            phastCons60wayPlacental.ggplot_boxplot.pdf \
            LINE \
            "" \
            "Conservation score" \
            0 \
            0.4 \
            0.2 \
            "LINE" \
            "enriched;nonenriched" \
            4 \
            3
    }
    phastCons60wayPlacental

    significance_log(){
        echo -e "
            [1] wilcox.test of LINE
            # A tibble: 1 × 8
                .y.    group1   group2          p p.adj p.format p.signif method
                <chr>  <chr>    <chr>       <dbl> <dbl> <chr>    <chr>    <chr>
            1 yValue enriched nonenriched     0     0 <2e-16   ****     Wilcoxon"
    }

    phastCons(){
        # D: Apr-18-2022 13:26 Mon
        echo -e "xValue\tyValue\tgroup" > phastCons.ggplot_basic.txt
        for group in enriched nonenriched
        do
            bash ${shPATH}/averageSignalInRegion.sh \
                subfml_${group}.bed \
                ${anPATH}/mm10/Conservation/mm10.60way.phastCons.bw \
                1 \
                ${group}.phastConsPlacental.txt
            awk -v GROUP=${group} 'BEGIN{FS=OFS="\t"}{print "LINE", $4, GROUP;}' \
                ${group}.phastConsPlacental.txt >> phastCons.ggplot_basic.txt
        done # for group end

        # Violin plot
        Rscript ${rPATH}/ggplot_violin.r \
            phastCons.ggplot_basic.txt \
            phastCons.ggplot_violin.pdf \
            LINE \
            "" \
            "Conservation score" \
            0 \
            0.4 \
            0.2 \
            "LINE" \
            "enriched;nonenriched" \
            4 \
            3

        # box plot
        Rscript ${rPATH}/ggplot_boxplot.r \
            phastCons.ggplot_basic.txt \
            phastCons.ggplot_boxplot.pdf \
            LINE \
            "" \
            "Conservation score" \
            0 \
            0.4 \
            0.2 \
            "LINE" \
            "enriched;nonenriched" \
            4 \
            3
    }
    # phastCons
}
LINE_F_sina_enrich_group

LTR_F_sina_enrich_group(){
    # # D: Apr-18-2022 13:15 Mon
    # # T:
    mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LTR/SubFamily;cd ${dirPATH}/CHMOrganization/Universal_specific/Features/RegionEnrichment/LTR/SubFamily
    # cat *.EnrichScore_*.txt | grep -v "EnrichScore" | \
    #     awk 'BEGIN{FS=OFS="\t";}{if($2!="-inf" && $2>=1){print $1;}}' | sort -u > subfml_enriched.txt #
    #
    # cat *.EnrichScore_*.txt | grep -v "EnrichScore" | \
    #     awk 'BEGIN{FS=OFS="\t";}{if($2=="-inf" || $2<1){print $1;}}' | sort -u | grep -wvf subfml_enriched.txt - > subfml_nonenriched.txt #
    #
    # # T:
    # for subfml in $(cut -f 1 subfml_enriched.txt)
    # do
    #     cat ${anPATH}/mm10/Repeats/LTR/*/mm10.${subfml}.bed | cut -f 1-3 >> subfml_enriched.bed # 238
    # done # for subfml
    #
    # for subfml in $(cut -f 1 subfml_nonenriched.txt)
    # do
    #     cat ${anPATH}/mm10/Repeats/LTR/*/mm10.${subfml}.bed | cut -f 1-3 >> subfml_nonenriched.bed # 413
    # done # for subfml

    enriched_manage_for_table(){
        # D: Apr-20-2022 00:23 Wed
        # T: executed in 269ms, finished 00:25:00 2022-04-20
        echo -e "Process\tFamily\tSubfamily" > enriched.LTR.txt
        for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific HeartSpecific
        do
            for fml in ERV1 "ERV1?" ERVK "ERVK?" ERVL ERVL-MaLR "ERVL?" Gypsy "Gypsy?"
            do
                tail -n +2 ${group}.EnrichScore_LTR_${fml}.txt | \
                    awk -v PROCESS=${group} -v FAMILY=${fml} 'BEGIN{FS=OFS="\t"}{\
                    if($2!="-inf" && $2>=1){print PROCESS, FAMILY, $1;}}' >> enriched.LTR.txt
            done # for fml end
        done # for group end
    }
    # enriched_manage_for_table

    enriched_overlap(){
        # D: Apr-20-2022 00:25 Wed
        Rscript UpSetR.r

    }
    # enriched_overlap

    phastConsPlacental(){
        # D:
        echo -e "xValue\tyValue\tgroup" > phastCons60wayPlacental.ggplot_basic.txt
        for group in enriched nonenriched
        do
            bash ${shPATH}/averageSignalInRegion.sh \
                subfml_${group}.bed \
                ${anPATH}/mm10/Conservation/mm10.60way.phastCons60wayPlacental.bw \
                1 \
                ${group}.phastConsPlacental.txt
            awk -v GROUP=${group} 'BEGIN{FS=OFS="\t"}{print "LTR", $4, GROUP;}' \
                ${group}.phastConsPlacental.txt >> phastCons60wayPlacental.ggplot_basic.txt
        done # for group end

        # Violin plot
        Rscript ${rPATH}/ggplot_violin.r \
            phastCons60wayPlacental.ggplot_basic.txt \
            phastCons60wayPlacental.ggplot_violin.pdf \
            LTR \
            "" \
            "Conservation score" \
            0 \
            0.6 \
            0.2 \
            "LTR" \
            "enriched;nonenriched" \
            4 \
            3

        # box plot
        Rscript ${rPATH}/ggplot_boxplot.r \
            phastCons60wayPlacental.ggplot_basic.txt \
            phastCons60wayPlacental.ggplot_boxplot.pdf \
            LTR \
            "" \
            "Conservation score" \
            0 \
            0.4 \
            0.2 \
            "LTR" \
            "enriched;nonenriched" \
            4 \
            3
    }
    # phastConsPlacental

    phastCons(){
        # D:
        echo -e "xValue\tyValue\tgroup" > phastCons.ggplot_basic.txt
        for group in enriched nonenriched
        do
            bash ${shPATH}/averageSignalInRegion.sh \
                subfml_${group}.bed \
                ${anPATH}/mm10/Conservation/mm10.60way.phastCons.bw \
                1 \
                ${group}.phastConsPlacental.txt
            awk -v GROUP=${group} 'BEGIN{FS=OFS="\t"}{print "LTR", $4, GROUP;}' \
                ${group}.phastConsPlacental.txt >> phastCons.ggplot_basic.txt
        done # for group end

        # Violin plot
        Rscript ${rPATH}/ggplot_violin.r \
            phastCons.ggplot_basic.txt \
            phastCons.ggplot_violin.pdf \
            LTR \
            "" \
            "Conservation score" \
            0 \
            0.6 \
            0.2 \
            "LTR" \
            "enriched;nonenriched" \
            4 \
            3

        # box plot
        Rscript ${rPATH}/ggplot_boxplot.r \
            phastCons.ggplot_basic.txt \
            phastCons.ggplot_boxplot.pdf \
            LTR \
            "" \
            "Conservation score" \
            0 \
            0.4 \
            0.2 \
            "LTR" \
            "enriched;nonenriched" \
            4 \
            3
    }
    # phastCons
}
# LTR_F_sina_enrich_group