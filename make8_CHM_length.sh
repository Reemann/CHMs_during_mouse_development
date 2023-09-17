#!/bin/bash
    length(){
        # D: Apr-19-2022 11:58 Tue
        mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features/Length;cd ${dirPATH}/CHMOrganization/Universal_specific/Features/Length
        manage_length(){
            echo -e "xValue\tyValue\tgroup" > length.ggplot_basic.txt
            for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific HeartSpecific LiverSpecific
            do
                awk -v XVALUE=${group} 'BEGIN{FS=OFS="\t"}{print XVALUE, $3-$2, "Length";}' ${dirPATH}/CHMOrganization/Universal_specific/${group}.bed >> length.ggplot_basic.txt
            done # for group end
        }
        # manage_length

        boxplot(){
            Rscript ${rPATH}/ggplot_boxplot.r \
                length.ggplot_basic.txt \
                length.ggplot_boxplot.pdf \
                Length \
                "" \
                "Length of stable CHMs" \
                0 \
                12000 \
                2000 \
                "Universal;EarlyEmbryoSpecific;PGCSpecific;SpermSpecific;RetinalSpecific;HeartSpecific;LiverSpecific" \
                Length \
                3.5 \
                3.5
        }
        # boxplot

        violinplot(){
            Rscript ${rPATH}/ggplot_violin.r \
                length.ggplot_basic.txt \
                length.ggplot_violinplot.pdf \
                Length \
                "" \
                "Length of stable CHMs" \
                0 \
                12000 \
                2000 \
                "Universal;EarlyEmbryoSpecific;PGCSpecific;SpermSpecific;RetinalSpecific;HeartSpecific;LiverSpecific" \
                Length \
                4.5 \
                3.5
        }
        violinplot
    }
    # length
