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


signal_boxplot(){
    prepare_CHM(){
        C_class=$1
        cd ${dirPATH}/CHMOrganization/Universal_specific
        ### universal + specific CHMs :
        cat Universal.${C_class}.bed EarlyEmbryoSpecific.${C_class}.bed PGCSpecific.${C_class}.bed SpermSpecific.${C_class}.bed RetinalSpecific.${C_class}.bed HeartSpecific.${C_class}.bed LiverSpecific.${C_class}.bed > ordered.${C_class}.bed
        intersectBed -wao -a ordered.${C_class}.bed -b ${dirPATH}/CHMOrganization/KMeansCluster/Union_refined.${C_class}.existence.txt -f 1.00 -r | \
            cut -f 4-11 > ordered.${C_class}.existence.txt
    }
    # prepare_CHM CHM
    # prepare_CHM CHnonM
    # prepare_CHM CMnonH

    prepare_signal(){
        cd ${dirPATH}/CHMOrganization/Universal_specific/Signals

        for stage in 2cell 8cell Morula ICM
        do
            ln -s ${dirPATH_Y}/EarlyEmbryogenesis/Merged/${stage}.methyl.sam.G.bed ${stage}.methyl.sam.G.bed
            ln -s ${dirPATH_Y}/EarlyEmbryogenesis/Merged/${stage}.H3K9me3.rmDup.bw ${stage}.H3K9me3.rmDup.bw
        done # for stage end

        for stage in E10.5 E13.5_female E13.5_male
        do
            ln -s ${dirPATH_Y}/PGCsDevelopment/Merged/${stage}.methyl.sam.G.bed PGC_${stage}.methyl.sam.G.bed
            ln -s ${dirPATH_Y}/PGCsDevelopment/Merged/${stage}.H3K9me3.rmDup.bw PGC_${stage}.H3K9me3.rmDup.bw
        done # for stage end

        for stage in US DS PS RS
        do
            ln -s ${dirPATH_Y}/Spermatogenesis_GSE137744/Merged/${stage}.methyl.sam.G.bed ${stage}.methyl.sam.G.bed
            ln -s ${dirPATH_Y}/Spermatogenesis_GSE137744/Merged/${stage}.H3K9me3.rmDup.bw ${stage}.H3K9me3.rmDup.bw
        done # for stage end

        for stage in E14.5 E17.5 P0 P3 P7 P10 P14 P21
        do
            ln -s ${dirPATH_Y}/RetinalDevelopment/GSE87064_Mouse/Merged/${stage}.methyl.sam.G.bed Retinal_${stage}.methyl.sam.G.bed
            ln -s ${dirPATH_Y}/RetinalDevelopment/GSE87064_Mouse/Merged/${stage}.H3K9me3.rmDup.bw Retinal_${stage}.H3K9me3.rmDup.bw
        done # for stage end

        for stage in E10.5 E11.5 E12.5 E13.5 E14.5 E15.5 E16.5 P0
        do
            ln -sf ${dirPATH}/HeartDevelopment/Merged/${stage}.methyl.sam.G.bed Heart_${stage}.methyl.sam.G.bed
            ln -sf ${dirPATH}/HeartDevelopment/Merged/${stage}.H3K9me3.rmDup.bw Heart_${stage}.H3K9me3.rmDup.bw
        done # for stage end

        for stage in E11.5 E12.5 E13.5 E14.5 E15.5 E16.5 P0
        do
            ln -sf ${dirPATH}/LiverDevelopment/PreparedBeforeCallCHM/${stage}.methyl.sam.G.bed Liver_${stage}.methyl.sam.G.bed
            ln -sf ${dirPATH}/LiverDevelopment/PreparedBeforeCallCHM/${stage}.H3K9me3.rmDup.bw Liver_${stage}.H3K9me3.rmDup.bw
        done # for stage end
    }
    # prepare_signal

    scan_signal(){
        methyl(){
            C_class=$1
            cd ${dirPATH}/CHMOrganization/Universal_specific/Signals
            echo -ne "#Chrom\tStart\tEnd" > title.txt
            # sort -k1,1 -k2,2n ${dirPATH}/CHMOrganization/Universal_specific/ordered.${C_class}.bed > ordered.${C_class}.methyl.txt
            for stage in 2cell 8cell Morula ICM PGC_E10.5 PGC_E13.5_female PGC_E13.5_male US DS PS RS Retinal_E14.5 Retinal_E17.5 Retinal_P0 Retinal_P3 Retinal_P7 Retinal_P10 Retinal_P14 Retinal_P21 Heart_E10.5 Heart_E11.5 Heart_E12.5 Heart_E13.5 Heart_E14.5 Heart_E15.5 Heart_E16.5 Heart_P0 Liver_E11.5 Liver_E12.5 Liver_E13.5 Liver_E14.5 Liver_E15.5 Liver_E16.5 Liver_P0
            do
                echo -ne "\t${stage}" >> title.txt
                bash ${shPATH}/averageMethylInRegionMultipleThreads.sh \
                    ${dirPATH}/CHMOrganization/Universal_specific/ordered.${C_class}.bed \
                    ${stage}.methyl.sam.G.bed \
                    ordered.${C_class}.methyl_${stage}.txt \
                    1
                sort -k1,1 -k2,2n ordered.${C_class}.methyl_${stage}.txt | cut -f 4 | paste ordered.${C_class}.methyl.txt - > tmp && mv tmp ordered.${C_class}.methyl.txt
            done # for stage end
            echo -ne "\n" >> title.txt
            intersectBed -wao -a ${dirPATH}/CHMOrganization/Universal_specific/ordered.${C_class}.bed -b ordered.${C_class}.methyl.txt -f 1.00 -r | \
                cut -f 4-40 >> title.txt
            mv title.txt ordered.${C_class}.methyl.final.txt
        }
        methyl CHM
        methyl CHnonM
        methyl CMnonH

        K9(){
            C_class=$1
            cd ${dirPATH}/CHMOrganization/Universal_specific/Signals
            echo -ne "#Chrom\tStart\tEnd" > title.txt
            cat ${dirPATH}/CHMOrganization/Universal_specific/ordered.bed > ordered.H3K9me3.txt
            for stage in 2cell 8cell Morula ICM PGC_E10.5 PGC_E13.5_female PGC_E13.5_male US DS PS RS Retinal_E14.5 Retinal_E17.5 Retinal_P0 Retinal_P3 Retinal_P7 Retinal_P10 Retinal_P14 Retinal_P21 Heart_E10.5 Heart_E11.5 Heart_E12.5 Heart_E13.5 Heart_E14.5 Heart_E15.5 Heart_E16.5 Heart_P0 Liver_E11.5 Liver_E12.5 Liver_E13.5 Liver_E14.5 Liver_E15.5 Liver_E16.5 Liver_P0
            # for stage in 2cell 8cell Morula ICM PGC_E10.5 PGC_E13.5_female PGC_E13.5_male US DS PS RS Retinal_E14.5 Retinal_E17.5 Retinal_P0 Retinal_P3 Retinal_P7 Retinal_P10 Retinal_P14 Retinal_P21 Heart_E10.5 Heart_E11.5 Heart_E12.5 Heart_E13.5 Heart_E14.5 Heart_E15.5 Heart_E16.5 Heart_P0
            do
                echo -ne "\t${stage}" >> title.txt
                bash ${shPATH}/averageSignalInRegionPrimaryOrder.sh \
                    ${dirPATH}/CHMOrganization/Universal_specific/ordered.bed \
                    ${stage}.H3K9me3.rmDup.bw \
                    1 \
                    ${stage}.H3K9me3_${stage}.txt
                paste ordered.H3K9me3.txt ${stage}.H3K9me3_${stage}.txt > tmp && mv tmp ordered.H3K9me3.txt
            done # for stage end
            echo -ne "\n" >> title.txt
            cat title.txt ordered.H3K9me3.txt > ordered.H3K9me3.final.txt
        }
    }
    scan_signal
}
signal_boxplot