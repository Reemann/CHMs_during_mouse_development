
###### ---------- sample correlation ----------
dirPATH=${HOME}/CHMsInOtherContexts/CellStateTransition
dirPATH_Y=/mnt/Storage/home/yanghui/imprinting/result.2021/Exploring/CHMsInOtherContexts/CellStateTransition
# shPATH_Y=/mnt/Storage/home/yanghui/scripts/shell
shPATH=/mnt/Storage/home/wangyiman/bin/utilities
# rPATH_Y=/mnt/Storage/home/yanghui/scripts/R
rPATH=/mnt/Storage/home/wangyiman/bin/utilities
anPATH=/mnt/Storage/home/yanghui/annotations

# mkdir -p ${HOME}/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/PreparedBeforeCallCHM
# cd ${HOME}/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/PreparedBeforeCallCHM
# scp -P 6666 wangyiman@10.10.196.114:~/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/PreparedBeforeCallCHM/*.methyl.sam.G.bed .
# scp -P 6666 wangyiman@10.10.196.114:~/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/PreparedBeforeCallCHM/*.H3K9me3.rmDup.bw .

process_cor(){
    name=${1} # Spermatogenesis_GSE137744
    CHM_STAGE=${2} # "US DS PS RS"
    CHM_STAGE_plot=${3} # "US;DS;PS;RS"

    cd ${dirPATH}/${name}/PreparedBeforeCallCHM
    basic_txt(){
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
                    "1"

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
        }
    # basic_txt

    plot_process(){
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
    # plot_process
}

# process_cor EarlyEmbryogenesis "2cell 8cell Morula ICM" "2cell;8cell;Morula;ICM"
# process_cor PGCsDevelopment "E9.5 E10.5 E13.5_female E13.5_male" "E9.5;E10.5;E13.5_female;E13.5_male"
# process_cor HeartDevelopment "E10.5 E11.5 E12.5 E13.5 E14.5 E15.5 E16.5 P0" "E10.5;E11.5;E12.5;E13.5;E14.5;E15.5;E16.5;P0"
# process_cor LiverDevelopment "E11.5 E12.5 E13.5 E14.5 E15.5 E16.5 P0" "E11.5;E12.5;E13.5;E14.5;E15.5;E16.5;P0"

# cd ~/CHMsInOtherContexts/figures
# ln -sf ~/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/PreparedBeforeCallCHM/cor.ggplot_barplot.pdf make4_barplot_methK9_corr_heart.pdf


###### ---------- process correlation ----------

merge_into_one(){
    # mkdir -p ${dirPATH}/Prepare
    cd ${dirPATH}/Prepare
    echo -e "xValue\tyValue\tgroup" > cor.ggplot_basic.txt
    tail -n +2 ${dirPATH}/EarlyEmbryogenesis/PreparedBeforeCallCHM/cor.ggplot_basic.txt >> cor.ggplot_basic.txt
    tail -n +2 ${dirPATH_Y}/PGCsDevelopment/Merged/cor.ggplot_basic.txt >> cor.ggplot_basic.txt
    tail -n +2 ${dirPATH_Y}/Spermatogenesis_GSE137744/Merged/cor.ggplot_basic.txt >> cor.ggplot_basic.txt
    tail -n +2 ${dirPATH_Y}/RetinalDevelopment/GSE87064_Mouse/Merged/cor.ggplot_basic.txt >> cor.ggplot_basic.txt
    tail -n +2 ${dirPATH}/HeartDevelopment/PreparedBeforeCallCHM/cor.ggplot_basic.txt | awk '{FS=OFS="\t"}{print "Heart_"$0}' >> cor.ggplot_basic.txt
    tail -n +2 ${dirPATH}/LiverDevelopment/PreparedBeforeCallCHM/cor.ggplot_basic.txt | awk '{FS=OFS="\t"}{print "Liver_"$0}' >> cor.ggplot_basic.txt

    Rscript ${rPATH}/ggplot_barplot.r \
        "cor.ggplot_basic.txt" \
        "cor.ggplot_barplot.pdf" \
        "Biological process" \
        "" \
        "Pearson's correlation coefficient__n(DNA methylation and H3K9me3)" \
        "-0.2;1;0.2" \
        "2cell;8cell;Morula;ICM;E10.5;E13.5_female;E13.5_male;US;DS;PS;RS;E14.5;E17.5;P0;P3;P7;P10;P14;P21;Heart_E10.5;Heart_E11.5;Heart_E12.5;Heart_E13.5;Heart_E14.5;Heart_E15.5;Heart_E16.5;Heart_P0;Liver_E11.5;Liver_E12.5;Liver_E13.5;Liver_E14.5;Liver_E15.5;Liver_E16.5;Liver_P0" \
        "low;intermediate;high" \
        "dodge" \
        "8" \
        "3.5" \
        "-007" \
        "F"
}
merge_into_one

cd ~/CHMsInOtherContexts/figures
ln -sf ~/CHMsInOtherContexts/CellStateTransition/Prepare/cor.ggplot_barplot.pdf make4_barplot_methK9_corr_all.pdf

