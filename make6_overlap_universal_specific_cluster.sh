#!/bin/bash

dirPATH=${HOME}/CHMsInOtherContexts/CellStateTransition
shPATH=/mnt/Storage/home/yanghui/scripts/shell
rPATH=/mnt/Storage/home/yanghui/scripts/R
dataPATH_CHM_earlyEmbryo=/mnt/Storage/home/yanghui/imprinting/result.2021/Filtering/CpGrichedHighMethylH3K9me3/Hybrid/ChromHMM/Manage
anPATH=/mnt/Storage/home/yanghui/annotations


cd ${dirPATH}/CHMOrganization

overlap(){
    # D: Apr-15-2022 22:52 Fri
    mkdir -p ${dirPATH}/CHMOrganization/Overlap;cd ${dirPATH}/CHMOrganization/Overlap
    prepare(){
        for C_class in CHM CHnonM CMnonH;do
            cd ${dirPATH}/CHMOrganization/Overlap
            # cat ${dirPATH}/CHMOrganization/*/stable.filteredBySignal.bed | sort -k1,1 -k2,2n | mergeBed -i - -d 2000 > Union.bed
            ### !!! modified to be the same with function `Cluster`
            # cat ${dirPATH}/CHMOrganization/*/stable.${C_class}s.bed | cut -f 1-3 | sort -k1,1 -k2,2n | mergeBed -i - -d 2000 > Union.${C_class}.bed
            cat ${dirPATH}/CHMOrganization/*/stable.${C_class}s.bed | cut -f 1-3 | sort -k1,1 -k2,2n | mergeBed -i - -d 2000 -c 2,3 -o median,median | cut -f 1,4-5 > Union.${C_class}.bed
            for process in EarlyEmbryogenesis PGCsDevelopment Spermatogenesis RetinalDevelopment HeartDevelopment LiverDevelopment
            do
                # intersectBed -u -a Union.bed -b ${dirPATH}/CHMOrganization/${process}/stable.filteredBySignal.bed | awk 'BEGIN{FS="\t";OFS="_"}{print $1, $2, $3;}' > ${process}.bed
                intersectBed -u -a Union.${C_class}.bed -b ${dirPATH}/CHMOrganization/${process}/stable.${C_class}s.bed | awk 'BEGIN{FS="\t";OFS="_"}{print $1, $2, $3;}' > ${process}.${C_class}.txt
            done # for process end
        done
    }
    prepare

    Rscript ${HOME}/CHMsInOtherContexts/bin/make6_UpSetR.r
    # Warning message:
    #     In Make_base_plot(Main_bar_plot = x$Main_bar, Matrix_plot = x$Matrix,:
    #     Plot might be out of range if ratio > 0.7 or < 0.3
}
# overlap

### 2023-04-27
overlap_process_specific_complementarySet(){
    cd ${dirPATH}/CHMOrganization/Overlap;cd ${dirPATH}/CHMOrganization/Overlap
    for process in EarlyEmbryogenesis PGCsDevelopment Spermatogenesis RetinalDevelopment HeartDevelopment LiverDevelopment;do
        # bedtools intersect -a Union.CHM.bed -b ${dirPATH}/CHMOrganization/${process}/stable.CHMs.bed -v | sort -k1,1 -k2,2n | mergeBed -d 2000 -i - > ${process}_complementarySet.CHM.bed
        awk 'BEGIN{FS="_";OFS="\t"}{print $1,$2,$3}' ${process}.CHM.txt | bedtools intersect -a Union.CHM.bed -b - -v | sort -k1,1 -k2,2n | mergeBed -d 2000 -i - > ${process}_complementarySet.CHM.bed
    done
}
overlap_process_specific_complementarySet

### 2023-04-27
overlap_process_specific_complementarySet_otherCommon(){
    cd ${dirPATH}/CHMOrganization/Overlap;cd ${dirPATH}/CHMOrganization/Overlap
    process=$1
    IFS=' ' read -r -a otherProcess <<< "$2"
    bedtools intersect -a Union.CHM.bed -b ${dirPATH}/CHMOrganization/${process}/stable.CHMs.bed | \
        bedtools intersect -a - -b ${dirPATH}/CHMOrganization/${otherProcess[0]}/stable.CHMs.bed -v |\
        bedtools intersect -a - -b ${dirPATH}/CHMOrganization/${otherProcess[1]}/stable.CHMs.bed -v |\
        bedtools intersect -a - -b ${dirPATH}/CHMOrganization/${otherProcess[2]}/stable.CHMs.bed -v |\
        bedtools intersect -a - -b ${dirPATH}/CHMOrganization/${otherProcess[3]}/stable.CHMs.bed -v |\
        bedtools intersect -a - -b ${dirPATH}/CHMOrganization/${otherProcess[4]}/stable.CHMs.bed -v |\
        sort -k1,1 -k2,2n | mergeBed -d 2000 -i - > ${process}_complementarySet_otherCommon.CHM.bed
}
# overlap_process_specific_complementarySet_otherCommon EarlyEmbryogenesis "PGCsDevelopment Spermatogenesis RetinalDevelopment HeartDevelopment LiverDevelopment"
# overlap_process_specific_complementarySet_otherCommon PGCsDevelopment "EarlyEmbryogenesis Spermatogenesis RetinalDevelopment HeartDevelopment LiverDevelopment"
# overlap_process_specific_complementarySet_otherCommon Spermatogenesis "EarlyEmbryogenesis PGCsDevelopment RetinalDevelopment HeartDevelopment LiverDevelopment"
# overlap_process_specific_complementarySet_otherCommon RetinalDevelopment "EarlyEmbryogenesis PGCsDevelopment Spermatogenesis HeartDevelopment LiverDevelopment"
# overlap_process_specific_complementarySet_otherCommon HeartDevelopment "EarlyEmbryogenesis PGCsDevelopment Spermatogenesis RetinalDevelopment LiverDevelopment"
# overlap_process_specific_complementarySet_otherCommon LiverDevelopment "EarlyEmbryogenesis PGCsDevelopment Spermatogenesis RetinalDevelopment HeartDevelopment"


###### ---------- group ----------

Cluster(){
    # D: Apr-16-2022 09:41 Sat
    mkdir -p ${dirPATH}/CHMOrganization/KMeansCluster;cd ${dirPATH}/CHMOrganization/KMeansCluster
    prepare(){
        # T: executed in 466ms, finished 09:54:24 2022-04-16
        C_class=$1
        # cat ${dirPATH}/CHMOrganization/*/stable.filteredBySignal.bed | cut -f 1-3 | sort -k1,1 -k2,2n | mergeBed -i - -d 2000 -c 2,3 -o median,median | cut -f 1,4-5 > Union_refined.bed # 13,919
        cat ${dirPATH}/CHMOrganization/*/stable.${C_class}s.bed | cut -f 1-3 | sort -k1,1 -k2,2n | mergeBed -i - -d 2000 -c 2,3 -o median,median | cut -f 1,4-5 > Union_refined.${C_class}.bed # 13,919
        echo -ne "#Chrom\tStart\tEnd" > title.txt
        cat Union_refined.${C_class}.bed > Union_refined.${C_class}.existence.txt
        for process in EarlyEmbryogenesis PGCsDevelopment Spermatogenesis RetinalDevelopment HeartDevelopment LiverDevelopment
        do
            echo -ne "\t${process}" >> title.txt
            # intersectBed -c -a Union_refined.${C_class}.bed -b ${dirPATH}/CHMOrganization/${process}/stable.filteredBySignal.bed | cut -f 4 | awk 'BEGIN{FS=OFS="\t"}{if($1==0){print 0;} else{print 1;}}' | paste Union_refined.${C_class}.existence.txt - > tmp && mv tmp Union_refined.${C_class}.existence.txt
            intersectBed -c -a Union_refined.${C_class}.bed -b ${dirPATH}/CHMOrganization/${process}/stable.${C_class}s.bed | cut -f 4 | awk 'BEGIN{FS=OFS="\t"}{if($1==0){print 0;} else{print 1;}}' | paste Union_refined.${C_class}.existence.txt - > tmp && mv tmp Union_refined.${C_class}.existence.txt
        done # for process end
        echo -ne "\n" >> title.txt
        cat title.txt Union_refined.${C_class}.existence.txt > tmp && mv tmp Union_refined.${C_class}.existence.txt
    }
    prepare CHM
    prepare CHnonM
    prepare CMnonH
}
# Cluster ### not used

gen_bed(){
    mkdir -p ${dirPATH}/CHMOrganization/Universal_specific;cd ${dirPATH}/CHMOrganization/Universal_specific
    for C_class in CHM CHnonM CMnonH;do
        awk 'BEGIN{FS=OFS="\t"}{if($4==1 && $5==1 && $6==1 && $7==1 && $8==1 && $9==1){print $1, $2, $3;}}' ${dirPATH}/CHMOrganization/KMeansCluster/Union_refined.${C_class}.existence.txt >           Universal.${C_class}.bed # 2677
        awk 'BEGIN{FS=OFS="\t"}{if($4==1 && $5==0 && $6==0 && $7==0 && $8==0 && $9==0){print $1, $2, $3;}}' ${dirPATH}/CHMOrganization/KMeansCluster/Union_refined.${C_class}.existence.txt > EarlyEmbryoSpecific.${C_class}.bed # 599
        awk 'BEGIN{FS=OFS="\t"}{if($4==0 && $5==1 && $6==0 && $7==0 && $8==0 && $9==0){print $1, $2, $3;}}' ${dirPATH}/CHMOrganization/KMeansCluster/Union_refined.${C_class}.existence.txt >         PGCSpecific.${C_class}.bed # 167
        awk 'BEGIN{FS=OFS="\t"}{if($4==0 && $5==0 && $6==1 && $7==0 && $8==0 && $9==0){print $1, $2, $3;}}' ${dirPATH}/CHMOrganization/KMeansCluster/Union_refined.${C_class}.existence.txt >       SpermSpecific.${C_class}.bed # 2006
        awk 'BEGIN{FS=OFS="\t"}{if($4==0 && $5==0 && $6==0 && $7==1 && $8==0 && $9==0){print $1, $2, $3;}}' ${dirPATH}/CHMOrganization/KMeansCluster/Union_refined.${C_class}.existence.txt >     RetinalSpecific.${C_class}.bed # 3121
        awk 'BEGIN{FS=OFS="\t"}{if($4==0 && $5==0 && $6==0 && $7==0 && $8==1 && $9==0){print $1, $2, $3;}}' ${dirPATH}/CHMOrganization/KMeansCluster/Union_refined.${C_class}.existence.txt >       HeartSpecific.${C_class}.bed # 314
        awk 'BEGIN{FS=OFS="\t"}{if($4==0 && $5==0 && $6==0 && $7==0 && $8==0 && $9==1){print $1, $2, $3;}}' ${dirPATH}/CHMOrganization/KMeansCluster/Union_refined.${C_class}.existence.txt >       LiverSpecific.${C_class}.bed # 74
    done
}
# gen_bed

### 2023-04-27
overlap_universal_complementarySet(){
    cd ${dirPATH}/CHMOrganization/Universal_specific
    bedtools intersect -a ../Overlap/Union.CHM.bed -b Universal.CHM.bed -v > Universal_complementarySet.CHM.bed
}
overlap_universal_complementarySet

rm_class_overlap(){
    cd ${dirPATH}/CHMOrganization/Universal_specific
    for process in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific HeartSpecific LiverSpecific;do
        for C_class in CHnonM CMnonH;do
            bedtools intersect -v -a ${process}.${C_class}.bed -b ${process}.CHM.bed > ${process}.${C_class}_final.bed
        done
    done
}
# rm_class_overlap