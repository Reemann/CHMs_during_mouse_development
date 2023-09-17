#!/bin/bash
dirPATH=${HOME}/CHMsInOtherContexts/CellStateTransition
anPATH=${HOME}/../yanghui/annotations


gen_region_CpGRiched(){
    cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Features/overlap_IAPEz
    awk 'BEGIN{FS=OFS="\t"}$4>=6{print $1,$2,$3,$4}' ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment/CHMs/b200bp.CpGNumber > b200bp_6CpG.bed
    for group in IAPEz_specific IAPEz_universalCHM_overlap universalCHM_IAPEz_overlap universalCHM_specific;do
        bedtools intersect -u -a b200bp_6CpG.bed -b ./${group}.bed > CpG_frequency/${group}.6CpG200bp.bed
    done
}
# gen_region_CpGRiched

# ---------- CHM ----------
CpG_freq(){
    cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Features/overlap_IAPEz/CpG_frequency

    CHM_200bp(){
        mkdir -p 200bp
        for group in universalCHM_IAPEz_overlap universalCHM_specific;do
            ### 200bp
            get_Kmer_frequency.py "CG" ${group}.6CpG200bp.bed 200bp/${group}_6CpG200bp_CpG.bed ${anPATH}/mm10/mm10.2bit
        done
    }
    # CHM_200bp

    total_region(){
        mkdir -p total_region
        for group in IAPEz_specific IAPEz_universalCHM_overlap universalCHM_IAPEz_overlap universalCHM_specific;do
            get_Kmer_frequency.py "CG" ../${group}.bed total_region/${group}_CpG.bed ${anPATH}/mm10/mm10.2bit
        done
    }
    total_region
}
CpG_freq

region_num(){
    cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Features/overlap_IAPEz/CpG_frequency
    for group in IAPEz_specific IAPEz_universalCHM_overlap universalCHM_IAPEz_overlap universalCHM_specific;do
        awk -v XVALUE=${group} 'BEGIN{FS=OFS="\t";i=0}{i+=1;print $1,$2,$3,XVALUE"-"i}' ../../${group}.bed > ../${group}.RegionNum.bed
    done
}
# region_num

assign_200bp_to_region(){
    cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Features/overlap_IAPEz/CpG_frequency/200bp
    for group in universalCHM_IAPEz_overlap universalCHM_specific;do
        bedtools intersect -b ../${group}.RegionNum.bed -a ${group}_6CpG200bp_CpG.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4,$5,$6,$15}' > ${group}_6CpG200bp_CpG.RegionNum.bed
    done
}
# assign_200bp_to_region
