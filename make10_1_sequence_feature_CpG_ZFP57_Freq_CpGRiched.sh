#!/bin/bash
dirPATH=${HOME}/CHMsInOtherContexts/CellStateTransition
anPATH=${HOME}/../yanghui/annotations

gen_nonCHM_old_version(){
    cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Sequence
    cat ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/*/merged.CHM.bed | sort -k1,1 -k2,2n > union_total_sorted.CHM.bed
    mergeBed -i union_total_sorted.CHM.bed > union_total_sorted_merged.CHM.bed
    bedtools subtract -a ~/CHMsInOtherContexts/CellStateTransition/Prepare/CpGs_high.bed -b union_total_sorted_merged.CHM.bed > nonCHM_CpG_riched.bed
    bedtools intersect -a ~/CHMsInOtherContexts/CellStateTransition/Prepare/CpGs_high.bed -b union_total_sorted_merged.CHM.bed -v > nonCHM_CpG_riched.1kb.bed
    for process in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific HeartSpecific LiverSpecific;do
        intersectBed -u -a ${anPATH}/mm10/Bins/mm10.b1kb.euchr.bed -b ${dirPATH}/CHMOrganization/Universal_specific/${process}.CHM.bed > ${process}.CHM.1kb.bed
    done

    # bedtools window -w 2000 -b ~/CHMsInOtherContexts/CellStateTransition/Prepare/CpGs_high.bed -a union_total_sorted_merged.CHM.bed -u | wc -l
}
# gen_nonCHM_old_version ##### abondoned


gen_nonCHM(){
    dirPATH=/mnt/Storage/home/yanghui/imprinting/result.2021/Exploring/CHMsInOtherContexts/CellStateTransition/CHMOrganization
    cd ${dirPATH}
    cat EarlyEmbryogenesis/*.CHM.bed \
    PGC_Development/E*.CHM.bed \
    Spermatogenesis/*S.CHM.bed \
    RetinalDevelopment/E*.CHM.bed \
    RetinalDevelopment/P*.CHM.bed \
    HeartDevelopment/*.CHM.bed \
    LiverDevelopment/*.CHM.bed | windowBed -w 2000 -v -a ~/CHMsInOtherContexts/CellStateTransition/Prepare/CpGs_high.bed -b - | sort -k1,1 -k2,2n -u > ${HOME}/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Sequence/NonCHMsCpGrich.bed # 22,069

}
# gen_nonCHM #####


gen_CHM_CpGRiched_1kb_universal(){
    cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific
    # cat *Specific.CHM.bed | sort -k1,1 -k2,2n -u | mergeBed -d 2000 -i - | sort -k1,1 -k2,2n -u > ProcessSpecific.CHM.bed
    cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Sequence
    for process in Universal;do
        bedtools intersect -u -f 1 -a ~/CHMsInOtherContexts/CellStateTransition/Prepare/CpGs_high.bed -b ${dirPATH}/CHMOrganization/Universal_specific/${process}.CHM.bed > ${process}.CHM.30CpG1kb.bed
    done
}
# gen_CHM_CpGRiched_1kb_universal
### to test the output file : bedtools intersect -a Universal.CHM.30CpG1kb.bed -b Universal.CHMnum.bed -wo | cut -f 8 | uniq

gen_CHM_CpGRiched_1kb_processComplement(){
    cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Sequence
    for process in EarlyEmbryogenesis PGCsDevelopment Spermatogenesis RetinalDevelopment HeartDevelopment LiverDevelopment;do
        bedtools intersect -u -f 1 -a ~/CHMsInOtherContexts/CellStateTransition/Prepare/CpGs_high.bed -b ${dirPATH}/CHMOrganization/Overlap/${process}_complementarySet.CHM.bed > ${process}.ComplementCHM.30CpG1kb.bed
    done
    bedtools intersect -u -f 1 -a ~/CHMsInOtherContexts/CellStateTransition/Prepare/CpGs_high.bed -b ${dirPATH}/CHMOrganization/Universal_specific/Universal_complementarySet.CHM.bed > Universal.ComplementCHM.30CpG1kb.bed
}
gen_CHM_CpGRiched_1kb_processComplement


# ---------- CHM ----------
CpG_freq(){
    universal_1kb(){
        cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Sequence/CpGFrequency_1kb
        for process in Universal;do 
            get_Kmer_frequency.py "CG" ${dirPATH}/CHMOrganization/Universal_specific/Sequence/${process}.CHM.30CpG1kb.bed ${process}_CHM30CpG1kb_CpG.bed ${anPATH}/mm10/mm10.2bit
            ### 1kb
        done
    }
    # universal_1kb

    processComplement(){
        cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Sequence/CpGFrequency_1kb
        for process in EarlyEmbryogenesis PGCsDevelopment Spermatogenesis RetinalDevelopment HeartDevelopment LiverDevelopment Universal;do 
            get_Kmer_frequency.py "CG" ${dirPATH}/CHMOrganization/Universal_specific/Sequence/${process}.ComplementCHM.30CpG1kb.bed ${process}Complement_CHM30CpG1kb_CpG.bed ${anPATH}/mm10/mm10.2bit
            ### 1kb
        done
    }
    processComplement

    NonCHMsCpGrich(){
        # get_Kmer_frequency.py "CG" ${dirPATH}/CHMOrganization/Universal_specific/Sequence/nonCHM_CpG_riched.1kb.bed nonCHM_CpG_riched1kb_CHM_CpG.bed ${anPATH}/mm10/mm10.2bit ##### old_version
        get_Kmer_frequency.py "CG" ${dirPATH}/CHMOrganization/Universal_specific/Sequence/NonCHMsCpGrich.bed nonCHM_CpG_riched1kb_CHM_CpG.bed ${anPATH}/mm10/mm10.2bit
    }
    # NonCHMsCpGrich
}
CpG_freq


ZFP57_freq(){
    universal_1kb(){
        mkdir -p ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Sequence/ZFP57Frequency_1kb
        cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Sequence/ZFP57Frequency_1kb
        for process in Universal;do
            get_Kmer_frequency.py "TGCCGC" ${dirPATH}/CHMOrganization/Universal_specific/Sequence/${process}.CHM.30CpG1kb.bed ${process}_CHM30CpG1kb_ZFP57.bed ${anPATH}/mm10/mm10.2bit
            ### 1kb
        done
    }
    # universal_1kb 

    processComplement(){
        cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Sequence/ZFP57Frequency_1kb
        for process in EarlyEmbryogenesis PGCsDevelopment Spermatogenesis RetinalDevelopment HeartDevelopment LiverDevelopment Universal;do 
            get_Kmer_frequency.py "TGCCGC" ${dirPATH}/CHMOrganization/Universal_specific/Sequence/${process}.ComplementCHM.30CpG1kb.bed ${process}Complement_CHM30CpG1kb_ZFP57.bed ${anPATH}/mm10/mm10.2bit
            ### 1kb
        done
    }
    processComplement

    NonCHMsCpGrich(){
        # get_Kmer_frequency.py "TGCCGC" ${dirPATH}/CHMOrganization/Universal_specific/Sequence/nonCHM_CpG_riched.1kb.bed nonCHM_CpG_riched1kb_CHM_ZFP57.bed ${anPATH}/mm10/mm10.2bit ##### old_version
        get_Kmer_frequency.py "TGCCGC" ${dirPATH}/CHMOrganization/Universal_specific/Sequence/NonCHMsCpGrich.bed nonCHM_CpG_riched1kb_CHM_ZFP57.bed ${anPATH}/mm10/mm10.2bit
    }
    # NonCHMsCpGrich

}
ZFP57_freq


assign_1kb_to_CHM(){
    kmer_freq(){
        cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Sequence/KmerFrequency_merged
        for process in Universal;do
            # awk -v XVALUE=${process} 'BEGIN{FS=OFS="\t";i=0}{i+=1;print $1,$2,$3,XVALUE"_"i}' ${dirPATH}/CHMOrganization/Universal_specific/${process}.CHM.bed > ../${process}.CHMnum.bed
            bedtools intersect -b ../${process}.CHMnum.bed -a ${process}_CHM30CpG1kb_kmer6.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4,$5,$6,$10}' > ${process}_CHM30CpG1kb_kmer6.CHMnum.bed
        done
    }
    # kmer_freq #### xxx_kmer6.bed are run on Zhanglab3

    CpG_freq_universal(){
        cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Sequence/CpGFrequency_1kb
        for process in Universal;do
            awk -v XVALUE=${process} 'BEGIN{FS=OFS="\t";i=0}{i+=1;print $1,$2,$3,XVALUE"_"i}' ${dirPATH}/CHMOrganization/Universal_specific/${process}.CHM.bed > ../${process}.CHMnum.bed
            bedtools intersect -b ../${process}.CHMnum.bed -a ${process}_CHM30CpG1kb_CpG.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4,$5,$6,$15}' > ${process}_CHM30CpG1kb_CpG.CHMnum.bed
        done
    }
    # CpG_freq_universal

    CpG_freq_processComplement(){
        cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Sequence/CpGFrequency_1kb
        for process in EarlyEmbryogenesis PGCsDevelopment Spermatogenesis RetinalDevelopment HeartDevelopment LiverDevelopment;do 
            awk -v XVALUE=${process}Complement 'BEGIN{FS=OFS="\t";i=0}{i+=1;print $1,$2,$3,XVALUE"_"i}' ${dirPATH}/CHMOrganization/Overlap/${process}_complementarySet.CHM.bed > ../${process}Complement.CHMnum.bed
            bedtools intersect -b ../${process}Complement.CHMnum.bed -a ${process}Complement_CHM30CpG1kb_CpG.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4,$5,$6,$15}' > ${process}Complement_CHM30CpG1kb_CpG.CHMnum.bed
        done

        for process in Universal;do 
            awk -v XVALUE=${process}Complement 'BEGIN{FS=OFS="\t";i=0}{i+=1;print $1,$2,$3,XVALUE"_"i}' ${dirPATH}/CHMOrganization/Universal_specific/${process}_complementarySet.CHM.bed > ../${process}Complement.CHMnum.bed
            bedtools intersect -b ../${process}Complement.CHMnum.bed -a ${process}Complement_CHM30CpG1kb_CpG.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4,$5,$6,$15}' > ${process}Complement_CHM30CpG1kb_CpG.CHMnum.bed
        done
    }
    CpG_freq_processComplement

    ZFP57_freq_universal(){
        cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Sequence/ZFP57Frequency_1kb
        for process in Universal;do
            # awk -v XVALUE=${process} 'BEGIN{FS=OFS="\t";i=0}{i+=1;print $1,$2,$3,XVALUE"_"i}' ${dirPATH}/CHMOrganization/Universal_specific/${process}.CHM.bed > ../${process}.CHMnum.bed
            bedtools intersect -b ../${process}.CHMnum.bed -a ${process}_CHM30CpG1kb_ZFP57.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4,$5,$6,$15}' >  ${process}_CHM30CpG1kb_ZFP57.CHMnum.bed
        done
    }
    # ZFP57_freq_universal

    ZFP57_freq_processComplement(){
        cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Sequence/ZFP57Frequency_1kb
        for process in EarlyEmbryogenesis PGCsDevelopment Spermatogenesis RetinalDevelopment HeartDevelopment LiverDevelopment Universal;do 
            bedtools intersect -b ../${process}Complement.CHMnum.bed -a ${process}Complement_CHM30CpG1kb_ZFP57.bed -wa -wb | awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4,$5,$6,$15}' >  ${process}Complement_CHM30CpG1kb_ZFP57.CHMnum.bed
        done
    }
    ZFP57_freq_processComplement    
}
assign_1kb_to_CHM 

