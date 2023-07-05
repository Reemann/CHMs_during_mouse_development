#!/bin/bash

###### ---------- merge replicates ----------
merge(){
    ### 1. WGBS
    mkdir -p ${HOME}/CHMsInOtherContexts/CellStateTransition/LiverDevelopment/PreparedBeforeCallCHM
    cd ${HOME}/CHMsInOtherContexts/CellStateTransition/LiverDevelopment/PreparedBeforeCallCHM
    ### merge SamG bed
    for stage in 'E11.5' 'E12.5' 'E13.5' 'E14.5' 'E15.5' 'E16.5' 'P0';do
        cat ${HOME}/CHMsInOtherContexts/CellStateTransition/LiverDevelopment/${stage}.WGBS.rep[1-6]/mSuite.G.bed | grep -v "#chrom" > ${stage}.WGBS.sam.G.bed
    done
    for stage in 'E11.5' 'E12.5' 'E13.5' 'E14.5' 'E15.5' 'E16.5' 'P0';do
        ln -s ${stage}.WGBS.sam.G.bed ${stage}.methyl.sam.G.bed
    done

    ### 2. H3K9me3
    cd ${HOME}/CHMsInOtherContexts/CellStateTransition/LiverDevelopment/PreparedBeforeCallCHM
    ### merge K9 bam; rm dup
    samtools view -H ${HOME}/CHMsInOtherContexts/CellStateTransition/LiverDevelopment/E11.5.H3K9me3.rep1/E11.5.H3K9me3.rep1.bam > header.sam
    for stage in 'E10.5' 'E11.5' 'E12.5' 'E13.5' 'E14.5' 'E15.5' 'E16.5' 'P0';do
        samtools merge -@ 5 -h header.sam ${stage}.H3K9me3.bam ${HOME}/CHMsInOtherContexts/CellStateTransition/LiverDevelopment/${stage}.H3K9me3.rep[12]/${stage}.H3K9me3.rep[12].bam
        samtools rmdup ${stage}.H3K9me3.bam ${stage}.H3K9me3.rmDup.bam &
    done
}
# merge

###### ---------- call CHM ----------
call_CHM(){
    ### 1. liver
    mkdir -p ${HOME}/CHMsInOtherContexts/CellStateTransition/LiverDevelopment/CHMs
    cd ${HOME}/CHMsInOtherContexts/CellStateTransition/LiverDevelopment/CHMs
    # for stage in 'E10.5' 'E11.5' 'E12.5' 'E13.5' 'E14.5' 'E15.5' 'E16.5' 'P0';do
    # for stage in 'E12.5' 'E16.5';do
    # for stage in 'E11.5';do
    # for stage in 'E15.5';do
    # for stage in 'E14.5';do
    for stage in 'E13.5' 'P0';do
        unset DISPLAY
        pcar -m callchm \
        -H ${HOME}/CHMsInOtherContexts/CellStateTransition/LiverDevelopment/PreparedBeforeCallCHM/${stage}.H3K9me3.rmDup.bam \
        -M ${HOME}/CHMsInOtherContexts/CellStateTransition/LiverDevelopment/PreparedBeforeCallCHM/${stage}.WGBS.sam.G.bed \
        -Z ${HOME}/../yanghui/annotations/mm10/mm10_euch.chrom.sizes \
        -Q ${HOME}/../yanghui/annotations/mm10/mm10.2bit \
        -B 200 \
        -N ${stage}_200bp \
        -T 10 > ${stage}_pcar.log 2>&1
    done
}
# call_CHM


###### ---------- generate bw ----------
gen_bw(){
    mkdir -p ${HOME}/CHMsInOtherContexts/CellStateTransition/LiverDevelopment/PreparedBeforeCallCHM
    cd ${HOME}/CHMsInOtherContexts/CellStateTransition/LiverDevelopment/PreparedBeforeCallCHM
    WGBS(){
        merge_reps(){
            ### 1. WGBS
            for stage in 'E11.5' 'E12.5' 'E13.5' 'E14.5' 'E15.5' 'E16.5' 'P0';do
                # ln -sf ${stage}.WGBS.sam.G.bed ${stage}.methyl.sam.G.bed
                cut -f 1-4 ${stage}.methyl.sam.G.bed | grep -v '#' | sort -k1,1 -k2,2n | bedtools merge -i - -c 4 -o mean > ${stage}.methyl.sam.G.merged.bed && \
                bedGraphToBigWig ${stage}.methyl.sam.G.merged.bed ${HOME}/../yanghui/annotations/mm10/mm10.chrom.sizes ${stage}.methyl.merged.bw &
            done
            wait;
        }
        merge_reps
        reps(){
            for stage in 'E10.5' 'E11.5' 'E12.5' 'E13.5' 'E14.5' 'E15.5' 'E16.5' 'P0';do
                for rep in 1 2 3 4;do
                    name=${stage}.WGBS.rep${rep}
                    ln -sf ../${name}/${name}.sam.G.bed .
                    grep -v '#' ${name}.sam.G.bed | cut -f 1-4 | sort -k1,1 -k2,2n > ${name}.sam.G.bed.tmp
                    bedGraphToBigWig ${name}.sam.G.bed.tmp ${HOME}/../yanghui/annotations/mm10/mm10.chrom.sizes ${name}.methyl.bw
                done
            done
        }
        # reps
    }
    WGBS

    k9(){
        ### 2. H3K9me3
        for stage in 'E11.5' 'E12.5' 'E13.5' 'E14.5' 'E15.5' 'E16.5' 'P0';do
            macs2 callpeak -g mm -n ${stage}.H3K9me3.rmDup -B -q 0.05 --nomodel --shift=73 --SPMR --broad -t ${stage}.H3K9me3.rmDup.bam && \
            intersectBed -a ${stage}.H3K9me3.rmDup_treat_pileup.bdg -b ${HOME}/../yanghui/annotations/mm10/mm10_main.chrom.limits -wa -f 1.00 | sort -k1,1 -k2,2n > ${stage}.H3K9me3.rmDup_treat_pileup.bdg.tmp && \
            bedGraphToBigWig ${stage}.H3K9me3.rmDup_treat_pileup.bdg.tmp ${HOME}/../yanghui/annotations/mm10/mm10_main.chrom.sizes ${stage}.H3K9me3.rmDup.bw && \
            rm ${stage}.H3K9me3.rmDup_treat_pileup.bdg.tmp &
        done
    }
    # k9
}
gen_bw
