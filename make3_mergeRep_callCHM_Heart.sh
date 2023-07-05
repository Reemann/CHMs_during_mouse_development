#!/bin/bash

###### ---------- merge replicates ----------
merge(){
    ### 1. WGBS - merge SamG bed (done)
    mkdir -p ${HOME}/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/PreparedBeforeCallCHM
    cd ${HOME}/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/PreparedBeforeCallCHM
    for stage in 'E10.5' 'E11.5' 'E12.5' 'E13.5' 'E14.5' 'E15.5' 'E16.5' 'P0';do
        cat ${HOME}/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/${stage}.WGBS.rep[1-4]/mSuite.G.bed | grep -v "#chrom" > ${stage}.WGBS.sam.G.bed
    done

    ### 2. H3K9me3 - merge K9 bam; rm dup (done)
    cd ${HOME}/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/PreparedBeforeCallCHM
    samtools view -H ${HOME}/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/E11.5.H3K9me3.rep1/E11.5.H3K9me3.rep1.bam > header.sam
    for stage in 'E10.5' 'E11.5' 'E12.5' 'E13.5' 'E14.5' 'E15.5' 'E16.5' 'P0';do
        samtools merge -@ 30 -h header.sam ${stage}.H3K9me3.bam ${HOME}/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/${stage}.H3K9me3.rep[12]/${stage}.H3K9me3.rep[12].bam
        samtools rmdup ${stage}.H3K9me3.bam ${stage}.H3K9me3.rmDup.bam &
    done

    cd ${HOME}/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/PreparedBeforeCallCHM
    samtools view -H ${HOME}/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/E12.5.input.rep1/E12.5.input.rep1.bam > header_input.sam
    # for stage in 'E10.5' 'E11.5' 'E12.5' 'E13.5' 'E14.5' 'E15.5' 'E16.5' 'P0';do
    # for stage in 'E10.5' 'E12.5' 'E13.5' 'E14.5' 'E15.5' 'E16.5' 'P0';do
    for stage in 'E11.5';do
        samtools merge -@ 30 -h header_input.sam ${stage}.input.bam ${HOME}/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/${stage}.input.rep[12]/${stage}.input.rep[12].bam
        samtools rmdup ${stage}.input.bam ${stage}.input.rmDup.bam &
    done
}
# merge

###### ---------- call CHM ----------
call_CHM(){
    mkdir -p ${HOME}/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/CHMs
    cd ${HOME}/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/CHMs
    # for stage in 'E10.5' 'E11.5' 'E12.5' 'E13.5' 'E14.5' 'E15.5' 'E16.5' 'P0';do
    # for stage in 'E15.5';do
    for stage in 'E10.5' 'E11.5' 'E12.5' 'E13.5' 'E14.5' 'E16.5' 'P0';do
        pcar -m callchm \
        -H ${HOME}/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/PreparedBeforeCallCHM/${stage}.H3K9me3.rmDup.bam \
        -C ${HOME}/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/PreparedBeforeCallCHM/${stage}.input.rmDup.bam \
        -M ${HOME}/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/PreparedBeforeCallCHM/${stage}.WGBS.sam.G.bed \
        -Z ${HOME}/annotations/mm10/mm10_euch.chrom.sizes \
        -Q ${HOME}/annotations/mm10/mm10.2bit \
        -B 200 \
        -N ${stage}_200bp \
        -T 1 > ${stage}_pcar.log
    done
}
# call_CHM



###### ---------- generate bw ----------
gen_bw(){
    mkdir -p ${HOME}/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/PreparedBeforeCallCHM
    cd ${HOME}/CHMsInOtherContexts/CellStateTransition/HeartDevelopment/PreparedBeforeCallCHM

    WGBS(){
        ### 1. WGBS (done)
        merge(){
            for stage in 'E10.5' 'E11.5' 'E12.5' 'E13.5' 'E14.5' 'E15.5' 'E16.5' 'P0';do
                ln -sf ${stage}.WGBS.sam.G.bed ${stage}.methyl.sam.G.bed
                ### (abandoned because of 'bedGraphToBigWig' can not operate on bed with Overlapping regions)
                cut -f 1-4 ${stage}.methyl.sam.G.bed | grep -v '#' | sort -k1,1 -k2,2n | bedtools merge -i - -c 4 -o mean > ${stage}.methyl.sam.G.merged.bed
                bedGraphToBigWig ${stage}.methyl.sam.G.merged.bed ${HOME}/annotations/mm10/mm10.chrom.sizes ${stage}.methyl.merged.bw
            done
        }
        merge

        reps(){
            for stage in 'E10.5' 'E11.5' 'E12.5' 'E13.5' 'E14.5' 'E15.5' 'E16.5' 'P0';do
                for rep in 1 2 3 4;do
                    name=${stage}.WGBS.rep${rep}
                    ln -sf ../${name}/${name}.sam.G.bed .
                    grep -v '#' ${name}.sam.G.bed | cut -f 1-4 | sort -k1,1 -k2,2n > ${name}.sam.G.bed.tmp
                    bedGraphToBigWig ${name}.sam.G.bed.tmp ${HOME}/annotations/mm10/mm10.chrom.sizes ${name}.methyl.bw
                done
            done
            }
        # reps
    }
    WGBS

    k9(){
        ### 2. H3K9me3 (done)
        for stage in 'E10.5' 'E11.5' 'E12.5' 'E13.5' 'E14.5' 'E15.5' 'E16.5' 'P0';do
            macs2 callpeak -g mm -n ${stage}.H3K9me3.rmDup -B -q 0.05 --nomodel --shift=73 --SPMR --broad -t ${stage}.H3K9me3.rmDup.bam && \
            intersectBed -a ${stage}.H3K9me3.rmDup_treat_pileup.bdg -b ${HOME}/annotations/mm10/mm10_main.chrom.limits -wa -f 1.00 | sort -k1,1 -k2,2n > ${stage}.H3K9me3.rmDup_treat_pileup.bdg.tmp && \
            bedGraphToBigWig ${stage}.H3K9me3.rmDup_treat_pileup.bdg.tmp ${HOME}/annotations/mm10/mm10_main.chrom.sizes ${stage}.H3K9me3.rmDup.bw && \
            rm ${stage}.H3K9me3.rmDup_treat_pileup.bdg.tmp &
        done

        for stage in 'E10.5' 'E11.5' 'E12.5' 'E13.5' 'E14.5' 'E15.5' 'E16.5' 'P0';do
            macs2 callpeak -g mm -n ${stage}.input.rmDup -B -q 0.05 --nomodel --shift=73 --SPMR --broad -t ${stage}.input.rmDup.bam && \
            intersectBed -a ${stage}.input.rmDup_treat_pileup.bdg -b ${HOME}/annotations/mm10/mm10_main.chrom.limits -wa -f 1.00 | sort -k1,1 -k2,2n > ${stage}.input.rmDup_treat_pileup.bdg.tmp && \
            bedGraphToBigWig ${stage}.input.rmDup_treat_pileup.bdg.tmp ${HOME}/annotations/mm10/mm10_main.chrom.sizes ${stage}.input.rmDup.bw && \
            rm ${stage}.input.rmDup_treat_pileup.bdg.tmp &
        done
    }
    # k9
}
gen_bw
