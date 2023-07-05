#!/bin/bash

dirPATH=/mnt/Storage/home/wangyiman/CHMsInOtherContexts/CellStateTransition
pyPATH=/mnt/Storage/home/yanghui/scripts/python
shPATH=/mnt/Storage/home/wangyiman/bin/utilities
rPATH=/mnt/Storage/home/yanghui/scripts/R
rscriptPATH=/mnt/Storage/home/yanghui/software/R-3.5.2/bin/Rscript
anPATH=/mnt/Storage/home/yanghui/annotations

DNAMethylationAmount(){
    mkdir -p ${dirPATH}/DNAMethylationAmount
    cd ${dirPATH}/DNAMethylationAmount
    PROCESS=${1}
    STAGE=${2}
    for stage in ${STAGE};do
        # step1. Obtain DNA methylation level in 1kb bins with 10bp sliding window
        cal_meth_level(){
            bash ${shPATH}/averageMethylInRegion_splitChrom.sh \
                "${anPATH}/mm10/Bins/mm10.b1kb.s10bp.euchr.bed" \
                "${dirPATH}/${PROCESS}/PreparedBeforeCallCHM/${stage}.methyl.sam.G.bed" \
                "${PROCESS}.${stage}.mm10.1kbB_10bpS.methyl_amount.txt" \
                "mm10"
        }
        cal_meth_level

        # step2. Calculate DNA methylation amount to the center
        meth2center(){
            [[ ! -e mm10.b1kb.s10bp.euchr.CpGnumber.sorted ]] && sort -k1,1 -k2,2n ${anPATH}/mm10/Bins/mm10.b1kb.s10bp.euchr.CpGnumber > mm10.b1kb.s10bp.euchr.CpGnumber.sorted
            paste ${PROCESS}.${stage}.mm10.1kbB_10bpS.methyl_amount.txt ${anPATH}/mm10/Bins/mm10.b1kb.s10bp.euchr.CpGnumber | grep -i -v na | \
                awk 'BEGIN{FS=OFS="\t"}{print $1,($2+$3)/2,($2+$3)/2+1,$4*$8}' > ${PROCESS}.${stage}.mm10.1kbB_10bpS.methyl_amount.bdg
        }
        meth2center

        # step3. Generate bw files
        gen_bw(){
            bedGraphToBigWig ${PROCESS}.${stage}.mm10.1kbB_10bpS.methyl_amount.bdg ${anPATH}/mm10/mm10_euch.chrom.sizes ${PROCESS}.${stage}.mm10.1kbB_10bpS.methyl_amount.bw
        }
        gen_bw

        # step4. gunzip temporatory files
        # rm ${PROCESS}.${stage}.mm10.1kbB_10bpS.methyl_amount.bdg
    done
}



# DNAMethylationAmount EarlyEmbryogenesis "2cell 8cell Morula ICM"
DNAMethylationAmount EarlyEmbryogenesis "2cell 8cell"
# DNAMethylationAmount PGCsDevelopment    "E10.5 E13.5_female E13.5_male"
# DNAMethylationAmount PGCsDevelopment    "E13.5_female"
# DNAMethylationAmount Spermatogenesis    "US DS RS PS"
# DNAMethylationAmount RetinalDevelopment "E14.5 E17.5 P0 P3 P7 P10 P14 P21"
# DNAMethylationAmount HeartDevelopment   "E10.5 E11.5 E12.5 E13.5 E14.5 E15.5 E16.5 P0"
# DNAMethylationAmount LiverDevelopment   "E11.5 E12.5 E13.5 E14.5 E15.5 E16.5 P0"
