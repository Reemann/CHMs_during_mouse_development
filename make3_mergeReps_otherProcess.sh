#!/bin/bash

gen_bw(){
    Process=${1}
    CHM_STAGE=${2}
    cd ${HOME}/CHMsInOtherContexts/CellStateTransition/${Process}/PreparedBeforeCallCHM
    ### 1. WGBS
    for stage in ${CHM_STAGE};do
        # ln -sf ${stage}.WGBS.sam.G.bed ${stage}.methyl.sam.G.bed
        cut -f 1-4 ${stage}.methyl.sam.G.bed | grep -v '#' | sort -k1,1 -k2,2n | bedtools merge -i - -c 4 -o mean > ${stage}.methyl.sam.G.merged.bed && \
        bedGraphToBigWig ${stage}.methyl.sam.G.merged.bed ${HOME}/../yanghui/annotations/mm10/mm10.chrom.sizes ${stage}.methyl.merged.bw &
    done
    wait;
}

gen_bw Spermatogenesis    "US DS RS PS"
gen_bw RetinalDevelopment "E14.5 E17.5 P0 P3 P7 P10 P14 P21"
gen_bw EarlyEmbryogenesis "2cell 8cell Morula ICM"
gen_bw PGCsDevelopment    "E10.5 E13.5_female E13.5_male"
