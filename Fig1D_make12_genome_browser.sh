#!/bin/bash

mkdir -p ~/CHMsInOtherContexts/CellStateTransition/MethK9Signal
cd ~/CHMsInOtherContexts/CellStateTransition/MethK9Signal

prepare_methAmount(){
    ln -s ~/CHMsInOtherContexts/CellStateTransition/DNAMethylationAmount/*.bw .
}
prepare_methAmount

prepare_signal(){
    PROCESS=${1}
    STAGE=${2}
    for stage in ${STAGE};do
        ln -s ~/CHMsInOtherContexts/CellStateTransition/${PROCESS}/PreparedBeforeCallCHM/${stage}.methyl.merged.bw ${PROCESS}_${stage}.methyl.merged.bw
        ln -s ~/CHMsInOtherContexts/CellStateTransition/${PROCESS}/PreparedBeforeCallCHM/${stage}.H3K9me3.rmDup.bw ${PROCESS}_${stage}.H3K9me3.rmDup.bw
    done
}
# prepare_signal EarlyEmbryogenesis "2cell 8cell Morula ICM"
# prepare_signal PGCsDevelopment    "E10.5 E13.5_female E13.5_male"
# prepare_signal Spermatogenesis    "US DS RS PS"
# prepare_signal RetinalDevelopment "E14.5 E17.5 P0 P3 P7 P10 P14 P21"
# prepare_signal HeartDevelopment   "E10.5 E11.5 E12.5 E13.5 E14.5 E15.5 E16.5 P0"
# prepare_signal LiverDevelopment   "E11.5 E12.5 E13.5 E14.5 E15.5 E16.5 P0"

prepare_CHM(){
    ln -s ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/*.CHM.bed .
}
prepare_CHM

prepare_CpG(){
    ln -s /mnt/Storage/home/yanghui/annotations/mm10/Bins/mm10.b1kb.s10bp.euchr.CpGnumber.bw .
}
prepare_CpG

