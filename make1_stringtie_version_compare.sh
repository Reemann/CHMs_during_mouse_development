#!/bin/bash

anPATH=/mnt/Storage/home/yanghui/annotations
pyPATH=/mnt/Storage/home/yanghui/scripts/python
pyPATH_w=/mnt/Storage/home/wangyiman/bin/utilities

cd

step2_quantify(){
    dirPATH=${1}
    sampleN=${2}
    genomeVersion=${3}

    mkdir -p ${dirPATH}/${sampleN}_stringtieV2
    cd ${dirPATH}/${sampleN}_stringtieV2
	/mnt/Storage/home/yanghui/anaconda3/bin/stringtie ${dirPATH}/${sampleN}/${sampleN}.bam -o ${sampleN}.gtf -p 30 -G ${anPATH}/${genomeVersion}/${genomeVersion}.refGene.withGeneSymbol.gtf -A ${sampleN}.gene_abund.tab -B -e
	mv t_data.ctab ${sampleN}.t_data.ctab
	mv e_data.ctab ${sampleN}.e_data.ctab
	mv i_data.ctab ${sampleN}.i_data.ctab
	mv e2t.ctab ${sampleN}.e2t.ctab
	mv i2t.ctab ${sampleN}.i2t.ctab
}

step2_quantify ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment E14.5.RNASeq.rep3 mm10


