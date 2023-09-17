#!/bin/bash

cd ~/CHMsInOtherContexts/CellStateTransition/pcar_test
pcar -m callchm -H E11.5.H3K9me3.rep1.bam -C E11.5.input.rep1.bam -M E11.5.WGBS.rep1.sam.G.bed -Z /mnt/Storage/home/yanghui/annotations/mm10/mm10_euch.chrom.sizes -Q /mnt/Storage/home/yanghui/annotations/mm10/mm10.2bit -P 0.005 -D 'CH-nonM' -N test_P0p005_CHnonM
echo "done"