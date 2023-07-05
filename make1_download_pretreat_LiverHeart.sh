#!/bin/bash

cd ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment


### ---------- Liver ----------

### download
cd ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment
# download_raw_seq.py -p PRJNA63471 -f SRA_WGBS_K9_liver.txt -o 0_raw_data
download_raw_seq.py -p PRJNA63471 -f /mnt/Storage/home/wangyiman/CHMsInOtherContexts/CellStateTransition/DataCollection/download_sra_list_liver.txt -o 0_raw_data
cd 0_raw_data
bash PRJNA63471.sh
# SRR3651382
# wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR365/002/SRR3651382/SRR3651382.fastq.gz
# SRR3651454
# SRR3651455
# wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR365/004/SRR3651454/SRR3651454.fastq.gz
# wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR365/005/SRR3651455/SRR3651455.fastq.gz
# SRA_WGBS_K9_liver.txt比download_sra_list_liver.txt少两个: 
# SRR3652945
# SRR3652946
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR365/005/SRR3652945/SRR3652945.fastq.gz &
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR365/006/SRR3652946/SRR3652946.fastq.gz &

### E14.5 input ChIP-seq
download_raw_seq.py -p PRJNA63475 -f input_liver_E14p5.txt -o 0_raw_data
cd 0_raw_data
bash PRJNA63475.sh
### 文章作者没有在ENA数据库上上传fastq file
### ！！！
### 这套数据不用了


# cd ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment
# mkdir E11.5.H3K9me3.rep1 && cd E11.5.H3K9me3.rep1
# ln -sf ../0_raw_data/SRR3653117.fastq.gz E11.5.H3K9me3.rep1_part1.fq.gz
# ln -sf ../0_raw_data/SRR3653118.fastq.gz E11.5.H3K9me3.rep1_part2.fq.gz && cd ../
# cat E11.5.H3K9me3.rep1/E11.5.H3K9me3.rep1_part*.fq.gz > E11.5.H3K9me3.rep1/E11.5.H3K9me3.rep1.fq.gz && rm E11.5.H3K9me3.rep1/E11.5.H3K9me3.rep1_part*.fq.gz
# mkdir E11.5.H3K9me3.rep2 && cd E11.5.H3K9me3.rep2
# ln -sf ../0_raw_data/SRR3653119.fastq.gz E11.5.H3K9me3.rep2_part1.fq.gz
# ln -sf ../0_raw_data/SRR3653120.fastq.gz E11.5.H3K9me3.rep2_part2.fq.gz && cd ../
# cat E11.5.H3K9me3.rep2/E11.5.H3K9me3.rep2_part*.fq.gz > E11.5.H3K9me3.rep2/E11.5.H3K9me3.rep2.fq.gz && rm E11.5.H3K9me3.rep2/E11.5.H3K9me3.rep2_part*.fq.gz
# mkdir E12.5.H3K9me3.rep1 && cd E12.5.H3K9me3.rep1
# ln -sf ../0_raw_data/SRR3651489.fastq.gz E12.5.H3K9me3.rep1.fq.gz && cd ../
# mkdir E12.5.H3K9me3.rep2 && cd E12.5.H3K9me3.rep2
# ln -sf ../0_raw_data/SRR3651490.fastq.gz E12.5.H3K9me3.rep2.fq.gz && cd ../
# mkdir E13.5.H3K9me3.rep1 && cd E13.5.H3K9me3.rep1
# ln -sf ../0_raw_data/SRR3651548.fastq.gz E13.5.H3K9me3.rep1.fq.gz && cd ../
# mkdir E13.5.H3K9me3.rep2 && cd E13.5.H3K9me3.rep2
# ln -sf ../0_raw_data/SRR3651549.fastq.gz E13.5.H3K9me3.rep2.fq.gz && cd ../
# mkdir E14.5.H3K9me3.rep1 && cd E14.5.H3K9me3.rep1
# ln -sf ../0_raw_data/SRR3651430.fastq.gz E14.5.H3K9me3.rep1_part1.fq.gz
# ln -sf ../0_raw_data/SRR3651431.fastq.gz E14.5.H3K9me3.rep1_part2.fq.gz
# ln -sf ../0_raw_data/SRR3651432.fastq.gz E14.5.H3K9me3.rep1_part3.fq.gz && cd ../
# cat E14.5.H3K9me3.rep1/E14.5.H3K9me3.rep1_part*.fq.gz > E14.5.H3K9me3.rep1/E14.5.H3K9me3.rep1.fq.gz && rm E14.5.H3K9me3.rep1/E14.5.H3K9me3.rep1_part*.fq.gz
# mkdir E14.5.H3K9me3.rep2 && cd E14.5.H3K9me3.rep2
# ln -sf ../0_raw_data/SRR3651433.fastq.gz E14.5.H3K9me3.rep2_part1.fq.gz
# ln -sf ../0_raw_data/SRR3651434.fastq.gz E14.5.H3K9me3.rep2_part2.fq.gz
# ln -sf ../0_raw_data/SRR3651435.fastq.gz E14.5.H3K9me3.rep2_part3.fq.gz && cd ../
# cat E14.5.H3K9me3.rep2/E14.5.H3K9me3.rep2_part*.fq.gz > E14.5.H3K9me3.rep2/E14.5.H3K9me3.rep2.fq.gz && rm E14.5.H3K9me3.rep2/E14.5.H3K9me3.rep2_part*.fq.gz
# mkdir E15.5.H3K9me3.rep1 && cd E15.5.H3K9me3.rep1
# ln -sf ../0_raw_data/SRR3653134.fastq.gz E15.5.H3K9me3.rep1.fq.gz && cd ../
# mkdir E15.5.H3K9me3.rep2 && cd E15.5.H3K9me3.rep2
# ln -sf ../0_raw_data/SRR3653135.fastq.gz E15.5.H3K9me3.rep2.fq.gz && cd ../
# mkdir E16.5.H3K9me3.rep1 && cd E16.5.H3K9me3.rep1
# ln -sf ../0_raw_data/SRR3651381.fastq.gz E16.5.H3K9me3.rep1.fq.gz && cd ../
# mkdir E16.5.H3K9me3.rep2 && cd E16.5.H3K9me3.rep2
# ln -sf ../0_raw_data/SRR3651382.fastq.gz E16.5.H3K9me3.rep2.fq.gz && cd ../
# mkdir P0.H3K9me3.rep1 && cd P0.H3K9me3.rep1
# ln -sf ../0_raw_data/SRR3651898.fastq.gz P0.H3K9me3.rep1_part1.fq.gz
# ln -sf ../0_raw_data/SRR3651899.fastq.gz P0.H3K9me3.rep1_part2.fq.gz && cd ../
# cat P0.H3K9me3.rep1/P0.H3K9me3.rep1_part*.fq.gz > P0.H3K9me3.rep1/P0.H3K9me3.rep1.fq.gz && rm P0.H3K9me3.rep1/P0.H3K9me3.rep1_part*.fq.gz
# mkdir P0.H3K9me3.rep2 && cd P0.H3K9me3.rep2
# ln -sf ../0_raw_data/SRR3651900.fastq.gz P0.H3K9me3.rep2_part1.fq.gz
# ln -sf ../0_raw_data/SRR3651901.fastq.gz P0.H3K9me3.rep2_part2.fq.gz
# ln -sf ../0_raw_data/SRR3651902.fastq.gz P0.H3K9me3.rep2_part3.fq.gz && cd ../
# cat P0.H3K9me3.rep2/P0.H3K9me3.rep2_part*.fq.gz > P0.H3K9me3.rep2/P0.H3K9me3.rep2.fq.gz && rm P0.H3K9me3.rep2/P0.H3K9me3.rep2_part*.fq.gz
# mkdir E11.5.input.rep1 && cd E11.5.input.rep1
# ln -sf ../0_raw_data/SRR3652454.fastq.gz E11.5.input.rep1.fq.gz && cd ../
# mkdir E11.5.input.rep2 && cd E11.5.input.rep2
# ln -sf ../0_raw_data/SRR3652455.fastq.gz E11.5.input.rep2.fq.gz && cd ../
# mkdir E12.5.input.rep1 && cd E12.5.input.rep1
# ln -sf ../0_raw_data/SRR3652762.fastq.gz E12.5.input.rep1_part1.fq.gz
# ln -sf ../0_raw_data/SRR3652763.fastq.gz E12.5.input.rep1_part2.fq.gz && cd ../
# cat E12.5.input.rep1/E12.5.input.rep1_part*.fq.gz > E12.5.input.rep1/E12.5.input.rep1.fq.gz && rm E12.5.input.rep1/E12.5.input.rep1_part*.fq.gz
# mkdir E12.5.input.rep2 && cd E12.5.input.rep2
# ln -sf ../0_raw_data/SRR3652764.fastq.gz E12.5.input.rep2_part1.fq.gz
# ln -sf ../0_raw_data/SRR3652765.fastq.gz E12.5.input.rep2_part2.fq.gz && cd ../
# cat E12.5.input.rep2/E12.5.input.rep2_part*.fq.gz > E12.5.input.rep2/E12.5.input.rep2.fq.gz && rm E12.5.input.rep2/E12.5.input.rep2_part*.fq.gz
# mkdir E13.5.input.rep1 && cd E13.5.input.rep1
# ln -sf ../0_raw_data/SRR3651584.fastq.gz E13.5.input.rep1_part1.fq.gz
# ln -sf ../0_raw_data/SRR3651585.fastq.gz E13.5.input.rep1_part2.fq.gz
# ln -sf ../0_raw_data/SRR3651586.fastq.gz E13.5.input.rep1_part3.fq.gz && cd ../
# cat E13.5.input.rep1/E13.5.input.rep1_part*.fq.gz > E13.5.input.rep1/E13.5.input.rep1.fq.gz && rm E13.5.input.rep1/E13.5.input.rep1_part*.fq.gz
# mkdir E13.5.input.rep2 && cd E13.5.input.rep2
# ln -sf ../0_raw_data/SRR3651587.fastq.gz E13.5.input.rep2_part1.fq.gz
# ln -sf ../0_raw_data/SRR3651588.fastq.gz E13.5.input.rep2_part2.fq.gz && cd ../
# cat E13.5.input.rep2/E13.5.input.rep2_part*.fq.gz > E13.5.input.rep2/E13.5.input.rep2.fq.gz && rm E13.5.input.rep2/E13.5.input.rep2_part*.fq.gz
mkdir E14.5.input.rep1 && cd E14.5.input.rep1
ln -sf ../0_raw_data/SRR3652945.fastq.gz E14.5.input.rep1.fq.gz && cd ../
mkdir E14.5.input.rep2 && cd E14.5.input.rep2
ln -sf ../0_raw_data/SRR3652946.fastq.gz E14.5.input.rep2.fq.gz && cd ../
# mkdir E15.5.input.rep1 && cd E15.5.input.rep1
# ln -sf ../0_raw_data/SRR3651754.fastq.gz E15.5.input.rep1_part1.fq.gz
# ln -sf ../0_raw_data/SRR3651755.fastq.gz E15.5.input.rep1_part2.fq.gz && cd ../
# cat E15.5.input.rep1/E15.5.input.rep1_part*.fq.gz > E15.5.input.rep1/E15.5.input.rep1.fq.gz && rm E15.5.input.rep1/E15.5.input.rep1_part*.fq.gz
# mkdir E15.5.input.rep2 && cd E15.5.input.rep2
# ln -sf ../0_raw_data/SRR3651756.fastq.gz E15.5.input.rep2_part1.fq.gz
# ln -sf ../0_raw_data/SRR3651757.fastq.gz E15.5.input.rep2_part2.fq.gz && cd ../
# cat E15.5.input.rep2/E15.5.input.rep2_part*.fq.gz > E15.5.input.rep2/E15.5.input.rep2.fq.gz && rm E15.5.input.rep2/E15.5.input.rep2_part*.fq.gz
# mkdir E16.5.input.rep1 && cd E16.5.input.rep1
# ln -sf ../0_raw_data/SRR3651335.fastq.gz E16.5.input.rep1.fq.gz && cd ../
# mkdir E16.5.input.rep2 && cd E16.5.input.rep2
# ln -sf ../0_raw_data/SRR3651336.fastq.gz E16.5.input.rep2.fq.gz && cd ../
# mkdir P0.input.rep1 && cd P0.input.rep1
# ln -sf ../0_raw_data/SRR3652553.fastq.gz P0.input.rep1.fq.gz && cd ../
# mkdir P0.input.rep2 && cd P0.input.rep2
# ln -sf ../0_raw_data/SRR3652554.fastq.gz P0.input.rep2.fq.gz && cd ../
# mkdir E11.5.WGBS.rep1 && cd E11.5.WGBS.rep1
# ln -sf ../0_raw_data/SRR3651222.fastq.gz E11.5.WGBS.rep1.fq.gz && cd ../
# mkdir E11.5.WGBS.rep2 && cd E11.5.WGBS.rep2
# ln -sf ../0_raw_data/SRR3651223.fastq.gz E11.5.WGBS.rep2.fq.gz && cd ../
# mkdir E11.5.WGBS.rep3 && cd E11.5.WGBS.rep3
# ln -sf ../0_raw_data/SRR3651224.fastq.gz E11.5.WGBS.rep3.fq.gz && cd ../
# mkdir E11.5.WGBS.rep4 && cd E11.5.WGBS.rep4
# ln -sf ../0_raw_data/SRR3651225.fastq.gz E11.5.WGBS.rep4.fq.gz && cd ../
# mkdir E12.5.WGBS.rep1 && cd E12.5.WGBS.rep1
# ln -sf ../0_raw_data/SRR3651878.fastq.gz E12.5.WGBS.rep1.fq.gz && cd ../
# mkdir E12.5.WGBS.rep2 && cd E12.5.WGBS.rep2
# ln -sf ../0_raw_data/SRR3651879.fastq.gz E12.5.WGBS.rep2.fq.gz && cd ../
# mkdir E12.5.WGBS.rep3 && cd E12.5.WGBS.rep3
# ln -sf ../0_raw_data/SRR3651880.fastq.gz E12.5.WGBS.rep3.fq.gz && cd ../
# mkdir E12.5.WGBS.rep4 && cd E12.5.WGBS.rep4
# ln -sf ../0_raw_data/SRR3651881.fastq.gz E12.5.WGBS.rep4.fq.gz && cd ../
# mkdir E13.5.WGBS.rep1 && cd E13.5.WGBS.rep1
# ln -sf ../0_raw_data/SRR3652697.fastq.gz E13.5.WGBS.rep1.fq.gz && cd ../
# mkdir E13.5.WGBS.rep2 && cd E13.5.WGBS.rep2
# ln -sf ../0_raw_data/SRR3652698.fastq.gz E13.5.WGBS.rep2.fq.gz && cd ../
# mkdir E13.5.WGBS.rep3 && cd E13.5.WGBS.rep3
# ln -sf ../0_raw_data/SRR3652699.fastq.gz E13.5.WGBS.rep3.fq.gz && cd ../
# mkdir E13.5.WGBS.rep4 && cd E13.5.WGBS.rep4
# ln -sf ../0_raw_data/SRR3652700.fastq.gz E13.5.WGBS.rep4.fq.gz && cd ../
# mkdir E14.5.WGBS.rep1 && cd E14.5.WGBS.rep1
# ln -sf ../0_raw_data/SRR3653025.fastq.gz E14.5.WGBS.rep1.fq.gz && cd ../
# mkdir E14.5.WGBS.rep2 && cd E14.5.WGBS.rep2
# ln -sf ../0_raw_data/SRR3653026.fastq.gz E14.5.WGBS.rep2.fq.gz && cd ../
# mkdir E14.5.WGBS.rep3 && cd E14.5.WGBS.rep3
# ln -sf ../0_raw_data/SRR3653027.fastq.gz E14.5.WGBS.rep3.fq.gz && cd ../
# mkdir E14.5.WGBS.rep4 && cd E14.5.WGBS.rep4
# ln -sf ../0_raw_data/SRR3653028.fastq.gz E14.5.WGBS.rep4.fq.gz && cd ../
# mkdir E14.5.WGBS.rep5 && cd E14.5.WGBS.rep5
# ln -sf ../0_raw_data/SRR3653029.fastq.gz E14.5.WGBS.rep5.fq.gz && cd ../
# mkdir E14.5.WGBS.rep6 && cd E14.5.WGBS.rep6
# ln -sf ../0_raw_data/SRR3653030.fastq.gz E14.5.WGBS.rep6.fq.gz && cd ../
# mkdir E15.5.WGBS.rep1 && cd E15.5.WGBS.rep1
# ln -sf ../0_raw_data/SRR3651847.fastq.gz E15.5.WGBS.rep1.fq.gz && cd ../
# mkdir E15.5.WGBS.rep2 && cd E15.5.WGBS.rep2
# ln -sf ../0_raw_data/SRR3651848.fastq.gz E15.5.WGBS.rep2.fq.gz && cd ../
# mkdir E15.5.WGBS.rep3 && cd E15.5.WGBS.rep3
# ln -sf ../0_raw_data/SRR3651849.fastq.gz E15.5.WGBS.rep3.fq.gz && cd ../
# mkdir E15.5.WGBS.rep4 && cd E15.5.WGBS.rep4
# ln -sf ../0_raw_data/SRR3651850.fastq.gz E15.5.WGBS.rep4.fq.gz && cd ../
# mkdir E16.5.WGBS.rep1 && cd E16.5.WGBS.rep1
# ln -sf ../0_raw_data/SRR3651452.fastq.gz E16.5.WGBS.rep1.fq.gz && cd ../
# mkdir E16.5.WGBS.rep2 && cd E16.5.WGBS.rep2
# ln -sf ../0_raw_data/SRR3651453.fastq.gz E16.5.WGBS.rep2.fq.gz && cd ../
# mkdir E16.5.WGBS.rep3 && cd E16.5.WGBS.rep3
# ln -sf ../0_raw_data/SRR3651454.fastq.gz E16.5.WGBS.rep3.fq.gz && cd ../
# mkdir E16.5.WGBS.rep4 && cd E16.5.WGBS.rep4
# ln -sf ../0_raw_data/SRR3651455.fastq.gz E16.5.WGBS.rep4.fq.gz && cd ../
# mkdir P0.WGBS.rep1 && cd P0.WGBS.rep1
# ln -sf ../0_raw_data/SRR3652380.fastq.gz P0.WGBS.rep1.fq.gz && cd ../
# mkdir P0.WGBS.rep2 && cd P0.WGBS.rep2
# ln -sf ../0_raw_data/SRR3652381.fastq.gz P0.WGBS.rep2.fq.gz && cd ../
# mkdir P0.WGBS.rep3 && cd P0.WGBS.rep3
# ln -sf ../0_raw_data/SRR3652382.fastq.gz P0.WGBS.rep3.fq.gz && cd ../
# mkdir P0.WGBS.rep4 && cd P0.WGBS.rep4
# ln -sf ../0_raw_data/SRR3652383.fastq.gz P0.WGBS.rep4.fq.gz && cd ../


ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192347 SINGLE E11.5.H3K9me3.rep1 mm10
ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192348 SINGLE E11.5.H3K9me3.rep2 mm10
ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191963 SINGLE E11.5.input.rep1 mm10
ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191964 SINGLE E11.5.input.rep2 mm10
# WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2190976 SINGLE E11.5.WGBS.rep1 mm10
# WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2190977 SINGLE E11.5.WGBS.rep2 mm10
# WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2190978 SINGLE E11.5.WGBS.rep3 mm10
# WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2190979 SINGLE E11.5.WGBS.rep4 mm10
# ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191248 SINGLE E12.5.H3K9me3.rep1 mm10
# ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191249 SINGLE E12.5.H3K9me3.rep2 mm10
ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192133 SINGLE E12.5.input.rep1 mm10
ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192134 SINGLE E12.5.input.rep2 mm10
# WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191548 SINGLE E12.5.WGBS.rep1 mm10
# WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191549 SINGLE E12.5.WGBS.rep2 mm10
WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191550 SINGLE E12.5.WGBS.rep3 mm10
WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191551 SINGLE E12.5.WGBS.rep4 mm10
# ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191298 SINGLE E13.5.H3K9me3.rep1 mm10
# ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191299 SINGLE E13.5.H3K9me3.rep2 mm10
# ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191328 SINGLE E13.5.input.rep1 mm10
# ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191329 SINGLE E13.5.input.rep2 mm10
WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192103 SINGLE E13.5.WGBS.rep1 mm10
WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192104 SINGLE E13.5.WGBS.rep2 mm10
WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192105 SINGLE E13.5.WGBS.rep3 mm10
WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192106 SINGLE E13.5.WGBS.rep4 mm10
# ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191198 SINGLE E14.5.H3K9me3.rep1 mm10
# ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191199 SINGLE E14.5.H3K9me3.rep2 mm10
ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192213 SINGLE E14.5.input.rep1 mm10
ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192214 SINGLE E14.5.input.rep2 mm10

WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192279 SINGLE E14.5.WGBS.rep1 mm10
WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192280 SINGLE E14.5.WGBS.rep2 mm10
WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192281 SINGLE E14.5.WGBS.rep3 mm10
WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192282 SINGLE E14.5.WGBS.rep4 mm10
WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192283 SINGLE E14.5.WGBS.rep5 mm10
WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192284 SINGLE E14.5.WGBS.rep6 mm10

ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192353 SINGLE E15.5.H3K9me3.rep1 mm10
ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192354 SINGLE E15.5.H3K9me3.rep2 mm10

ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191450 SINGLE E15.5.input.rep1 mm10
ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191451 SINGLE E15.5.input.rep2 mm10
# WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191530 SINGLE E15.5.WGBS.rep1 mm10
# WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191531 SINGLE E15.5.WGBS.rep2 mm10
WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191532 SINGLE E15.5.WGBS.rep3 mm10
WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191533 SINGLE E15.5.WGBS.rep4 mm10
# ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191141 SINGLE E16.5.H3K9me3.rep1 mm10
ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191142 SINGLE E16.5.H3K9me3.rep2 mm10
# ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191074 SINGLE E16.5.input.rep1 mm10
# ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191075 SINGLE E16.5.input.rep2 mm10
# WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191212 SINGLE E16.5.WGBS.rep1 mm10
# WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191213 SINGLE E16.5.WGBS.rep2 mm10
WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191214 SINGLE E16.5.WGBS.rep3 mm10
WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191215 SINGLE E16.5.WGBS.rep4 mm10

ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191568 SINGLE P0.H3K9me3.rep1 mm10
ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191569 SINGLE P0.H3K9me3.rep2 mm10

ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192021 SINGLE P0.input.rep1 mm10
ChIPseq_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2192022 SINGLE P0.input.rep2 mm10
WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191921 SINGLE P0.WGBS.rep1 mm10
WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191922 SINGLE P0.WGBS.rep2 mm10
WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191923 SINGLE P0.WGBS.rep3 mm10
WGBS_processing_publicData.sh ~/CHMsInOtherContexts/CellStateTransition/LiverDevelopment GEO GSM2191924 SINGLE P0.WGBS.rep4 mm10







### ---------- Heart ----------

### download
cd ~/CHMsInOtherContexts/CellStateTransition/HeartDevelopment
download_raw_seq.py -p PRJNA63471 -f /mnt/Storage/home/wangyiman/CHMsInOtherContexts/CellStateTransition/DataCollection/download_sra_list_heart.txt -o 0_raw_data
cd 0_raw_data
bash PRJNA63471.sh

### md5 wrong files:
# SRR3651163
# re-download on zhanglab3
# SRR3651255.fastq.gz
# wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR365/005/SRR3651255/SRR3651255.fastq.gz
# SRR3651737
# wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR365/007/SRR3651737/SRR3651737.fastq.gz
