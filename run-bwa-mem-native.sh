#!/bin/bash -ue
## bwa mem -M
bwa mem -M -t 30 /home/ubuntu/scratch/genomes/bwa_index/hs37d5-viral-prok.fa NIST7035_TAAGGCGA_L001_R1_001.fastq.gz NIST7035_TAAGGCGA_L001_R2_001.fastq.gz | samblaster --addMateTags --excludeDups -d NA12878_NIST7035_TAAGGCGA.disc.sam -s NA12878_NIST7035_TAAGGCGA.split.sam -u NA12878_NIST7035_TAAGGCGA.unmapped.fastq | sambamba view -t 30 -S -f bam /dev/stdin | sambamba sort -t 30 -m 8GB --tmpdir=/home/ubuntu/scratch/alignment_files/tmp -o NA12878_NIST7035_TAAGGCGA.dupemk.bam /dev/stdin
## index
sambamba index NA12878_NIST7035_TAAGGCGA.dupemk.bam
