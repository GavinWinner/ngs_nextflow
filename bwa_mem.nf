#!/usr/bin/env nextflow


params.genome_fasta = "/home/ubuntu/scratch/genomes/bwa_index/hs37d5-viral-prok.fa"
params.threads = 30
params.fastq_r1 =
params.fastq_r2 =
params.bam_out =


process bwa_mem {

"""
bwa mem -M -t 28 ${GENOMEINDEX} ${FQ1} ${FQ2} | \
  samblaster --addMateTags --excludeDups \
  --discordantFile ${SOUTDocker}/alignments/${BAM_PREFIX}.discordant.sam \
  --splitterFile ${SOUTDocker}/alignments/${BAM_PREFIX}.splitread.sam \
  --unmappedFile ${SOUTDocker}/alignments/${BAM_PREFIX}.unmapped.fastq | \
  sambamba view -t ${NCPU} -S -f bam /dev/stdin | \
  sambamba sort -t ${NCPU} -m 2GB --tmpdir=${SOUTDocker}/tmp -o ${SOUTDocker}/alignments/${BAM_PREFIX}.dupemk.bam /dev/stdin && \
  sambamba index ${SOUTDocker}/alignments/${BAM_PREFIX}.dupemk.bam && \
  sambamba flagstat -t ${NCPU} ${SOUTDocker}/alignments/${BAM_PREFIX}.dupemk.bam > ${SOUTDocker}/reports/${BAM_PREFIX}.dupemk.bam.flagstat && \
  bedtools bamtobed -i ${SOUTDocker}/alignments/${BAM_PREFIX}.dupemk.bam | bedtools merge > ${SOUTDocker}/reports/${BAM_PREFIX}.dupemk.bed && \
  sambamba view -t ${NCPU} -S -f bam ${SOUTDocker}/alignments/${BAM_PREFIX}.discordant.sam | \
  sambamba sort -t ${NCPU} -m 2GB --tmpdir=${SOUTDocker}/tmp -o ${SOUTDocker}/alignments/${BAM_PREFIX}.discordant.bam /dev/stdin && \
  sambamba index ${SOUTDocker}/alignments/${BAM_PREFIX}.discordant.bam && \
  sambamba view -t ${NCPU} -S -f bam ${SOUTDocker}/alignments/${BAM_PREFIX}.splitread.sam | \
  sambamba sort -t ${NCPU} -m 2GB --tmpdir=${SOUTDocker}/tmp -o ${SOUTDocker}/alignments/${BAM_PREFIX}.splitread.bam /dev/stdin && \
  sambamba index ${SOUTDocker}/alignments/${BAM_PREFIX}.splitread.bam && \
  rm -v ${SOUTDocker}/alignments/${BAM_PREFIX}.discordant.sam && \
  rm -v ${SOUTDocker}/alignments/${BAM_PREFIX}.splitread.sam
"""
}
