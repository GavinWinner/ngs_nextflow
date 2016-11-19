#!/usr/bin/env nextflow

/*
 * Defines some parameters in order to specify the reference genomes
 * read pairs, threads and output by using the command line options
 */
params.threads = 30
params.genome_fasta = "/home/ubuntu/scratch/genomes/bwa_index/hs37d5-viral-prok.fa"
params.fastq_r1 = "/home/ubuntu/scratch/fastq/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz"
params.fastq_r2 = "/home/ubuntu/scratch/fastq/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz"
params.outdir = "/home/ubuntu/scratch/alignment_files"
params.bam_prefix = "NA12878_NIST7035_TAAGGCGA"


/*
 * ref genome and fastq R1 and R2 files
 */

fastq1 = file(params.fastq_r1)
fastq2 = file(params.fastq_r2)
genome = file(params.genome_fasta)
threads = params.threads

/*
 * log.info
 */
log.info "BWA MEM PIPELINE                 "
log.info "================================="
log.info "genome             : ${params.genome_fasta}"
log.info "R1                 : ${params.fastq_r1}"
log.info "R2                 : ${params.fastq_r2}"
log.info "bam_prefix         : ${params.bam_prefix}"
log.info "outdir             : ${params.outdir}"
log.info "threads            : ${params.threads}"

/*
 * bwa alignment
 */

process bwa_mem {

  input:
    file genome
    file fastq1
    file fastq2
    val n_cpu from params.threads

  output:
    set bam_out from params.bam_prefix
    set out from params.outdir

"""
mkdir $out
mkdir $out/tmp

bwa mem -M -t $n_cpu $genome $fastq1 $fastq2 | \
  sambamba view -t $n_cpu -S -f bam /dev/stdin | \
  sambamba sort -t $n_cpu -m 2GB --tmpdir=${out}/tmp -o ${out}/${bam_out}.dupemk.bam /dev/stdin && \
  sambamba index ${out}/${bam_out}.dupemk.bam
"""
}
