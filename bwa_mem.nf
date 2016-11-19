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
params.outdir_tmp = "/home/ubuntu/scratch/alignment_files/tmp"
params.bam_prefix = "NA12878_NIST7035_TAAGGCGA"


/*
 * ref genome and fastq R1 and R2 files
 */

fastq1 = file(params.fastq_r1)
fastq2 = file(params.fastq_r2)

/*
 * log.info
 */
log.info "BWA MEM PIPELINE : NGSeasy-ish>"
log.info "================================="
log.info "genome             : ${params.genome_fasta}"
log.info "R1                 : ${params.fastq_r1}"
log.info "R2                 : ${params.fastq_r2}"
log.info "bam_prefix         : ${params.bam_prefix}"
log.info "outdir             : ${params.outdir}"
log.info "temp_dir           : ${params.outdir_tmp}"
log.info "threads            : ${params.threads}"

/*
 * bwa alignment
 */

process bwa_mem {
  echo true
  tag { params.bam_prefix }
  cpus params.threads
  publishDir params.outdir, mode: 'copy'

  input:
    file fastq1
    file fastq2

"""
## bwa mem -M
bwa mem -M \
-t ${task.cpus} \
$params.genome_fasta \
$fastq1 $fastq2 | \
samblaster --addMateTags --excludeDups \
-d ${params.bam_prefix}.disc.sam \
-s ${params.bam_prefix}.split.sam \
-u ${params.bam_prefix}.unmapped.fastq | \
sambamba view -t ${task.cpus} -S -f bam /dev/stdin | \
sambamba sort -t ${task.cpus} -m 8GB \
--tmpdir=${params.outdir_tmp} \
-o ${params.bam_prefix}.dupemk.bam /dev/stdin
## index
sambamba index ${params.bam_prefix}.dupemk.bam
"""
}

/*
 * using nf synatx and setting genome = params.genome_fasta
 * breaks nf. using file genome in input
 * nf can't find the index files and breaks
 */

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
