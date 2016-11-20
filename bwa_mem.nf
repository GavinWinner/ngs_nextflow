#!/usr/bin/env nextflow

/*
 * Defines some parameters in order to specify the reference genomes
 * read pairs, threads and output by using the command line options
 */
params.threads = 28
params.genome_fasta = "/home/ubuntu/scratch/genomes/bwa_index/hs37d5-viral-prok.fa"
params.fastq_r1 = "./test/chr20.R1.fq.gz"
params.fastq_r2 = "./test/chr20.R2.fq.gz"
params.outdir = "/home/ubuntu/scratch/alignment_files"
params.outdir_tmp = "/home/ubuntu/scratch/alignment_files/tmp"
params.bam_prefix = "NA12878_NIST7035_TAAGGCGA_Chr20_test_01"


/*
 * fastq R1 and R2 files
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
 * Step 1.0 : bwa alignment: samblaster and sambamba
 */

process bwa_mem {
  echo true
  tag { params.bam_prefix }
  cpus params.threads
  publishDir params.outdir,  mode: 'copy', overwrite: false

  input:
    file fastq1
    file fastq2

  output:
    file '*.dupemk.bam' into dupemk_bam
    file '*.disc.sam' into disc_sam
    file '*.split.sam' into split_sam
    file '*.unmapped.fastq' into unmapped_fastq

"""
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
"""
}

/*
 * Step 1.1 index bwa aligned bam file
 */

process index_dupemk_bam{
  echo true
  tag { params.bam_prefix }
  cpus params.threads
  publishDir params.outdir,  mode: 'copy', overwrite: false

  input:
  file bam from dupemk_bam

  output:
    file '*.bai' into bai_file

"""
sambamba index --nthreads=${task.cpus} $bam
"""
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
