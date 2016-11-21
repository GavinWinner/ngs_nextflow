#!/usr/bin/env nextflow

/*
 * Defines some parameters in order to specify the reference genomes
 */

params.genome_fasta = "/home/ubuntu/scratch/genomes/bwa_index/hs37d5-viral-prok.fa"
params.outdir = "/home/ubuntu/scratch/bwa_index"

/*
 * fastq R1 and R2 files
 */

fasta = file(params.genome_fasta)
params.outdir = "/home/ubuntu/scratch/bwa_index"

/*
 * log.info: bwa Version: 0.7.15-r1140
 */
log.info "BWA Index > samtools faidx > pircard seq dict "
log.info "=============================================="
log.info "genome             : ${params.genome_fasta}   "
log.info "out dir            : ${params.outdir}         "

/*
 * Step 1.0 : Generate the BWA index
 */

process bwa_index {
  echo true
  tag { fasta }
  publishDir params.outdir,  mode: 'copy', overwrite: false

  input:
    file fasta

  output:
    file '*.amb' into bwa_index_amb
    file '*.ann' into bwa_index_ann
    file '*.bwt' into bwa_index_bwt
    file '*.pac' into bwa_index_pac
    file '*.sa'  into bwa_index_sa

"""
bwa index -a bwtsw $fasta
"""
}

/*
 * Step 2.0 : Generate the fasta file index
 * https://github.com/nextflow-io/faq
 */

process samtools_faidx {
  echo true
  tag { fasta }
  publishDir params.outdir,  mode: 'copy', overwrite: false

input:
  file fasta

output:
  file '*.fai' into samtools_fai

"""
samtools faidx $fasta
"""
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
