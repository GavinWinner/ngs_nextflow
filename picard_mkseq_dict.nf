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
 * Step 3.0 : Generate the sequence dictionary
 * https://github.com/nextflow-io/faq
 */


process picard_sequence_dictionary {
  echo true
  tag { fasta }
  publishDir params.outdir,  mode: 'copy', overwrite: false

input:
  file fasta

output:
  file "${fasta}.dict" into sequence_dict

shell:
"""
picard CreateSequenceDictionary \
REFERENCE=$fasta \
OUTPUT="${fasta}.dict"
"""
}
