#!/usr/bin/env nextflow

/*
 * Defines some parameters in order to specify the reference genomes
 * read pairs, threads and output by using the command line options
 */

//number of threads
params.threads = 28

//genome version
params.genome_fasta = "${HOME}/Dropbox/ngs_nextflow/test/chr20.fa"

//sample information
params.fastq_r1 = "${HOME}/Dropbox/ngs_nextflow/test/chr20.R1.fq.gz"
params.fastq_r2 = "${HOME}/Dropbox/ngs_nextflow/test/chr20.R2.fq.gz"
params.bam_prefix = "chr20"
params.sample_name="NA12878"
params.platform="ILLUMINA"
params.platform_unit="PU1"
params.library="LB1"
params.rundate="2016"

//output and tempory directories
params.outdir = "${HOME}/Dropbox/ngs_nextflow/test"
params.outdir_tmp = "/home/ubuntu/scratch/alignment_files/tmp"

/*
 * fastq R1 and R2 files
 */

fastq1 = file(params.fastq_r1)
fastq2 = file(params.fastq_r2)

/*
 * log.info
 */
log.info "===================================================================="
log.info "BWA MEM PIPELINE : Map, mark duplicates, sort and index             "
log.info "===================================================================="
log.info "genome             : ${params.genome_fasta}"
log.info "R1                 : ${params.fastq_r1}"
log.info "R2                 : ${params.fastq_r2}"
log.info "RG:SM              : ${params.sample_name}"
log.info "RG:PL              : ${params.platform}"
log.info "RG:PU              : ${params.platform_unit}"
log.info "RG:LB              : ${params.library}"
log.info "RG:DT              : ${params.rundate}"
log.info "bam_prefix         : ${params.bam_prefix}"
log.info "outdir             : ${params.outdir}"
log.info "temp_dir           : ${params.outdir_tmp}"
log.info "threads            : ${params.threads}"

/*
 * Step 1.0 : bwa alignment: samblaster and sambamba
 * https://github.com/nextflow-io/faq
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
    file '*.dupemk.bam' into dupemk_bam, dupemk_bam_to_index, dupemk_bam_for_flagstat
    file '*.disc.sam' into disc_sam
    file '*.split.sam' into split_sam
    file '*.unmapped.fastq' into unmapped_fastq
    file '*.err' into error_file

"""
bwa mem \
-M \
-t ${task.cpus} \
-R '@RG\tID:${params.bam_prefix}\tSM:${params.sample_name}\tPU:${params.platform_unit}\tPL:${params.platform}\tLB:${params.library}\tDT:${params.rundate}' \
$params.genome_fasta \
$fastq1 $fastq2 | \
samblaster --addMateTags --excludeDups \
-d ${params.bam_prefix}.disc.sam \
-s ${params.bam_prefix}.split.sam \
-u ${params.bam_prefix}.unmapped.fastq | \
sambamba view -t ${task.cpus} -S -f bam /dev/stdin | \
sambamba sort -t ${task.cpus} -m 8GB \
--tmpdir=${params.outdir_tmp} \
-o ${params.bam_prefix}.dupemk.bam /dev/stdin 2>${params.bam_prefix}.bwa_mem.err
"""
}

/*
 * Step 1.1 index bwa aligned bam file
 */

process index_dupemk_bam {
  echo true
  tag { params.bam_prefix }
  cpus params.threads
  publishDir params.outdir,  mode: 'copy', overwrite: false

  input:
    file bam_to_index from dupemk_bam_to_index

  output:
    file '*.bai' into bai_file

"""
sambamba index --nthreads=${task.cpus} $bam_to_index
"""
}

/*
 * Step 1.2 flagstats
 */

process sambamba_flagstat {
  echo true
  tag { params.bam_prefix }
  cpus params.threads
  publishDir params.outdir,  mode: 'copy', overwrite: false

  input:
    file bam_to_flagstat from dupemk_bam_for_flagstat
    file bai from bai_file
    
  output:
    file '*.flagstat.txt' into flagstat_file

"""
sambamba flagstat -t ${task.cpus} $bam_to_flagstat > ${params.bam_prefix}.dupemk.bam.flagstat.txt
"""
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
