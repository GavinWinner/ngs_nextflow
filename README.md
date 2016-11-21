# NGS Pipelines - Nextflow

Playing around with nextflow for basic NGS alignment and variant calling pipelines.

## Requirements

To do....set up and installation and dockerfiles

- Linux, 64bit
- Anaconda  
- using bioconda packages  
- register gatk with bioconda  
- full tool list coming soon
- `bwa`, `samtools`, `samblaster`, `picard`, `freebayes`  

## Test fatsq

chr20 from `NIST7035_H7AP8ADXX_TAAGGCGA NA1287`

To do....how i got the data

## NF Test Pipelines

- [bwa_index.nf](https://github.com/snewhouse/ngs_nextflow/blob/master/bwa_index.nf). NOT TESTED  
- [bwa_mem.nf](https://github.com/snewhouse/ngs_nextflow/blob/master/bwa_mem.nf)  

## Running a test pipelines

```bash
./run-bwa-nextflow.sh
```

```
N E X T F L O W  ~  version 0.22.5
Launching `./bwa_mem.nf` [amazing_wright] - revision: a1b8a4ea96
BWA MEM PIPELINE : NGSeasy-ish>
=================================
genome             : /home/ubuntu/scratch/genomes/bwa_index/hs37d5-viral-prok.fa
R1                 : ./test/chr20.R1.fq.gz
R2                 : ./test/chr20.R2.fq.gz
bam_prefix         : NA12878_NIST7035_TAAGGCGA_Chr20_test_01
outdir             : /home/ubuntu/scratch/alignment_files
temp_dir           : /home/ubuntu/scratch/alignment_files/tmp
threads            : 28
[warm up] executor > local
[05/31e956] Submitted process > bwa_mem (NA12878_NIST7035_TAAGGCGA_Chr20_test_01)
[99/275f80] Submitted process > index_dupemk_bam (NA12878_NIST7035_TAAGGCGA_Chr20_test_01)
[d0/58f8ad] Submitted process > call_variants_freebayes (NA12878_NIST7035_TAAGGCGA_Chr20_test_01)
Done!
```

****

## Future...

- GATK best practices pipeline
