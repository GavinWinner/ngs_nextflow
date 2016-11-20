#!/usr/bin/env bash
nextflow run ./bwa_mem.nf \
-with-trace /home/ubuntu/Dropbox/nextflow-work/trace.txt \
-with-timeline /home/ubuntu/Dropbox/nextflow-work/timeline.html \
-work-dir /home/ubuntu/scratch/nextflow-work
