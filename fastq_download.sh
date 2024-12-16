#!/bin/bash
module load sra-tools/3.0.3-gcc-13.2.0
cd /scratch_tmp/users/k24052661/alignment/raw_fastq
fastq-dump --split-files SRR7059704
fastq-dump --split-files SRR7059705
fastq-dump --split-files SRR7059706
fastq-dump --split-files SRR7059707
fastq-dump --split-files SRR7059708
fastq-dump --split-files SRR7059709


