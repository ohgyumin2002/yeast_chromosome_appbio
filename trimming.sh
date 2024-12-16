#!/bin/bash

cd /scratch_tmp/users/k24052661/alignment/raw_fastq

module load fastqc/0.12.1-gcc-13.2.0
module load trimgalore/0.6.6-gcc-13.2.0-python-3.11.6

trim_galore --fastqc --output_dir /scratch_tmp/users/k24052661/alignment/trimmed_fastq *.fastq
