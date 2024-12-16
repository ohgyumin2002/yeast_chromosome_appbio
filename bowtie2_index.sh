#!/bin/bash

module load bowtie2/2.5.1-gcc-13.2.0-python-3.11.6
bowtie2-build S288C.fasta yeast_S288C
bowtie2-inspect -n yeast_S288C
export BOWTIE2_INDEXES=/scratch_tmp/users/k24052661/alignment/bowtie_index
cp *.bt2 $BOWTIE2_INDEXES
