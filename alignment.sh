#!/bin/bash
module load bowtie2/2.5.1-gcc-13.2.0-python-3.11.6
bowtie2 -x bowtie_index/yeast_S288C SRR7059706_1_trimmed.fq --no-unal -S wt_r2_1.sam
bowtie2 -x bowtie_index/yeast_S288C SRR7059706_2_trimmed.fq --no-unal -S wt_r2_2.sam
