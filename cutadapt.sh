#!/bin/bash

#Remove Illumina universal adapter and UMIs (first 5bp)

for file in *.fq
do
	cutadapt -b AGATCGGAAGAG -j 12 -o $file-trim.fq $file.fq
	cutadapt -j 24 -u 5 -m 21 -o $file-trim-noUMI.fq $file-trim.fq
done
