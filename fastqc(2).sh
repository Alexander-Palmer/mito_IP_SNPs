#!/bin/bash

#Perform fastQC quality control on trimmed sequencing data

for file in trim_no_UMI*.fq
do
	fastqc $file
done
