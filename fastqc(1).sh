#!/bin/bash

#Perform fastQC quality control on raw sequencing data

for file in *.fq
do
	fastqc $file
done
