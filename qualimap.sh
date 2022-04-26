#!/bin/bash

#Find sequencing metrics for non-duplicated & duplicated reads

for file in mrkdup*.bam
do
	cd packages/qualimap_v2.2.1
	./qualimap bamqc -sd -bam ../../$file
	cd ../..
done

