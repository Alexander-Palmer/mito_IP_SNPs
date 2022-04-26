#!/bin/bash

#Remove empty reads using picard

for file in mrkdup*.bam
do
	java -jar java_dir/picard.jar FilterSamReads \
	I=$file \
	O=noempty_$file \
	READ_LIST_FILE=namesempty.txt FILTER=excludeReadList
done

################################################################

#Mark duplicates using picard

for file in *sort.bam
do
	java -jar java_dir/picard.jar MarkDuplicates \
	I=$file \
	O=mrkdup_$file \
	M=mrkdup_metric_$file.txt \
	TAG_DUPLICATE_SET_MEMBERS=true \
	BARCODE_TAG=MI
done

################################################################

#Input for mutect2.sh
#Input for qualimap.sh
