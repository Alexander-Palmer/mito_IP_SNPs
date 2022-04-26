#!/bin/bash

#Use mutect2 to identify low-frequency SNVs & indels in the data

################################################################

#Create reference index (dictionary) file
if [! -f GRCm38_68.dict ]; then
	./packages/gatk-4.2.4.1/gatk CreateSequenceDictionary \
		-R GRCm38_68.fa
fi

################################################################

#Add read groups for Mutect2
for file in mrkdup*.bam
do
	java -jar java_dir/picard.jar AddOrReplaceReadGroups \
		I=$file \
		O=RG_$file \
		RGID=1 \
		RGLB=lib1 \
		RGPL=illumina \
		RGPU=barcode1 \
		RGSM=20
done

################################################################

#Index .bam files for mutect2
for file in RG*.bam
do
	samtools index $file
done

################################################################

#Validate bam/sam file
for file in RG*.bam
do
	java -jar java_dir/picard.jar ValidateSamFile \
		I=$file \
		MODE=VERBOSE
done

################################################################

#Run mutect2 program
for file in RG*.bam
do
	./packages/gatk-4.2.4.1/gatk Mutect2 \
		-R reference_genome/GRCm38_68.fa \
		-I $file \
		-O mutect2_$file.vcf.gz \
		--mitochondria-mode \
		-L MT \
		--native-pair-hmm-threads 8
done

################################################################

#Input for NCBI genome workbench