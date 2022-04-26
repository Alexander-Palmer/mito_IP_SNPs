#!/bin/bash

#Create index for BWA

#bwa index GRCm38_68.fa

#Use BWA to align reads to reference genome

bwa mem -t 48 -M reference_genome/GRCm38_68.fa S1825-LDK0749-12_L3_1.fq S1825-LDK0749-12_L3_2.fq > S1825.bam
bwa mem -t 48 -M reference_genome/GRCm38_68.fa S1829-LDK0749-4_L3_1.fq S1829-LDK0749-4_L3_2.fq > S1829.bam
bwa mem -t 48 -M reference_genome/GRCm38_68.fa S3444-LDK0749-9_L3_1.fq S3444-LDK0749-9_L3_2.fq > S3444.bam
bwa mem -t 48 -M reference_genome/GRCm38_68.fa S3446-LDK0749-3_L3_1.fq S3446-LDK0749-3_L3_2.fq > S3446.bam
bwa mem -t 48 -M reference_genome/GRCm38_68.fa S3450-LDK0749-2_L3_1.fq S3450-LDK0749-2_L3_2.fq > S3450.bam
bwa mem -t 48 -M reference_genome/GRCm38_68.fa S3451-LDK0749-1_L3_1.fq S3451-LDK0749-1_L3_2.fq > S3451.bam

bwa mem -t 48 -M reference_genome/GRCm38_68.fa S3453-LDK0749-7_L3_1.fq S3453-LDK0749-7_L3_2.fq > S3453.bam
bwa mem -t 48 -M reference_genome/GRCm38_68.fa S3454-LDK0749-11_L3_1.fq S3454-LDK0749-11_L3_2.fq > S3454.bam
bwa mem -t 48 -M reference_genome/GRCm38_68.fa S3458-LDK0749-5_L3_1.fq S3458-LDK0749-5_L3_2.fq > S3458.bam
bwa mem -t 48 -M reference_genome/GRCm38_68.fa S3460-LDK0749-8_L3_1.fq S3460-LDK0749-8_L3_2.fq > S3460.bam
bwa mem -t 48 -M reference_genome/GRCm38_68.fa S3461-LDK0749-6_L3_1.fq S3461-LDK0749-6_L3_2.fq > S3461.bam
bwa mem -t 48 -M reference_genome/GRCm38_68.fa S3476-LDK0749-10_L3_1.fq S3476-LDK0749-10_L3_2.fq > S3476.bam

#Used as input for samtools.sh
