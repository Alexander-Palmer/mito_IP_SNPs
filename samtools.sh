#!/bin/bash

#Use samtools to sort reads by coordinate for MarkDuplicates in Picard

samtools sort -@ 24 -o S1825_samtools-sort.bam S1825.bam
samtools sort -@ 24 -o S1829_samtools-sort.bam S1829.bam
samtools sort -@ 24 -o S3444_samtools-sort.bam S3444.bam
samtools sort -@ 24 -o S3446_samtools-sort.bam S3446.bam
samtools sort -@ 24 -o S3450_samtools-sort.bam S3450.bam
samtools sort -@ 24 -o S3451_samtools-sort.bam S3451.bam

samtools sort -@ 24 -o S3453_samtools-sort.bam S3453.bam
samtools sort -@ 24 -o S3454_samtools-sort.bam S3454.bam
samtools sort -@ 24 -o S3458_samtools-sort.bam S3458.bam
samtools sort -@ 24 -o S3460_samtools-sort.bam S3460.bam
samtools sort -@ 24 -o S3461_samtools-sort.bam S3461.bam
samtools sort -@ 24 -o S3476_samtools-sort.bam S3476.bam

#Used for picard.sh
