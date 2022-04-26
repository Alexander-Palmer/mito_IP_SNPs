# mito_IP_SNPs
Quantifying the incidence of single nucleotide polymorphisms (SNPs) of neuron-specific mitochondrial from mouse brain tissue \
A collaboration between Ina Kirmes (Zuryn Lab) and the Bartlett Lab, Queensland Brain Institute (QBI)

This repository describes our workflow for SNP analysis, from the filtering and processing of raw sequencing data to the creation of count tables for biological intepretation

## Dependencies
* fastQC (v0.10.1)
* cutadapt (v1.3)
* BBMap (v38.95)
* bwa mem (v0.7.17)
* Samtools (v0.1.19)
* Picard (v2.27.1)
* qualimap (v2.2.1)
* mutect2 (GATK, v4.2.4.1)
* NCBI genome workbench (v3.7.1)

## Pipeline
_M. musculus_ reference genome: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.23/ \
_M. musculus_ mtDNA annotation: **mtDNA_features.csv**

1. **fastQC.sh**: _initial assessment of read quality_
2. **cutadapt.sh**: _removal of sequencing adapters_
3. **fastQC.sh**: _assessment of read quality folowing adapter removal_
4. **bbtools.sh**: _removal of reads with unmatched names_
5. **bwa.sh**: _index and alignment of reads to the_ M. musculus _genome_
6. **samtools.sh**: _sort reads by coordinate for duplicate marking_
7. **picard.sh**: _mark duplicates and remove empty reads_
8. **qualimap.sh**: _assess quality of processed reads_
9. **mutect2.sh**: _create genome reference index, add arbitrary read groups, validate and index input, and run mutect2 program_
10. **NCBI genome workbench**: _export count table of SNPs including relevant information_
11. **statistical_analysis_and_visualisation.R**: _investigate SNPs significantly enriched between libraries_ 

## Citation
Manuscript currently pre-submission
Written by Alexander Palmer for Ina Kirmes (Zuryn Lab) and Bartlett Lab, Queensland Brain Institute (QBI), Australia (2022)
