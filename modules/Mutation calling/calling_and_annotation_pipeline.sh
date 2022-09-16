#!/bin/bash

## Tell BASH to fail on error, pipefail, undefined vars
## and to print commands.
set -e
set -o pipefail
set -u
set -x

gatk_path=~/gatk-4.2.6.1/

## Add GATK to path
PATH=${PATH}:${gatk_path}

## Input tumor FASTQs
inputTumor_FASTQ_1=chr22.TCRBOA6-Tumor_1.fastq.gz
inputTumor_FASTQ_2=chr22.TCRBOA6-Tumor_2.fastq.gz

## Sample Name for Tumor
inputTumorSampleName=TCRBOA6-Tumor


## INPUT Normal FASTQs
inputNormal_FASTQ_1=chr22.TCRBOA6-Normal_1.fastq.gz
inputNormal_FASTQ_2=chr22.TCRBOA6-Normal_2.fastq.gz

## Sample name for normal
inputNormalSampleName=TCRBOA6-Normal

## Prefix the output with the following prefix:
## you can leave this as "" if you don't want a prefix
## or make it a path to use a path (but make sure the path exists).
outputPrefix="backup.chr22_"

threads=8
memory=30
readsPerBatch=100000000

## Align the normal sample
time bwa mem -Y -t ${threads} \
    -K ${readsPerBatch} \
    -R "@RG\tID:${inputNormalSampleName}-RG1\tLB:lib1\tPL:Illumina\tSM:${inputNormalSampleName}\tPU:${inputNormalSampleName}-RG1" \
    references/Homo_sapiens_assembly38.fasta ${inputNormal_FASTQ_1} ${inputNormal_FASTQ_2} \
    | samtools sort \
    -O BAM \
    -@ 2 \
    -o ${outputPrefix}${inputNormalSampleName}.bam

## Align the tumor sample
time bwa mem -Y -t ${threads} \
    -K ${readsPerBatch} \
    -R "@RG\tID:${inputTumorSampleName}-RG1\tLB:lib1\tPL:Illumina\tSM:${inputTumorSampleName}\tPU:${inputTumorSampleName}-RG1" \
    references/Homo_sapiens_assembly38.fasta ${inputTumor_FASTQ_1} ${inputTumor_FASTQ_2} \
    | samtools sort \
    -O BAM \
    -@ 2 \
    -o ${outputPrefix}${inputTumorSampleName}.bam 

## Mark duplicates in the normal sample
time gatk MarkDuplicates \
    --java-options -Xmx${memory}g \
    -I  ${outputPrefix}${inputNormalSampleName}.bam \
    -O  ${outputPrefix}${inputNormalSampleName}.markdups.bam \
    -M  ${outputPrefix}${inputNormalSampleName}.markdups.metrics.txt

## Mark duplicates in the tumor sample
time gatk MarkDuplicates \
    --java-options -Xmx${memory}g \
    -I ${outputPrefix}${inputTumorSampleName}.bam \
    -O ${outputPrefix}${inputTumorSampleName}.markdups.bam \
    -M ${outputPrefix}${inputTumorSampleName}.markdups.metrics.txt


time gatk BaseRecalibrator \
    --java-options -Xmx${memory}g \
    --input ${outputPrefix}${inputNormalSampleName}.markdups.bam \
    --output ${outputPrefix}${inputNormalSampleName}.markdups.BQSR-REPORT.txt \
    --known-sites references/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --reference references/Homo_sapiens_assembly38.fasta


time gatk BaseRecalibrator \
    --java-options -Xmx${memory}g \
    --input ${outputPrefix}${inputTumorSampleName}.bam \
    --output ${outputPrefix}${inputTumorSampleName}.markdups.BQSR-REPORT.txt \
    --known-sites references/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --reference references/Homo_sapiens_assembly38.fasta

## Apply BQSR to normal BAM
time gatk ApplyBQSR \
    --java-options -Xmx${memory}g \
    -R references/Homo_sapiens_assembly38.fasta \
    -I ${outputPrefix}${inputNormalSampleName}.markdups.bam \
    --bqsr-recal-file ${outputPrefix}${inputNormalSampleName}.markdups.BQSR-REPORT.txt \
    -O ${outputPrefix}${inputNormalSampleName}.markdups.baseRecal.bam

## Apply BQSR to tumor BAM
time gatk ApplyBQSR \
    --java-options -Xmx${memory}g \
    -R references/Homo_sapiens_assembly38.fasta \
    -I ${outputPrefix}${inputTumorSampleName}.markdups.bam \
    --bqsr-recal-file ${outputPrefix}${inputTumorSampleName}.markdups.BQSR-REPORT.txt \
    -O ${outputPrefix}${inputTumorSampleName}.markdups.baseRecal.bam

## Index our normal BAM
time samtools index ${outputPrefix}${inputNormalSampleName}.markdups.baseRecal.bam

## Index our tumor BAM
time samtools index ${outputPrefix}${inputTumorSampleName}.markdups.baseRecal.bam


time gatk Mutect2 \
    -R ~/references/Homo_sapiens_assembly38.fasta \
    --input ${outputPrefix}${inputTumorSampleName}.markdups.baseRecal.bam \
    --tumor-sample ${inputTumorSampleName} \
    --input ${outputPrefix}${inputNormalSampleName}.markdups.baseRecal.bam \
    --normal-sample ${inputNormalSampleName} \
    -L chr22 \
    --output ${outputPrefix}${inputTumorSampleName}.${inputNormalSampleName}.vcf

time gatk Funcotator \
     --variant ${outputPrefix}${inputTumorSampleName}.${inputNormalSampleName}.vcf \
     --reference references/Homo_sapiens_assembly38.fasta \
     --ref-version hg38 \
     --data-sources-path funcotator_dataSources.v1.7.20200521s \
     --output-file-format MAF \
     --output ${outputPrefix}${inputTumorSampleName}.${inputNormalSampleName}.funcotated.maf
