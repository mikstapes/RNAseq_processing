#! /bin/bash

set -eo pipefail 

build="mm39 galGal6"
HERE=$(realpath ..)
REF=$HERE/references

build="mm39 galGal6"

for B in $build; do
    if [ $B == "mm39" ]; then
        GTF="$REF"/gencode.vM32.primary_assembly.annotation.gtf
        for LIB in $(cat $HERE/_fastq/mm39); do
            STAR --genomeDir "$REF"/STARv2.7.9a/"$B" --readFilesCommand zcat \
            --readFilesIn $HERE/_fastq/"$LIB".fastq.gz \
            --runThreadN 8 \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            --sjdbGTFfile "$GTF" \
            --outFileNamePrefix "$HERE"/_bam/"$LIB"
        done
    else 
        GTF="$REF"/galGal6.GRCg7b.109.gtf
        for LIB in $(cat $HERE/_fastq/galGal6); do
            STAR --genomeDir "$REF"/STARv2.7.9a/"$B" --readFilesCommand zcat \
            --readFilesIn "$HERE"/_fastq/"$LIB".fastq.gz   \
            --runThreadN 8 \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            --sjdbGTFfile "$GTF" \
            --outFileNamePrefix "$HERE"/_bam/"$LIB"
        done
    fi

done