#! /bin/bash
## STAR v2.7.9a ref indexing with annotation ##

set -euo pipefail 


REF="/path/to/references"

mmGenome="$REF"/GRCm39.primary_assembly.genome.fa
mmAnno="$REF"/gencode.vM32.primary_assembly.annotation.gtf

ggGenome="$REF"/galGal6.GRCg7b.dna.toplevel.fa
ggAnno="$REF"/galGal6.GRCg7b.109.gtf

STAR --runMode genomeGenerate \
    --runThreadN 16 \
    --genomeDir "$REF"/STARv2.7.9a/mm39 \
    --genomeFastaFiles "$mmGenome" --sjdbGTFfile "$mmAnno"

STAR --runMode genomeGenerate \
    --runThreadN 16 \
    --genomeDir "$REF"/STARv2.7.9a/galGal6 \
    --genomeFastaFiles "$ggGenome" --sjdbGTFfile "$ggAnno"