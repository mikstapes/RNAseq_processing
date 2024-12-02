import os
import configparser

# Load configuration
configfile: "config.yaml"

# Define global variables
HERE = os.path.realpath(".")  # where HERE = working dir containing all ref_dirs and sample files

# Gather sample information from configuration
builds = config["builds"]

# Generate list of libraries based on build
def get_samples(build_name):
    samples_file = builds[build_name]["fastq"]
    with open(os.path.join(HERE, samples_file)) as f:
        return [line.strip() for line in f]

# Define all targets
all_samples_mm39 = expand(
    os.path.join(HERE, "_bam/{sample}Aligned.sortedByCoord.out.bam"),
    sample=get_samples("mm39")
)
all_samples_galGal6 = expand(
    os.path.join(HERE, "_bam/{sample}Aligned.sortedByCoord.out.bam"),
    sample=get_samples("galGal6")
)

star_idx = expand(
    os.path.join(config["ref_dir"], 'STARv2.7.9a/{build}/exonInfo.tab'),
    build=builds
)

rule all:
    input:
        star_idx + all_samples_mm39 + all_samples_galGal6 

rule build:
    input:
        fasta=lambda wc: os.path.join(config["ref_dir"], builds[wc.build]["fasta"])
        anno=lambda wc: os.path.join(config["ref_dir"], builds[wc.build]["GTF"])
    output: lambda wc: directory('%s/STARv2.7.9a/{wc.build}' % config["ref_dir"])
    shell:
    """
    STAR --runMode genomeGenerate --runThreadN 16 \
         --genomeDir {output} \
         --genomeFastaFiles {input.fasta} \
         --sjdbGTFfile {input.anno} 
    """

rule align:
    input:
        fastq=os.path.join(HERE, "_fastq/{sample}.fastq.gz"),
        genomeDir=lambda wildcards: os.path.join(config["ref_dir"], "STARv2.7.9a", wildcards.build),
        gtf=lambda wildcards: os.path.join(config["ref_dir"], builds[wildcards.build]["GTF"])
    output:
        bam=os.path.join(HERE, "_bam/{sample}Aligned.sortedByCoord.out.bam")
    params:
        pref_dirix=os.path.join(HERE, "_bam/{sample}")
    threads: 8
    shell:
        """
        STAR --genomeDir {input.genomeDir} --readFilesCommand zcat \
             --readFilesIn {input.fastq} \
             --runThreadN {threads} \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode GeneCounts \
             --sjdbGTFfile {input.gtf} \
             --outFileNamePref_dirix {params.pref_dirix}
        """


