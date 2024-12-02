# RNAseq_processing
Currently repo consists of (predominantly) codes for RNAseq processing for thesis project, as part of this [preprint](https://www.biorxiv.org/content/10.1101/2024.05.13.590087v1) (See Figure 1). 
Will update for more generalized workflows in the future.

---


1. Codes for data processing from fastq are found under `workflow/`. Where there are:
   - Several shell scripts for running STAR alignment with input gene annotations for generating counts in `workflow/shell_scripts`
   - A generalized Snakemake file for the STAR workflow, from building genome indices to alignment & gene counts `STARalign.smk`

2. Downstream analysis in R are found under the appropriately named dir `R`, consisting of R codes for tasks such as:
   - Building custom TxDb object to be used as referenced annotations: `buildTxDb.R`
   - Generating TPM counts from STAR output: `getSTARquantsTPM.R`
   - Characterizing the orthology relationship in genes expressed in the heart of chicken and mouse (counts obtained from 2.). Orthology annotations are obtained from the public Ensembl database and provided via biomaRt: `GeneOrthology.ipynb`
   - Comparing differentially expressed genes, specifically heart-specific genes (which are one-to-one orthologs, as determined from 3.), between chicken and mouse using DESeq2. Differentially expressed genes are annoted by relevant Gene Ontolgy: `OGDiffExpAnalysis.ipynb`
