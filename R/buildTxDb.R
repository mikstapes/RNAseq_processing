#########################################
### Build custom TxDb objects
#########################################


libs <- c("biomaRt", "GenomicFeatures", "SummarizedExperiment")
invisible(suppressMessages(lapply(libs, library, character.only = TRUE)))

## GFF from GENCODE
mm39_gtf <-'/project/ibr_procs/data/_processing/rna/references/references/gencode.vM32.primary_assembly.annotation.gtf'
txdb.mmus <- makeTxDbFromGFF(mm39_gtf, dataSource = "GENCODE vM32", organism = "Mus musculus")

## GRCg7b
txdb.ggal <- makeTxDbFromEnsembl(organism = "Gallus gallus", server = "ensembldb.ensembl.org")

## Save for future use
saveDb(txdb.mmus, file = "/project/ibr_procs/data/_processing/rna/references/GENCODE_vM32_mm39_Txdb.sqlite")
saveDb(txdb.ggal, file = "/project/ibr_procs/data/_processing/rna/references/galGal6_Ensembl_GRCg7b.sqlite")