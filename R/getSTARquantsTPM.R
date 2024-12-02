#! /pkg/R-4.2.1-0/bin/R
#! /pkg/R-4.2.1-0/bin/Rscript

##################################################################

## Take STAR quants & convert to TPM, save as RDS

##################################################################


args <- commandArgs(trailingOnly=T)
if (length(args) != 3) stop('Usage: ./getSTARquantsTPM.R  <path_to_counts> <ref_genome> <out_path>')

##--- Load libs

libs <- c("data.table", "biomaRt", "GenomicFeatures", "SummarizedExperiment", "stringr")
invisible(suppressMessages(lapply(libs, library, character.only = TRUE)))

##--- Parsing inputs ---##

tab.path <-  args[1]
    ## get sample info from path
    names(tab.path) <- gsub('.ReadsPerGene.out.tab$', '', basename(tab.path))

ref_genome <-  args[2]

outpath  <- args[3]

    ## load txdb from ref ##
    if (ref_genome =="mm10") {
        pkg <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
        suppressMessages(require(pkg, character.only = TRUE))
        assign("txdb", eval(parse(text = pkg)))
    } else if (ref_genome=="galGal6") {
        txdb <- loadDb("/path/to/references/galGal6_Ensembl_GRCg7b.sqlite")
    } else if (ref_genome=="mm39") {
        txdb <- loadDb("/path/to/references/GENCODE_vM32_mm39_Txdb.sqlite")
    } else { 
       stop('genome not supported') 
        } 


###########

##--- Functions ---##


### compute gene length from sum exon length
getGeneLength <- function(txdb, gene.ids){
    ## get all exons
    exons <- GenomicFeatures::exonsBy(txdb, by = 'gene')
    ## remove versions from gene ids
    names(exons) <-  gsub("\\.[0-9]+$", "", names(exons))
    ## calculate gene length = sum exon length (in kb)
    len <- sapply(gene.ids, function(gene){
        if (!is.null(exons[[gene]])) {
            return(sum(width(exons[[gene]])/1000))
        }
    })
    
    return(len)
}

### build TPM count matrices
computeTPM <- function(count_mat, txdb){
    row.names(count_mat) <- gsub("\\.[0-9]+$", "", row.names(count_mat))
    
    ## estimate gene length = sum exon lengths in kb
    tx.len <- getGeneLength(txdb=txdb, gene.ids = row.names(count_mat))
    
    ## filter count mats for genes with estimated length
    count_mat <- count_mat[names(tx.len), ]
    
    ## compute TPM from RPK
    rpk <- count_mat/tx.len
    tpm <- t(t(rpk)*1e6/sum(rpk))
    return(tpm)
}

###########

##--- Load unstranded count matrix ---##

df <- data.table::fread(tab.path, select=c(1,2), skip = 4)
# convert into matrix where gene names = rownames
    mat <- matrix(df$V2, ncol=1, dimnames=list(df$V1, NULL))

##--- Convert count mat to SE obj ---##

sample_pattern <- ifelse(ref_genome=='galGal6', 'HH[:digit:]{2}', 'E[:digit:]{3}')

# se <- SummarizedExperiment(assays = list(counts = mat), 
#                             colData = DataFrame(row.names = names(tab.path),
#                                                stage = str_extract(names(tab.path), sample_pattern),
#                                                rep = sub('.*Rep', 'Rep', names(tab.path))))

##--- Compute TPM & store as SE object ---##
se_tpm <- SummarizedExperiment(assays = list(counts = computeTPM(count_mat = mat ,txdb = txdb)),
                               colData = DataFrame(row.names = names(tab.path),
                                                   stage = str_extract(names(tab.path), sample_pattern),
                                                   rep = sub('.*Rep', 'Rep', names(tab.path))))



saveRDS(se_tpm, file.path(outpath, paste0('SE_TPM_', names(tab.path), '.rds')))


