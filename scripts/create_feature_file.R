#!/usr/bin/env Rscript
# Author: M.J. Bonder (Adapted from 'https://github.com/sc-eQTLgen-consortium/WG3-pipeline-QTL/blob/master/scripts/createFeatureAnnotation.R' by M. Vochteloo)

# options parser
.libPaths("/usr/local/lib/R/site-library")
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list <- list(
  make_option(c("--in_gtf"), action="store", default=NA, type='character', 
              help="Input GTF file. (required)"),
  make_option(c("--feature_name"), action="store", default="GeneName", 
              help="The feature_id used in the expression matrix."), 
  make_option(c("--biotype_flag"), action="store", default="gene_biotype", 
              help="The feature_id used in the expression matrix."), 
  make_option(c("--autosomes_only"), action="store_true", default=FALSE, type='character', 
              help="Flag to only keep autosomal genes, default false. "), 
  make_option(c("--out_dir"), action="store", default=NA, type='character', 
              help="Output main directory. (required)"))
opt <- parse_args(OptionParser(option_list=option_list))

print("Options in effect:")
for (name in names(opt)) {
	print(paste0("  --", name, " ", opt[[name]]))
}
print("")

# check date is provided
if ((is.na(opt$in_gtf) || is.na(opt$out_dir))) {
  stop("required parameters (in_gtf & out_dir) must be provided. ")
}

#### Read gene annotation
gtfInfo <- read.delim(opt$in_gtf, as.is=T, comment="#", header=F)

##S1. make annotation file.

##Select only genes.
gtfInfo <- gtfInfo[which(gtfInfo$V3=="gene"),]
##Make sure genes are mapped to numeric
gtfInfo$V1 <- gsub("chr", "", gtfInfo$V1)
##Select only main chromosome mappings.
gtfInfo <- gtfInfo[which(gtfInfo$V1 %in% c(1:22, "X", "MT", "M", "Y")),]

##Extend matrix with relevant entries.
partialMatrix <- strsplit(gtfInfo$V9, split="; ")

ensgIdLoc <- which(startsWith(partialMatrix[[1]], "gene_id"))
gtfInfo["ENSG"] <- gsub( "gene_id ", "", unlist(lapply(partialMatrix, "[[", ensgIdLoc)))
nameLoc <- which(startsWith(partialMatrix[[1]], "gene_name"))
gtfInfo["Gene_Name"] <- gsub( "gene_name ", "", unlist(lapply(partialMatrix, "[[", nameLoc)))
if(opt$biotype_flag!="None"){
	btLoc <- which(startsWith(partialMatrix[[1]], opt$biotype_flag))
	gtfInfo["biotype"] <- gsub(paste(opt$biotype_flag, " ", sep=""), "", unlist(lapply(partialMatrix, "[[", btLoc)))
} else {
	gtfInfo["biotype"] <- NA
}


if(opt$feature_name=="ENSG"){
  geneInfo <- gtfInfo[, c(10, 1, 4, 5, 11, 12)]
  colnames(geneInfo) <- c("feature_id", "chromosome", "start", "end", "GeneName", "biotype")
} else {
  geneInfo <- gtfInfo[, c(11, 1, 4, 5, 10, 12)]
  colnames(geneInfo) <- c("feature_id", "chromosome", "start", "end", "ENSG", "biotype")
}

if(opt$autosomes_only){
  geneInfo <- geneInfo[which(geneInfo$chromosome %in% c(1:22)),]
}

##drop all genes that have the same names, to avoid having issues with matching and mapping after.
toDrop <- geneInfo$feature_id[which(duplicated(geneInfo$feature_id))]
geneInfo <- geneInfo[which(!(geneInfo$feature_id %in% toDrop)),]

write.table(geneInfo, paste0(opt$out_dir, "/LimixAnnotationFile.txt"), quote=F, sep="\t", row.names=F)