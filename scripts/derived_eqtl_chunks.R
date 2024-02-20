#!/usr/bin/env Rscript
# Author: M.J. Bonder (Adapted from 'https://github.com/sc-eQTLgen-consortium/WG3-pipeline-QTL/blob/master/scripts/createFeatureAnnotation.R' by M. Vochteloo)

# options parser
.libPaths("/usr/local/lib/R/site-library")
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list <- list(
  make_option(c("--feature_file"), action="store", default=NA, type='character',
              help="Input GTF file. (required)"),
  make_option(c("--n_genes"), action="store", default=100, type='integer',
              help="Number of genes to test in one QTL map job. "),
  make_option(c("--out_dir"), action="store", default=NA, type='character',
              help="Output main directory. (required)"))
opt <- parse_args(OptionParser(option_list=option_list))

print("Options in effect:")
for (name in names(opt)) {
	print(paste0("  --", name, " ", opt[[name]]))
}
print("")

# Very important to prevent scientific notation
options(scipen = 999)

# check date is provided
if (is.na(opt$feature_file)) {
  stop("required parameters feature_file must be provided. ")
}

print("Running eQTL window script")
print("reading and preparing gene features file..")

#### Read gene annotation
geneInfo <- read.delim(opt$feature_file, as.is=T, header=T)

##S2. make chunk file.
testCombinations <- NULL
#
nGenes <- opt$n_genes
startPos <- 0
endOffset <- 1000000000

lines <- NULL
for(chr in unique(geneInfo$chromosome)){
  print(paste0("  Creating eQTL windows for chromosome ", chr,".."))
  annotationRel <- geneInfo[which(geneInfo$chromosome == chr),]
  annotationRel <- annotationRel[order(annotationRel$start, annotationRel$end),]
  ##First go through the list to fix genes to they all touch.
  annotationRel$start[1] <- startPos
  for(i in 2:nrow(annotationRel)){
    if(i == nrow(annotationRel)){
      annotationRel$end[i] <- annotationRel$end[i]+endOffset
    }
    #If "overlapping" than we don't need to do anything.
    if((annotationRel$start[i] > annotationRel$end[i-1])){
      #print(i)
      distance <- (annotationRel$start[i]-annotationRel$end[i-1])/2
      annotationRel$start[i] <- ceiling(annotationRel$start[i]-distance)
      annotationRel$end[i-1] <- floor(annotationRel$end[i-1]+distance)
    }
  }
  chunks <- seq(1, nrow(annotationRel), nGenes)
  #Need to add the last one as a separate entry.
  if(chunks[length(chunks)] < nrow(annotationRel)){
    chunks <- c(chunks,(nrow(annotationRel)+1))
  }
  for(i in 1:(length(chunks)-1)){
    lines <- c(lines,(paste(chr, ":", annotationRel$start[chunks[i]], "-", annotationRel$end[(chunks[i + 1] - 1)], sep="")))
  }
}

print(paste0("Created ", length(lines), " chunks"))
write.table(lines,paste0(opt$out_dir,"eQTLChunkingFile.txt"), quote=F, sep="\t", row.names=F, col.names=F)