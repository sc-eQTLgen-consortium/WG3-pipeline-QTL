#!/usr/bin/env Rscript

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--in_gtf"), action="store", default=NA, type='character',
              help="Input GTF file. (required)"),
  make_option(c("--n_genes"), action="store", default=100, type='integer',
              help="Number of genes to test in one QTL map job. "),
  make_option(c("--feature_name"), action="store", default="GeneName", type='character',
              help="Number of genes to test in one QTL map job. "),
  make_option(c("--out_dir"), action="store", default=NA, type='character',
              help="Output main directory. (required)"))
opt = parse_args(OptionParser(option_list=option_list))

# check date is provided
if ((is.na(opt$in_gtf) || is.na(opt$out_dir))) {
  stop("required parameters (in_gtf & out_dir) must be provided. ")
}

#### Read gene annotation
gtfInfo = read.delim(opt$in_gtf,as.is=T,comment="#",header=F)

##S1. make annotation file.

##Select only genes.
gtfInfo = gtfInfo[which(gtfInfo$V3=="gene"),]
##Make sure genes are mapped to numeric
gtfInfo$V1 = gsub("chr","",gtfInfo$V1)
##Select only main chromosome mappings.
gtfInfo = gtfInfo[which(gtfInfo$V1 %in% c(1:22,"X","MT", "M","Y")),]

##Extend.
gtfInfo["ENSG"] = gsub("gene_id ","",unlist(lapply(strsplit(gtfInfo$V9,split = ";"),"[[",1)))
gtfInfo["Gene_Name"] = trimws(gsub("gene_name","",unlist(lapply(strsplit(gtfInfo$V9,split = ";"),"[[",3))))
gtfInfo["biotype"] = trimws(gsub("gene_biotype","",unlist(lapply(strsplit(gtfInfo$V9,split = ";"),"[[",5))))

if(opt$feature_name=="ENSG"){
  geneInfo = gtfInfo[,c(10,1,4,5,11,12)]
  colnames(geneInfo)=c("feature_id","chromosome","start","end","GeneName","biotype")
} else {
  geneInfo = gtfInfo[,c(11,1,4,5,10,12)]
  colnames(geneInfo)=c("feature_id","chromosome","start","end","ENSG","biotype")
}

##drop all genes that have the same names.
toDrop = geneInfo$feature_id[which(duplicated(geneInfo$feature_id))]
geneInfo = geneInfo[which(!(geneInfo$feature_id %in% toDrop)),]

write.table(geneInfo,paste0(opt$out_dir,"/LimixAnnotationFile.txt"),quote = F,sep="\t",row.names = F)

##S2. make chunk file.
testCombinations = NULL
#
nGenes = opt$n_genes
startPos = 0
endOffset = 1000000000

lines = NULL
for(chr in unique(geneInfo$chromosome)){
  #print(chr)
  annotationRel = geneInfo[which(geneInfo$chromosome==chr),]
  annotationRel = annotationRel[order(annotationRel$start,annotationRel$end),]
  ##First go through the list to fix genes to they all touch.
  annotationRel$start[1] = startPos
  for(i in 2:nrow(annotationRel)){
    if(i == nrow(annotationRel)){
      annotationRel$end[i] = annotationRel$end[i]+endOffset
    }
    #If "overlapping" than we don't need to do anything.
    if((annotationRel$start[i]>annotationRel$end[i-1])){
      #print(i)
      distance = (annotationRel$start[i]-annotationRel$end[i-1])/2
      annotationRel$start[i] = ceiling(annotationRel$start[i]-distance)
      annotationRel$end[i-1] = floor(annotationRel$end[i-1]+distance)
    }
  }
  chunks = seq(1, nrow(annotationRel),nGenes)
  #Need to add the last one as a separate entry.
  if(chunks[length(chunks)] < nrow(annotationRel)){
    chunks = c(chunks,(nrow(annotationRel)+1))
  }
  for(i in 1:(length(chunks)-1)){
    lines = c(lines,(paste(chr,":",annotationRel$start[chunks[i]],"-",annotationRel$end[(chunks[i+1]-1)],sep="")))
  }
}
write.table(lines,paste0(opt$out_dir,"/ChunkingFile.txt"),quote = F,sep="\t",row.names = F,col.names=F)
