#!/usr/bin/env Rscript

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--wg1_data"), action="store", default=NA, type='character',
              help="WG1 seurat object out (singlets and filtered)."),
  make_option(c("--metadata"), action="store", default=NA, type='character',
              help="MAD filter and donor annotation information."),
  make_option(c("--cell_level"), action="store", default="l1", type='character',
              help="Azimuth l1 or l2."),
  make_option(c("--aggregate_fun"), action="store", default='mean', type='character',
              help="sum or mean --> if mean, we need to normalize before."),
  make_option(c("--out_dir"), action="store", default='wg1-preprocessing/00.pseudobulk_by_celltype/', type='character',
              help="Output main directory"))
opt = parse_args(OptionParser(option_list=option_list))

# Packages
shhh(library(data.table))
shhh(library(Seurat))
shhh(library(matrixStats))
shhh(library(textTinyR))
shhh(library(pbapply))

#################### Set Variables and load Data #################### 

# Report
print(paste0('input1: ', opt$wg1_data))
print(paste0('input2: ', opt$metadata))
print(paste0('Azimuth level: ', opt$cell_level))
print(paste0('Aggegate function: ', opt$aggregate_fun))
cat('\n\n')

#################### Read data ####################
# Seurat object
pbmc_fn <- opt$wg1_data
print(paste0('Reading pbmc seurat object file in: ',pbmc_fn))
system.time(pbmc <- readRDS(pbmc_fn))

# metadata object
pbmcMD_fd <- opt$metadata
print(paste0('Reading pbmc metadata object file in: ',pbmcMD_fd))
system.time(pbmcMetaD <- readRDS(pbmcMD_fd))

#Filtering to cells to keep.
pbmcMetaD  = pbmcMetaD[which(pbmcMetaD$tag=="NotOutlier"),]
pbmcMetaD  = pbmcMetaD[which(pbmcMetaD$predicted.celltype.l1 %in% c("CD4T","Mono","CD8T","NK","B","T_other","DC")),]

#Make sure the entries match.
pbmcMetaD  = pbmcMetaD[which(pbmcMetaD$Barcode %in% colnames(pbmc)),]
pbmc = pbmc[,which(colnames(pbmc) %in% pbmcMetaD$Barcode)]

pbmc = pbmc[,order(colnames(pbmc))]
pbmcMetaD = pbmcMetaD[match(colnames(pbmc),pbmcMetaD$Barcode),]
pbmcMetaD = cbind(pbmcMetaD,paste0(pbmcMetaD$predicted.celltype.l1,";;",pbmcMetaD$Assignment))
colnames(pbmcMetaD)[ncol(pbmcMetaD)] ="CT_Donor"

aggregate_fun = opt$aggregate_fun

if(aggregate_fun!='mean'){
  print('Not implemented yet...')
  stop()
}
print(all(colnames(pbmc)==pbmcMetaD$Barcode))
if(all(colnames(pbmc)==pbmcMetaD$Barcode)){
  pbmc@meta.data = pbmcMetaD
  
  ## grab the countmatrix
  countMatrix <- GetAssayData(pbmc, slot = "counts")
  countMatrix = countMatrix[which(rowSums(countMatrix)!=0),]
  IDs = pbmc@meta.data$CT_Donor
  unique_ID_list = unique(pbmc@meta.data$CT_Donor)
  
  # Info
  n_cells <- ncol(countMatrix)
  n_genes <- nrow(countMatrix)
  print(paste0('Calculating pseudo-bulk expression matrix using: ', aggregate_fun))
  print(paste0('   ### n_cells: ', n_cells))
  print(paste0('   ### n_genes (sc-data): ', n_genes))
  cat('\n')

  print('Normalizing the sc-data..')
  cat('\n')
  
  # PF
  sampleSumInfo = colSums(countMatrix)
  meanSampleSum = mean(sampleSumInfo)
  sampleScale = sampleSumInfo/meanSampleSum
  
  # Scale to match meanRowSums
  system.time(countMatrix <- sweep(countMatrix, 2, sampleScale, FUN="/")) #divide each column by sampleScale
  
  # Log
  countMatrix = log(countMatrix+1)
  
  # PF
  sampleSumInfo = colSums(countMatrix)
  meanSampleSum = mean(sampleSumInfo)
  sampleScale = sampleSumInfo/meanSampleSum
  
  # Scale to match meanRowSums
  system.time(countMatrix <- sweep(countMatrix, 2, sampleScale, FUN="/")) #divide each column by sampleScale
  
  print('Aggregating count matrix using lapply + textTinyR::sparse_Means() ...')
  countMatrix <- as(countMatrix, "dgCMatrix")
  
  system.time(aggregate_countMatrix <- as.data.frame(pblapply(unique_ID_list, FUN = function(x){sparse_Means(countMatrix[,IDs == x, drop = FALSE], rowMeans = TRUE)})))
  cellcount <- pblapply(unique_ID_list, FUN = function(x){ncol(countMatrix[,IDs == x, drop = FALSE])})
  colnames(aggregate_countMatrix) <- names(cellcount) <- unique_ID_list
  rownames(aggregate_countMatrix) <- rownames(countMatrix)
  
  ##Filter to 10 cells minimum.
  cellCount = unlist(cellcount)
  write.table(cellCount,paste0(opt$out_dir,"/cellCounts.txt"),quote=F,sep="\t",col.names=NA)
  aggregate_countMatrix = aggregate_countMatrix[,which(colnames(aggregate_countMatrix) %in% names(which(cellCount>9)))]
  ##Split per L1 cell type
  ctList = unique(unlist(lapply(strsplit(colnames(aggregate_countMatrix),split=";;"),'[[',1)))
  for(ct in ctList){
    ##select relevant columns
    selMatrix = aggregate_countMatrix[,grep(ct,colnames(aggregate_countMatrix))]
    colnames(selMatrix) = unlist(lapply(strsplit(colnames(selMatrix),split=";;"),'[[',2))
    if(ncol(selMatrix)>15){
      print(paste0("Writing: ",ct))
      ##Drop rows that are not varying, which will include genes that are only zero.
      rowVarInfo = rowVars(as.matrix(selMatrix))
      selMatrix = selMatrix[which(rowVarInfo!=0),]
      ## Do inverse normal transform per gene.
      for(rN in 1:nrow(selMatrix)){
         selMatrix[rN,] = qnorm((rank(selMatrix[rN,],na.last="keep")-0.5)/sum(!is.na(selMatrix[rN,])))
      }
      ##Do PCA, and select first 10 components.
      pcOut = prcomp(t(selMatrix))
      covOut = pcOut$x[,1:10]
      ##Write out PCs and input matrix for QTL.
      write.table(selMatrix,paste0(opt$out_dir,"/",ct,".inputExpression.txt"),quote=F,sep="\t",col.names=NA)
      write.table(covOut,paste0(opt$out_dir,"/",ct,".covariates.txt"),quote=F,sep="\t",col.names=NA)
    }
  }
}
