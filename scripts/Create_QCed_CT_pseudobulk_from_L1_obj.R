shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(Seurat))
shhh(library(matrixStats))
shhh(library(textTinyR))
shhh(library(pbapply))
shhh(library(optparse))

option_list = list(
  make_option(c("--in_rds"), action="store", default=NA, type='character',
              help="Qced L1 object to make L2 data for."),
  make_option(c("--out_dir"), action="store", default=NA, type='character',
              help="Output main directory"),
  make_option(c("--aggregate_fun"), action="store", default='mean', type='character',
              help="sum or mean --> if mean, we need to normalize before."))
opt = parse_args(OptionParser(option_list=option_list))

pbmc = readRDS(opt$in_rds)
## ##Here I assume that the RDS file is in the right main directory.
##expressionCovs <- read.delim(gsub(".Qced.Normalized.SCs.Rds",".qtlInput.Pcs.txt","/scratch/hb-sceqtlgen/oneK1K/input/L1/CD4_T.Qced.Normalized.SCs.Rds"),as.is=T)
## samplesToKeep = expressionCovs[,1]


out_dir = opt$out_dir
aggregate_fun = opt$aggregate_fun

cellCts = pbmc@meta.data$predicted.celltype.l2
ctList = unique(na.omit(cellCts))

for(ct in ctList){

  ##Grab meta.data
  meta.d.full = as.data.frame(pbmc@meta.data[which(cellCts == ct),])
  meta.d = meta.d.full
  rownames(meta.d) = NULL
  meta.d = unique(meta.d)
  ##Grab the countmatrix
  countMatrixFull <- GetAssayData(pbmc[,which(cellCts == ct)], slot = "counts")
  countMatrix = countMatrixFull[which(rowSums(countMatrixFull)!=0),]
  
  rm(countMatrixFull)
  
  IDs = meta.d.full$Donor_Pool
  unique_ID_list = unique(IDs)
  
  normCountMatrix = countMatrix ##To store in the new object in.
  
  # Info
  print(paste0('Calculating pseudo-bulk for ',ct,' expression matrix using: ', aggregate_fun))
  print(paste0('   ### n_cells: ', ncol(normCountMatrix)))
  print(paste0('   ### n_genes (sc-data): ', nrow(normCountMatrix)))
  cat('\n')
  
  if(sum(cellCts == ct) <= 1){
    print(paste(ct,'has 1 or 0 cells, skip'))
    next
  }
  # PF
  sampleSumInfo = colSums(normCountMatrix)
  meanSampleSum = mean(sampleSumInfo)
  sampleScale = sampleSumInfo/meanSampleSum
  
  # Scale to match meanRowSums
  system.time(normCountMatrix@x <- normCountMatrix@x / rep.int(sampleScale, diff(normCountMatrix@p))) #divide each column by sampleScale
  
  # Log
  normCountMatrix@x  = log(normCountMatrix@x +1)
  gc()
  # PF
  sampleSumInfo = colSums(normCountMatrix)
  meanSampleSum = mean(sampleSumInfo)
  sampleScale = sampleSumInfo/meanSampleSum
  
  # Scale to match meanRowSums
  system.time(normCountMatrix@x <- normCountMatrix@x / rep.int(sampleScale, diff(normCountMatrix@p))) #divide each column by sampleScale
  
  ## Single cell object
  print(paste0("Writing QC-ed and Normalized SC data: ",ct))
  rownames(meta.d.full)=meta.d.full$Barcode
  gc()
  seuratObj <- CreateSeuratObject(countMatrix, assay = "RNA",meta.data = meta.d.full)
  seuratObj[["data"]] = CreateAssayObject(data=normCountMatrix)
  
  saveRDS(seuratObj,paste0(out_dir,"/",gsub(" ","_",ct),".Qced.Normalized.SCs.Rds"))
  rm(seuratObj,meta.d.full,countMatrix)
  
  print('Aggregating count matrix using lapply + textTinyR::sparse_Means() ...')
  ##normCountMatrix <- as(normCountMatrix, "dgCMatrix")
  
  system.time(aggregate_normCountMatrix <- as.data.frame(pblapply(unique_ID_list, FUN = function(x){sparse_Means(normCountMatrix[,IDs == x, drop = FALSE], rowMeans = TRUE)})))
  cellCount <- pblapply(unique_ID_list, FUN = function(x){ncol(normCountMatrix[,IDs == x, drop = FALSE])})
  colnames(aggregate_normCountMatrix) <- names(cellCount) <- unique_ID_list
  rownames(aggregate_normCountMatrix) <- rownames(normCountMatrix)
  
  cellCount = unlist(cellCount)
  
  print(paste0("Writing pseudo bulk data: ",ct))
  ##Write out psuedo-bulk.
  write.table(aggregate_normCountMatrix,paste0(out_dir,"/",gsub(" ","_",ct),".Exp.txt"),quote=F,sep="\t",col.names=NA)
  
  ##Write covariates global covariates.
  meta.d = cbind(meta.d,cellCount[match(meta.d$Donor_Pool,names(cellCount))]) ##Add cell numbers.
  colnames(meta.d)[ncol(meta.d)]="CellCount" ##Fix nameing.
  write.table(meta.d,paste0(out_dir,"/",gsub(" ","_",ct),".covariates.txt"),quote=F,sep="\t",row.names=F)
  
  ##Write QTL input.
  
  ##For QTL first filter, take only obersvations from 5 cells, and rows that are varying.
  aggregate_normCountMatrix = aggregate_normCountMatrix[,which(colnames(aggregate_normCountMatrix) %in% names(which(cellCount>4)))]
  ##Only write-out files when there are more then 15 observations.
  if(!is.null(dim(aggregate_normCountMatrix))){
    if(ncol(aggregate_normCountMatrix)>15){
      ##Filter based on gene variance
      rowVarInfo = rowVars(as.matrix(aggregate_normCountMatrix))
      aggregate_normCountMatrix = aggregate_normCountMatrix[which(rowVarInfo!=0),]
      
      ## Do inverse normal transform per gene.
      for(rN in 1:nrow(aggregate_normCountMatrix)){
        aggregate_normCountMatrix[rN,] = qnorm((rank(aggregate_normCountMatrix[rN,],na.last="keep")-0.5)/sum(!is.na(aggregate_normCountMatrix[rN,])))
      }
      ##Do PCA, and select first 10 components.
      pcOut = prcomp(t(aggregate_normCountMatrix))
      covOut = pcOut$x[,1:10]
      ##Write out PCs and input matrix for QTL.
      write.table(aggregate_normCountMatrix,paste0(out_dir,"/",gsub(" ","_",ct),".qtlInput.txt"),quote=F,sep="\t",col.names=NA)
      write.table(covOut,paste0(out_dir,"/",gsub(" ","_",ct),".qtlInput.Pcs.txt"),quote=F,sep="\t",col.names=NA)
    }
  }
  print("done.")
  rm(normCountMatrix,aggregate_normCountMatrix)
  gc()
}
