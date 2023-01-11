#!/usr/bin/env Rscript

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--wg1_data"), action="store", default=NA, type='character',
              help="WG1 seurat object out (singlets and filtered)."),
  make_option(c("--metadata"), action="store", default=NA, type='character',
              help="MAD filter and donor annotation information."),
    make_option(c("--wg1_psam"), action="store", default=NA, type='character',
              help="psam file created in WG1."),
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
print(paste0('input2: ', opt$wg1_psam))
print(paste0('input3: ', opt$metadata))
print(paste0('Azimuth level: ', opt$cell_level))
print(paste0('Aggegate function: ', opt$aggregate_fun))
cat('\n\n')

#################### Read data ####################
# Seurat object
pbmc_fn <- opt$wg1_data
print(paste0('Reading seurat object: ',pbmc_fn))
system.time(pbmc <- readRDS(pbmc_fn))

# Seurat object
pbmc_psam_fn <- opt$wg1_psam
print(paste0('Reading psam info: ',pbmc_psam_fn))
system.time(psamMetaD <- read.delim(pbmc_psam_fn,as.is=T,check.names=F))

# metadata object
pbmc_metadata_fn <- opt$metadata
print(paste0('Reading pbmc WG1-WG2 metadata: ',pbmc_metadata_fn))
system.time(pbmcMetaD <- readRDS(pbmc_metadata_fn))
pbmcMetaD = as.data.frame(pbmcMetaD)
rm(pbmc_fn,pbmc_metadata_fn,pbmc_psam_fn)

##Make combined L1 annotation.
pbmcMetaD["celltype.l1"] = pbmcMetaD$predicted.celltype.l1
pbmcMetaD$celltype.l1[which(pbmcMetaD$celltype.l1!=pbmcMetaD$scpred.l1.prediction)] = NA

#Filtering to cells to keep.
pbmcMetaD  = pbmcMetaD[which(pbmcMetaD$tag=="NotOutlier"),]
pbmcMetaD  = pbmcMetaD[which(pbmcMetaD$celltype.l1 != "Doublet"),]
#pbmcMetaD  = pbmcMetaD[which(pbmcMetaD$celltype.l1 %in% c("CD4_T","Mono","CD8_T","NK","B","DC","other_T")),]

#Make sure the entries match.
pbmcMetaD  = pbmcMetaD[which(pbmcMetaD$Barcode %in% colnames(pbmc)),]
pbmc = pbmc[,which(colnames(pbmc) %in% pbmcMetaD$Barcode)]

pbmc = pbmc[,order(colnames(pbmc))]
pbmcMetaD = pbmcMetaD[match(colnames(pbmc),pbmcMetaD$Barcode),]
pbmcMetaD = cbind(pbmcMetaD,paste0(pbmcMetaD$celltype.l1,";;",pbmcMetaD$Pool,";;",pbmcMetaD$Assignment))
colnames(pbmcMetaD)[ncol(pbmcMetaD)] ="CT_Pool_Donor"

pbmcMetaD = cbind(pbmcMetaD,paste0(pbmcMetaD$Assignment,";;",pbmcMetaD$Pool))
colnames(pbmcMetaD)[ncol(pbmcMetaD)] ="Donor_Pool"

##Extend meta_data with info from psam.
pbmcMetaD = cbind(pbmcMetaD,psamMetaD[match(pbmcMetaD$Assignment,psamMetaD$IID),])
relCol = c("Donor_Pool",colnames(psamMetaD))
aggregate_fun = opt$aggregate_fun

if(aggregate_fun!='mean'){
  print('Not implemented yet...')
  stop()
}
print(all(colnames(pbmc)==pbmcMetaD$Barcode))

write.table(cbind(pbmc@meta.data,pbmcMetaD),paste0(opt$out_dir,"/AllMetaData.debug.txt"),quote=F,sep="\t",row.names=F)

if(all(colnames(pbmc)==pbmcMetaD$Barcode)){
  pbmc@meta.data = pbmcMetaD
  rm(pbmcMetaD,psamMetaD)
  cellCts = pbmc$celltype.l1
  ctList = unique(na.omit(cellCts))
  
  ##Temporary save Seurate object.
  saveRDS(pbmc,paste0(opt$out_dir,"/tmpFiltered.Seurat.Rds"))
  for(ct in ctList){
    if(is.null(pbmc)){
      pbmc = readRDS(paste0(opt$out_dir,"/tmpFiltered.Seurat.Rds"))
    }
    ##Grab meta.data
    meta.d.full = as.data.frame(pbmc@meta.data[which(cellCts == ct),])
    meta.d = unique(meta.d.full[,which(colnames(meta.d.full) %in% relCol)])
    ##Grab the countmatrix
    countMatrixFull <- GetAssayData(pbmc[,which(cellCts == ct)], slot = "counts")
    countMatrix = countMatrixFull[which(rowSums(countMatrixFull)!=0),]
    ##Clean
    pbmc=NULL
    rm(countMatrixFull)
    gc();
    
    IDs = meta.d.full$Donor_Pool
    unique_ID_list = unique(IDs)
    
    normCountMatrix = countMatrix ##To store in the new object in.
    
    # Info
    print(paste0('Calculating pseudo-bulk for ',ct,' expression matrix using: ', aggregate_fun))
    print(paste0('   ### n_cells: ', ncol(normCountMatrix)))
    print(paste0('   ### n_genes (sc-data): ', nrow(normCountMatrix)))
    cat('\n')
    
    # PF
    sampleSumInfo = colSums(normCountMatrix)
    meanSampleSum = mean(sampleSumInfo)
    sampleScale = sampleSumInfo/meanSampleSum
    
    # Scale to match meanRowSums
    system.time(normCountMatrix@x <- normCountMatrix@x / rep.int(sampleScale, diff(normCountMatrix@p))) #divide each column by sampleScale
    
    # Log
    normCountMatrix@x  = log(normCountMatrix@x +1)
    
    # PF
    sampleSumInfo = colSums(normCountMatrix)
    meanSampleSum = mean(sampleSumInfo)
    sampleScale = sampleSumInfo/meanSampleSum
  
    # Scale to match meanRowSums
    system.time(normCountMatrix@x <- normCountMatrix@x / rep.int(sampleScale, diff(normCountMatrix@p))) #divide each column by sampleScale
    
    ## Single cell object
    print(paste0("Writing QC-ed and Normalized SC data: ",ct))
    rownames(meta.d.full)=meta.d.full$Barcode
    seuratObj <- CreateSeuratObject(countMatrix, assay = "RNA",meta.data = meta.d.full)
    seuratObj[["data"]] = CreateAssayObject(data=normCountMatrix)
    saveRDS(seuratObj,paste0(opt$out_dir,"/",ct,".Qced.Normalized.SCs.Rds"))
    rm(seuratObj,meta.d.full,countMatrix)
    gc();
    
    print('Aggregating count matrix using lapply + textTinyR::sparse_Means() ...')
    ##normCountMatrix <- as(normCountMatrix, "dgCMatrix")
    
    system.time(aggregate_normCountMatrix <- as.data.frame(pblapply(unique_ID_list, FUN = function(x){sparse_Means(normCountMatrix[,IDs == x, drop = FALSE], rowMeans = TRUE)})))
    cellCount <- pblapply(unique_ID_list, FUN = function(x){ncol(normCountMatrix[,IDs == x, drop = FALSE])})
    colnames(aggregate_normCountMatrix) <- names(cellCount) <- unique_ID_list
    rownames(aggregate_normCountMatrix) <- rownames(normCountMatrix)
    
    cellCount = unlist(cellCount)
    
    print(paste0("Writing pseudo bulk data: ",ct))
    ##Write out psuedo-bulk.
    write.table(aggregate_normCountMatrix,paste0(opt$out_dir,"/",ct,".Exp.txt"),quote=F,sep="\t",col.names=NA)
    
    ##Write covariates global covariates.
    meta.d = cbind(meta.d,cellCount[match(meta.d$Donor_Pool,names(cellCount))]) ##Add cell numbers.
    colnames(meta.d)[ncol(meta.d)]="CellCount" ##Fix nameing.
    write.table(meta.d,paste0(opt$out_dir,"/",ct,".covariates.txt"),quote=F,sep="\t",row.names=F)
    
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
        write.table(aggregate_normCountMatrix,paste0(opt$out_dir,"/",ct,".qtlInput.txt"),quote=F,sep="\t",col.names=NA)
        write.table(covOut,paste0(opt$out_dir,"/",ct,".qtlInput.Pcs.txt"),quote=F,sep="\t",col.names=NA)
      }
    }
  }
}
