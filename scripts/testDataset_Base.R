setwd("/")

##Ye EUR (on E:) 
datasetFolder = "./input/L1/"

threshold = 2.5
toCheckC = list.files(datasetFolder,pattern = ".covariates.",full.names = T)

##Need to fix the cov file (make sure there are no rownames)

##Lazily first calculate total number of cells per donor.
donorCells = NULL
for(cts in toCheckC){
  ##print(cts)
  covData <- read.delim(cts,as.is=T)
  
  covData = unique(covData)
  row.names(covData) = covData$Donor_Pool
  
  
  if(is.null(donorCells)){
    donorCells = covData[,c("Donor_Pool","CellCount")]
  } else {
    if(length(which(!(covData$Donor_Pool %in% donorCells$Donor_Pool)))>0){
      toAdd = cbind(covData$Donor_Pool[which(!(covData$Donor_Pool %in% donorCells$Donor_Pool))],0)
      colnames(toAdd) = colnames(donorCells)
      toAdd[,2] = as.numeric(toAdd[,2])
      donorCells = rbind(donorCells,toAdd)
    }
    donorCells[,2] = as.numeric(donorCells[,2])
    additionVector= covData[match(donorCells$Donor_Pool,covData$Donor_Pool),"CellCount"]
    additionVector[which(is.na(additionVector))] = 0
    donorCells[,2] = donorCells[,2]+additionVector
  }
}

hist(donorCells[,2])
cCount = median(donorCells[,2])
cCount = cCount - (median(c(abs(donorCells[,2]-cCount)))*threshold)
abline(v=cCount)

donorCells[which(donorCells[,2]<cCount),]

## Test fraction of cell-types.
threshold = 5
for(cts in toCheckC){
  print(cts)
  covData <- read.delim(cts)

  covData = unique(covData)
  row.names(covData) = covData$Donor_Pool
  
  
  # cCount = median(covData$CellCount)
  # cCount = cCount - (median(c(abs(covData$CellCount-cCount)))*threshold)
  # 
  # hist(covData$CellCount)
  # abline(v=cCount)
  # which(covData$CellCount<cCount)
  
  ##Fraction
  cellFraction = covData$CellCount/donorCells[match(covData$Donor_Pool,donorCells[,1]),2]
  fractionThreshold = median(cellFraction)
  madFractThresh = median(abs(fractionThreshold-cellFraction))
  print(paste(fractionThreshold,madFractThresh))
  hist(cellFraction,main = cts)
  abline(v=(fractionThreshold-(madFractThresh*threshold)))
  abline(v=(fractionThreshold+(madFractThresh*threshold)))
  toDrop = c(which(cellFraction<(fractionThreshold-(madFractThresh*threshold))),which(cellFraction>(fractionThreshold+(madFractThresh*threshold))))
  print(rownames(covData)[toDrop])
  print(cellFraction[toDrop])
  
}


## Test within cell type, PCs
toCheck = list.files(datasetFolder,pattern = ".Pcs.",full.names = T)
overall=NA
for(cts in toCheck){
  print(cts)
  pcData <- read.delim(cts,row.names = 1)
  
  #plot(pcData[,c(1:2)])
  
  pc1 = median(pcData[,1])
  pc2 = median(pcData[,2])
  madPc1 = median(abs(pcData[,1]-pc1))
  madPc2 = median(abs(pcData[,2]-pc2))
  
  toDrop = unique(c(which(pcData[,1]<(pc1-(madPc1*threshold))),which(pcData[,1]>(pc1+(madPc1*threshold))),
             which(pcData[,2]<(pc2-(madPc2*threshold))),which(pcData[,2]>(pc2+(madPc2*threshold)))))
  print(rownames(pcData)[toDrop])
  overall = c(overall,rownames(pcData)[toDrop])
  colors = rep("black",nrow(pcData))
  colors[toDrop] = "red"
  plot(pcData[,c(1,2)], main=cts,col=colors)
  View(pcData[,c(1:2)])
  print(toDrop)
}
table(overall)/length(toCheck)


##Genotypes.
sampleCov <- read.delim("./input/sample.kinship",as.is=T,row.names = 1,check.names=F)
out = prcomp(sampleCov)

## Look for samples with 
for(colId in c(1:ncol(sampleCov))){
  sameSample = which((sampleCov[colId,]>0.8))
  if(length(sameSample)>1){
    print(colnames(sampleCov)[sameSample])
  }
}

####
plot(out$x[,c(1,2)])

##Filter to donors seen in the exprssion data (already did european filter.)
samplesToKeep = NULL
toCheck = list.files(datasetFolder,pattern = ".Exp.txt",full.names = T)
for(cts in toCheck){
  expInfo = read.delim(cts,as.is=T,row.names=1,check.names=F)
  samplesToKeep = c(samplesToKeep,colnames(expInfo))
}
samplesToKeep = unique(samplesToKeep)

samplesToKeep = unlist(lapply(strsplit(samplesToKeep,split = ";;"),"[[",1))

sampleCov = sampleCov[(which(rownames(sampleCov) %in% samplesToKeep)),(which(colnames(sampleCov) %in% samplesToKeep))]
out = prcomp(sampleCov)
plot(out$x[,c(1,2)])

##Back to original for plotting outliers
sampleCov <- read.delim("./input/sample.kinship",as.is=T,row.names = 1,check.names=F)
out = prcomp(sampleCov)

plot(out$x[,c(1,2)])
colVec = rep("black",ncol(sampleCov))
colVec[which(colnames(sampleCov) %in% samplesToKeep)] = "red"
plot(out$x[,c(1,2)],col=colVec)

##Add PC1 as extra covariate.

## Test within cell type, PCs
toCheck = list.files(datasetFolder,pattern = ".Pcs.",full.names = T)
#Drop current originals. (the ones I already changed)
toCheck = toCheck[grep(".org",toCheck,invert = T,ignore.case = T)]
overall=NA
for(cts in toCheck){
  print(cts)
  pcData <- read.delim(cts,row.names = 1)
  if(!file.exists(gsub(pattern = ".Pcs.",replacement = ".Pcs.Org.",cts))){
    write.table(pcData,gsub(pattern = ".Pcs.",replacement = ".Pcs.Org.",cts),quote=F,sep="\t",col.names = NA)
  }
  sampleNames = unlist(lapply(strsplit(rownames(pcData),split = ";;"),"[[",1))
  #all(rownames(out$x)[match(sampleNames,rownames(out$x))] == sampleNames)
  
  pcData = cbind(pcData,out$x[match(sampleNames,rownames(out$x)),1])
  colnames(pcData)[11] = "genoPC1"
  write.table(pcData,cts,quote=F,sep="\t",col.names = NA)
}

##Update chunking file,

geneInfo <-  read.delim("./input/LimixAnnotationFile.txt", sep="\t")

##S2. make chunk file [Optional]
testCombinations = NULL
#
nGenes = 50
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
write.table(lines,"./input/ChunkingFile.txt",quote = F,sep="\t",row.names = F,col.names=F)

