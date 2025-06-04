#####
# DEFINING INPUT AND OUTPUT FILES
#####

cat("Running script to generate input files for cell-type composition GWAS\n")

args <- commandArgs(trailingOnly = TRUE)
print("Input arguments: RDS input file, PSAM input file, Output path")
print(args)

wg2_input_rds_file <- args[1] 
psam_file <- args[2] 

output_path <- args[3]

output_ct_prop <- paste0(output_path,"ct_prop_per_donor.txt")
output_dummy_var <- paste0(output_path,"dummy_annotation.txt")
output_cells_per_donor <- paste0(output_path,"cells_per_donor.txt")
output_ct_chunk_file <- paste0(output_path,"ct_ChunkingFile.txt")

n_cells_threshold <- 350

#####
# SETTING UP ANNOTATION FILES
#####

cat(paste("Loading RDS file:",wg2_input_rds_file,"\n"))
file <- readRDS(wg2_input_rds_file)

file_meta <- file@meta.data

## Locate and load cell type reference file
azimuth_pairing <- "/tools/wg1-qc_filtering/azimuth_l1_l2.csv"
pairs.df <- read.csv2(azimuth_pairing)
colnames(pairs.df) <- c('l1','l2')

## Make some changes to azimuth and input files to better define cell-types we are interested in
pairs.df$l2[pairs.df$l2=="NK"]="NK_L2"
file_meta$predicted.celltype.l2[file_meta$predicted.celltype.l2=="NK"] <- "NK_L2"
file_meta$scpred_prediction[file_meta$scpred_prediction=="NK"] <- "NK_L2"

## Define cell-types we are not interested to remove at the end
# L2 cell-types that are identical to L1 cell-types
# other_T due to it being a poorly defined cell-type
# myeloid due to lymphocyte implicitly containing this information by itself
dup_ct_to_rm <- c(pairs.df$l2[pairs.df$l1 == pairs.df$l2])
other_ct_to_rm <- c("other_T","myeloid")

## From L2 celltypes calculate L1 celltypes for both Azimuth and scpred
file_meta["predicted.celltype.l1"] <- pairs.df$l1[match(file_meta$predicted.celltype.l2, pairs.df$l2)]

colnames(file_meta)[which(colnames(file_meta)=="scpred_prediction")]="scpred.l2.prediction"

file_meta["scpred.l1.prediction"] <- pairs.df$l1[match(file_meta$scpred.l2.prediction,pairs.df$l2)]

## Remove celltype labels that do not match between scpred and azimuth
file_meta["celltype.l1"] <- file_meta$predicted.celltype.l1
file_meta$celltype.l1[which(file_meta$celltype.l1!=file_meta$scpred.l1.prediction)] = NA

#Define lists of myeloid and lymphocyte l1 cell-types
myeloid <- c("Mono","DC")
lymphocyte <- c("CD4_T","CD8_T","NK","B","other_T")

#Add l0 cell-types to pairs.df
l0 <- pairs.df$l1
l0[l0 %in% myeloid] <- "myeloid"
l0[l0 %in% lymphocyte] <- "lymphocyte"
l0[!(l0 %in% c("myeloid","lymphocyte"))] <- "other"
pairs.df$l0 <- l0

#Add annotation to file_meta
file_meta$celltype.l0 <- file_meta$celltype.l1
file_meta$celltype.l0[file_meta$celltype.l0 %in% pairs.df[pairs.df$l0=="myeloid","l1"]] <- "myeloid"
file_meta$celltype.l0[file_meta$celltype.l0 %in% pairs.df[pairs.df$l0=="lymphocyte","l1"]] <- "lymphocyte"
file_meta$celltype.l0[file_meta$celltype.l0 %in% pairs.df[pairs.df$l0=="other","l1"]] <- "other"


#####
#L0 CELL-TYPES
#####
cat("Calculating proportions of myeloid vs lymphocyte cell-types\n")

#Making filtered dataframes containing information on assigned cell-type and donor
# Dropping 'other' as we are only interested in myeloid and lymphocyte cells
filt.df.ct <- file_meta$celltype.l0[!is.na(file_meta$celltype.l0) & file_meta$celltype.l0 !="Doublet" & file_meta$celltype.l0 !="other"]
filt.df.as <- file_meta$Assignment[!is.na(file_meta$celltype.l0) & file_meta$celltype.l0 !="Doublet" & file_meta$celltype.l0 !="other"]

## Calculating fraction of cell type per donor
cat("Initialize a new dataframe to store the results\n")
new.df.l0 <- data.frame(matrix(nrow =length(unique(filt.df.ct)),ncol =0))
rownames(new.df.l0) <- c(unique(filt.df.ct))

cat("Making dataframe for fractions of total cells belonging to that celltype\n")
for (donor in unique(filt.df.as)){
  new.df.l0[[donor]] <- lapply(rownames(new.df.l0), function(celltype) sum(filt.df.ct[filt.df.as==donor]==celltype)/sum(filt.df.as==donor))
}
cat("remove rows that are only 0\n")
new.mat.l0 <- matrix(unlist(new.df.l0),ncol=length(new.df.l0))
colnames(new.mat.l0)<-colnames(new.df.l0)
row.names(new.mat.l0)<-row.names(new.df.l0)
new.df.l0 <- new.mat.l0[rowSums(new.mat.l0)>0,]


#####
#L1 CELL-TYPES
#####
cat("Calculating proportions of L1 cell-types\n")

#Making filtered dataframes containing information on assigned cell-type and donor
filt.df.ct <- file_meta$celltype.l1[!is.na(file_meta$celltype.l1) & file_meta$celltype.l1 !="Doublet"]
filt.df.as <- file_meta$Assignment[!is.na(file_meta$celltype.l1) & file_meta$celltype.l1 !="Doublet"]

## Calculating fraction of cell type per donor
cat("Initialize a new dataframe to store the results\n")
new.df.l1 <- data.frame(matrix(nrow =length(unique(filt.df.ct)),ncol =0))
rownames(new.df.l1) <- c(unique(filt.df.ct))

cat("Making dataframe for fractions of total cells belonging to that celltype\n")
for (donor in unique(filt.df.as)){
  new.df.l1[[donor]] <- lapply(rownames(new.df.l1), function(celltype) sum(filt.df.ct[filt.df.as==donor]==celltype)/sum(filt.df.as==donor))
}
cat("remove rows that are only 0\n")
new.mat.l1 <- matrix(unlist(new.df.l1),ncol=length(new.df.l1))
colnames(new.mat.l1)<-colnames(new.df.l1)
row.names(new.mat.l1)<-row.names(new.df.l1)
new.df.l1 <- new.mat.l1[rowSums(new.mat.l1)>0,]


#####
#MAKE COVARIATE FILE
#####
cat("Generating covariates file with age and sex as covariates\n")

#Using PSAM file to access age and sex information
psam <- read.csv(psam_file,sep='\t')

covariates <- data.frame(samples=unique(filt.df.as))
n_cells <- lapply(covariates$samples, function(x) sum(filt.df.as==x))
age <- as.matrix(lapply(covariates$samples, function(x) psam$age[psam$IID==x]))
sex <- as.matrix(lapply(covariates$samples, function(x) psam$SEX[psam$IID==x]))

# Check if there are missing values in age and/or sex columns
if (any(is.na(sex),is.na(age))){
  stop("There are missing values in age and/or sex columns")
}

covariates$n_cells <- n_cells
covariates$age <- age
covariates$sex <- sex

n_cells_df <- covariates

#We will remove donors based on a minimum number of cells
donors_to_rm <- n_cells_df[n_cells_df$n_cells < n_cells_threshold,]$samples

cat(paste0("Number of donors that pass ",n_cells_threshold," cells filter: ",length(n_cells_df[n_cells_df$n_cells >= n_cells_threshold,]$samples),"\n"))


#####
#L2 CELL-TYPES
#####
cat("Calculating proportions of L2 cell-types\n")

## Create L2 df filtered on L1 celltype
#Making filtered dataframes containing information on assigned cell-type and donor
filt.df.ct <- file_meta$predicted.celltype.l2[!is.na(file_meta$celltype.l1) & file_meta$celltype.l1 !="Doublet"]
filt.df.as <- file_meta$Assignment[!is.na(file_meta$celltype.l1) & file_meta$celltype.l1 !="Doublet"]

## Calculating fraction of cell type per donor
# Initialize a new dataframe to store the results
new.df.l2 <- data.frame(matrix(nrow =length(unique(filt.df.ct)),ncol =0))
rownames(new.df.l2) <- c(unique(filt.df.ct))

#Making dataframe for fractions of total cells belonging to that celltype
for (donor in unique(filt.df.as)){
  new.df.l2[[donor]] <- lapply(rownames(new.df.l2), function(celltype) sum(filt.df.ct[filt.df.as==donor]==celltype)/sum(filt.df.as==donor))
}
# remove rows that are only 0
new.mat.l2 <- matrix(unlist(new.df.l2),ncol=length(new.df.l2))
colnames(new.mat.l2) <- colnames(new.df.l2)
row.names(new.mat.l2) <- row.names(new.df.l2)
new.df.l2 <- new.mat.l2[rowSums(new.mat.l2)>0,]

# Dropping cell-types that are identical between l1 and l2
new.df.l2 <- new.df.l2[!(rownames(new.df.l2) %in% dup_ct_to_rm),]

#####
#FINAL OUTPUT FILE
#####
### The ouput file which has the fraction of each cell-type per donor
final.df <- rbind(new.df.l0, new.df.l1, new.df.l2)


## Adding extra comparisons
#Proportion between L1 cell-types
ct_l1 <- c("CD4_T","CD8_T","B","Mono","NK","DC")
for (i in seq(length(ct_l1)-1)){
  for (j in seq(i+1, length(ct_l1))){
    test <- rbind(final.df, final.df[rownames(final.df)==ct_l1[i],] / (final.df[rownames(final.df)==ct_l1[j],] + final.df[rownames(final.df)==ct_l1[i],]))
    rownames(test) <- c(rownames(final.df), paste0(ct_l1[i],"_vs_",ct_l1[j]))
    final.df <- test
  }
}


#Fraction of L2 cell-type in relevant L1 cell-type
#selecting L2 cell-types that are present in final.df
ct_l2 <- rownames(final.df[rownames(final.df) %in% pairs.df$l2,])
ct_l2 <- ct_l2[!(ct_l2 %in% dup_ct_to_rm)]
for (celltype in ct_l2){
  test <- rbind(final.df, final.df[rownames(final.df)==celltype,] / final.df[rownames(final.df)==pairs.df[pairs.df$l2==celltype,]$l1,])
  rownames(test) <- c(rownames(final.df), paste0(celltype,"_in_",pairs.df[pairs.df$l2==celltype,]$l1))
  final.df <- test
}

#Remove Infinite values in case there were divisions by 0 above
final.df[final.df==Inf] <- NA

# Dropping some cell-types as outlined above and in readme
final.df <- final.df[!(rownames(final.df) %in% other_ct_to_rm),]

# Dropping donors with less than required number of cells
final.df <- final.df[,!(colnames(final.df) %in% donors_to_rm)]

# Replace NA with 0
final.df[is.na(final.df)] <- 0

#####
#CREATE ANNOTATION FILE AND CHUNKING FILE
#####
len.col <- length(rownames(final.df))
start.vals <- seq(10, len.col*10, by=10)
end.vals <- start.vals+9
# this is a dummy annotation file; it does not contain real information but is required for the pipeline to run
annot.df <- data.frame(feature_id=rownames(final.df), chromosome=rep(1,len.col), start=start.vals+3, end=start.vals+4)

chunk.df <- data.frame(chromosome=rep(1,len.col), window=paste(start.vals,end.vals,sep='-'))

cat("Writing output\n")
write.table(final.df, file=output_ct_prop, sep="\t", quote=F)
write.table(as.matrix(annot.df), file=output_dummy_var, sep="\t", row.names=FALSE, quote=F)
write.table(as.matrix(n_cells_df), file=output_cells_per_donor, sep="\t", row.names=FALSE, quote=F)
write.table(chunk.df, file=output_ct_chunk_file, sep=":", row.names=FALSE, col.names=FALSE, quote=F)

