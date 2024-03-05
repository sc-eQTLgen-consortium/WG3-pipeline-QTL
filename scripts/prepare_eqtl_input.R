#!/usr/bin/env Rscript
# Author: M.J. Bonder (Adapted from 'https://github.com/sc-eQTLgen-consortium/WG3-pipeline-QTL/blob/master/scripts/Create_QCed_CT_pseudobulk_Azi_L1.R' by M. Vochteloo)

# options parser
.libPaths("/usr/local/lib/R/site-library")
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list <- list(
  make_option(c("--seurat"), action="store", default=NA, type='character',
              help="the QC filtered PF-log-PF seurat object."),
  make_option(c("--wg1_psam"), action="store", default=NA, type='character',
              help="psam file created in WG1."),
  make_option(c("--cell_annotation"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--aggregate_fun"), action="store", default='mean', type='character',
              help="sum or mean --> if mean, we need to normalize before."),
  make_option(c("--min_n_cells"), action="store", default=5, type='integer',
              help=""),
  make_option(c("--min_n_samples"), action="store", default=15, type='integer',
              help=""),
  make_option(c("--n_pcs"), action="store", default=10, type='integer',
              help=""),
  make_option(c("--out_dir"), action="store", default=NA, type='character',
              help="Output main directory"),
  make_option(c("--out_file"), action="store", default=NA, type='character',
              help="Output main directory"))
opt <- parse_args(OptionParser(option_list=option_list))

print("Options in effect:")
for (name in names(opt)) {
	print(paste0("  --", name, " ", opt[[name]]))
}
print("")

if(opt$aggregate_fun != 'mean'){
  print('Not implemented yet...')
  stop()
}

# check date is provided
if (any(unlist(lapply(c(opt$wg1_psam, opt$cell_annotation, opt$out_dir, opt$out_file), is.na)))) {
  stop("required parameters (wg1_psam & cell_annotation & out_dir & out_file) must be provided.")
}

shhh(library(readr))
shhh(library(Seurat))
shhh(library(matrixStats))
shhh(library(textTinyR))
shhh(library(pbapply))

print("Loading data")
seurat <- readRDS(file = opt$seurat)

wg1_psam <- read_delim(opt$wg1_psam, delim = "\t")
cell_annotation <- read_delim(opt$cell_annotation, delim = "\t")
covariate_cols <- c("Donor_Pool",
                    gsub(x = colnames(wg1_psam), pattern = "\\#", replacement = ""),
                    colnames(cell_annotation)[!colnames(cell_annotation) %in% c("Barcode")])
rm(wg1_psam, cell_annotation)

print("Writing smf file")
smf <- seurat@meta.data[, c("IID", "Donor_Pool")]
smf <- smf[!duplicated(smf), ]
rownames(smf) <- NULL
write_delim(smf, paste0(opt$out_dir, opt$out_file, ".smf.txt"), delim = "\t", col_names = FALSE)
rm(smf)

# Create the Pseudobulk matrix.
normCountMatrix <- GetAssayData(object = seurat, assay = "data", slot = "data")
##normCountMatrix <- as(normCountMatrix, "dgCMatrix")

IDs <- seurat@meta.data$Donor_Pool
unique_ID_list <- unique(IDs)

nindividuals <- length(unique_ID_list)
ncells <- ncol(normCountMatrix)
ngenes <- nrow(normCountMatrix)

# Info
print(paste0('Calculating pseudo-bulk for expression matrix using: ', opt$aggregate_fun))
print(paste0('   ### n_individuals: ', nindividuals))
print(paste0('   ### n_cells: ', ncells))
print(paste0('   ### n_genes: ', ngenes))
cat('\n')

if (ncells <= 1) {
  stop('Cell type has 1 or 0 cells, skip')
}

print('Aggregating count matrix using lapply + textTinyR::sparse_Means() ...')
system.time(aggregate_normCountMatrix <- as.data.frame(pblapply(unique_ID_list, FUN = function(x) { sparse_Means(normCountMatrix[, IDs == x, drop = FALSE], rowMeans = TRUE) })))
cellCount <- pblapply(unique_ID_list, FUN = function(x) { ncol(normCountMatrix[, IDs == x, drop = FALSE]) })
colnames(aggregate_normCountMatrix) <- names(cellCount) <- unique_ID_list
rownames(aggregate_normCountMatrix) <- rownames(normCountMatrix)

cellCount <- unlist(cellCount)

##Write out psuedo-bulk.
print("Writing pseudo bulk data")
write.table(aggregate_normCountMatrix, paste0(opt$out_dir, opt$out_file, ".Exp.txt"), quote = F, sep = "\t", col.names = NA)

##Write covariates global covariates.
print("Writing covariates")
meta.d <- seurat@meta.data[, which(colnames(seurat@meta.data) %in% covariate_cols)]
rownames(meta.d) <- NULL
meta.d <- meta.d[!duplicated(meta.d), ]
meta.d <- cbind(meta.d, cellCount[match(meta.d$Donor_Pool, names(cellCount))]) ##Add cell numbers.
colnames(meta.d)[ncol(meta.d)] <- "CellCount" ##Fix nameing.
write.table(meta.d, paste0(opt$out_dir, opt$out_file, ".covariates.txt"), quote = F, sep = "\t", row.names = F)
rm(meta.d)

##Write QTL input.
print("Writing QTL input")

##For QTL first filter, take only obersvations from min_n_cells cells, and rows that are varying.
aggregate_normCountMatrix <- aggregate_normCountMatrix[, which(colnames(aggregate_normCountMatrix) %in% names(which(cellCount >= opt$min_n_cells)))]

##Only write-out files when there are more then min_n_sample observations.
if (is.null(dim(aggregate_normCountMatrix)) || ncol(aggregate_normCountMatrix) <= opt$min_n_sample) {
  stop("  No samples with minimum number of cells")
}

##Filter based on gene variance
rowVarInfo <- rowVars(as.matrix(aggregate_normCountMatrix))
aggregate_normCountMatrix <- aggregate_normCountMatrix[which(rowVarInfo != 0),]

## Do inverse normal transform per gene.
for (rN in 1:nrow(aggregate_normCountMatrix)) {
  aggregate_normCountMatrix[rN,] <- qnorm((rank(aggregate_normCountMatrix[rN,], na.last = "keep") - 0.5) / sum(!is.na(aggregate_normCountMatrix[rN,])))
}

##Write out input matrix for QTL.
write.table(aggregate_normCountMatrix, paste0(opt$out_dir, opt$out_file, ".qtlInput.txt"), quote = F, sep = "\t", col.names = NA)

##Do PCA, and select first N components.
print("Writing PCA")
pcOut <- prcomp(t(aggregate_normCountMatrix))
var_explained <- as.data.frame(pcOut$sdev^2 / sum(pcOut$sdev^2))
colnames(var_explained) <- c("VarianceExplained")
write.table(var_explained, paste0(opt$out_dir, opt$out_file, ".qtlInput.Pcs.var.txt"), quote = F, sep = "\t", col.names = NA)

n_pcs <- ncol(pcOut$x)
if (n_pcs < opt$n_pcs) {
  print(paste0("  Warning, only ", n_pcs, " PCs found."))
  opt$n_pcs <- n_pcs
}
covOut <- pcOut$x[, 1:opt$n_pcs]

##Write out PCs for QTL.
write.table(covOut, paste0(opt$out_dir, opt$out_file, ".qtlInput.Pcs.txt"), quote = F, sep = "\t", col.names = NA)
