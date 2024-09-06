#!/usr/bin/env Rscript
# Author: M.J. Bonder (Adapted from 'https://github.com/sc-eQTLgen-consortium/WG3-pipeline-QTL/blob/master/scripts/Create_QCed_CT_pseudobulk_Azi_L1.R' by M. Vochteloo)

# options parser
.libPaths("/usr/local/lib/R/site-library")
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list <- list(
  make_option(c("--seurat"), action="store", default=NA, type='character',
              help="the QC filtered PF-log-PF seurat object."),
  make_option(c("--aggregate_fun"), action="store", default='mean', type='character',
              help="sum or mean --> if mean, we need to normalize before."),
  make_option(c("--min_n_cells"), action="store", default=5, type='integer',
              help=""),
  make_option(c("--min_n_samples"), action="store", default=15, type='integer',
              help=""),
  make_option(c("--individual_aggregate"), action="store", default="Assignment", type='character',
              help="Metadata on which to aggregate individuals."),
  make_option(c("--sample_aggregate"), action="store", default="Assignment_Run_Lane", type='character',
              help="Metadata on which to aggregate samples."),
  make_option(c("--n_pcs"), action="store", default=5, type='integer',
              help=""),
  make_option(c("--data_out"), action="store", default=NA, type='character',
              help="Output main directory"),
  make_option(c("--plot_out"), action="store", default=NA, type='character',
              help="Output plot directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

print("Options in effect:")
for (name in names(opt)) {
	print(paste0("  --", name, " ", opt[[name]]))
}
print("")

if(opt$aggregate_fun != 'mean' & opt$aggregate_fun != 'sum'){
  stop("Aggregation function not implemented...")
}

# check date is provided
if (any(unlist(lapply(c(opt$seurat, opt$data_out, opt$plot_out), is.na)))) {
  stop("required parameters (seurat, data_out and plot_out) must be provided.")
}

shhh(library(plyr))
shhh(library(readr))
shhh(library(Seurat))
shhh(library(matrixStats))
shhh(library(textTinyR))
shhh(library(pbapply))
shhh(library(reshape2))
shhh(library(ggplot2))

print("Loading data")
seurat <- readRDS(file = opt$seurat)

# Extract the data we want to use.
metadataMatrix <- seurat@meta.data
normCountMatrix <- GetAssayData(object = seurat, assay = "data", slot = "data")
##normCountMatrix <- as(normCountMatrix, "dgCMatrix")
rm(seurat)
gc()

print("Writing smf file")
smf <- metadataMatrix[, c(opt$individual_aggregate, opt$sample_aggregate)]
smf <- smf[!duplicated(smf), ]
rownames(smf) <- NULL
write_delim(smf, paste0(opt$data_out, "smf.txt"), delim = "\t", col_names = FALSE)

# Create the Pseudobulk matrix.
IDs <- metadataMatrix[, opt$sample_aggregate]
unique_ID_list <- unique(IDs)

nindividual_samples <- nrow(smf)
nindividuals <- length(unique(smf[, opt$individual_aggregate]))
nsamples <- length(unique(smf[, opt$sample_aggregate]))
ncells <- ncol(normCountMatrix)
ngenes <- nrow(normCountMatrix)

# Info
print(paste0('Calculating pseudo-bulk for expression matrix using: ', opt$aggregate_fun))
print(paste0('   ### n_individuals - n_samples: ', nindividual_samples))
print(paste0('   ### n_individuals: ', nindividuals))
print(paste0('   ### n_samples: ', nsamples))
print(paste0('   ### n_cells: ', ncells))
print(paste0('   ### n_genes: ', ngenes))
cat('\n')

if (ncells <= 1) {
  stop('Cell type has 1 or 0 cells, skip')
}

if (opt$aggregate_fun == 'mean') {
  print('Aggregating count matrix using lapply + textTinyR::sparse_Means() ...')
  system.time(aggregate_normCountMatrix <- as.data.frame(pblapply(unique_ID_list, FUN = function(x) { sparse_Means(normCountMatrix[, IDs == x, drop = FALSE], rowMeans = TRUE) })))
} else if (opt$aggregate_fun == 'sum') {
  print('Aggregating count matrix using lapply + textTinyR::sparse_Sums() ...')
  system.time(aggregate_normCountMatrix <- as.data.frame(pblapply(unique_ID_list, FUN = function(x) { sparse_Sums(normCountMatrix[, IDs == x, drop = FALSE], rowSums = TRUE) })))
} else {
  print('Aggregation function not implemented...')
  stop()
}

cellCount <- pblapply(unique_ID_list, FUN = function(x) { ncol(normCountMatrix[, IDs == x, drop = FALSE]) })
colnames(aggregate_normCountMatrix) <- names(cellCount) <- unique_ID_list
rownames(aggregate_normCountMatrix) <- rownames(normCountMatrix)
rm(normCountMatrix, smf, IDs, unique_ID_list)
gc()

cellCount <- unlist(cellCount)
cell.count.d <- data.frame(cellCount)
colnames(cell.count.d) <- c("CellCount")
cell.count.d[opt$sample_aggregate] <- rownames(cell.count.d)
rownames(cell.count.d) <- NULL

##Write out psuedo-bulk.
print("Writing pseudo bulk data")
write.table(aggregate_normCountMatrix, paste0(opt$data_out, "Exp.txt"), quote = F, sep = "\t", col.names = NA)

##Write covariates global covariates.
print("Writing covariates")
covariate_cols <- c(c(opt$individual_aggregate, opt$sample_aggregate),
                    c("FID", "IID", "PAT", "MAT", "SEX", "Provided_Ancestry", "genotyping_platform", "array_available",
                     "wgs_available", "wes_available", "age", "age_range", "Study", "smoking_status",
                     "hormonal_contraception_use_currently", "menopause", "pregnancy_status"),
                    c("sequencing_platform", "sequencing_run", "sequencing_lane", "scrna_platform", "plate_based", "umi_based", "biomaterial", "sorting", "cell_treatment", "sample_condition"))
meta.d <- metadataMatrix[, which(colnames(metadataMatrix) %in% covariate_cols)]
rownames(meta.d) <- NULL
meta.d <- meta.d[!duplicated(meta.d), ]
meta.d <- join(x = meta.d, y = cell.count.d, by = opt$sample_aggregate, type = "left", match = "first")
write.table(meta.d, paste0(opt$data_out, "covariates.txt"), quote = F, sep = "\t", row.names = F)
rm(meta.d)

##Write QTL input.
print("Writing QTL input")

##For QTL first filter, take only obsevations from min_n_cells cells, and rows that are varying.
include_samples = which(colnames(aggregate_normCountMatrix) %in% names(which(cellCount >= opt$min_n_cells)))
print(paste0('   ### excluded ', nsamples - length(include_samples) ,' samples due to number of cells <', opt$min_n_cells))
aggregate_normCountMatrix <- aggregate_normCountMatrix[, include_samples]

##Only write-out files when there are more then min_n_sample observations.
if (is.null(dim(aggregate_normCountMatrix)) || ncol(aggregate_normCountMatrix) <= opt$min_n_sample) {
  stop("  Not sufficient number of samples with minimum number of cells")
}

##Filter based on gene variance
rowVarInfo <- rowVars(as.matrix(aggregate_normCountMatrix))
aggregate_normCountMatrix <- aggregate_normCountMatrix[which(rowVarInfo != 0),]

## Do inverse normal transform per gene.
# TODO: implement improvement of Roy
for (rN in 1:nrow(aggregate_normCountMatrix)) {
  aggregate_normCountMatrix[rN,] <- qnorm((rank(aggregate_normCountMatrix[rN,], na.last = "keep") - 0.5) / sum(!is.na(aggregate_normCountMatrix[rN,])))
}

##Write out input matrix for QTL.
write.table(aggregate_normCountMatrix, paste0(opt$data_out, "qtlInput.txt"), quote = F, sep = "\t", col.names = NA)

##Do PCA
print("Writing PCA")
pcOut <- prcomp(t(aggregate_normCountMatrix))

##Write out PCs for QTL.
write.table(pcOut$x, paste0(opt$data_out, "qtlInput.Pcs.txt"), quote = F, sep = "\t", col.names = NA)

var_explained <- as.data.frame(pcOut$sdev^2 / sum(pcOut$sdev^2))
colnames(var_explained) <- c("VarianceExplained")
write.table(var_explained, paste0(opt$data_out, "qtlInput.Pcs.var.txt"), quote = F, sep = "\t", col.names = NA)

print("Plotting PCA")
pca <- as.data.frame(pcOut$x)[,1:(opt$n_pcs + 1)]
colnames(pca) <- paste0(colnames(pca), " ", round(var_explained[1:(opt$n_pcs + 1), "VarianceExplained"] * 100, 2), "%")

plot_df1 <- melt(as.matrix(pca))
plot_df1$var_index <- as.numeric(gsub("PC", "", gsub(" .*", "", plot_df1$Var2)))
plot_df2 <- data.frame(plot_df1)

colnames(plot_df1) <- c("SID", "variable1", "value1", "var_index1")
colnames(plot_df2) <- c("SID", "variable2", "value2", "var_index2")

plot_df <- merge(x = plot_df1, y = plot_df2, by = "SID", all.x = TRUE, all.y = TRUE)
plot_df <- plot_df[plot_df$var_index1 < plot_df$var_index2, ]
rm(plot_df1, plot_df2)

pca_plot <- ggplot(plot_df, aes(value1, value2)) +
  geom_point() +
  facet_grid(rows=vars(variable2), cols=vars(variable1), scales="free") +
  xlab("") +
  ylab("") +
  ggtitle("PCA Plot")

ggsave(pca_plot, filename = paste0(opt$plot_out, "qtlInput.Pcs.png"), width = 29.7, height = 21, units = c("cm"))

var_explained$x <- as.numeric(rownames(var_explained))
scree_plot <- ggplot(var_explained, aes(x=x, y=VarianceExplained, group=1)) +
  geom_line() +
  geom_bar(stat="identity", color="black") +
  xlab("Principal Component") +
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  xlim(0, min(nrow(var_explained), 100))
ggsave(scree_plot, filename = paste0(opt$plot_out, "qtlInput.Pcs.var.png"), width = 29.7, height = 21, units = c("cm"))
