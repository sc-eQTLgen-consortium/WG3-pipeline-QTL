#!/usr/bin/env Rscript
# Author: M.J. Bonder (Adapted from 'https://github.com/sc-eQTLgen-consortium/WG3-pipeline-QTL/blob/master/scripts/Create_QCed_CT_pseudobulk_Azi_L1.R' by M. Vochteloo)

# options parser
.libPaths("/usr/local/lib/R/site-library")
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list <- list(
  make_option(c("--poolsheet"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--input_dir"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--metadata"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--ancestry"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--cell_level"), action="store", default="celltype.l1", type='character',
              help="Azimuth l1 or l2."),
  make_option(c("--cell_type"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--exclude"), action="store", default=NULL, type='character',
              help=""),
  make_option(c("--individual_aggregate"), action="store", default="Assignment", type='character',
              help="Metadata on which to aggregate individuals."),
  make_option(c("--sample_aggregate"), action="store", default="Assignment_Run_Lane", type='character',
              help="Metadata on which to aggregate samples."),
  make_option(c("--out_dir"), action="store", default=NA, type='character',
              help="Output main directory"))
opt <- parse_args(OptionParser(option_list=option_list))

print("Options in effect:")
for (name in names(opt)) {
	print(paste0("  --", name, " ", opt[[name]]))
}
print("")

# check date is provided
if (any(unlist(lapply(c(opt$poolsheet, opt$input_dir, opt$metadata, opt$cell_type, opt$out_dir), is.na)))) {
  stop("required parameters (poolsheet, input_dir, metadata, cell_type & out_dir) must be provided.")
}

if (is.na(opt$ancestry)) {
    opt$ancestry <- "ALL"
}

shhh(library(readr))
shhh(library(Seurat))
shhh(library(matrixStats))
shhh(library(textTinyR))
shhh(library(pbapply))
shhh(library(dplyr))

print(paste0('Combining pools and filter on QC passing ', opt$cell_type, ' cells'))
pools <- read_delim(opt$poolsheet, delim = "\t")
metadata <- as.data.frame(read_delim(opt$metadata, delim="\t"))
rownames(metadata) <- metadata$Barcode

exclude <- NULL
if (!is.null(opt$exclude)) {
  print("Loading exclude file.")
  exclude_smf <- read_delim(opt$exclude, delim = "\t", col_names = FALSE)
  exclude <- paste0(exclude_smf$X1, "_", exclude_smf$X2)
  rm(exclude_smf)
}

filter_seurat_object <- function(row){
  print(paste0("  ", row[["Pool"]]))

  # Load the filter stats.NotOutlier
  filter.stats <- read_delim(paste0(opt$input_dir, row[["Pool"]], ".filter.stats.tsv"), delim = "\t")
  filter.stats.list <- unlist(filter.stats[, row[["Pool"]]], use.names=F)
  names(filter.stats.list) <- unlist(filter.stats[, "filter"], use.names=F)
  rm(filter.stats)

  # Load the seurat pool.
  counts <- tryCatch({
    seurat <- readRDS(paste0(opt$input_dir, row[["Pool"]], ".rds"))
  },error = function(e){
    print(e)
    return(list("pool" = row[["Pool"]], "seurat" = NULL, "filter" = filter.stats.list))
  })

  # Add the QC tag.
  seurat <- AddMetaData(seurat, metadata)
  # print(paste0("    Input object size: ", format(object.size(seurat), units = "Mb")))

  n_barcodes <- nrow(seurat@meta.data)
  if (n_barcodes != filter.stats.list[["DropletType - Singlet"]]) {
    print("Error, filter.stats.tsv and seurat file do not match.")
    quit()
  }
  if (sum(seurat@meta.data["DropletType"] != "singlet")) {
    print("    Error, not all 'DropletType' are equal to 'Singlet'. These barcodes should have been removed in 'create_seurat.R''")
    quit()
  }

  # THIS HAS TO BE IDENTICAL TO qtl_sample_qc.py filter except for the cell type filter!!!!!!!!!!!
  mask1 <- seurat@meta.data["DropletType"] == "singlet" # Should already have been filtered on.
  mask2 <- seurat@meta.data["tag"] == "NotOutlier"
  mask3 <- seurat@meta.data[opt$cell_level] == opt$cell_type
  mask4 <- opt$ancestry == "ALL" | seurat@meta.data["Provided_Ancestry"] == opt$ancestry
  mask4[is.na(mask4)] <- FALSE # this info might be missing
  mask5 <- seurat@meta.data["cell_treatment"] == "UT"

  mask6 <- rep(TRUE, n_barcodes)
  if (!is.null(exclude)) {
    mask6 <- !(paste0(seurat@meta.data[,opt$individual_aggregate], "_", seurat@meta.data[,opt$sample_aggregate]) %in% exclude)
  }

  mask <- mask1 & mask2 & mask3 & mask4 & mask5 & mask6
  mask[is.na(mask)] <- FALSE

  # Save and print the filter stats.
  filter.stats.list["DropletType - Singlet"] = sum(mask1, na.rm=T)
  filter.stats.list["tag - NotOutlier"] = sum(mask2, na.rm=T)
  filter.stats.list[paste0(opt$cell_level, " - ", opt$cell_type)] = sum(mask3, na.rm=T)
  filter.stats.list[paste0("Provided_Ancestry - ", opt$ancestry)] = sum(mask4, na.rm=T)
  filter.stats.list["cell_treatment - UT"] = sum(mask5, na.rm=T)
  filter.stats.list[paste0(opt$individual_aggregate, " - ", opt$sample_aggregate, " - include")] = sum(mask6, na.rm=T)
  filter.stats.list["Cells"] = sum(mask, na.rm=T)

  print(paste0("    Filtering barcodes:"))
  for (name in names(filter.stats.list)) {
    print(paste0("      ", name, ": N = ", filter.stats.list[[name]]))
  }

  if (sum(mask) > 0) {
    seurat <- seurat[, mask]
  } else {
    seurat <- NULL
  }

  # print(paste0("    Filtered object size: ", format(object.size(seurat), units = "Mb")))
  cat('\n')
  return(list("pool" = row[["Pool"]], "seurat" = seurat, "filter" = filter.stats.list))
}
combined_data <- apply(pools, 1, filter_seurat_object)

# Combine seurat objects into one big object.
seurat_objects <- unlist(lapply(combined_data, function(pool) pool$seurat))
seurat <- merge(seurat_objects[[1]], y = seurat_objects[2:length(seurat_objects)])
# print(paste0("    Combined object size: ", format(object.size(seurat), units = "Mb")))
nsamples <- length(seurat_objects)

# Combine filter stats into one data frame.
filter.stats <- bind_cols(lapply(combined_data, function(pool) as.data.frame(pool$filter)))
colnames(filter.stats) <- lapply(combined_data, function(pool) pool$pool)
print(filter.stats)
write.table(filter.stats, gzfile(paste0(opt$out_dir, opt$cell_type, ".filter.stats.tsv.gz")), sep="\t")

##Clean
rm(pools, metadata, exclude, combined_data, seurat_objects, filter.stats)
gc()

metadataMatrixFull <- seurat@meta.data
countMatrixFull <- GetAssayData(seurat, slot = "counts")
countMatrix <- countMatrixFull[which(rowSums(countMatrixFull) != 0), ] # TODO: if Dan removes this then the results differ

##Clean
rm(seurat, countMatrixFull)
gc()

ncells <- ncol(countMatrix)
ngenes <- nrow(countMatrix)

# Info
print(paste0('PFlogPF normalising expression matrix'))
print(paste0('   ### n_samples: ', nsamples))
print(paste0('   ### n_cells: ', ncells))
print(paste0('   ### n_genes: ', ngenes))
cat('\n')

if (ncells <= 1) {
  stop('Cell type has 1 or 0 cells, skip')
}

normCountMatrix <- countMatrix ##To store in the new object in.

# PF, scale to match meanRowSums
sampleSumInfo <- colSums(normCountMatrix)
meanSampleSum <- mean(sampleSumInfo)
sampleScale <- sampleSumInfo / meanSampleSum
system.time(normCountMatrix@x <- normCountMatrix@x / rep.int(sampleScale, diff(normCountMatrix@p))) #divide each column by sampleScale

# Log
normCountMatrix@x <- log(normCountMatrix@x + 1)

# PF, scale to match meanRowSums
sampleSumInfo <- colSums(normCountMatrix)
meanSampleSum <- mean(sampleSumInfo)
sampleScale <- sampleSumInfo / meanSampleSum
system.time(normCountMatrix@x <- normCountMatrix@x / rep.int(sampleScale, diff(normCountMatrix@p))) #divide each column by sampleScale

## Single cell object
print("Writing QC-ed and Normalized SC data")
seuratObj <- CreateSeuratObject(countMatrix, assay = "RNA", meta.data = metadataMatrixFull)
seuratObj[["data"]] <- CreateAssayObject(data = normCountMatrix)
saveRDS(seuratObj, paste0(opt$out_dir, opt$cell_type, ".Qced.Normalized.SCs.Rds"))
write_delim(seuratObj@meta.data, gzfile(paste0(opt$out_dir, opt$cell_type, ".metadata.tsv.gz")), "\t")
