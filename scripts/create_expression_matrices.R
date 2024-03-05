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

print(paste0('Combining pools and filter on QC passing ', opt$cell_type, ' cells'))
pools <- read_delim(opt$poolsheet, delim = "\t")
metadata <- as.data.frame(read_delim(opt$metadata, delim="\t"))
rownames(metadata) <- metadata$Barcode

exclude <- NULL
if (!is.null(opt$exclude)) {
  print("Filtering on assignment_pool from include.")
  exclude_smf <- read_delim(opt$exclude, delim = "\t", col_names = FALSE)
  exclude <- paste0(exclude_smf$X1, "_", exclude_smf$X2)
  rm(exclude_smf)
}

filter_seurat_object <- function(row){
  print(paste0("  ", row[["Pool"]]))

  # Load the seurat pool.
  seurat <- readRDS(paste0(opt$input_dir, row[["Pool"]], ".rds"))

  # Add the QC tag.
  seurat <- AddMetaData(seurat, metadata)

  print(paste0("    Barcodes: N = ", length(colnames(seurat))))

  if (sum(seurat@meta.data["DropletType"] != "singlet")) {
    print("    Error, not all 'DropletType' are equal to 'Singlet'. These barcodes should have been removed in 'create_seurat.R''")
    quit()
  }

  mask1 <- seurat@meta.data["DropletType"] == "singlet" # Should already have been filtered on.
  mask2 <- seurat@meta.data["tag"] == "NotOutlier"
  mask3 <- seurat@meta.data[opt$cell_level] == opt$cell_type
  mask4 <- opt$ancestry == "ALL" | seurat@meta.data["Provided_Ancestry"] == opt$ancestry
  mask4[is.na(mask4)] <- FALSE # this info might be missing
  mask5 <- seurat@meta.data["cell_treatment"] == "UT"

  mask6 <- TRUE
  if (!is.null(exclude)) {
    mask6 <- !(paste0(seurat@meta.data$IID, "_", seurat@meta.data$Donor_Pool) %in% exclude)
  }

  mask <- mask1 & mask2 & mask3 & mask4 & mask5 & mask6

  print(paste0("    DropletType - Singlet: N = ", sum(mask1)))
  print(paste0("    tag - NotOutlier: N = ", sum(mask2)))
  print(paste0("    ", opt$cell_level, " - ", opt$cell_type, ": N = ", sum(mask3)))
  print(paste0("    Provided_Ancestry - ", opt$ancestry, ": N = ", sum(mask4)))
  print(paste0("    cell_treatment - UT: N = ", sum(mask5)))
  print(paste0("    IID_SID - include: N = ", sum(mask6)))
  print(paste0("    pass filter: N = ", sum(mask)))
  cat('\n')

  if (sum(mask) > 0 ) {
    return(seurat[, mask])
  } else {
    return()
  }
}
seurat_objects <- apply(pools, 1, filter_seurat_object)
seurat_objects[sapply(seurat_objects, is.null)] <- NULL
nindividuals <- length(seurat_objects)
seurat <- merge(seurat_objects[[1]], y = seurat_objects[2:length(seurat_objects)])

countMatrixFull <- GetAssayData(seurat, slot = "counts")
countMatrix <- countMatrixFull[which(rowSums(countMatrixFull) != 0), ]

##Clean
rm(countMatrixFull)
gc();

ncells <- ncol(countMatrix)
ngenes <- nrow(countMatrix)

# Info
print(paste0('PFlogPF normalising expression matrix'))
print(paste0('   ### n_individuals: ', nindividuals))
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
seuratObj <- CreateSeuratObject(countMatrix, assay = "RNA", meta.data = seurat@meta.data)
seuratObj[["data"]] <- CreateAssayObject(data = normCountMatrix)
saveRDS(seuratObj, paste0(opt$out_dir, opt$cell_type, ".Qced.Normalized.SCs.Rds"))
write_delim(seurat@meta.data, gzfile(paste0(opt$out_dir, opt$cell_type, ".metadata.tsv.gz")), "\t")
