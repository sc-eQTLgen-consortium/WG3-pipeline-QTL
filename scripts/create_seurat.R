#!/usr/bin/env Rscript
# Author: M. Vochteloo

# options parser
.libPaths("/usr/local/lib/R/site-library")
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list <- list(
  make_option(c("--poolsheet"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--pool"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--wg1_metadata"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--wg2_metadata"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--wg2_pairing"), action="store", default=NULL, type='character',
              help=""),
  make_option(c("--cell_annotation"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--wg1_psam"), action="store", default=NA, type='character',
              help="psam file created in WG1."),
  make_option(c("--rb_genes"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--mt_genes"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--out_dir"), action="store", default=NA, type='character',
              help=""))
opt <- parse_args(OptionParser(option_list=option_list))

print("Options in effect:")
for (name in names(opt)) {
	print(paste0("  --", name, " ", opt[[name]]))
}
print("")

# check date is provided
if (any(unlist(lapply(c(opt$poolsheet, opt$pool, opt$wg1_metadata, opt$wg2_metadata, opt$cell_annotation, opt$wg1_psam, opt$rb_genes, opt$mt_genes, opt$out_dir), is.na)))) {
  stop("all parameters except --wg2_pairing must be provided.")
}

shhh(library(tidyverse))
shhh(library(plyr))
shhh(library(Seurat))
shhh(library(scCustomize))

# Read the poolsheet file.
pools <- read_delim(opt$poolsheet, delim = "\t")
pools <- pools %>% filter(Pool == opt$pool)
counts_file <- pools %>% getElement("Counts")
barcodes_file <- pools %>% getElement("Barcodes")

# Read the counts matrix.
print(paste0("Creating Seurat object for pool ", opt$pool))
counts <- tryCatch({
  Read10X_h5(counts_file)
}, error = function(e) {
  Read_CellBender_h5_Mat(counts_file)
})
print(paste0("Loaded matrix with shape (", counts@Dim[1], ", ", counts@Dim[2], ")"))
rm(counts_file)

print("Constructing metadata table")
meta.data <- read_delim(barcodes_file, delim = "\t", col_names = c("Barcode"))
print(paste0("  Loaded metadata with shape (", nrow(meta.data), ", ", ncol(meta.data), ")"))
nbarcodes = nrow(meta.data)
if (!identical(colnames(counts), meta.data$Barcode)) {
  print("Error, barcodes do not match count matrix.")
  quit()
}

# Add the pool and make the barcodes unique for merging later.
meta.data$Pool <- opt$pool
meta.data$Barcode <- gsub("-1", "", meta.data$Barcode)
meta.data$Barcode <- paste0(meta.data$Barcode, "_", opt$pool)
meta.data$Order <- 1:dim(meta.data)[1] # This enables merging file that are not the same order.
colnames(counts) <- meta.data$Barcode

# Function that merges two meta.data files and checks if the number of rows are still identical.
merge_metadata <- function(meta.data, filepath, delim="\t", col_names=TRUE, remap_colnames=NULL, copy_colnames=NULL, by=NULL, type="left") {
  meta.data2 <- read_delim(filepath, delim = delim, col_names = col_names, col_types = cols(.default = "c")) # make sure to read as character to prevent typing issues in pools with just numbers
  print(paste0("  Loaded ", basename(filepath), " with shape (", dim(meta.data2)[1], ", ", dim(meta.data2)[2], ")"))

  colnames(meta.data2) <- gsub(x = colnames(meta.data2), pattern = "\\#", replacement = "")
  if (!is.null(remap_colnames)) {
    meta.data2 <- rename(meta.data2, remap_colnames)
  }
  if (!is.null(copy_colnames)) {
    for (name in names(copy_colnames)) {
      meta.data2[[name]] <- meta.data2[[copy_colnames[[name]]]]
    }
    meta.data2 <- rename(meta.data2, remap_colnames)
  }

  if (is.numeric(by)) {
    by <- colnames(meta.data2)[by]
  }

  meta.data <- join(x = meta.data, y = meta.data2, by = by, type = type, match = "first")
  meta.data <- meta.data[order(meta.data$Order), ] # This enables merging file that are not the same order.
  if (!identical(colnames(counts), meta.data$Barcode)) {
    print("Error, metadata files do not match.")
    quit()
  }
  print("  Added metadata file")

  rm(meta.data2)
  return(meta.data)
}

# Add extra cell annotation.
meta.data <- merge_metadata(meta.data = meta.data,
                            filepath = opt$cell_annotation,
                            by = c("Barcode"),
                            type = "inner")

# Add the WG1 metadata.
meta.data <- merge_metadata(meta.data = meta.data,
                            filepath = opt$wg1_metadata,
                            by = c("Pool", "Barcode"),
                            copy_colnames = list("IID" = "Assignment"),
                            type = "inner")

# Add the WG2 metadata. This file may miss some Barcodes.
meta.data <- merge_metadata(meta.data = meta.data,
                            filepath = opt$wg2_metadata,
                            by = c("Pool", "Barcode"))

# Add extra WG2 cell annotation if applicable.
if (!is.null(opt$wg2_pairing)) {
  meta.data <- merge_metadata(meta.data = meta.data,
                              filepath = opt$wg2_pairing,
                              delim = ";",
                              by = 2)
}

# Add the WG1 PSAM.
meta.data <- merge_metadata(meta.data = meta.data,
                            filepath = opt$wg1_psam,
                            by = c("IID"))

meta.data$Assignment_Run_Lane <- paste0(meta.data$Assignment, ";;", meta.data$sequencing_run, "_", meta.data$sequencing_lane)

print(paste0("  Combined full metadata with shape (", nrow(meta.data), ", ", ncol(meta.data), ")"))
if (nrow(meta.data) != nbarcodes) {
    print("Error, metadata has a different number of barcodes before and after merging.")
}

print("Saving full metadata")
write_delim(meta.data, gzfile(paste0(opt$out_dir, opt$pool, ".full.metadata.tsv.gz")), "\t")

# Create the Seurat object.
print("Creating Seurat object")
rownames(meta.data) <- meta.data$Barcode
min.cells <- 1
min.features <- 1
seurat <- CreateSeuratObject(counts, min.cells = min.cells, min.features = min.features, meta.data = meta.data)
nbarcodes.post.cells.features.filter = nrow(seurat@meta.data)
print(paste0("  CreateSeuratObject with min.cells = ", min.cells, " and min.features = ", min.features ," removed ", nbarcodes - nbarcodes.post.cells.features.filter, " cells."))

# Add the complexity meta data column.
print("  Adding complexity")
seurat@meta.data$complexity <- seurat@meta.data$nFeature_RNA / seurat@meta.data$nCount_RNA

# Read in the mt and rb gene lists
MT_genes <- read_delim(opt$mt_genes, delim = "\t")
RB_genes <- read_delim(opt$rb_genes, delim = "\t")

# Get the mitochondiral and ribosomal percentage QC metrics for each cell
print("  Adding mitochondiral %")
if ((sum(which(MT_genes$GeneID %in% rownames(seurat))) > sum(which(MT_genes$ENSG %in% rownames(seurat)))) & (sum(which(MT_genes$GeneID %in% rownames(seurat))) > length(grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))))){
    mt_features <- MT_genes$GeneID[MT_genes$GeneID %in% rownames(seurat)]
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features = mt_features)
} else if ((sum(which(MT_genes$ENSG %in% rownames(seurat))) > sum(which(MT_genes$GeneID %in% rownames(seurat)))) & (sum(which(MT_genes$ENSG %in% rownames(seurat))) > length(grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))))){
    mt_features <- MT_genes$ENSG[MT_genes$ENSG %in% rownames(seurat)]
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features = mt_features)
} else if ((length(grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))) > sum(which(MT_genes$GeneID %in% rownames(seurat)))) & (length(grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))) > sum(which(MT_genes$ENSG %in% rownames(seurat))))){
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features = rownames(seurat)[grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))])
} else {
    message("Either you do not have mitochondrial genes in your dataset or they are not labeled with ENSG IDs or Gene IDs")
}

print("  Adding ribosomal %")
if ((sum(which(RB_genes$GeneID %in% rownames(seurat))) > sum(which(RB_genes$ENSG %in% rownames(seurat)))) & (sum(which(RB_genes$GeneID %in% rownames(seurat))) > length(grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))))){
    rb_features <- RB_genes$GeneID[RB_genes$GeneID %in% rownames(seurat)]
    seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, features = rb_features)
} else if ((sum(which(RB_genes$ENSG %in% rownames(seurat))) > sum(which(RB_genes$GeneID %in% rownames(seurat)))) & (sum(which(RB_genes$ENSG %in% rownames(seurat))) > length(grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))))){
    rb_features <- RB_genes$ENSG[RB_genes$ENSG %in% rownames(seurat)]
    seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, features = rb_features)
} else if ((length(grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))) > sum(which(RB_genes$GeneID %in% rownames(seurat)))) & (length(grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))) > sum(which(RB_genes$ENSG %in% rownames(seurat))))){
    seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, features = rownames(seurat)[grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))])
} else {
    message("Either you do not have ribosomal genes in your dataset or they are not labeled with ENSG IDs or Gene IDs")
}

print("  Adding MALAT1")
if ("MALAT1" %in% rownames(seurat)) {
    seurat[["MALAT1"]] <- FetchData(object=seurat, vars = "MALAT1", layer = "counts")
} else {
    seurat[["MALAT1"]] <- 0
}

print("Saving QC metrics metadata")
write_delim(seurat@meta.data[, c("Barcode", "nCount_RNA", "nFeature_RNA", "complexity", "percent.mt", "percent.rb", "MALAT1")], gzfile(paste0(opt$out_dir, opt$pool, ".qc_metrics.tsv.gz")), "\t")

print("  Filtering for just singlets")
nbarcodes.post.singlet.filter <- 0
counts <- tryCatch({
    seurat <- subset(seurat, DropletType == "singlet")
    nbarcodes.post.singlet.filter <- nrow(seurat@meta.data)
},error = function(e){
    print(e)
})
print(paste0("  subset on DropletType 'singlet' removed ", nbarcodes.post.cells.features.filter - nbarcodes.post.singlet.filter, " cells."))

# Create a data frame over the filter stats.
filter.stats <- data.frame(c("Barcodes", paste0("cells >", min.cells, " - features >", min.features), "DropletType - Singlet"), c(nbarcodes, nbarcodes.post.cells.features.filter, nbarcodes.post.singlet.filter))
colnames(filter.stats) <- c("filter", opt$pool)

print("Saving filter stats")
write_delim(filter.stats, paste0(opt$out_dir, opt$pool, ".filter.stats.tsv"), "\t")

if (nbarcodes.post.singlet.filter == 0) {
    print("  Final Seurat object:")
    print("NA")
    quit()
}

print("  Final Seurat object:")
print(str(seurat))

# Save the seurat object and meta data
print("Saving output")
write_delim(seurat@meta.data, gzfile(paste0(opt$out_dir, opt$pool, ".metadata.tsv.gz")), "\t")
saveRDS(seurat, paste0(opt$out_dir, opt$pool, ".rds"))

print("Done")