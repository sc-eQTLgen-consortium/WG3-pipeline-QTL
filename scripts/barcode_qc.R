#!/usr/bin/env Rscript
# Author: M.J. Bonder (Adapted from 'QC_statistics_final.R' by M. Vochteloo)

# options parser
.libPaths("/usr/local/lib/R/site-library")
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list <- list(
  make_option(c("--functions_fn"), action="store", default='/tools/wg1-qc_filtering/scripts/QC_functions.R', type='character',
              help=""),
  make_option(c("--poolsheet"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--input_dir"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--level"), action="store", default="dataset", type='character',
              help="Level at which the QC filtering is peformed. By default, it will be performed at the whole dataset. You can chose to peform it at the level of a specific metadata variable (such as, lane/pool)."),
  make_option(c("--qc_mad"), action="store", default="/tools/wg1-qc_filtering/qc_mad_final.tab", type='character',
              help="QC metrics file. By default, the QC metrics-thresholds will be: lower nCount_RNA and upper percent.mt, from MAD 1 to 5."),
  make_option(c("--out_dir"), action="store", default=NA, type='character',
              help=""))
opt <- parse_args(OptionParser(option_list=option_list))

print("Options in effect:")
for (name in names(opt)) {
	print(paste0("  --", name, " ", opt[[name]]))
}
print("")

# check date is provided
if ((is.na(opt$poolsheet) || is.na(opt$input_dir) || is.na(opt$out_dir))) {
  stop("required parameters (poolsheet & input_dir & out_dir) must be provided.")
}

shhh(library(readr))

# Loading functions
print(paste0('Loading functions from: ', opt$functions_fn))
source(opt$functions_fn)

# Read QC-bound-MAD file
print(paste0('Reading QC metrics and MAD combinations file in: ',opt$qc_mad))
qc_mad.fn <- paste0(opt$qc_mad)
qc.mad_df <- read.table(qc_mad.fn, header = T)
qc.mad.temp_list <- split(qc.mad_df, qc.mad_df$QC_metric)
qc.mad_list <- lapply(qc.mad.temp_list, function(i) qc_mad_format(i))

# Setting MAD combination order
mad_label.list <- lapply(qc.mad_list, function(x) paste0('MAD_',x$mad))
mad.order <- unlist(lapply(mad_label.list[[1]], function(i) paste0(i,':',mad_label.list[[2]])))

# Print report
print(paste0('Considering ', nrow(qc.mad_df), ' QC metrics: ', paste(qc.mad_df$QC_metric, collapse = ', ')))
info <- lapply(names(qc.mad_list), function(qc){
  print(qc)
  print(paste0('MAD min: ', min(qc.mad_list[[qc]][["mad"]])))
  print(paste0('MAD max: ', max(qc.mad_list[[qc]][["mad"]])))
  cat('\n')
  return(NULL)
})

#################################### Summarize ####################################

# Read the poolsheet file.
pools <- read_delim(opt$poolsheet, delim = "\t")

# Merge the metadata files.
load_metadata <- function(row){
  message("Reading ", row[["Pool"]])
  metadata <- read_delim(paste0(opt$input_dir, row[["Pool"]], ".metadata.tsv.gz"), delim="\t")
  return(metadata)
}
metadata <- apply(pools, 1, load_metadata)
metadata <- as.data.frame(do.call(rbind,metadata))
metadata$index <- rownames(metadata)

# Calculate MADs and assigning Outlier/NotOutlier tags per cell barcode
print('Calculating MADs and assigning Outlier/NotOutlier tags per cell barcode...')
metadata <- qc_tag(so_md = metadata,
                   qc_mad_list = qc.mad_list,
                   filter_level = opt$level)
metadata <- metadata[order(as.numeric(metadata$index)),]
metadata$index <- NULL

# Save the full meta data
write_delim(metadata, gzfile(paste0(opt$out_dir, "metadata.tsv.gz")), "\t")

# Only keep those columns we just created new.
metadata <- metadata[, c("Pool", "Barcode", "nCount_RNA.tag", "nCount_RNA.mad", "percent.mt.tag", "percent.mt.mad", "mad_comb", "tag")]

# Save the reduced meta data
write_delim(metadata, gzfile(paste0(opt$out_dir, "metadata.tagged.tsv.gz")), "\t")