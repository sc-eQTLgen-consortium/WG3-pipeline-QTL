#!/usr/bin/env Rscript
# Author: M.J. Bonder (Adapted from 'https://github.com/sc-eQTLgen-consortium/WG3-pipeline-QTL/blob/master/scripts/process_output.R' by M. Vochteloo)


# options parser
.libPaths("/usr/local/lib/R/site-library")
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list <- list(
  make_option(c("--input"), action="store", default=NULL, type='character',
              help=""),
  make_option(c("--maf"), action="store", default=0.01, type='double',
              help=""),
  make_option(c("--r2"), action="store", default=0.6, type='double',
              help=""),
  make_option(c("--call"), action="store", default=0.99, type='double',
              help=""),
  make_option(c("--out_dir"), action="store", default="QC_statistics", type='character',
              help="Output directory.")
)
opt <- parse_args(OptionParser(option_list=option_list))

print("Options in effect:")
for (name in names(opt)) {
	print(paste0("  --", name, " ", opt[[name]]))
}
print("")

##Read information
variants <- data.table::fread(opt$input)

##Filter on MAF
variants <- variants[which(variants$MAF>=opt$maf),]

##Filter on imputation score
variants <- variants[which(variants$MACH_R2>=opt$r2),]

##Filter on call rate
variants <- variants[which(variants$CALL>=opt$call),]

#To upload to central site.
metadata_fn <- paste0(opt$out, gsub("_stats", "_stats_filtered", basename(opt$input)))
print(paste0('Saving output in: ', metadata_fn))
write.table(variants, file=gzfile(metadata_fn), quote=F, sep="\t", row.names=F)

for(chr in 1:22){
	#To filter in genotype harmonizer
	metadata_chr_fn <- paste0(opt$out, gsub("_stats.vars.gz", paste0("_inclusion_", chr, ".vars"), basename(opt$input)))
	print(paste0('Saving chromosome output in: ', metadata_chr_fn))
	write.table(variants$ID[which(variants$CHR==chr)], file=metadata_chr_fn, quote=F, sep="\t", row.names=F)
}

##End R

