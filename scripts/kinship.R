#!/usr/bin/env Rscript
# Author: M.J. Bonder (Adapted from 'https://github.com/sc-eQTLgen-consortium/WG3-pipeline-QTL/blob/master/scripts/kinship.R' by M. Vochteloo)


# options parser
.libPaths("/usr/local/lib/R/site-library")
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list <- list(
  make_option(c("--king"), action="store", type='character',
              help=""),
  make_option(c("--king_id"), action="store", type='character',
              help=""),
  make_option(c("--out"), action="store", type='character',
              help="")
)
opt <- parse_args(OptionParser(option_list=option_list))

print("Options in effect:")
for (name in names(opt)) {
	print(paste0("  --", name, " ", opt[[name]]))
}
print("")

kin <- read.delim(opt$king,header=F)
kinIds <- read.delim(opt$king_id)

kin <- kin * 2

rownames(kin) <- colnames(kin) <- kinIds[,1]

write.table(kin, opt$out, sep="\t", quote=F, col.names=NA)
