#!/usr/bin/env Rscript
# Author: M. Vochteloo

# options parser
.libPaths("/usr/local/lib/R/site-library")
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list <- list(
  make_option(c("--input"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--pvalue"), action="store", default="pvalue", type='character',
              help=""),
  make_option(c("--qvalue"), action="store", default="qvalue", type='character',
              help=""),
  make_option(c("--alpha"), action="store", default=0.05, type='numeric',
              help=""),
  make_option(c("--data_out"), action="store", default=NA, type='character',
              help="Output main directory"),
  make_option(c("--plot_out"), action="store", default=NA, type='character',
              help="Output plot directory"),
  make_option(c("--suffix"), action="store", default="_qvalues", type='character',
              help="Output suffix")
)
opt <- parse_args(OptionParser(option_list=option_list))

print("Options in effect:")
for (name in names(opt)) {
	print(paste0("  --", name, " ", opt[[name]]))
}
print("")

# check date is provided
if (any(unlist(lapply(c(opt$input, opt$data_out, opt$plot_out), is.na)))) {
  stop("required parameters (input, data_out & plot_out) must be provided.")
}

shhh(library(readr))
shhh(library(ggplot2))
shhh(library(qvalue))

print("Loading data")
data <- as.data.frame(read_delim(opt$input, delim = "\t"))
data <- data[!is.na(data[[opt$pvalue]]), ]
data <- data[order(data[[opt$pvalue]]), ]
print(paste0("  ", nrow(data), " genes loaded"))

print("Counting pvalues")
print(paste0("  significant at pvalue<", opt$alpha, ": ", sum(data[[opt$pvalue]] < opt$alpha), " genes"))

print("Calculating qvalues")
qobj <- qvalue(p = data[[opt$pvalue]])
data[[opt$qvalue]] <- qobj$qvalues
print(paste0("  significant at qvalue<", opt$alpha, ": ", sum(data[[opt$qvalue]] < opt$alpha), " genes"))

print("Saving data")
write.table(data, paste0(opt$data_out, opt$suffix, ".txt"), quote = F, sep = "\t", row.names = FALSE)

dir.create(dirname(opt$plot_out), recursive = TRUE, showWarnings = FALSE)
setwd(dirname(opt$plot_out))

print("Creating figures")
png(paste0(opt$plot_out, "_overview.png"))
plot1 <- plot(qobj)
dev.off()

plot2 <- hist(qobj)
ggsave(plot2, filename = paste0(opt$plot_out, "_hist.png"), width = 29.7, height = 21 ,units = c("cm"))

print("Done")