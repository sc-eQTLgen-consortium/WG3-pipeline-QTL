#!/usr/bin/env Rscript
# Author: M. Vochteloo


# options parser
.libPaths("/usr/local/lib/R/site-library")
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list <- list(
  make_option(c("--kinship"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--wg1_psam"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--smf"), action="store", default=NA, type='character',
              help=""),
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

# check date is provided
if (any(unlist(lapply(c(opt$kinship, opt$wg1_psam, opt$data_out, opt$plot_out), is.na)))) {
  stop("required parameters (kinship, wg1_psam, data_out & plot_out) must be provided.")
}

shhh(library(reshape2))
shhh(library(ggplot2))

print("Loading data")
kin <- read.delim(opt$kinship, as.is = T, row.names = 1, check.names = F)
if (!identical(rownames(kin), colnames(kin))) {
   print("Rownames and colnames are not identical in kinship file.")
   quit()
 }

if (!is.na(opt$smf)) {
  print("Subsetting kinship matrix.")
  include_smf <- read.delim(opt$smf, header = F)
  samples <- unique(include_smf$V1)
  samples <- samples[samples %in% rownames(kin)]
  kin <- kin[samples, samples]
  if (!identical(rownames(kin), colnames(kin)) | !identical(rownames(kin), samples)) {
    print("Error in subsetting")
    quit()
  }
  print(paste0("  ", nrow(kin), " individuals remain."))
}

print("Calculating PCA")
pcOut <- prcomp(kin)
write.table(pcOut$x, paste0(opt$data_out, "Pcs.txt"), quote = F, sep = "\t", col.names = NA)

var_explained <- as.data.frame(pcOut$sdev^2 / sum(pcOut$sdev^2))
colnames(var_explained) <- c("VarianceExplained")
write.table(var_explained, paste0(opt$data_out, "Pcs.var.txt"), quote = F, sep = "\t", col.names = NA)

print("Plotting PCA")
pca <- as.data.frame(pcOut$x)[,1:(opt$n_pcs + 1)]
colnames(pca) <- paste0(colnames(pca), " ", round(var_explained[1:(opt$n_pcs + 1), "VarianceExplained"] * 100, 2), "%")

psam <- read.delim(opt$wg1_psam, header=T)
plot_df1 <- melt(merge(x = pca, y = psam[, c("IID", "Provided_Ancestry")], by.x = "row.names", by.y = "IID", all.x = TRUE), id = c("Row.names", "Provided_Ancestry"))
plot_df1$var_index <- as.numeric(gsub("PC", "", gsub(" .*", "", plot_df1$variable)))
plot_df2 <- data.frame(plot_df1)

colnames(plot_df1) <- c("IID", "Provided_Ancestry", "variable1", "value1", "var_index1")
colnames(plot_df2) <- c("IID", "Provided_Ancestry", "variable2", "value2", "var_index2")

plot_df <- merge(x = plot_df1, y = plot_df2, by = c("IID", "Provided_Ancestry"), all.x = TRUE, all.y = TRUE)
plot_df <- plot_df[plot_df$var_index1 < plot_df$var_index2, ]
rm(plot_df1, plot_df2)

pca_plot <- ggplot(plot_df, aes(value1, value2, colour=Provided_Ancestry)) +
  geom_point() +
  facet_grid(rows=vars(variable2), cols=vars(variable1), scales="free") +
  scale_color_manual(breaks = c("AFR", "AMR", "SAS", "EAS", "EUR"),
                     values=c("#d95f02", "#7570b3", "#e7298a", "#1b9e77", "#66a61e")) +
  xlab("") +
  ylab("") +
  ggtitle("PCA Plot")

ggsave(pca_plot, filename = paste0(opt$plot_out, "Pcs.png"), width = 29.7, height = 21, units = c("cm"))

var_explained$x <- as.numeric(rownames(var_explained))
scree_plot <- ggplot(var_explained, aes(x=x, y=VarianceExplained, group=1)) +
  geom_line() +
  geom_bar(stat="identity", color="black") +
  xlab("Principal Component") +
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  xlim(0, min(nrow(var_explained), 100))
ggsave(scree_plot, filename = paste0(opt$plot_out, "Pcs.var.png"), width = 29.7, height = 21, units = c("cm"))

