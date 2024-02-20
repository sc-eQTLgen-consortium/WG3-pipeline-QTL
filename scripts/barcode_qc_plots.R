#!/usr/bin/env Rscript
# Adapted from drneavin https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/blob/master/Demultiplexing/scripts/Singlet_QC_Figures.R by M. Vochteloo

# options parser
.libPaths("/usr/local/lib/R/site-library")
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list <- list(
  make_option(c("--metadata"), action="store", default=NA, type='character',
              help=""),
  make_option(c("--out"), action="store", default=NA, type='character',
              help=""))
opt <- parse_args(OptionParser(option_list=option_list))

print("Options in effect:")
for (name in names(opt)) {
	print(paste0("  --", name, " ", opt[[name]]))
}
print("")

# check date is provided
if ((is.na(opt$metadata) || is.na(opt$out))) {
  stop("required parameters (metadata & out) must be provided. ")
}

shhh(library(tidyr))
shhh(library(tidyverse))
shhh(library(ggplot2))
shhh(library(ggforce))
shhh(library(ggnewscale))

metadata <- as.data.frame(read_delim(opt$metadata, delim="\t"))
metadata$Pool <- as.character(metadata$Pool)
pools_list <- unique(metadata$Pool)

##### Create a dataframe of MADs #####
MAD_df_list <- lapply(pools_list, function(x){
    as.data.frame(matrix(ncol = 2, nrow = length((which(c("percent.rb","percent.mt","nCount_RNA","nFeature_RNA") %in% colnames(metadata))))))
})
names(MAD_df_list) <- pools_list

MAD_df_list <- lapply(names(MAD_df_list), function(x){
    rownames(MAD_df_list[[x]]) <- colnames(metadata)[which(colnames(metadata) %in% c("percent.rb","percent.mt","nCount_RNA","nFeature_RNA"))]
    colnames(MAD_df_list[[x]]) <- c("Median","MAD")
    for (QC in rownames(MAD_df_list[[x]])){
        MAD_df_list[[x]][QC, "Median"] <- median(metadata[which(metadata$Pool == x),QC], na.rm = TRUE)
        MAD_df_list[[x]][QC, "MAD"] <- mad(metadata[which(metadata$Pool == x),QC], center = MAD_df_list[[x]][QC, "Median"],  constant = 1.4826, na.rm = TRUE,low = FALSE, high = FALSE)
    }
    MAD_df_list[[x]]$Pool <- x
    MAD_df_list[[x]]$QC_Metric <- rownames(MAD_df_list[[x]])
    return(MAD_df_list[[x]])
})
MAD_df <- do.call(rbind,MAD_df_list)
MAD_df$Pool <- as.factor(MAD_df$Pool)

MAD_df_All <- as.data.frame(matrix(ncol = 2, nrow = length((which(c("percent.rb","percent.mt","nCount_RNA","nFeature_RNA") %in% colnames(metadata))))))
rownames(MAD_df_All) <- colnames(metadata)[which(colnames(metadata) %in% c("percent.rb","percent.mt","nCount_RNA","nFeature_RNA"))]
colnames(MAD_df_All) <- c("Median","MAD")
for (QC in rownames(MAD_df_All)){
    MAD_df_All[QC, "Median"] <- median(metadata[,QC], na.rm = TRUE)
    MAD_df_All[QC, "MAD"] <- mad(metadata[,QC], center = MAD_df_All[QC, "Median"],  constant = 1.4826, na.rm = TRUE,low = FALSE, high = FALSE)
}
MAD_df_All$QC_Metric <- rownames(MAD_df_All)

color_gradient <- c("grey", "blue3", "darkviolet", "firebrick1", "darkorange1", "gold1")
# color_gradient <- c("grey0", "grey10", "grey30", "grey50", "grey70", "grey90")
color_values <- c("Median"=color_gradient[1],
                  "1 MAD"=color_gradient[2],
                  "2 MAD"=color_gradient[3],
                  "3 MAD"=color_gradient[4],
                  "4 MAD"=color_gradient[5],
                  "5 MAD"=color_gradient[6],
                  "MAD1"=color_gradient[2],
                  "MAD2"=color_gradient[3],
                  "MAD3"=color_gradient[4],
                  "MAD4"=color_gradient[5],
                  "MAD5"=color_gradient[6])


# MAD_df_long <- pivot_longer(MAD_df, names_to = "QC_metric", cols = c("mt.percent","nUMI","nFeature"), values_to = "value")

##### Make figures #####
### Mt % ###
## Add multiple MAD lines ##
violins_MADperPOOL <- list()
violins_MAD_ALL <- list()
violins <- list()

for (QC in unique(MAD_df$QC_Metric)){
    if (QC == "percent.mt" | QC == "percent.rb"){
        violins_MADperPOOL[[QC]] <- ggplot(metadata, aes( x = Pool, y = metadata[,QC])) +
                    geom_violin() +
                    geom_sina(size = 1, alpha = 0.6) +
                    theme_classic() +
                    ylim(0, min(max(metadata[,QC]),100)) +
                    labs(title = paste0(QC, " Violin Plot"), y = QC, x = "Pool") +
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median, yend=Median, col = "Median"), size = 1) +
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+MAD, yend=Median+MAD,col = "1 MAD"), size = 1, linetype = "longdash") +
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-MAD, yend=Median-MAD, col = "1 MAD"), size = 1, linetype = "longdash") +
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+2*MAD, yend=Median+2*MAD, col = "2 MAD"), size = 1, linetype = "longdash") +
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-2*MAD, yend=Median-2*MAD, col = "2 MAD"), size = 1, linetype = "longdash") +
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+3*MAD, yend=Median+3*MAD, col = "3 MAD"), size = 1, linetype = "longdash") +
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-3*MAD, yend=Median-3*MAD, col = "3 MAD"), size = 1, linetype = "longdash") +
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+4*MAD, yend=Median+4*MAD, col = "4 MAD"), size = 1, linetype = "longdash") +
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-4*MAD, yend=Median-4*MAD, col = "4 MAD"), size = 1, linetype = "longdash") +
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+5*MAD, yend=Median+5*MAD, col = "5 MAD"), size = 1, linetype = "longdash") +
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-5*MAD, yend=Median-5*MAD, col = "5 MAD"), size = 1, linetype = "longdash") +
                    scale_color_manual("MAD", values=color_values) +
                    theme(text = element_text(size=14),
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        plot.title = element_text(hjust = 0.5),
                        plot.subtitle = element_text(hjust = 0.5))

        violins_MAD_ALL[[QC]] <- ggplot(metadata, aes( x = Pool, y = metadata[,QC])) +
                    geom_violin() +
                    geom_sina(size = 1, alpha = 0.6) +
                    theme_classic() +
                    ylim(0, min(max(metadata[,QC]),100)) +
                    labs(title = paste0(QC, " Violin Plot"), y = QC, x = "Pool") +
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median, yend=Median, col = "Median"), size = 1) +
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median+MAD, yend=Median+MAD,col = "1 MAD"), size = 1, linetype = "longdash") +
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median-MAD, yend=Median-MAD, col = "1 MAD"), size = 1, linetype = "longdash") +
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median+2*MAD, yend=Median+2*MAD, col = "2 MAD"), size = 1, linetype = "longdash") +
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median-2*MAD, yend=Median-2*MAD, col = "2 MAD"), size = 1, linetype = "longdash") +
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median+3*MAD, yend=Median+3*MAD, col = "3 MAD"), size = 1, linetype = "longdash") +
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median-3*MAD, yend=Median-3*MAD, col = "3 MAD"), size = 1, linetype = "longdash") +
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median+4*MAD, yend=Median+4*MAD, col = "4 MAD"), size = 1, linetype = "longdash") +
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median-4*MAD, yend=Median-4*MAD, col = "4 MAD"), size = 1, linetype = "longdash") +
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median+5*MAD, yend=Median+5*MAD, col = "5 MAD"), size = 1, linetype = "longdash") +
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median-5*MAD, yend=Median-5*MAD, col = "5 MAD"), size = 1, linetype = "longdash") +
                    scale_color_manual("MAD", values=color_values) +
                    theme(text = element_text(size=14),
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        plot.title = element_text(hjust = 0.5),
                        plot.subtitle = element_text(hjust = 0.5))
    } else {
        violins_MADperPOOL[[QC]] <- ggplot(metadata, aes( x = Pool, y = metadata[,QC])) +
                            geom_violin() +
                            geom_sina(size = 1, alpha = 0.6) +
                            theme_classic() +
                            ylim(0, NA) +
                            labs(title = paste0(QC, " Violin Plot"), y = QC, x = "Pool") +
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median, yend=Median, col = "Median"), size = 1) +
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+MAD, yend=Median+MAD,col = "1 MAD"), size = 1, linetype = "longdash") +
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-MAD, yend=Median-MAD, col = "1 MAD"), size = 1, linetype = "longdash") +
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+2*MAD, yend=Median+2*MAD, col = "2 MAD"), size = 1, linetype = "longdash") +
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-2*MAD, yend=Median-2*MAD, col = "2 MAD"), size = 1, linetype = "longdash") +
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+3*MAD, yend=Median+3*MAD, col = "3 MAD"), size = 1, linetype = "longdash") +
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-3*MAD, yend=Median-3*MAD, col = "3 MAD"), size = 1, linetype = "longdash") +
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+4*MAD, yend=Median+4*MAD, col = "4 MAD"), size = 1, linetype = "longdash") +
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-4*MAD, yend=Median-4*MAD, col = "4 MAD"), size = 1, linetype = "longdash") +
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+5*MAD, yend=Median+5*MAD, col = "5 MAD"), size = 1, linetype = "longdash") +
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-5*MAD, yend=Median-5*MAD, col = "5 MAD"), size = 1, linetype = "longdash") +
                            scale_color_manual("MAD", values=color_values) +
                            theme(text = element_text(size=14),
                                axis.text.x = element_text(angle = 45, hjust = 1),
                                plot.title = element_text(hjust = 0.5),
                                plot.subtitle = element_text(hjust = 0.5))

        violins_MAD_ALL[[QC]] <- ggplot(metadata, aes( x = Pool, y = metadata[,QC])) +
                            geom_violin() +
                            geom_sina(size = 1, alpha = 0.6) +
                            theme_classic() +
                            ylim(0, NA) +
                            labs(title = paste0(QC, " Violin Plot"), y = QC, x = "Pool") +
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median, yend=Median, col = "Median"), size = 1) +
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median+MAD, yend=Median+MAD,col = "1 MAD"), size = 1, linetype = "longdash") +
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median-MAD, yend=Median-MAD, col = "1 MAD"), size = 1, linetype = "longdash") +
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median+2*MAD, yend=Median+2*MAD, col = "2 MAD"), size = 1, linetype = "longdash") +
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median-2*MAD, yend=Median-2*MAD, col = "2 MAD"), size = 1, linetype = "longdash") +
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median+3*MAD, yend=Median+3*MAD, col = "3 MAD"), size = 1, linetype = "longdash") +
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median-3*MAD, yend=Median-3*MAD, col = "3 MAD"), size = 1, linetype = "longdash") +
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median+4*MAD, yend=Median+4*MAD, col = "4 MAD"), size = 1, linetype = "longdash") +
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median-4*MAD, yend=Median-4*MAD, col = "4 MAD"), size = 1, linetype = "longdash") +
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median+5*MAD, yend=Median+5*MAD, col = "5 MAD"), size = 1, linetype = "longdash") +
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median-5*MAD, yend=Median-5*MAD, col = "5 MAD"), size = 1, linetype = "longdash") +
                            scale_color_manual("MAD", values=color_values) +
                            theme(text = element_text(size=14),
                                axis.text.x = element_text(angle = 45, hjust = 1),
                                plot.title = element_text(hjust = 0.5),
                                plot.subtitle = element_text(hjust = 0.5))
    }
    violins[[QC]] <- ggplot(metadata, aes( x = Pool, y = metadata[,QC])) +
                    geom_violin() +
                    geom_sina(size = 1, alpha = 0.6) +
                    theme_classic() +
                    ylim(0, NA) +
                    labs(title = paste0(QC, " Violin Plot"), y = QC, x = "Pool") +
                    theme(text = element_text(size=14),
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        plot.title = element_text(hjust = 0.5),
                        plot.subtitle = element_text(hjust = 0.5))
}

for (QC in names(violins_MAD_ALL)){
    ggsave(violins_MADperPOOL[[QC]], filename = paste0(opt$out, QC, "_violin_MADper_Pool.png"), width = 29.7, height = 21 ,units = c("cm"))
    ggsave(violins_MAD_ALL[[QC]], filename = paste0(opt$out, QC, "_violin_MAD_All.png"), width = 29.7, height = 21 ,units = c("cm"))
    ggsave(violins[[QC]], filename = paste0(opt$out, QC, "_violin_noMADlines.png"), width = 29.7, height = 21 ,units = c("cm"))
}

MAD_df_All$MAD1_up <- MAD_df_All$Median + MAD_df_All$MAD
MAD_df_All$MAD1_down <- MAD_df_All$Median - MAD_df_All$MAD
MAD_df_All$MAD2_up <- MAD_df_All$Median + 2*MAD_df_All$MAD
MAD_df_All$MAD2_down <- MAD_df_All$Median - 2*MAD_df_All$MAD
MAD_df_All$MAD3_up <- MAD_df_All$Median + 3*MAD_df_All$MAD
MAD_df_All$MAD3_down <- MAD_df_All$Median - 3*MAD_df_All$MAD
MAD_df_All$MAD4_up <- MAD_df_All$Median + 4*MAD_df_All$MAD
MAD_df_All$MAD4_down <- MAD_df_All$Median - 4*MAD_df_All$MAD
MAD_df_All$MAD5_up <- MAD_df_All$Median + 5*MAD_df_All$MAD
MAD_df_All$MAD5_down <- MAD_df_All$Median -5*MAD_df_All$MAD

MAD_df_All_long <- pivot_longer(MAD_df_All, cols = c("Median", "MAD1_up", "MAD1_down", "MAD2_up","MAD2_down", "MAD3_up", "MAD3_down", "MAD4_up", "MAD4_down","MAD5_up", "MAD5_down"), names_to = "Metric")
MAD_df_All_long <- separate(MAD_df_All_long, col = Metric, sep = "_", into = c("MAD_metric", "direction"))


### Number Genes vs Number UMI per Cell ###
## Color by Individual ##
## Add multiple MAD lines ##
if ("percent.mt" %in% unique(MAD_df_All$QC_Metric) & ("nCount_RNA" %in% MAD_df_All$QC_Metric) & ("nFeature_RNA" %in% MAD_df_All$QC_Metric)){
    pUMI_MTscatter_MAD <- ggplot(metadata, mapping = aes(x = nCount_RNA, y = percent.mt)) +
        geom_point(size = 0.5, alpha = 0.5) +
        theme_classic() +
        labs(x = "Number UMI", y = "Percent Mitochondial Genes") +
        xlim(0, NA) +
        ylim(0, NA) +
        geom_segment(data=MAD_df_All_long[which(MAD_df_All_long$QC_Metric == "nCount_RNA"),], aes(x=value, xend = value, color = MAD_metric), y = min(metadata$percent.mt), yend = max(metadata$percent.mt), linetype = "longdash") +
        geom_segment(data=MAD_df_All_long[which(MAD_df_All_long$QC_Metric == "percent.mt"),], aes(y=value, yend = value, color = MAD_metric), x = min(metadata$nCount_RNA), xend = max(metadata$nCount_RNA), linetype = "longdash") +
        scale_color_manual(values = color_values)
    ggsave(pUMI_MTscatter_MAD, filename = paste0(opt$out, "UMI_vs_percentMT_QC_scatter_w_MADlines.png"))

    pNfeatures_MTscatter_MAD <- ggplot(metadata, mapping = aes(x = nFeature_RNA, y = percent.mt)) +
        geom_point(size = 0.5, alpha = 0.5) +
        theme_classic() +
        labs(x = "Number Features", y = "Percent Mitochondial Genes") +
        xlim(0, NA) +
        ylim(0, NA) +
        geom_segment(data=MAD_df_All_long[which(MAD_df_All_long$QC_Metric == "nFeature_RNA"),], aes(x=value, xend = value, color = MAD_metric), y = min(metadata$percent.mt), yend = max(metadata$percent.mt), linetype = "longdash") +
        geom_segment(data=MAD_df_All_long[which(MAD_df_All_long$QC_Metric == "percent.mt"),], aes(y=value, yend = value, color = MAD_metric), x = min(metadata$nFeature_RNA), xend = max(metadata$nFeature_RNA), linetype = "longdash") +
        scale_color_manual(values = color_values)
    ggsave(pNfeatures_MTscatter_MAD, filename = paste0(opt$out, "nFeatures_vs_percentMT_QC_scatter_w_MADlines.png"))

    pUMI_MTscatter <- ggplot(metadata, aes(x = nCount_RNA, y = percent.mt, color = Pool)) +
        geom_point(alpha = 0.5, size = 0.5) +
        theme_classic() +
        scale_color_manual(values = rainbow(length(unique(metadata$Pool))), name = "Pool")+
        labs(x = "Number UMI", y = "Percent Mitochondial Genes")
    ggsave(pUMI_MTscatter, filename = paste0(opt$out, "/UMI_vs_percentMT_QC_scatter_colorPool.png"))

    pNfeature_MTscatter <- ggplot(metadata, aes(x = nFeature_RNA, y = percent.mt, color = Pool)) +
        geom_point(alpha = 0.5, size = 0.5) +
        theme_classic() +
        scale_color_manual(values = rainbow(length(unique(metadata$Pool))), name = "Pool")+
        labs(x = "Number Features", y = "Percent Mitochondial Genes")
    ggsave(pNfeature_MTscatter, filename = paste0(opt$out, "nFeatures_vs_percentMT_QC_scatter_colorPool.png"))

}

if ("nFeature_RNA" %in% unique(MAD_df_All$QC_Metric) & ("nCount_RNA" %in% MAD_df_All$QC_Metric)){
    pUMI_Genes_scatter_MAD <- ggplot() +
        geom_point(data = metadata, aes(x = nCount_RNA, y = nFeature_RNA), alpha = 0.25, size = 0.5) +
        theme_classic() +
        labs(x = "Number UMI", y = "Number Genes")  +
        xlim(0, NA) +
        ylim(0, NA) +
        geom_segment(data=MAD_df_All_long[which(MAD_df_All_long$QC_Metric == "nCount_RNA"),], aes(x=value, xend = value, color = MAD_metric), y = min(metadata$nFeature_RNA), yend = max(metadata$nFeature_RNA), linetype = "longdash") +
        geom_segment(data=MAD_df_All_long[which(MAD_df_All_long$QC_Metric == "nFeature_RNA"),], aes(y=value, yend = value, color = MAD_metric), x = min(metadata$nCount_RNA), xend = max(metadata$nCount_RNA), linetype = "longdash") +
        scale_color_manual(values = color_values)
    ggsave(pUMI_Genes_scatter_MAD, filename = paste0(opt$out, "UMI_vs_Genes_QC_scatter_w_MADlines.png"))

    pUMI_Genes_scatter<- ggplot() +
        geom_point(data = metadata, aes(x = nCount_RNA, y = nFeature_RNA, color = as.factor(Pool)), alpha = 0.25, size = 0.5) +
        theme_classic() +
        scale_color_manual(values = rainbow(length(unique(metadata$Pool))), name = "Pool")+
        labs(x = "Number UMI", y = "Number Genes")

    ggsave(pUMI_Genes_scatter, filename = paste0(opt$out, "UMI_vs_Genes_QC_scatter.png"))
}
