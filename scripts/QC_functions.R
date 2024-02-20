#!/usr/bin/env Rscript

################################## Load R packages ################################
print('Loading R packages...')
shhh <- suppressPackageStartupMessages
shhh(library(Seurat))
shhh(library(SeuratObject))
shhh(library(stringr))
shhh(library(stringi))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(reshape2))
shhh(library(data.table))
shhh(library(ggplot2))
shhh(library(ggpubr))
shhh(library(scales))
shhh(library(RColorBrewer))
shhh(library(pheatmap))
shhh(library(grid))

########################### Functions used in QC_statistics.R #####################
# 1. Format QC and MAD
qc_mad_format <- function(df){
  th_vec <- unlist(str_split(df$bound,','))
  mad_vec <- seq(df$MAD_min, df$MAD_max, 0.5)
  out_list <- list(th = th_vec,
                   mad = mad_vec)
  return(out_list)
}

# 2. Add Azimuth l1 classification to the metadata
add_azimuth_l1 <- function(so_md, pairing){
  # dicctionary 
  ct.vec <- pairing$l1
  names(ct.vec) <- pairing$l2
  
  # check
  setdiff(unique(so_md$predicted.celltype.l2),names(ct.vec)) #it should be 0
  
  # add the new column using the dictionary 'ct.vec'
  so_md$predicted.celltype.l1 <- ct.vec[so_md$predicted.celltype.l2]
  
  return(so_md)
}

# 3. Downsample by metadata (optional)
downsample_by_metadata <- function(so_md, n_cells, md_var){
  ds_md <- so_md
  md_var.level <- unique(so_md[[md_var]])
  n_cells.level <- nrow(so_md)
  if(n_cells<n_cells.level){
    set.seed(123)
    n_row <- sample(nrow(so_md),n_cells)
    ds_md <- droplevels(so_md[n_row,])
  }else{
    n_cells <- n_cells.level
  }
  print(paste0('Downsampling ', md_var, ' (', md_var.level, ') to: ', n_cells))
  return(ds_md)
}

# 4. MAD function (called in qc_mad.func)
mad.func <- function(number_mad, qc, th, so_md){
  print(paste0('MAD: ', number_mad))
  mad <- mad(so_md[,qc])
  low <- median(so_md[,qc]) - number_mad*mad
  high <- median(so_md[,qc]) + number_mad*mad
  if(all(c("lower","upper")%in%th)){
    print('lower and upper boundaries...')
    print("The lower bound is:")
    print(low)
    print("The upper bound is:")
    print(high)
    so_md$tag <- ifelse((so_md[,qc] > low & so_md[,qc] < high),"NotOutlier", "Outlier")
  }else if(th=="upper"){
    print('only upper boundary...')
    print("The upper bound is:")
    print(high)
    so_md$tag <- ifelse(so_md[,qc] < high, "NotOutlier", "Outlier")
  }else if(th=="lower"){
    print('only lower boundary...')
    print("The lower bound is:")
    print(low)
    so_md$tag <- ifelse(so_md[,qc] > low, "NotOutlier", "Outlier")
  }else{
    print(paste0('Provide a valid bound option in the qc bound in the QC_metrics.tab: lower,upper / lower / upper.'))
  }
  so_md$qc_metric <- qc
  so_md$mad <- paste0('MAD_',number_mad)
  cat('\n')
  return(so_md)
}

# 5. MAD function by QC metric
qc_mad.func <- function(qc, qc_mad_list, so_md, filter_level){
  print(qc)
  qc_th <- qc_mad_list[[qc]][['th']]
  mad_vec <- qc_mad_list[[qc]][['mad']]
  print(qc_th)
  if(filter_level!='dataset'){
    print(paste0('Peforming the QC filtering by ', filter_level, '...'))
    cat('\n')
    so_level.list <- split(so_md, so_md[[filter_level]])
    mad_list <- lapply(mad_vec, function(n){
      tag_by_level <- do.call("rbind",lapply(so_level.list, function(so_level) mad.func(number_mad = n, 
                                                                                        qc = qc, 
                                                                                        th = qc_th, 
                                                                                        so_md = so_level)))
    })
  }else{
    print(paste0('Peforming the QC filtering across the whole dataset...'))
    cat('\n')
    mad_list <- lapply(mad_vec, function(n) mad.func(number_mad = n, 
                                                     qc = qc, 
                                                     th = qc_th, 
                                                     so_md = so_md))
  }
  names(mad_list) <- paste0('MAD_',mad_vec)
  return(mad_list)
}

# 6. Join tag (Outlier/NotOutlier)
join_tag.func <- function(qc_list, qc_mad_list){
  system.time(qc.list_temp <- sapply(names(qc_list), function(qc){
    qc.df <- qc_list[[qc]]
    cnames <- c('tag','mad')
    colnames(qc.df)[colnames(qc.df)%in%cnames] <- paste0(qc,'.',cnames)
    qc.df <- qc.df[,-which(colnames(qc.df)=='qc_metric')]
    return(qc.df)
  }, simplify = F))
  system.time(cnames_int <- Reduce(intersect, lapply(qc.list_temp, colnames)))
  system.time(qc.list_temp <- lapply(qc.list_temp, function(x) data.table(x)))
  system.time(qc.merged <- Reduce(function(x, y) data.table::merge.data.table(x, y, by = cnames_int, allow.cartesian=TRUE), qc.list_temp)) #revisit --> 0.862 (10K)
  system.time(qc_mad.cols <- paste0(names(qc_mad_list),'.mad'))
  system.time(qc.merged$mad_comb <- apply(qc.merged[,..qc_mad.cols],1,paste,collapse = ":")) #revisit (sub) --> 3.181 (10K)
  system.time(qc_tag.cols <- paste0(names(qc_mad_list),'.tag'))
  system.time(qc.merged$tag <- ifelse(qc.merged[[qc_tag.cols[1]]]=='NotOutlier' & qc.merged[[qc_tag.cols[2]]]=='NotOutlier', 'NotOutlier', 'Outlier'))
  # system.time(qc.outliers <- qc.merged[apply(qc.merged[,..qc_tag.cols]=='Outlier', 1, all),])
  # system.time(qc.outliers$tag <- 'Outlier')
  # system.time(qc.merged <- data.table::merge.data.table(qc.merged, qc.outliers, by = colnames(qc.merged), all.x = TRUE)) #revisit --> 0.821 (10K)
  # system.time(qc.merged[is.na(qc.merged$tag),]$tag <- 'NotOutlier')
  return(qc.merged)
}

# 7. QC and tag function
qc_tag <- function(so_md, qc_mad_list, filter_level){
  ## apply qc_mad.func
  print('Applying qc_mad.func function to each QC metric: assigning Outlier/NotOutlier tag to each cell depending on different MADs... ')
  system.time(qc_mad.list <- sapply(names(qc_mad_list), function(i) qc_mad.func(qc =i, 
                                                                                qc_mad_list = qc_mad_list, 
                                                                                so_md = so_md,
                                                                                filter_level = filter_level), simplify = F))
  
  ## rbind MADs dataframes inside the QC list
  print('Joining Outlier/NotOutlier for each QC metric...')
  system.time(qc.list <- lapply(qc_mad.list, function(qc){
    mad.df <- do.call("rbind",qc)
    rownames(mad.df) <- NULL
    return(mad.df)
  }))
  system.time(qc.df <- do.call("rbind", qc.list))
  system.time(rownames(qc.df) <- NULL)
  
  ## join tag (Outlier/NotOutlier)
  print('Join tag (Outlier/NotOutlier)...')
  system.time(join_tag.df <- join_tag.func(qc_list = qc.list, 
                                           qc_mad_list = qc_mad_list))
  return(join_tag.df)
}

# 8. Add missing tags (Outlier/NotOutlier) (called in summarize_by_md.func)
add_missing.func <- function(mad_c, df, md_var=NULL, md_level=NULL, tags = c('NotOutlier','Outlier')){
  if(!is.null(md_var)){
    df.mad_c <- droplevels(df[df$mad_comb==mad_c & df[[md_var]]==md_level,])
  }else{
    df.mad_c <- droplevels(df[df$mad_comb==mad_c,])
  }
  df.out <- df.mad_c
  if(nrow(df.mad_c)==1){
    tag_missing <- setdiff(tags, df.mad_c$tag)
    df.add <- data.frame(mad_comb = mad_c,
                         tag = tag_missing,
                         n = 0,
                         prop = 0,
                         pct = 0)
    if(!is.null(md_var)){
      print(paste0('Adding missing tag in: ', mad_c, ' - ', md_level, ' (', tag_missing, ')'))
      df.add[[md_var]] <- md_level
    }else{
      print(paste0('Adding missing tag in: ', mad_c, ' (', tag_missing, ')'))
    }
    df.out <- rbind(df.out, df.add)
  }
  return(df.out)
}

# 9. Summarize count and add BY metadata variable (inside metadata type)
summarize.func <- function(df, md_var = NULL){
  # Count
  print('Counting Outlier/NotOutlier tags...')
  vars_i <- c(md_var,'mad_comb')
  vars_all <- c(vars_i,'tag')
  
  df %>%
    group_by_at(vars(vars_all), .drop=FALSE) %>%
    count() -> df.c
  
  df.c %>%
    group_by_at(vars(vars_i)) %>%
    mutate(prop = n / sum(n),
           pct = (n / sum(n))*100) %>% 
    as.data.frame() -> tag.c
  
  # Add missing mad comb
  print('Adding missing Outlier/NotOutlier tags...')
  if(is.null(md_var)){
    md_levels <- NULL
    print('Across the whole dataset...')
    add_missing.res <- lapply(mad.order, function(i) add_missing.func(mad_c = i, 
                                                                      df = tag.c))
    tag.out <- do.call("rbind",add_missing.res)
  }else{
    md_levels <- md_vec.order[[md_var]]
    print('By metadata variable...')
    add_missing.res <- lapply(mad.order, function(i)
      lapply(md_levels, function(j) add_missing.func(mad_c = i,
                                                     df = tag.c,
                                                     md_var = md_var,
                                                     md_level = j)))
    
    tag.out <- do.call("rbind",lapply(add_missing.res, function(x) do.call("rbind",x)))
    tag.out[[md_var]] <- factor(tag.out[[md_var]],
                                levels = md_levels)
  }
  tag.out$mad_comb <- factor(tag.out$mad_comb,
                             levels = mad.order)
  return(tag.out)
}

# 10. Main MAD function (by dataset)
count_by_dataset <- function(qc_tag, md_type=NULL){
  print('Counting nCells in each MAD combination')
  if(is.null(md_type)){
    print('Summarizing by dataset...')
    summarize.res <- summarize.func(df = qc_tag)
  }else{
    print(paste0('Summarizing by metadata variable: ', md_type, '...'))
    md_vars <- md_list[[md_type]]
    summarize.res <- sapply(md_vars, function(i) summarize.func(df = qc_tag,
                                                                md_var = i), simplify = F)
  }
  return(summarize.res)
}

# 11. Save pheatmap in pdf
save_pheatmap_pdf <- function(x, filename, width, height) {
  pdf(filename, width = width, height = height)
  # vp = viewport(height=unit(50, "inches"), width=unit(50, "inches"))
  grid::grid.newpage()
  # grid.rect(vp=vp,gp=gpar(col="white"))
  grid::grid.draw(x$gtable)
  dev.off()
}

# 10. Accessory heatmap function (by tags, values and metadata variable)
pheatmap.func <- function(md_var, tag, value, display_numbers = "none", cluster_mads = F, cluster_md = F, tag_md_list, mat_cols = colorRampPalette(brewer.pal(n = 9, name = "Reds"))(100), out_dir, width_list){
  # Convert long to wide fromat (data fram --> matrix)
  tag.df <- tag_md_list[[md_var]][[tag]]
  colnames(tag.df)[1] <- 'md_level'
  tag.mat <- reshape2::dcast(tag.df, mad_comb~md_level, value.var=value)
  tag.mat[is.na(tag.mat)] <- 0
  rownames(tag.mat) <- tag.mat[,1]
  if(ncol(tag.mat)>2){
    tag.mat <- tag.mat[,-1]
    mat <- as.matrix(tag.mat)
  }else{
    rnames <- rownames(tag.mat)
    cnames <- unique(tag.df$md_level)
    mat <- matrix(tag.mat[,-1], dimnames = list(rnames, cnames))
    cluster_md <- FALSE
  }
  
  # Pheatmap
  print(md_var)
  print(tag)
  print(value)
  print(display_numbers)
  print(cluster_mads)
  print(cluster_md)
  ## Annotation
  ### annotation rows: MAD combinations
  tag.df[,c('nCount_RNA','percent.mt')] <- str_split_fixed(tag.df$mad_comb, ':', 2)[,c(1,2)]
  mad_row <- tag.df[,c('mad_comb','nCount_RNA','percent.mt')]
  mad_row <- unique(mad_row)
  rownames(mad_row) <- mad_row$mad_comb
  mad_row <- mad_row[-1]
  
  ## List of annotation colors
  nCount_RNA.vec <- rev(brewer.pal(9,'Greens'))
  names(nCount_RNA.vec) <- unique(mad_row$nCount_RNA)
  percent.mt.vec <- rev(brewer.pal(9,'Purples'))
  names(percent.mt.vec) <- unique(mad_row$percent.mt)
  annot_cols <- list(nCount_RNA = nCount_RNA.vec,
                     percent.mt = percent.mt.vec)
  
  
  ## Pheatmap
  main.var <- paste0(dataset, '\n', md_var, ' - ', tag, ' (', value, ')')
  p.fn <- paste0(out_dir, md_var, '_', tag, '_', value)
  if(display_numbers!="none"){
    p.fn <- paste0(p.fn, '.label')
    if(display_numbers=="raw"){
      mat_text <- round(mat,2)
      mat_text <- matrix(as.character(mat_text), ncol=ncol(mat_text))
      rownames(mat_text) <- rownames(mat)
      colnames(mat_text) <- colnames(mat)
      p.fn <- paste0(p.fn, '_raw')
    }
    if(display_numbers=="n"){
      mat_text <- reshape2::dcast(tag.df, mad_comb~md_level, value.var='n')
      mat_text[is.na(mat_text)] <- 0
      rownames(mat_text) <- mat_text[,1]
      mat_text <- mat_text[,-1]
      mat_text <- as.matrix(mat_text)
      p.fn <- paste0(p.fn, '_n')
    }
    mat_text[mat_text=="0"] <- ""
    p <- pheatmap(mat,
                  color = mat_cols,
                  main = main.var,
                  annotation_row = mad_row,
                  annotation_colors = annot_cols,
                  cluster_rows = cluster_mads,
                  cluster_cols = cluster_md,
                  fontsize_col = 12,
                  display_numbers = mat_text,
                  number_color = "black",
                  fontsize_number = 10, 
                  silent = T)
  }else{
    p <- pheatmap(mat,
                  color = mat_cols,
                  main = main.var,
                  annotation_row = mad_row,
                  annotation_colors = annot_cols,
                  cluster_rows = cluster_mads,
                  cluster_cols = cluster_md,
                  fontsize_col = 12,
                  silent = T)
  }
  
  ## Save plot
  if(cluster_mads){
    p.fn <- paste0(p.fn, '.cluster_mads')
  }
  if(cluster_md){
    p.fn <- paste0(p.fn, '.cluster_md')
  }
  p.fn <- paste0(p.fn,'.pdf')
  print(paste0('Saving pheatmap in: ',p.fn))
  save_pheatmap_pdf(x = p,
                    filename = p.fn,
                    width = width_list[[md_var]], height = 14.5)
  return(p)
}

# 12. Main heatmap function (by metadata type)
pheatmap_main <- function(md_type, tag_list, out_dir, tags = c('NotOutlier','Outlier'), values = c('n','prop','pct'), display_numbers_vec = c('none', 'raw', 'n'), logical_vec = c(TRUE,FALSE), width_list){
  print(md_type)
  
  # Subdirectory
  out.sdir <- paste0(out_dir,md_type, '/')
  if(!dir.exists(out.sdir)){dir.create(out.sdir, recursive = T)}
  md_vars <- md_list[[md_type]]
  
  # apply pheatmap to all md_vars from a md_type
  tag_md.list <- tag_list[[md_type]] 
  pheatmap.res <- lapply(md_vars, function(md)
    lapply(tags, function(t)
      lapply(values, function(v)
        lapply(display_numbers_vec, function(d)
          lapply(logical_vec, function(c_mads) 
            lapply(logical_vec, function(c_md) pheatmap.func(md_var = md,
                                                             tag = t, 
                                                             value = v,
                                                             display_numbers = d,
                                                             cluster_mads = c_mads,
                                                             cluster_md = c_md,
                                                             tag_md_list = tag_md.list,
                                                             out_dir = out.sdir,
                                                             width_list = width_list)))))))
  return(pheatmap.res)
}

######################### Functions used in QC_extract_files.R ####################
# 1. Accessory extract function
extract.func <- function(fn, out_dir){
  # File and directories (mymic)
  fn_split <- unlist(str_split(fn, '/'))
  out_temp.idx <- length(fn_split)
  out_fn <- fn_split[out_temp.idx]
  out_sdir <- paste(fn_split[-length(fn_split)],collapse='/')
  out_fdir <- paste0(out_dir,  out_sdir, '/') 
  if(!dir.exists(out_fdir)){dir.create(out_fdir, recursive = T)}
  
  # Copy file
  out_fn.label <- paste0(out_fdir, out_fn)
  print(paste0('Copying file into: ', out_fn.label))
  cp_cmd <- paste0('cp ', fn, ' ', out_fdir)
  system(cp_cmd)
  return(NULL)
}

# 2. Main extract function
extract_files <- function(type, outdirs, files){
  out_dir.type <- outdirs[[type]]
  fn_list.type <- files[[type]]
  print(paste0('Extracting files in: ', out_dir.type))
  res <- lapply(fn_list.type, function(i) extract.func(fn = i, out_dir = out_dir.type))
}

########################### Functions used in QC_heatmaps.R #####################
# 1. Read outputs from QC_statistics.R
read_in <- function(dataset, in_dir = in.dir, qc_label = qc.label, filter_level = opt$level){
  print(dataset)
  in_fn <- paste0(in_dir, dataset, '/', qc_label, '/by_dataset/', filter_level, '/tag.rds')
  in_df <- readRDS(in_fn)
  in_df$dataset <- dataset
  return(in_df)
}

# 2. Heatmap function across datasets (by tags and values)
pheatmap_across_datasets <- function(tag, value, display_numbers = "none", cluster_mads = F, cluster_datasets = T, tag_list, mat_cols = colorRampPalette(brewer.pal(n = 9, name = "Reds"))(100), out_dir, width_var){
  # Convert long to wide fromat (data fram --> matrix)
  tag.df <- tag_list[[tag]]
  tag.mat <- reshape2::dcast(tag.df, mad_comb~dataset, value.var=value)
  tag.mat[is.na(tag.mat)] <- 0
  rownames(tag.mat) <- tag.mat[,1]
  if(ncol(tag.mat)>2){
    tag.mat <- tag.mat[,-1]
    mat <- as.matrix(tag.mat)
  }else{
    rnames <- rownames(tag.mat)
    cnames <- unique(tag.df$dataset)
    mat <- matrix(tag.mat[,-1], dimnames = list(rnames, cnames))
    cluster_datasets <- FALSE
  }
  
  # Pheatmap
  print(tag)
  print(value)
  print(display_numbers)
  print(cluster_mads)
  print(cluster_datasets)
  
  ## Annotation
  ### annotation rows: MAD combinations
  tag.df[,c('nCount_RNA','percent.mt')] <- str_split_fixed(tag.df$mad_comb, ':', 2)[,c(1,2)]
  mad_row <- tag.df[,c('mad_comb','nCount_RNA','percent.mt')]
  mad_row <- unique(mad_row)
  rownames(mad_row) <- mad_row$mad_comb
  mad_row <- mad_row[-1]
  
  ## List of annotation colors
  nCount_RNA.vec <- rev(brewer.pal(9,'Greens'))
  names(nCount_RNA.vec) <- unique(mad_row$nCount_RNA)
  percent.mt.vec <- rev(brewer.pal(9,'Purples'))
  names(percent.mt.vec) <- unique(mad_row$percent.mt)
  annot_cols <- list(nCount_RNA = nCount_RNA.vec,
                     percent.mt = percent.mt.vec)
  
  
  ## Pheatmap
  main.var <- paste0(tag, ' (', value, ')')
  p.fn <- paste0(out_dir, tag, '_', value)
  if(display_numbers!="none"){
    p.fn <- paste0(p.fn, '.label')
    if(display_numbers=="raw"){
      mat_text <- round(mat,2)
      mat_text <- matrix(as.character(mat_text), ncol=ncol(mat_text))
      rownames(mat_text) <- rownames(mat)
      colnames(mat_text) <- colnames(mat)
      p.fn <- paste0(p.fn, '_raw')
    }
    if(display_numbers=="n"){
      mat_text <- reshape2::dcast(tag.df, mad_comb~dataset, value.var='n')
      mat_text[is.na(mat_text)] <- 0
      rownames(mat_text) <- mat_text[,1]
      mat_text <- mat_text[,-1]
      mat_text <- as.matrix(mat_text)
      p.fn <- paste0(p.fn, '_n')
    }
    mat_text[mat_text=="0"] <- ""
    p <- pheatmap(mat,
                  color = mat_cols,
                  main = main.var,
                  annotation_row = mad_row,
                  annotation_colors = annot_cols,
                  cluster_rows = cluster_mads,
                  cluster_cols = cluster_datasets,
                  fontsize_col = 12,
                  display_numbers = mat_text,
                  number_color = "black",
                  fontsize_number = 11,
                  silent = T)
    
  }else{
    p <- pheatmap(mat,
                  color = mat_cols,
                  main = main.var,
                  annotation_row = mad_row,
                  annotation_colors = annot_cols,
                  cluster_rows = cluster_mads,
                  cluster_cols = cluster_datasets,
                  fontsize_col = 12,
                  silent = T)
  }
  
  ## Save plot
  if(cluster_mads){
    p.fn <- paste0(p.fn, '.cluster_mads')
  }
  p.fn <- paste0(p.fn,'.pdf')
  print(paste0('Saving pheatmap in: ',p.fn))
  save_pheatmap_pdf(x = p,
                    filename = p.fn,
                    width = width_var, height = 14.5)
  return(p)
}