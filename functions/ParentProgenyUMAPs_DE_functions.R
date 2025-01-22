##############################################
#### Create seurat object subset containing only fusion trios (parent-parent-fusion) with at least cell_min cells each
subset_fusion_trios_min_cell <- function(dat, parent_conds = c('presort','initial','ctrl'), cell_min=10){
  # Find high confidence fusion clones and related parental barcodes
  fusion_pairs <- dat@meta.data %>% filter(isConfidentFusion==TRUE) %>% select(gbarid, rbarid, barpairid) %>% distinct()
  
  # subset data to only contain high confidence fusion cells and their parents
  dat_fusionset <- subset(dat, subset = 
                            (gbarid %in% fusion_pairs$gbarid & barpairid %notin% fusion_pairs$barpairid & condition %in% parent_conds & fusion != 'fusion' & exp_color!='gfp_mch')|
                            (rbarid %in% fusion_pairs$rbarid & barpairid %notin% fusion_pairs$barpairid & condition %in% parent_conds & fusion != 'fusion' & exp_color!='gfp_mch')|
                            (barpairid %in% fusion_pairs$barpairid & condition=='fusion'  & fusion == 'fusion') )
  
  # classify cells as fusion parents or fusion progeny
  dat_fusionset@meta.data <- dat_fusionset@meta.data %>% 
    mutate(cell_class = case_when(        condition=='fusion' & fusion=='fusion' ~ 'progeny',
                                  condition %in% parent_conds & fusion!='fusion' ~ 'parent',
                                  .default='nd'))
  
  ### find fusion parents and progeny that have at least n cells per clone for the set ('n' parent 1, 'n' parent 2, and 'n' fusion progeny)
  # Identify clones which pass min cells
  min_cell_pass_clones <- dat_fusionset@meta.data %>% 
    dplyr::select(barpairid) %>% 
    pivot_longer(c(barpairid), values_to = 'barid') %>% 
    drop_na() %>% group_by(barid) %>% 
    summarize(cells_per_clone = n()) %>% 
    filter(cells_per_clone >= cell_min) %>% pull(barid) %>% unique()
  
  fusion_pairs_keep <- fusion_pairs %>% 
    filter(gbarid %in% min_cell_pass_clones & rbarid %in% min_cell_pass_clones & barpairid %in% min_cell_pass_clones) 
  
  keep_clones <- fusion_pairs_keep  %>% pivot_longer(everything(), values_to = 'barpairid') %>% pull(barpairid) #%>% pivot_longer(-min_clone_all, values_to = 'barpairid') %>% pull(barpairid) 
  
  # subset to only contain high confidence fusion cells with >= 10 cells per clone in trio
  dat_fusionset_mincells <- subset(dat_fusionset, subset = barpairid %in% keep_clones)
  
  return(dat_fusionset_mincells)
}



##############################################
## function to save figures
fig_path <- "./plots/ParentsAndProgeny/"

SaveFigure <- function(plots, name, type = "png", width, height, res=500){
  if(type == "png") {
    png(paste0(fig_path, name, ".", type),
        width = width, height = height, units = "in", res = res)
  } else {
    pdf(paste0(fig_path, name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}


##############################################
## function to plot fusion trio big on top of all grey cells
parent_progeny_umap_onFull <- function(embedding, pair=fusepairs[i,], legend.title='CMDuo barcode'){
  embedding <- embedding %>% 
    mutate(trace_fusion = case_when(gbarid==pair$gbarid & barpairid==pair$gbarid ~'green_parent', 
                                    rbarid==pair$rbarid & barpairid==pair$rbarid ~'red_parent', 
                                    barpairid==pair$barpairid ~ 'fusion_progeny', .default='x')) 
  
  if('umap_1' %in% colnames(embedding)==TRUE){
    embedding <- embedding %>% dplyr::rename('UMAP_1'='umap_1', 'UMAP_2'='umap_2')
  }
  
  ggplot() +
    geom_point(alpha=0.2, data=embedding %>% filter(trace_fusion == 'x'), 
               aes(x=UMAP_1, y=UMAP_2), color='grey', size=0.1) + 
    geom_point(size=1, alpha=0.9, data=embedding %>% filter(trace_fusion == 'red_parent'| trace_fusion=='green_parent'), 
               aes(x=UMAP_1, y=UMAP_2, color=trace_fusion, shape=trace_fusion)) + 
    geom_point(size=2, alpha=0.8, data=embedding %>% filter(trace_fusion=='fusion_progeny'), 
               aes(x=UMAP_1, y=UMAP_2, color=trace_fusion, shape=trace_fusion)) + 
    theme_void() + 
    theme(aspect.ratio = 0.8) + 
    xlab('UMAP 1') + 
    ylab('UMAP 2') +
    theme(text=element_text(size=14)) + 
    scale_color_manual(values = c('goldenrod','forestgreen','red3','NA'), 
                       breaks = c("fusion_progeny", "green_parent", "red_parent"),
                       labels = c(pair$barpairid, pair$gbarid, pair$rbarid),
                       name=legend.title) + 
    scale_shape_manual(values = c(19,19,19,19), 
                       breaks = c("fusion_progeny", "green_parent", "red_parent"),
                       labels = c(pair$barpairid, pair$gbarid, pair$rbarid),
                       name=legend.title) 
}


##############################################
## function to plot fusion trio big on top of all grey cells but colored by cluster
parent_progeny_umap_onFull_clustcolor <- function(embedding, pair=fusepairs[i,], legend.title='CMDuo barcode'){
  embedding <- embedding %>% 
    mutate(trace_fusion = case_when(gbarid==pair$gbarid & barpairid==pair$gbarid ~'green_parent', 
                                    rbarid==pair$rbarid & barpairid==pair$rbarid ~'red_parent', 
                                    barpairid==pair$barpairid ~ 'fusion_progeny', .default='x')) 
  
  if('umap_1' %in% colnames(embedding)==TRUE){
    embedding <- embedding %>% dplyr::rename('UMAP_1'='umap_1', 'UMAP_2'='umap_2')
  }
  
  ggplot() +
    geom_point(alpha=0.2, data=embedding %>% filter(trace_fusion == 'x'), 
               aes(x=UMAP_1, y=UMAP_2), color='grey', size=0.1) + 
    geom_point(size=1, alpha=0.9, data=embedding %>% filter(trace_fusion == 'red_parent'| trace_fusion=='green_parent'), 
               aes(x=UMAP_1, y=UMAP_2, color=clust, shape=trace_fusion)) + 
    geom_point(size=2, alpha=0.8, data=embedding %>% filter(trace_fusion=='fusion_progeny'), 
               aes(x=UMAP_1, y=UMAP_2, color=clust, shape=trace_fusion)) + 
    theme_void() + 
    theme(aspect.ratio = 0.8) + 
    xlab('UMAP 1') + 
    ylab('UMAP 2') +
    theme(text=element_text(size=14)) + 
    scale_color_manual(values = clust_colors, name=legend.title)  +
    scale_shape_manual(values = c(17,15,16), 
                       breaks = c("fusion_progeny", "green_parent", "red_parent"),
                       labels = c(pair$barpairid, pair$gbarid, pair$rbarid),
                       name=legend.title)
}



##############################################
### function to plot grey UMAP with parents and fusion cells marked for a specific trio
# Here you choose one parent clone and this function will print UMAPs for all fusion clones with that parent
plot_ppf_trios <- function(obj, subset_obj, pick_clone, save_fig=FALSE, legend.pos = 'right',legend.title='CMDuo barcode'){
  fusion_pairs_mincells <- subset_obj@meta.data %>% 
    filter(fusion=='fusion') %>% 
    select(gbarid, rbarid, barpairid) %>% 
    distinct()
  
  plot_pairs <- fusion_pairs_mincells %>% filter(gbarid %in% c(pick_clone) | rbarid %in% c(pick_clone))
  umapcell <- obj@reductions$umap@cell.embeddings %>% data.frame()
  embedding <- AddMetaData(obj, umapcell)@meta.data
  
  for(i in seq(1,length(plot_pairs$barpairid))){ #seq(1,2)){# length(plot_pairs$barpairid)
    pair      <- plot_pairs[i,]
    embedding <- embedding %>% 
      mutate(trace_fusion = case_when(gbarid==pair$gbarid & barpairid==pair$gbarid ~ 'green_parent', 
                                      rbarid==pair$rbarid & barpairid==pair$rbarid ~ 'red_parent', 
                                      barpairid==pair$barpairid ~ 'fusion_progeny', .default='x')) 
    fuse_umap <- parent_progeny_umap_onFull(embedding, pair, legend.title=legend.title) + theme(legend.position=legend.pos)
    print(fuse_umap)
    if(save_fig ==TRUE){
      SaveFigure(fuse_umap, paste0(cellline,"_trace_umap_fusionpair_",pair$barpairid,sep=''), width = 6, height = 4, res=300)
    }
  }
}


##############################################
### function to plot grey UMAP with parents and fusion cells marked for a specific trio
# Here you choose one parent clone and this function will print UMAPs for all fusion clones with that parent
plot_ppf_trios_cluster_color <- function(obj, subset_obj, pick_clone, save_fig=FALSE, legend.pos = 'right',legend.title='CMDuo barcode'){
  fusion_pairs_mincells <- subset_obj@meta.data %>% 
    filter(fusion=='fusion') %>% 
    select(gbarid, rbarid, barpairid) %>% 
    distinct()
  
  plot_pairs <- fusion_pairs_mincells %>% filter(gbarid %in% c(pick_clone) | rbarid %in% c(pick_clone))
  umapcell <- obj@reductions$umap@cell.embeddings %>% data.frame()
  embedding <- AddMetaData(obj, umapcell)@meta.data
  
  for(i in seq(1,length(plot_pairs$barpairid))){ #seq(1,2)){# length(plot_pairs$barpairid)
    pair      <- plot_pairs[i,]
    embedding <- embedding %>% 
      mutate(trace_fusion = case_when(gbarid==pair$gbarid & barpairid==pair$gbarid ~ 'green_parent', 
                                      rbarid==pair$rbarid & barpairid==pair$rbarid ~ 'red_parent', 
                                      barpairid==pair$barpairid ~ 'fusion_progeny', .default='x')) 
    fuse_umap <- parent_progeny_umap_onFull_clustcolor(embedding, pair, legend.title=legend.title) + theme(legend.position=legend.pos)
    print(fuse_umap)
    if(save_fig ==TRUE){
      SaveFigure(fuse_umap, paste0(cellline,"_trace_umap_fusionpair_",pair$barpairid,sep=''), width = 6, height = 4, res=300)
    }
  }
}


##############################################
#### subset full data to include any clones matching the parent of interest or its progeny and progeny-matched parent
subset_on_parent_clone_with_progeny <- function(obj, subset_obj, pick_clone){
  obj_marked <- obj
  
  fusion_pairs_mincells <- subset_obj@meta.data %>% 
    filter(fusion=='fusion') %>% 
    select(gbarid, rbarid, barpairid) %>% 
    distinct()
  
  plot_pairs <- fusion_pairs_mincells %>% filter(gbarid %in% c(pick_clone) | rbarid %in% c(pick_clone))
  
  obj_marked@meta.data <- obj_marked@meta.data  %>% 
    mutate(trace_fusion = case_when(
      gbarid %in% plot_pairs$gbarid & barpairid %in% plot_pairs$gbarid ~ 'green_parent', 
      rbarid %in% plot_pairs$rbarid & barpairid %in% plot_pairs$rbarid ~ 'red_parent', 
      barpairid %in% plot_pairs$barpairid ~ 'fusion_progeny', .default='x'))
  
  parent_clone_subset <- subset(obj_marked, subset = trace_fusion != 'x')
}



##############################################
#### subset full data to contain a specific parent-parent-fusion trio
subset_parent_parent_fusion <- function(obj, pick_clones, parent_pairs){
  obj_marked <- obj
  
  plot_pairs <- parent_pairs %>% filter(barpairid %in% pick_clones$gbarid | 
                                        barpairid %in% pick_clones$rbarid |
                                        barpairid %in% pick_clones$barpairid)
  
  obj_marked@meta.data <- obj_marked@meta.data  %>% 
    mutate(trace_fusion = case_when(
      condition %in% c('initial') & gbarid %in% plot_pairs$gbarid & barpairid %in% plot_pairs$gbarid ~ 'green_parent', 
      condition %in% c('initial') & rbarid %in% plot_pairs$rbarid & barpairid %in% plot_pairs$rbarid ~ 'red_parent', 
      condition %in% c('fusion') & barpairid %in% plot_pairs$barpairid ~ 'fusion_progeny', .default='x'))
  
  ppf_subset <- subset(obj_marked, subset = trace_fusion != 'x')
  
  return(ppf_subset)
}





##############################################
## get top de genes from FindMarkers wilcox test
pull_de_genes <- function(de_out, n_genes, pos_only=FALSE, padj_thresh){
  ordered <- de_out %>% rownames_to_column('gene') %>% arrange(desc(avg_log2FC))
  top <- ordered %>% filter(gene %notin% c('EGFP','MCHERRY'), avg_log2FC > 0, p_val_adj <= padj_thresh) %>% 
    slice_head(n = n_genes)  %>% pull(gene)
  if(pos_only==TRUE){
    return(top)
  }
  if(pos_only==FALSE){
    bot <- ordered %>% filter(gene %notin% c('EGFP','MCHERRY', avg_log2FC < 0, p_val_adj <= padj_thresh)) %>% 
      slice_tail(n = n_genes)  %>% arrange(avg_log2FC) %>% pull(gene)
    return(c(top, bot))
  }
}



##############################################
## get top de genes with clone_color from FindMarkers wilcox test
pull_de_genes_df <- function(de_out, n_genes, pos_only=FALSE, padj_thresh, clone_color = c('fusion')){
  ordered <- de_out %>% rownames_to_column('gene') %>% arrange(desc(avg_log2FC))
  top <- ordered %>% filter(gene %notin% c('EGFP','MCHERRY'), avg_log2FC > 0, p_val_adj <= padj_thresh) %>% 
    slice_head(n = n_genes)  %>% pull(gene)
  tdf <- data.frame(gene = top, clone_color=clone_color[1])
  if(pos_only==TRUE){
    return(tdf)
  }
  if(pos_only==FALSE){
    bot <- ordered %>% filter(gene %notin% c('EGFP','MCHERRY', avg_log2FC < 0, p_val_adj <= padj_thresh)) %>% 
      slice_tail(n = n_genes)  %>% arrange(avg_log2FC) %>% pull(gene)
    bdf <- data.frame(gene = bot, clone_color=clone_color[2])
    df <- rbind(tdf, bdf)
    return(df)
  }
}





##############################################
### heatmap for clone average gene expression

heatmap_genes_clone_average <- function(obj, plotGenes, color_scale, row_clust=TRUE){
  plot_meta <-obj@meta.data %>% select(barpairid, trace_fusion, rbarid) %>% arrange(barpairid, trace_fusion)
  plot_annotation <- plot_meta %>% remove_rownames() %>% distinct() %>% column_to_rownames('barpairid')
  plot_annotation$rbarid[is.na(plot_annotation$rbarid)==TRUE] <- 'x'
  plot_annotation$trace_fusion <- factor(plot_annotation$trace_fusion, levels = c("green_parent", "fusion_progeny", "red_parent"))
  plot_annotation <- plot_annotation %>% arrange(trace_fusion)
  
  exp_scaled <- data.frame(obj@assays$RNA$scale.data, check.names = FALSE)
  exp_scaled_plotgenes <- exp_scaled[row.names(exp_scaled) %in% plotGenes==TRUE,]
  
  barmap <- plot_meta %>% distinct()
  
  mean_exp <- exp_scaled_plotgenes %>% rownames_to_column('gene') %>% 
    pivot_longer(-gene, names_to = 'cellid', values_to='norm_expression') %>% 
    left_join(plot_meta %>% rownames_to_column('cellid')) %>% 
    group_by(barpairid, gene) %>% select(-cellid) %>% 
    summarize(mean_exp = mean(norm_expression)) %>% 
    pivot_wider(names_from = 'barpairid', values_from = 'mean_exp') %>% 
    column_to_rownames('gene') %>% 
    select(barmap$barpairid[barmap$trace_fusion=='green_parent'],
           barmap$barpairid[barmap$trace_fusion=='fusion_progeny'],
           barmap$barpairid[barmap$trace_fusion=='red_parent'])
  
  mycolors <- list(`clone type` = c('green_parent' = 'forestgreen','fusion_progeny' = 'goldenrod', 'red_parent' = 'red3'))
  
  Heatmap(as.matrix((mean_exp)), 
          col=colorRamp2(color_scale, c("magenta", "black", "yellow")),
          cluster_columns = FALSE, 
          cluster_rows = row_clust,
          show_column_names = T,
          column_names_rot = 30,
          show_row_names = T, #clustering_distance_rows = 'spearman',
          row_names_side = "left", 
          column_title = paste(parent_pairs[i,]$barpairid), 
          column_title_gp = gpar(fontsize = 10),
          name = 'scaled expression',
          border = TRUE,
          column_names_gp = gpar(fontsize = 8), 
          row_names_gp = gpar(fontsize = 6), #column_split = plot_annotation$trace_fusion,
          top_annotation = HeatmapAnnotation( `clone type` = plot_annotation$trace_fusion,
                                              col = mycolors,
                                              border = TRUE, show_annotation_name = FALSE))
}





##############################################
### heatmap for clone average gene expression
heatmap_genes_clone_allcell <- function(obj, plotGenes, color_scale = c(-2,0,3), cluster_alg = 'euclidean', split_col=TRUE, col_clust = TRUE,
                                        row_clust=T){
  plot_meta <-obj@meta.data %>% select(barpairid, trace_fusion, rbarid) %>% arrange(barpairid, trace_fusion)
  plot_annotation <- plot_meta %>% remove_rownames() %>% arrange(trace_fusion)
  plot_annotation$rbarid[is.na(plot_annotation$rbarid)==TRUE] <- 'x'
  plot_annotation$trace_fusion <- factor(plot_annotation$trace_fusion, levels = c("green_parent", "fusion_progeny", "red_parent"))
  plot_annotation <- plot_annotation %>% arrange(trace_fusion, barpairid)
  
  exp_scaled <- data.frame(obj@assays$RNA$scale.data, check.names = FALSE)
  exp_scaled_plotgenes <- exp_scaled[row.names(exp_scaled) %in% plotGenes==TRUE,]
  
  barmap <- plot_meta
  
  mean_exp <- exp_scaled_plotgenes %>% 
    select(rownames(barmap)[barmap$trace_fusion=='green_parent'],
           rownames(barmap)[barmap$trace_fusion=='fusion_progeny'], 
           rownames(barmap)[barmap$trace_fusion=='red_parent'])
  
  mycolors <- list(`clone type` = c('green_parent' = 'forestgreen', 'fusion_progeny' = 'goldenrod', 'red_parent' = 'red3'))
  
  if(split_col==FALSE){
    h <- Heatmap(as.matrix((mean_exp)), col=colorRamp2(color_scale, c("magenta", "black", "yellow")),
            cluster_columns = col_clust, 
            cluster_rows = row_clust,
            show_column_names = FALSE,
            show_row_names = TRUE, 
            clustering_distance_columns = cluster_alg,
            row_names_side = "left", 
            column_title = paste(parent_pairs[i,]$barpairid),
            border = TRUE,
            name = 'scaled expression',
            column_title_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 6), 
            row_names_gp = gpar(fontsize = 5), #column_split = plot_annotation$trace_fusion,
            top_annotation = HeatmapAnnotation( `clone type` = plot_annotation$trace_fusion,
                                                col = mycolors,
                                                border = TRUE, show_annotation_name = TRUE))
  }
  if(split_col==TRUE){
    h <- Heatmap(as.matrix((mean_exp)), col=colorRamp2(color_scale, c("magenta", "black", "yellow")),
            cluster_columns = col_clust, 
            cluster_rows = row_clust,
            show_column_names = FALSE,
            show_row_names = TRUE, 
            clustering_distance_columns = cluster_alg,
            row_names_side = "left", 
            column_title = paste(parent_pairs[i,]$barpairid),
            border = TRUE,
            name = 'scaled expression',
            column_title_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 6), 
            row_names_gp = gpar(fontsize = 5), 
            column_split = plot_annotation$trace_fusion,
            top_annotation = HeatmapAnnotation( `clone type` = plot_annotation$trace_fusion,
                                                col = mycolors,
                                                border = TRUE, show_annotation_name = TRUE))
  }
  return(h)
}

