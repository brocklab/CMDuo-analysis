plot_cluster_marker_genes <- function(obj_trio_min, 
                                      fusepairs, 
                                      markers_full = h1806_markers, 
                                      n_genes = 5,
                                      include_fusion_markers = FALSE,
                                      include_control_markers = FALSE,
                                      max_cell_per_clone_plot = 20,
                                      minmax_scale = TRUE,
                                      scale_exp_min = -1,
                                      scale_exp_max = 1,
                                      scale_exp_min_sc = -1,
                                      scale_exp_max_sc = 1,
                                      min_counts = 10,
                                      title_size = 6,
                                      color_box_width = 0.2,
                                      clust_genes=FALSE,
                                      xlabel_size = 12,
                                      ylabel_size = 6, 
                                      heatmap_color_low = NA,
                                      heatmap_color_mid = 'grey90',
                                      heatmap_color_hi = 'grey10'
                                      ){
       plots <- list()
  cell_plots <- list()

  all <- markers_full %>% filter(gene %notin% c('EGFP','MCHERRY')) %>% group_by(cluster) %>% arrange(desc(avg_log2FC))
  
  markergenes_df <- data.frame()
  for(i in seq(1, length(unique(markers_full$cluster)))){
    new <- all %>% filter(cluster == LETTERS[i], gene %notin% markergenes_df$gene) %>% head(n_genes) %>% pull(gene)
    df <- data.frame(gene = new, cluster = LETTERS[i])
    markergenes_df <- rbind(markergenes_df, df)
  }

  
  markergenes <- unique(markergenes_df$gene)

  for(i in seq(1,nrow(fusepairs))){
    # subset to only contain parent-parent-fusion set of interest
    sub <- subset_parent_parent_fusion(obj_trio_min, pick_clones = fusepairs[i,], parent_pairs = fusepairs)
    
    # get marker genes for set
    markergenes_sub <- markergenes_df
    pick_genes <- markergenes_sub %>% pull(gene) %>% unique()
    
    # # # rescale data on markergenes for plotting
    #sub@assays$RNA$scale.data <- NULL  #### edit
    #sub <- ScaleData(sub, features = pick_genes)  #### edit
    #sub <- ScaleData(sub, features = rownames(sub))
    
    # get counts to calculate % gene expression
    counts <- FetchData(sub, vars = pick_genes, layer='counts') %>% rownames_to_column('cellid') 
    meta <- sub@meta.data %>% rownames_to_column('cellid') %>% select(cellid, barpairid, cell_class, clone_color)
    
    # filter out cells with low expression of all genes
    pick_gene_counts <- rowSums(counts %>% column_to_rownames('cellid'))
    keep_cells <- names(pick_gene_counts[pick_gene_counts >= min_counts])
    
    # calculate % expression of each gene
    pct_out<- counts %>% 
      left_join(meta, by='cellid') %>% 
      filter(cellid %in% keep_cells) %>% 
      group_by(barpairid) %>% 
      mutate(n_clone=n()) %>% 
      ungroup() %>% 
      pivot_longer(-c(cellid, barpairid, cell_class, clone_color, n_clone), names_to = 'gene', values_to='log_exp') %>% 
      mutate(is_exp = case_when(log_exp > 0 ~ 1, .default = 0)) %>% 
      group_by(barpairid, gene) %>% mutate(pct_exp = round(100*sum(is_exp)/n_clone)) %>% 
      select(barpairid, cell_class, clone_color, gene, pct_exp) %>% 
      distinct()
    
    # get scaled data for plotting
    scaled_out_cells <- FetchData(sub, vars = pick_genes, layer='scale.data') %>% 
      rownames_to_column('cellid')  %>% 
      left_join(meta, by='cellid') %>% 
      filter(cellid %in% keep_cells) %>% 
      group_by(barpairid) %>% mutate(n_clone=n()) %>% 
      ungroup() %>% 
      pivot_longer(-c(cellid, barpairid, cell_class, clone_color, n_clone), names_to = 'gene', values_to='scaled_exp') 
    
    scaled_out <- scaled_out_cells %>% 
      group_by(barpairid, gene) %>% 
      mutate(mean_exp = mean(scaled_exp)) %>% 
      select(barpairid, gene, mean_exp) %>% 
      distinct()
    
    # join scaled expression and pct expressed data
    plot_data <- full_join(pct_out, scaled_out, by=c('barpairid', 'gene'))
    
    # factor and reorder data
    order_barcodes <- plot_data %>% arrange(cell_class, barpairid) %>% pull(barpairid) %>% unique()
    plot_data$barpairid <- factor(plot_data$barpairid, levels = order_barcodes)
    
    
    ## data for making same plot, but showing each cell 
    plot_data_cells <- scaled_out_cells
    
    if(minmax_scale == TRUE){
      # set upper and lower normalized expression values for plotting
      plot_data$scaled_minmax = MinMax(data = plot_data$mean_exp, min = scale_exp_min, max = scale_exp_max)
      plot_data_cells$scaled_minmax = MinMax(data = plot_data_cells$scaled_exp, min = scale_exp_min_sc, max = scale_exp_max_sc)
    }else if(minmax_scale == FALSE){
      # set upper and lower normalized expression values for plotting
      plot_data$scaled_minmax = plot_data$mean_exp
      plot_data_cells$scaled_minmax = plot_data_cells$scaled_exp
    }
    # 
    # 
    # ## data for making same plot, but showing each cell 
    # plot_data_cells <- scaled_out_cells
    # 
    # # set upper and lower normalized expression values for plotting
    #       plot_data$scaled_minmax = MinMax(data = plot_data$mean_exp, min = scale_exp_min, max = scale_exp_max)
    # plot_data_cells$scaled_minmax = MinMax(data = plot_data_cells$scaled_exp, min = scale_exp_min, max = scale_exp_max)
    # 
  if(clust_genes == TRUE){
    ### function to order genes by similarity of expression within all cells
              cluster_genes_within_group <- function(cell_scaled_exp = scaled_out_cells, pick_group = 'A' ){
                group_genes <- markergenes_sub$gene[markergenes_sub$cluster == pick_group]

                mat <- cell_scaled_exp %>%
                  filter(gene %in% group_genes) %>%
                  dplyr::select(cellid, gene, scaled_minmax) %>%
                  pivot_wider(names_from = cellid, values_from = scaled_minmax) %>%
                  column_to_rownames('gene')

                distance_mat <- dist(mat, method = 'euclidean')

                set.seed(123)
                hclust <- hclust(distance_mat, method = 'ward.D2')
                ordered_genes<- rownames(mat)[hclust$order]

                return(ordered_genes)
              }

    ## find best gene order
    gene_order <- c()
    for(group in markergenes_sub$cluster %>% unique()){
      add <- cluster_genes_within_group(plot_data_cells, pick_group = group)
      gene_order <- c(gene_order, add)
    }

    #set gene order
    plot_data$gene <- factor(plot_data$gene, levels = rev(gene_order %>%  unique()))
    plot_data_cells$gene <- factor(plot_data_cells$gene, levels = rev(gene_order %>%  unique()))
    plot_data_cells$gene <- factor(plot_data_cells$gene, levels = rev(gene_order %>%  unique()))
  }else{
    plot_data$gene <- factor(plot_data$gene, levels = rev(markergenes_sub$gene %>% unique()))
    plot_data_cells$gene <- factor(plot_data_cells$gene, levels = rev(markergenes_sub$gene %>% unique()))
  }
    
    
    
    # scale function for control dot size on plot
    scale.func <- switch(EXPR = 'radius', 'size' = scale_size, 'radius' = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
    
    # plotting parameters
    ls = 3
    line_alpha = 0.2
    cols = c("lightgrey", "black")
    n_sections = length(unique(markergenes_sub$cluster))-1
    
    # make plot
    p <- 
      ggplot(plot_data, aes(y = gene, x = barpairid)) +
      geom_hline(yintercept = 0.5+n_genes*seq(1,6), 
                 color = 'black',     size=0.05, alpha=1) +
      theme_classic() + 
      theme(strip.background = element_blank()) +
      geom_point(aes(size=pct_exp, color=scaled_minmax))  +
      scale.func(range = c(0, 3), limits = c(NA, NA)) + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      scale_color_gradient2(low=heatmap_color_low, mid = heatmap_color_mid, high = heatmap_color_hi) + #scale_color_gradient(low = cols[1], high = cols[2]) + 
      labs(color=paste0('Scaled','\n','expression'), size=paste0('% expressed')) + 
      xlab('') + ylab('') + #ggtitle(paste0(fusepairs[i,]$barpairid)) + 
      ggtitle(paste0(fusepairs[i,]$gbarid,' + ',fusepairs[i,]$rbarid)) +
      theme(plot.title = element_text(hjust = 0.5, size = title_size))+ 
      theme(axis.text.y = element_text(size = ylabel_size),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
            legend.title = element_text( size=6), 
            legend.text=element_text(size=6), 
            legend.key.size = unit(0.07, 'in'),
            legend.position = 'none')   

      
      for(k in seq(1, length(unique(markers_full$cluster)))){
        p <- p +
          geom_tile(data = markergenes_sub %>% filter(cluster == LETTERS[k]) %>% select(gene), 
                    aes(x=0.35, y=gene), 
                    fill=clust_colors[k], color=clust_colors[k], width = color_box_width, height = 1, inherit.aes = FALSE) 
      }
    
    plots[[i]] <- p
    
    # function to order plotting of cells by similarity of expression of the chosen marker genes
    cluster_cells_within_clone <- function(cell_scaled_exp, bc_color='green'){
      mat <- cell_scaled_exp %>% 
        filter(clone_color==bc_color) %>% 
        dplyr::select(cellid, gene, scaled_minmax) %>% 
        pivot_wider(names_from = gene, values_from = scaled_minmax) %>% 
        column_to_rownames('cellid')
    
      distance_mat <- dist(mat, method = 'euclidean')
      
      set.seed(123) 
      hclust <- hclust(distance_mat, method = 'ward.D2')
      ordered_cells <- rownames(mat)[hclust$order]
      return(ordered_cells)
    }
    
    grn_order <- cluster_cells_within_clone(plot_data_cells, bc_color = 'green')
    red_order <- cluster_cells_within_clone(plot_data_cells, bc_color = 'red')
    fus_order <- cluster_cells_within_clone(plot_data_cells, bc_color = 'fusion')
    
    plot_data_cells$cellid <- factor(plot_data_cells$cellid, levels = c(grn_order, red_order, fus_order))
    
    #####
    
    #####
   
    # Down sample cells for plotting
      if(length(grn_order) > max_cell_per_clone_plot){
        set.seed(123) 
        gpick <- sample(grn_order, max_cell_per_clone_plot, replace=FALSE)
      } else{
        gpick <- grn_order
      }
      
      if(length(red_order) > max_cell_per_clone_plot){
        set.seed(123) 
        rpick <- sample(red_order, max_cell_per_clone_plot, replace=FALSE)
      } else{
        rpick <- red_order
      } 
      
      if(length(fus_order) > max_cell_per_clone_plot){
        set.seed(123) 
        fpick <- sample(fus_order, max_cell_per_clone_plot, replace=FALSE)
      } else{
        fpick <- fus_order
      }  
    
    plot_data_cells_downsampled <- plot_data_cells[plot_data_cells$cellid %in% c(gpick, rpick, fpick),]
    
    
    ###########################
    # MAKE SC HEATMAPS 
    line_alpha =0.5
    hm <-
      ggplot(plot_data_cells_downsampled, aes(y = gene, x = cellid)) +
      theme(strip.background = element_blank()) +
      geom_tile(aes(fill=scaled_minmax)) +
      scale_fill_gradient2(low=heatmap_color_low, mid = heatmap_color_mid, high = heatmap_color_hi) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      theme_classic() + 
      labs(color=paste0('Scaled','\n','expression'), size=paste0('% expressed')) + 
      xlab('') + 
      ylab('') + 
      ggtitle(paste0(fusepairs[i,]$gbarid,' + ',fusepairs[i,]$rbarid)) + 
      theme(plot.title = element_text(hjust = 0.5, size = title_size),
            axis.text.y = element_text(size = ylabel_size), 
            axis.ticks.x=element_blank(),
            legend.position = 'none',
            axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size=xlabel_size),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=0.75)) +
      geom_vline(size=0.25, 
                 xintercept = c(0.5+length(gpick), 0.5+length(c(gpick,rpick)) )) + #scale_color_gradient(low = NA, high = cols[2]) + geom_hline(yintercept = (c(n_genes + 0.5, (2*n_genes)+0.5)), color = 'black',     size=0.25, alpha=1) +
      geom_hline(yintercept = (n_genes*seq(1,n_sections)+0.5), 
                 color = 'black',     size=0.25, alpha=1) +
      scale_x_discrete(breaks = 
                         c(sort(unique(plot_data_cells_downsampled$cellid[plot_data_cells_downsampled$clone_color=='green']))[length(gpick)/2], 
                           sort(unique(plot_data_cells_downsampled$cellid[plot_data_cells_downsampled$clone_color=='red']))[length(rpick)/2], 
                           sort(unique(plot_data_cells_downsampled$cellid[plot_data_cells_downsampled$clone_color=='fusion']))[length(fpick)/2]), 
                       labels = c(fusepairs[i,]$gbarid, fusepairs[i,]$rbarid, fusepairs[i,]$barpairid))
    
    for(k in seq(1, length(unique(markers_full$cluster)))){
      hm <- hm +
        geom_tile(data = markergenes_sub %>% filter(cluster == LETTERS[k]) %>% select(gene), 
                  aes(x=-0.5, y=gene), 
                  fill=clust_colors[[k]], color=clust_colors[[k]], width = 1, height = 1, inherit.aes = FALSE)
    }
    cell_plots[[i]] <- hm
  }
  
  p_legend <- ggplot(plot_data, aes(y = gene, x = barpairid)) +
    theme_classic() + 
    theme(strip.background = element_blank()) +
    geom_point(aes(size=pct_exp, color=scaled_minmax))  +
    scale.func(range = c(0, 3), limits = c(NA, NA)) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    scale_color_gradient2(low=heatmap_color_low, mid = heatmap_color_mid, high = heatmap_color_hi) + #scale_color_gradient(low = cols[1], high = cols[2]) + 
    labs(color=paste0('Scaled','\n','expression'), size=paste0('% expressed')) + 
    theme(axis.text.y = element_text(size = 6)) + 
    theme(legend.title = element_text( size=6), legend.text=element_text(size=6), legend.key.size = unit(0.09, 'in'))
  
  return(list(plots, cell_plots, p_legend))
}


