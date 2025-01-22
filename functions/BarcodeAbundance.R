########################################################
# plot umaps by condition
## requires ('functions/ggumap.R')
condition_plots <- function(obj, strip.font.size=12){
  emb <- get_embedding(obj)
  emb$condition[emb$condition=='ctrl'] <- 'control'
  emb$condition[emb$condition=='presort'] <- 'pre-sort'
  emb$condition[emb$condition=='fusion'] <- 'fusion-\nenriched'
  emb$condition <- factor(emb$condition,levels=c('initial','pre-sort','control','fusion-\nenriched'))
  
  p1 <- ggumap(emb) + 
    geom_point(aes(color=condition),alpha=0.8,size=0.5, shape=16) + 
    scale_color_manual(values=condition_colors) + 
    labs(color='Condition') + ggtitle('') + theme_void() + 
    theme(aspect.ratio = 0.8,
          legend.title = element_text(size=strip.font.size), 
          legend.text=element_text(size=strip.font.size-2)) 
  
  p2 <- ggumap(emb) + 
    geom_point(aes(color=condition),alpha=0.8,size=0.5, shape=16) + 
    scale_color_manual(values=condition_colors) + 
    labs(color='Condition') + ggtitle('') + theme_void() + 
    theme(aspect.ratio = 0.8) + facet_wrap(~condition) + theme(legend.position = 'none')
  
  p3 <- ggumap(emb %>% filter(replicate !=0)) + 
    geom_point(aes(color=condition),alpha=0.8,size=0.5, shape=16) + 
    scale_color_manual(values=condition_colors) + 
    labs(color='Condition') + ggtitle('') + theme_void() + 
    theme(aspect.ratio = 0.8, strip.text.x = element_text(size = strip.font.size)) + 
    facet_wrap(~replicate, ncol=2) + theme(legend.position = 'none')
  
  p4 <- ggumap(emb %>% filter(replicate !=0)) + 
    geom_point(aes(color=condition),alpha=0.8,size=0.5, shape=16) + 
    scale_color_manual(values=condition_colors) + 
    labs(color='Condition') + ggtitle('') + theme_void() + 
    theme(aspect.ratio = 0.8) + 
    facet_wrap(~condition+replicate, ncol=4) + theme(legend.position = 'none')
  
  return(list(p1, p2, p3, p4))
}

########################################################
# calculate % cells with n barcodes per cell
find_perc_nbars <- function(obj, show_unassigned = FALSE){
  meta <- obj@meta.data %>% dplyr::select(condition, barpairid, replicate, cellLine, n_bar_cell) 
  if(show_unassigned == FALSE){
    meta <- meta %>% filter(n_bar_cell >= 1)
  }
  if(show_unassigned == TRUE){
    meta$n_bar_cell[is.na(meta$n_bar_cell)==TRUE] <- 'none assigned'
  }
  perc_nbars <- meta %>% 
    group_by(condition, replicate, cellLine) %>% 
    mutate(barcount_total = n()) %>% 
    group_by(condition, replicate, cellLine, n_bar_cell) %>% 
    mutate(barcount_each = n()) %>% ungroup() %>% 
    mutate(abundance = barcount_each/barcount_total) %>% 
    arrange(desc(abundance)) %>% 
    dplyr::select(-barpairid) %>% distinct()
  return(perc_nbars)
}

########################################################
# plot % cells with n barcodes per cell
plot_perc_nbars <- function(perc_nbars, font.size=10, strip.font.size=12, legend.box.scale=2){
  perc_nbars$condition[perc_nbars$condition=='ctrl'] <- 'control'
  perc_nbars$condition[perc_nbars$condition=='fusion'] <- 'fusion-\nenriched'
  perc_nbars$condition[perc_nbars$condition=='presort'] <- 'pre-sort'
  
  ggplot(perc_nbars, aes(fill=reorder(as.factor(n_bar_cell),n_bar_cell), y=abundance, x=as.factor(replicate))) + 
    geom_bar(position="stack", stat="identity", color='black', size=0.1) + 
    theme_classic() + 
    theme(strip.background = element_blank(), 
          legend.title = element_text(size=font.size), 
          legend.text=element_text(size=font.size),
          strip.text.x = element_text(size = strip.font.size)) +
    scale_fill_manual(values=c(wes_palette("GrandBudapest2"), rep('grey',20)), name = "bars per cell") + 
    ylab('Proportion of Cells') + xlab('') + 
    facet_wrap(~factor(condition, levels=c('initial','pre-sort','control','fusion-\nenriched')), scales='free_x', ncol=4) +
    guides(fill = guide_legend(override.aes = list(size = legend.box.scale))) 
}

########################################################
# plot % cells with n barcodes per cell but black and white
plot_perc_nbars_bw <- function(perc_nbars, 
                               strip.font.size=12,
                               axis.title.font.size = 12,
                               leg.font.size=10,
                               axis.text.font.size = 12,
                               legend.box.scale=2,
                               add_color_bars = TRUE,
                               show.legend=FALSE
){
  
  perc_nbars$condition[perc_nbars$condition=='ctrl'] <- 'control'
  perc_nbars$condition[perc_nbars$condition=='fusion'] <- 'fusion-\nenriched'
  perc_nbars$condition[perc_nbars$condition=='presort'] <- 'pre-sort'
  
  perc_nbars$n_bar_cell[perc_nbars$n_bar_cell=='none assigned'] <- 0
  
  p <- ggplot(perc_nbars, aes(fill=reorder(as.factor(n_bar_cell),n_bar_cell), y=abundance, x=as.factor(replicate))) + 
    geom_bar(position="stack", stat="identity", color='black', size=0.1) + 
    theme_classic() + 
    theme(strip.background = element_blank(), 
          legend.title = element_text(size=leg.font.size), 
          legend.text=element_text(size=leg.font.size),
          strip.text.x = element_text(size = strip.font.size),
          axis.title = element_text(size = axis.title.font.size),
          axis.text= element_text(size = axis.text.font.size)
    ) +
    scale_fill_manual(values=c('0'='grey95','1'='grey90','2'='#baa99b', '3'='#736254', '4'='#3d3025', condition_colors), #name = "barcodes per cell",
                      breaks = c('0','1','2','3') ) +
    labs(fill=paste0('CMDuo\nbarcodes\nper cell')) +
    ylab('Proportion of Cells') + 
    xlab('Replicate') + 
    facet_wrap(~factor(condition, levels=c('initial','pre-sort','control','fusion-\nenriched')), scales='free_x', ncol=4) +
    guides(fill = guide_legend(override.aes = list(size = legend.box.scale))) 
  
  if(add_color_bars==TRUE){
    p <- p + 
      geom_tile(data = perc_nbars %>% group_by(condition, replicate) %>% summarize(n_clones = n()) %>% ungroup(), 
                aes(x=replicate, y=1.1, fill = condition), 
                width = 1, height = 0.05, inherit.aes = FALSE) +
      scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1) ) 
  }
  if(show.legend==FALSE){
    p <- p + theme(legend.position='none')
  }
  return(p)
}


########################################################
# calculate abundance of each clone by condition
find_clonal_abund <- function(obj, n_bar_cell_min=1){
  abund <- obj@meta.data %>% 
    dplyr::select(condition, barpairid, replicate, cellLine, n_bar_cell) %>% 
    filter(n_bar_cell >= n_bar_cell_min) %>% 
    group_by(condition, replicate, cellLine) %>% 
    mutate(barcount_total = n()) %>% 
    group_by(condition, replicate, barpairid, cellLine) %>% 
    mutate(barcount_each = n()) %>% rowwise() %>% 
    mutate(abundance = barcount_each/barcount_total) %>% ungroup() %>% distinct() %>% arrange(desc(abundance))
  return(abund)
}

########################################################
# plot abundance of each clone by condition
plot_clonal_abund <- function(obj, 
                              abund, 
                              add_clone_counts = TRUE, 
                              add_color_bars = TRUE,
                              clone_count_size = 3,
                              strip.font.size = 14,
                              axis.title.font.size = 14,
                              axis.text.font.size = 12
){
  
  colors_df <- obj@meta.data %>% filter(n_bar_cell >= 1) %>% 
    dplyr::select(barpairid, barpairid_color) %>% 
    distinct()
  color_map <- colors_df$barpairid_color
  names(color_map) <- colors_df$barpairid
  
  abund$condition[abund$condition=='ctrl'] <- 'control'
  abund$condition[abund$condition=='fusion'] <- 'fusion-\nenriched'
  abund$condition[abund$condition=='presort'] <- 'pre-sort'
  
  p <- ggplot(abund, aes(x = as.factor(replicate), 
                         y = abundance)) + 
    geom_bar(position="stack", stat="identity", color='black', size=0.1,
             aes( fill = reorder(as.factor(barpairid),dplyr::desc(abundance)) )) + 
    theme_classic() + 
    theme(legend.position = 'none',
          strip.background = element_blank(), 
          strip.text.x = element_text(size = strip.font.size),
          axis.title = element_text(size = axis.title.font.size),
          axis.text= element_text(size = axis.text.font.size)
    ) +
    scale_fill_manual(values=c(color_map, condition_colors))  +
    ylab('Clonal Abundance') + 
    xlab('Replicate') + 
    facet_wrap(~factor(condition, levels=c('initial','pre-sort','control','fusion-\nenriched')), scales='free_x', ncol=4)
  
  if(add_clone_counts==TRUE){
    p <- p + 
      geom_text(data = abund %>% group_by(condition, replicate) %>% summarize(n_clones = n()) %>% ungroup(),
                position = position_dodge(width=0.9),
                vjust = 0.5,
                size=clone_count_size,
                aes(x=as.factor(replicate), y=1.05, 
                    label=paste0(n_clones))) 
  }
  if(add_color_bars==TRUE){
    p <- p + 
      geom_tile(data = abund %>% group_by(condition, replicate) %>% summarize(n_clones = n()) %>% ungroup(), 
                aes(x=replicate, y=1.15, fill = condition), 
                width = 1, height = 0.05, inherit.aes = FALSE) +
      scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1) ) 
  }
  return(p)
}

########################################################
# find abundance of each clone by cluster
find_clonal_abund_cluster <- function(obj, n_bar_cell_min=1){
  abund <- obj@meta.data %>% 
    dplyr::select(condition, barpairid, clust, cellLine, n_bar_cell) %>% 
    filter(n_bar_cell >= n_bar_cell_min) %>% 
    group_by(condition, clust, cellLine) %>% 
    mutate(barcount_total = n()) %>% 
    group_by(condition, clust, barpairid, cellLine) %>% 
    mutate(barcount_each = n()) %>% rowwise() %>% 
    mutate(abundance = barcount_each/barcount_total) %>% ungroup() %>% distinct() %>% arrange(desc(abundance))
  return(abund)
}

########################################################
# plot abundance of each clone by cluster
plot_clonal_abund_cluster <- function(obj, abund, font.size=10, strip.font.size=12){
  colors_df <- obj@meta.data %>% filter(n_bar_cell >= 1) %>% 
    dplyr::select(barpairid, barpairid_color) %>% 
    distinct()
  color_map <- colors_df$barpairid_color
  names(color_map) <- colors_df$barpairid
  
  abund$condition[abund$condition=='ctrl'] <- 'control'
  abund$condition[abund$condition=='fusion'] <- 'fusion-\nenriched'
  abund$condition[abund$condition=='presort'] <- 'pre-sort'
  
  ggplot(abund, aes(x = as.factor(clust), 
                    y = abundance, 
                    fill = reorder(as.factor(barpairid),dplyr::desc(abundance)))) + 
    geom_bar(position="stack", stat="identity", color='black', size=0.1) + 
    theme_classic() + 
    theme(legend.position = 'none',
          strip.background = element_blank(), 
          legend.title = element_text(size=font.size), 
          legend.text=element_text(size=font.size),
          strip.text.x = element_text(size = strip.font.size)
    ) +
    scale_fill_manual(values=color_map)  +
    ylab('Clonal Abundance') + xlab('') + 
    facet_wrap(~factor(condition, levels=c('initial','pre-sort','control','fusion-\nenriched')), scales='free_x', ncol=4)
}

########################################################
# save plots
fig_path <- "./plots/barcode_abundance/"

SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {png(paste0(fig_path, name, ".", type),width = width, height = height, units = "in", res = res)} 
  else {pdf(paste0(fig_path, name, ".", type),width = width, height = height)}
  print(plots)
  dev.off()
}