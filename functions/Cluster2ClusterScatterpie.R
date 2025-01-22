###############################################################
### save figure function
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {png(paste0(fig_path, name, ".", type),width = width, height = height, units = "in", res = res)} 
  else {pdf(paste0(fig_path, name, ".", type),width = width, height = height)}
  print(plots)
  dev.off()
}


###############################################################
### subset on fusion trios with at least `cell_min` cells each
fusionset_mincell_subset <- function(obj, fusion_pairs, cell_min=10){
  # classify cells as fusion parents or fusion progeny
  obj@meta.data <- obj@meta.data %>% mutate(cell_class = case_when(fusion=='fusion' ~ 'progeny', .default='parent'))
  
  # find fusion parents and progeny that have at least n cells per clone for the set ('n' parent 1, 'n' parent 2, and 'n' fusion progeny)
  cell_min = 10
  
  # Identify clones which pass min cells
  min_cell_pass_clones <- dat_fusionset@meta.data %>% 
    dplyr::select(barpairid) %>% 
    pivot_longer(c(barpairid), values_to = 'barid') %>% 
    drop_na() %>% group_by(barid) %>% 
    summarize(cells_per_clone = n()) %>% 
    filter(cells_per_clone >= cell_min) %>% 
    pull(barid) %>% 
    unique()
  
  fusion_pairs_keep <- fusion_pairs %>% 
    mutate(min_clone_all = case_when( 
      (gbarid %in% min_cell_pass_clones & rbarid %in% min_cell_pass_clones & barpairid %in% min_cell_pass_clones) ~ TRUE, 
      .default=FALSE)) %>% 
    filter(min_clone_all==TRUE) 
  
  keep_clones <- fusion_pairs_keep %>% pivot_longer(-min_clone_all, values_to = 'barpairid') %>% pull(barpairid) 
  
  # subset to only contain high confidence fusion cells with >= 10 cells per clone in trio
  obj_fusionset_mincells <- subset(dat_fusionset, subset = barpairid %in% keep_clones)
  
  return(obj_fusionset_mincells)
}

###############################################################
### format data for Alluvial plot
format_flow_data <- function(assigned_metaclusters, format='alluvial'){
  #assigned_metaclusters <- test_assign
  greens <- assigned_metaclusters %>% 
    filter(clone_color=='green') %>% 
    dplyr::select(-clone_color) %>% 
    dplyr::rename('gbarid'='barpairid', 'g_parent'='clone_metacluster')
  
  reds <- assigned_metaclusters %>% 
    filter(clone_color=='red') %>% 
    dplyr::select(-clone_color) %>% 
    dplyr::rename('rbarid'='barpairid', 'r_parent'='clone_metacluster')
  
  fusions <- assigned_metaclusters %>% 
    filter(clone_color=='fusion') %>% 
    select(-clone_color) %>% 
    dplyr::rename('fusion'='clone_metacluster')
  
  map_set <- dat_fusionset_mincells@meta.data %>% 
    dplyr::select(gbarid, rbarid, barpairid) %>% 
    distinct() %>% 
    filter(barpairid %in% fusions$barpairid)
  
  map_set_annotated <- left_join(left_join(left_join(map_set, greens, by='gbarid', relationship = "many-to-many"), 
                                           reds, by='rbarid', relationship = "many-to-many"), 
                                 fusions, by='barpairid', relationship = "many-to-many")
  
  map_set_alluvial <-  map_set_annotated %>% 
    mutate(`green parent cluster` = g_parent) %>% 
    mutate(`red parent cluster` = r_parent)
  
  map_set_alluvial <- map_set_alluvial %>% dplyr::rename('fusion cluster'='fusion')
  
  if(format=='alluvial'){
    map_set <- map_set_alluvial
  }
  if(format=='sankey'){
    map_set <- map_set_alluvial %>% 
      drop_na() %>% 
      make_long(`green parent cluster`, `fusion cluster`, `red parent cluster`)
  }
  return(map_set)
}

###############################################################
#### alluvial plot
plot_alluvial <- function(map_set_sankey_wide, color_map=hue_pal()(10)){
  
  alluvial_wide( select(map_set_sankey_wide, `green parent cluster`, `red parent cluster`,`fusion cluster`), 
                 fill_by = 'last_variable', 
                 order_levels = LETTERS[1:7],
                 col_vector_flow = color_map,
                 col_vector_value = color_map,
                 stratum_width = 1/4) +
    theme_void() + 
    theme(axis.text.x= element_text())
}


###############################################################
### find percent of each clone assigned to each metacluster (clust)
per_clone_cluster <- function(dat_obj, min_pct_in_clust = 0){
  bar_clust_count_all <- dat_obj@meta.data %>% 
    select(n_bar_cell, condition, barpairid, clust, cell_class, clone_color) %>%  
    group_by(cell_class, barpairid, clust, clone_color) %>% 
    summarize(n_bar_in_cluster=n()) %>% 
    ungroup()
  
  # count total number of instances of a clone
  total_bar_count <- bar_clust_count_all %>% 
    group_by(barpairid) %>% 
    summarize(total_bar_count = sum(n_bar_in_cluster)) %>% 
    ungroup()
  
  # calculate percentage of each clone in each cluster
  pct_bar_clust <- left_join(bar_clust_count_all[bar_clust_count_all$barpairid %in% total_bar_count$barpairid,], total_bar_count, by='barpairid') %>% 
    mutate(pct_bar_in_clust = n_bar_in_cluster/total_bar_count) %>%
    filter(pct_bar_in_clust >= min_pct_in_clust)
  
  return(pct_bar_clust)
}


###############################################################
### assign a row for each clone metacluster
assign_metacluster_row_pct <- function(pct_bar_clust, scale_factor = 100){
  scaled_count <- pct_bar_clust %>% 
    dplyr::select(barpairid, clust, pct_bar_in_clust, clone_color) %>% 
    mutate(scaled_count = round(scale_factor*pct_bar_in_clust)) 
  
  uniq_bars <- unique(scaled_count$barpairid)
  
  row_pct_clone <- data.frame()
  
  for(i in seq(1,length(uniq_bars))){
    sub <- scaled_count %>% filter(barpairid == uniq_bars[i])
    
    for(clust in unique(sub$clust)){
      if(sub$scaled_count[sub$clust == clust] > 0){
        new <- data.frame(clone_color = sub$clone_color[1], 
                          barpairid = uniq_bars[i] , 
                          clone_metacluster = rep(clust, sub$scaled_count[sub$clust == clust]))
        row_pct_clone <- rbind(row_pct_clone, new)
      }
    }
  }
  row_pct_clone$clone_metacluster <- as.factor(row_pct_clone$clone_metacluster)
  
  return(row_pct_clone)
}


##################################################################################
## scatter pie function
make_scatterpie <- function(dat, pct_bar_clust_allgroups, clust_df, pick_sim_clust, flip_coord=FALSE){
  
  # plotting parameters
  rad = 1.2
  pospie = 4*rad

  ####  order by most abundant cluster
  # find most represented cluster in this subset
  max_clust <- pct_bar_clust_allgroups %>% 
    filter(pct_clust %in% c(pick_sim_clust)) %>% 
    filter(clone_color=='fusion') %>% group_by(clust) %>%
    summarize(mean_clust_pct = mean(pct_bar_in_clust)) %>% 
    filter(mean_clust_pct == max(mean_clust_pct)) %>% pull(clust)
  
  if(flip_coord==FALSE){
    # order piecharts based on descending abundance in that cluster
    fuse_meta <- pct_bar_clust_allgroups %>% ungroup() %>% 
      filter(pct_clust %in% c(pick_sim_clust), clone_color=='fusion', clust == max_clust) %>%
      arrange(desc(pct_bar_in_clust)) %>% 
      dplyr::select(groupid) %>% distinct()
  }
  if(flip_coord==TRUE){
    # order piecharts based on ascending abundance in that cluster
    fuse_meta <- pct_bar_clust_allgroups %>% ungroup() %>% 
      filter(pct_clust %in% c(pick_sim_clust), clone_color=='fusion', clust == max_clust) %>%
      arrange(pct_bar_in_clust) %>%
      dplyr::select(groupid) %>% distinct()
  }

  # and make sure all clones are accounted for
  fuse_extra <- pct_bar_clust_allgroups %>% filter(pct_clust %in% c(pick_sim_clust), clone_color=='fusion') %>%
    filter(groupid %notin% fuse_meta$groupid) %>% 
    arrange(desc(pct_bar_in_clust)) %>% 
    dplyr::select(groupid) %>% distinct()
  
  if(flip_coord==FALSE){
    fuse_meta <- rbind(fuse_meta, fuse_extra) %>% ungroup() %>% mutate(groupid_num=seq(1,3*n(), by=3))
  }
  
  if(flip_coord==TRUE){
    fuse_meta <- rbind(fuse_extra, fuse_meta) %>% ungroup() %>% mutate(groupid_num=seq(1,3*n(), by=3))
  }
  
  ######################  ######################  ######################
  
  
  fuse_order <- fuse_meta$groupid
  
  pct_bar_clust_allgroups <- left_join(pct_bar_clust_allgroups %>% ungroup() %>% filter(pct_clust %in% c(pick_sim_clust)), 
                                       fuse_meta
                                       )
  
  pct_bar_clust_allgroups$groupid <- factor(pct_bar_clust_allgroups$groupid, levels = fuse_order)
  
  pct_bar_clust_wide <- pct_bar_clust_allgroups %>% 
    pivot_wider(names_from = clust, values_from = pct_bar_in_clust) %>% 
    mutate(across(everything(), ~ replace_na(.x, 0))) 
  
  clust_cols <- colnames(pct_bar_clust_wide[colnames(pct_bar_clust_wide) %in% names(clust_colors)])
 
  map_cols <- clust_colors[names(clust_colors) %in% clust_cols]
  
  pct_bar_clust_wide <- pct_bar_clust_wide %>% ungroup() %>% 
    mutate(x.pos = groupid_num,
           y.pos = case_when(clone_color == 'green' ~ pospie, 
                             clone_color == 'red' ~ 0, 
                             clone_color =='fusion' ~ -pospie)
           ) %>% 
    mutate(n_cells_trios = sum(n_bar_in_cluster)) %>% 
    group_by(barpairid) %>% 
    mutate(trio_size = total_bar_count/n_cells_trios) %>% 
    ungroup()
  
  barmapping <- pct_bar_clust_wide %>% 
    dplyr::select(clone_color, pct_clust, barpairid, x.pos, y.pos) %>% 
    distinct() %>% 
    ungroup()
  
  pct_bar_clust_wide <- pct_bar_clust_wide %>% 
    mutate(class_color = case_when(clone_color == 'green'     ~ 'green parent',
                                   clone_color == 'red' ~ 'red parent',
                                   clone_color == 'fusion'  ~ 'fusion')
           )
  
  if(flip_coord==FALSE){
  p <- ggplot() + 
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin =  2*rad+0.6, ymax =  2*rad-0.6+pospie), fill = "forestgreen", alpha = 0.1) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -2*rad+0.6, ymax =  2*rad-0.6), fill = "red3",        alpha = 0.1) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -2*rad+0.6-pospie      , ymax = -2*rad-0.6), fill = "goldenrod",   alpha = 0.1) +
    geom_segment(aes(x = seq(1,max(pct_bar_clust_wide$groupid_num),by=3),         # v plus
                     xend = seq(1,max(pct_bar_clust_wide$groupid_num),by=3),
                     y = 2*rad+0.4,
                     yend = 2*rad-0.4), size=0.2)+
    geom_segment(aes(x = seq(1,max(pct_bar_clust_wide$groupid_num),by=3)-0.4,   # h plus
                     xend = seq(1,max(pct_bar_clust_wide$groupid_num),by=3)+0.4,
                     y = 2*rad, 
                     yend = 2*rad), size=0.2) +
    geom_segment(arrow = arrow(type = "closed", length = unit(0.025, "inches")),  # arrow
                 aes(x = seq(1,max(pct_bar_clust_wide$groupid_num),by=3),
                     xend = seq(1,max(pct_bar_clust_wide$groupid_num),by=3),
                     y = -rad-0.3,
                     yend = -pospie+rad+0.3), size=0.2) +
    geom_scatterpie(aes(x=x.pos, y=y.pos, group=barpairid, r=rad),
                    data=pct_bar_clust_wide,
                    size=0.2,
                    cols=clust_cols,
                    legend_name = 'cluster') + 
    scale_fill_manual(values=clust_colors) +
    theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #theme(legend.position = 'none') +
    coord_equal() +
    scale_x_continuous(breaks = pct_bar_clust_wide$groupid_num, labels=pct_bar_clust_wide$groupid) +
    scale_y_continuous(breaks = pct_bar_clust_wide$y.pos, labels=pct_bar_clust_wide$class_color) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust=1, size = 8),
          axis.text.y = element_text(size = 9),
          plot.margin = margin(10, 10, 10, 10)) +
    xlab('') + ylab('')
  }
  
  if(flip_coord==TRUE){
    p <- ggplot() + 
      geom_rect(aes(ymin = -Inf, ymax = Inf, xmin =  2*rad+0.6, xmax =  2*rad-0.6+pospie), fill = "goldenrod", alpha = 0.1) +
      geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = -2*rad+0.6, xmax =  2*rad-0.6), fill = "red3",        alpha = 0.1) +
      geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = -2*rad+0.6-pospie      , xmax = -2*rad-0.6), fill = "forestgreen",   alpha = 0.1) +
      geom_segment(aes(y = seq(1,max(pct_bar_clust_wide$groupid_num),by=3),         # v plus
                       yend = seq(1,max(pct_bar_clust_wide$groupid_num),by=3),
                       x = -2*rad+0.4,
                       xend = -2*rad-0.4), size=0.2) +
      geom_segment(aes(y = seq(1,max(pct_bar_clust_wide$groupid_num),by=3)-0.4,   # h plus
                       yend = seq(1,max(pct_bar_clust_wide$groupid_num),by=3)+0.4,
                       x = -2*rad, 
                       xend = -2*rad), size=0.2) +
      geom_segment(arrow = arrow(type = "closed", length = unit(0.025, "inches")),  # arrow
                   aes(   y = seq(1,max(pct_bar_clust_wide$groupid_num),by=3),
                          yend = seq(1,max(pct_bar_clust_wide$groupid_num),by=3),
                          xend = pospie-rad-0.3,
                          x = rad+0.3), size=0.2) +
      geom_scatterpie(aes(y=x.pos, x=-y.pos, group=barpairid, r=rad),
                      data=pct_bar_clust_wide,
                      size=0.1,
                      cols=clust_cols,
                      legend_name = 'cluster') + 
      scale_fill_manual(values=clust_colors) +
      theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #theme(legend.position = 'none') +
      coord_equal() +
      scale_y_continuous(breaks = pct_bar_clust_wide$groupid_num, labels=pct_bar_clust_wide$groupid) +
      scale_x_continuous(breaks = pct_bar_clust_wide$y.pos, labels=pct_bar_clust_wide$class_color) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust=1, size = 8),
            axis.text.y = element_text(size = 12),
            plot.margin = margin(0, 10, 0, 0)) + # top right bottom left
      xlab('') + ylab('')
  }
  
  return(list(p, pct_bar_clust_wide))
}


