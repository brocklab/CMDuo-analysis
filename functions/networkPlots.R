##### CREATE GGPLOT COMPATIBLE NETWORK PLOT LAYOUT USING IGRAPH FUNCTIONS
layout_network <- function(cellline='mb231', 
                           obj = mb231, 
                           min_cell_per_fusion_clone = 10, 
                           min_cell_per_clone = 1, 
                           layout_type='nicely',
                           plot_scaling = 10){
  ## MAKE NETWORK OBJECT
  # add color information to metadata
  obj@meta.data <- obj@meta.data %>% 
    mutate(color=case_when(is.na(gbarid)==TRUE  & is.na(rbarid)==FALSE ~ 'mcherry', 
                           is.na(gbarid)==FALSE & is.na(rbarid)==TRUE  ~ 'gfp',
                           is.na(gbarid)==FALSE & is.na(rbarid)==FALSE ~ 'fusion',
                           .default = 'other')) %>% 
    mutate(replicate = as.numeric(replicate))
  
  # count number of cells per fusion clone per replicate
  # find barcode groups which are present in the minimum number of cells required within each replicate
  df <- obj@meta.data %>% filter(condition == 'fusion', color == 'fusion') %>% 
    rownames_to_column('cellid') %>% 
    group_by(barpairid, replicate) %>% mutate(f_cell_count = n()) %>% 
    filter(f_cell_count >= min_cell_per_fusion_clone) %>% 
    dplyr::select(replicate, gbarid, rbarid, barpairid, color, cellid) %>% distinct() %>% arrange(barpairid)
  
  # single cell ids for all fusion cells which passed minimum clone number filtering
  keeper_fusions <- df$cellid
  
  # some fusion barpairs have two green or two red barcodes, need to represent this
  df_quant <- df %>% dplyr::select(-cellid) %>% dplyr::select(gbarid, rbarid, barpairid) %>% distinct() %>% 
    mutate(n_gfp_bars = 1+str_count(gbarid, "_"), n_mch_bars = 1+str_count(rbarid, "_"))
  rownames(df_quant) <- NULL
  
  # identify fusion clones with only one red and one green barcode
  nodes_single <- df_quant %>% filter(n_gfp_bars <=1 & n_mch_bars <= 1) %>% dplyr::select(-n_gfp_bars, -n_mch_bars) %>% 
    dplyr::rename('from'='gbarid', 'to'='rbarid') %>% distinct()
  
  # identify fusion clones with more than one red or green barcode
  nodes_multi <-  df_quant %>% filter(n_gfp_bars > 1 | n_mch_bars > 1)
  
  # if there are not any fusion clones with multiple red or green barcodes, move on
  if(dim(nodes_multi)[1] == 0){
    edges_all <- nodes_single
  }
  
  # if there are fusion clones with more than 1 green or red barcode (e.g. gNNN_gMMM_rXXX or gXXX_rNNN_rMMM)
  # then split the multi-barcoded parent node into two nodes                                        
  ####### TO-DO: build in a way to test if just a double barcoded cell from viral integration
  if(dim(nodes_multi)[1] != 0){
    multimaps <- data.frame()
    for(i in seq(1, nrow(nodes_multi))){
      # split nodes for cells with two barcodes in green or red and pull all nodes
      m_nodes <- nodes_multi[i,] %>% 
        separate(gbarid, into = c("g1", "g2"), sep = "_") %>% 
        separate(rbarid, into = c("r1", "r2"), sep = "_") %>% 
        dplyr::select(-n_gfp_bars, -n_mch_bars) %>% 
        pivot_longer(-c(barpairid, replicate), values_to = 'nodes') %>% 
        dplyr::select(-name) %>% drop_na() %>% pull(nodes)
      
      ## for each barcode set, create a line representing a network pairing between each node within the set
      # find all node combinations within a set of barcodes
      bar_combos <- as.data.frame(t(combn(m_nodes, 2)), stringsAsFactors = FALSE)
      colnames(bar_combos) <- c("from", "to")
      bar_combos <- bar_combos %>% mutate(barpairid = nodes_multi[i,]$barpairid, replicate = nodes_multi[i,]$replicate)
      multimaps <- rbind(multimaps, bar_combos)
      multimaps <- multimaps %>% dplyr::select(replicate, from, to, barpairid)
    }
    # add to data frame with fusion cells with single-bardoded parents
    edges_all <- rbind(nodes_single, multimaps)
  }
  
  
  # count clones 
  clone_counts <- obj@meta.data %>% 
    rownames_to_column('cellid') %>% 
    filter( (condition %in% c('fusion') & cellid %in% keeper_fusions) ) %>% 
    group_by(barpairid, replicate) %>% 
    summarize(cell_count = n())
  
  # filter non-fusion clones which do not meet them minimum count (`min_cell_per_clone`)
  barcellmin_all <- clone_counts %>% 
    mutate(n_bars = 1+str_count(barpairid, "_")) %>% 
    filter( (n_bars >= 2) | (n_bars >= 2 & barpairid %in% unique(multimaps$barpairid)) | (cell_count >= min_cell_per_clone) ) %>% 
    dplyr::select(barpairid, replicate)
  
  # find meta data for cells which passed filtering
  dfx <- left_join(barcellmin_all, obj@meta.data, by=c('barpairid','replicate')) %>% 
    dplyr::select(replicate, gbarid, rbarid, barpairid, color) %>% 
    distinct()
  
  # add meta data to edges data frame
  edges_all_meta <- left_join(dfx %>% dplyr::select(barpairid, replicate), edges_all, by=c('barpairid','replicate'))
  
  # count number of replicates a barcode set appears in
  nreps <- dfx %>% dplyr::select(barpairid, replicate, color) %>% 
    group_by(barpairid) %>% filter(replicate %in% c(1,2,3,4)) %>% distinct() %>%  
    mutate(n_reps_pair = n()) %>% ungroup() %>% 
    dplyr::select(barpairid, n_reps_pair, color) %>% distinct()
  
  # add replicate occurence to edges data frame
  edges_all_meta <- left_join(edges_all_meta, nreps, by='barpairid') %>% drop_na()
  
  ## use igraph functions to create point layouts
  # generate igraph object from edges data frame
  g <- graph_from_data_frame(edges_all_meta %>% dplyr::select(from, to, barpairid) %>% drop_na() %>% distinct() %>% ungroup())
  
  # generate layout  
  set.seed(1212) # set seed for consistency
  # can layout with different igraph network algorithms if desired. kk and fr tend to look the best 
  if(layout_type == 'kk'){
    layout <- layout_with_kk(g)
  }
  if(layout_type == 'fr'){
    layout <- layout_with_fr(g)
  }
  if(layout_type == 'tree'){
    layout <- layout_as_tree(g)
  }
  if(layout_type == 'nicely'){
    layout <- layout_nicely(g)
  }
  if(layout_type == 'lgl'){
    layout <- layout_with_lgl(g)
  }
  
  # format node position information and scale the node positions to exist between 0 and 10 for plotting consistency
  layoutpoints <- data.frame(x.coord = rescale(layout[, 1], to=c(0,plot_scaling)),
                             y.coord = rescale(layout[, 2], to=c(0,plot_scaling)),
                             vertex.names = V(g)$name) %>% 
    mutate(barcolor=case_when(startsWith(vertex.names,'g')==TRUE ~ 'gfp', 
                              startsWith(vertex.names,'r')==TRUE ~ 'mcherry'))
  
  
  # make edges data from long and add node positions to edges data frame
  edge_meta_long_xy <- edges_all_meta %>% 
    pivot_longer(c(from,to), names_to = 'start', values_to = 'vertex.names') %>% 
    dplyr::select(-start) %>% 
    drop_na() %>% 
    distinct() %>% 
    left_join(layoutpoints, by="vertex.names") %>% 
    rowwise() %>% 
    mutate(barnum = list(as.numeric(str_extract_all(vertex.names, "\\d+")[[1]]))) %>% 
    mutate(barnum = paste(unlist(barnum), collapse='_'))
  
  # calculate abundance of each fusion clone in the population
  prop_barpair <- left_join(barcellmin_all, obj@meta.data, by=c('barpairid', 'replicate')) %>% 
    filter(condition=='fusion') %>% 
    group_by(replicate) %>% mutate(ncell_sample = n()) %>% 
    group_by(replicate, barpairid) %>% mutate(n_barpair = n()) %>% ungroup() %>% 
    rowwise() %>% mutate(prop_bar = n_barpair/ncell_sample) %>% 
    dplyr::select(replicate, barpairid, prop_bar) %>% distinct() %>% arrange(desc(prop_bar))
  
  # add abundance information meta data to edges data frame
  edge_meta_long_xy <- left_join(edge_meta_long_xy, prop_barpair, by=c('barpairid', 'replicate'))
  
  
  
  ###### ADD ALL OTHER CLONES NOT INVOLVED IN FUSION
  
  # to be able to plot all clones in the population regardless of whether they participated in fusion, 
  # we will identify all single barcoded clones in the early populations and map an edges to themselves
  all_n1_clones <- obj@meta.data %>% 
    filter(condition %in% c('presort','initial')) %>% 
    dplyr::select(n_bar_cell, barpairid) %>% 
    filter(n_bar_cell == 1) %>% filter(barpairid %notin% unique(c(df_quant$rbarid,df_quant$gbarid, df_quant$barpairid))) %>% 
    group_by(barpairid) %>% summarize(n_cells = n()) %>% 
    filter(n_cells >= min_cell_per_clone) %>% 
    pull(barpairid) %>% unique() %>% sort()
  
  ## map each node to itself in every replicate so that the background plot consistently shows all potential clones
  other_edges <-  rbind(data.frame(replicate = 0, from = all_n1_clones, to = all_n1_clones, barpairid = all_n1_clones))
  
  
  # add replicate occurence to edges data frame
  other_edges_meta <- other_edges %>% mutate(n_reps_pair = 4) %>% drop_na()
  
  ## use igraph functions to create point layouts
  # generate igraph object from edges data frame
  g_other <- graph_from_data_frame(other_edges_meta %>% dplyr::select(from, to, barpairid) %>% drop_na() %>% distinct() %>% ungroup())
  
  # generate layout  
  layout_other <- layout_nicely(g_other)
  
  if(layout_type == 'kk' | layout_type == 'fr'){
    layoutpoints_other <- data.frame(x.coord = rescale(layout_other[, 1], to=c(-2,plot_scaling+0.5)),
                                     y.coord = rescale(layout_other[, 2], to=c(-2,plot_scaling+0.5)),
                                     vertex.names = V(g_other)$name) %>%
      mutate(barcolor=case_when(startsWith(vertex.names,'g')==TRUE ~ 'gfp',
                                startsWith(vertex.names,'r')==TRUE ~ 'mcherry'))
  }
  else{
    layoutpoints_other <- data.frame(x.coord = rescale(layout_other[, 1], to=c(-0.5,plot_scaling+0.5)),
                                     y.coord = rescale(layout_other[, 2], to=c(-0.5,plot_scaling+0.5)),
                                     vertex.names = V(g_other)$name) %>%
      mutate(barcolor=case_when(startsWith(vertex.names,'g')==TRUE ~ 'gfp',
                                startsWith(vertex.names,'r')==TRUE ~ 'mcherry'))
  }

  # make edges data from long and add node positions to edges data frame
  other_edges_meta_long_xy <- other_edges_meta %>% 
    pivot_longer(c(from,to), names_to = 'start', values_to = 'vertex.names') %>% 
    dplyr::select(-start) %>% 
    drop_na() %>% 
    distinct() %>% 
    left_join(layoutpoints_other, by="vertex.names") %>% 
    rowwise() %>% 
    mutate(barnum = list(as.numeric(str_extract_all(vertex.names, "\\d+")[[1]]))) %>% 
    mutate(barnum = paste(unlist(barnum), collapse='_'))
  
  
  # Add network data for all clones to network data for fusion clones
  plot_all <- rbind(edge_meta_long_xy, other_edges_meta_long_xy %>% mutate(prop_bar = 0, color=barcolor) %>% dplyr::select(colnames(edge_meta_long_xy))) %>% drop_na()
  
  return(list(plot_all, edges_all_meta, df_quant, nodes_multi))
}



############################
############################
############################
###########################


##### NETWORK PLOT
plot_network <- function(edge_meta_long_xy, which_reps = c(1,2,3,4), nchar_bar = 1, show_labels = TRUE){
  text_nudge = 0.4
  text_size = 3
  axmin = -1
  axmax = 11
  
  p <- ggplot(edge_meta_long_xy %>% filter(replicate %in% which_reps) %>% filter(nchar(barpairid) > nchar_bar), 
              aes(x=x.coord, y=y.coord)) +
    theme_void() +
    theme(legend.position='none', aspect.ratio = 1) +
    scale_color_manual(values=c('fusion'='goldenrod','gfp'='forestgreen','mcherry'='red3')) + 
    ylim(axmin,axmax) + xlim(axmin,axmax)
  
  if(show_labels==TRUE){
    p <- p +
      geom_shape(aes(group=barpairid), color='goldenrod', alpha=0.8, size=2, fill=NA) + 
      geom_text(data= net  %>% filter(nchar(barpairid) > 5) %>% 
                  dplyr::select(barnum, barcolor, x.coord, y.coord) %>% distinct(), 
                aes(x=x.coord, y=y.coord, label = barnum, color=barcolor), color='white', size = text_size, fontface ='bold', nudge_x = 0, nudge_y = 0) + #geom_text(aes(label = barnum), color='white', size = text_size, fontface ='bold', nudge_x = 0, nudge_y = 0) +
      geom_text(data= net  %>% filter(nchar(barpairid) > 5) %>% 
                  dplyr::select(barnum, barcolor, x.coord, y.coord) %>% distinct(), 
                aes(x=x.coord, y=y.coord, label = barnum, color=barcolor), 
                size = text_size, fontface ='bold', nudge_x = 0, nudge_y = 0) #geom_text(aes(label = barnum, color=barcolor), size = text_size, fontface='bold', nudge_x = 0, nudge_y = 0) 
    return(p)
  }
  
  if(show_labels==FALSE){
    p <- p + geom_shape(aes(group=barpairid), color='goldenrod', alpha=0.8, size=2, fill=NA) + #, size=as.factor(n_reps_pair)
      geom_point(aes(color= barcolor), size=2)   + 
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
    return(p)
  }
}

### PRINT ALL PLOTS
library(ggforce)
make_all_plots <- function(net, cellline, scale='none'){
  axmin = -1
  axmax = 11
  
  plots <- list()
  plots[[1]] <- plot_network(net, which_reps = c(0,1,2,3,4), show_labels = FALSE, nchar_bar=1) + ggtitle(paste0(cellline), subtitle='')
  
  for(i in seq(1,4)){
    #print(plot_network(net, which_reps = c(i), nchar_bar = 5) + ggtitle(paste0(cellline), subtitle=paste0('Replicate ', i)) )
    if(scale=='none'){
      p <- ggplot(net %>% filter(replicate %in% c(i), nchar(barpairid) > 5), aes(x=x.coord, y=y.coord)) +
        geom_shape(aes(group=barpairid), color='black', alpha=0.8, fill=NA) 
    }
    
    if(scale=='abundance'){
      p <- ggplot(net %>% filter(replicate %in% c(i), nchar(barpairid) > 5), aes(x=x.coord, y=y.coord)) +
        geom_shape(aes(group=barpairid, size=prop_bar), color='black', alpha=0.8, fill=NA)
    }
    
    if(scale=='reps'){
      p <- ggplot(net %>% filter(replicate %in% c(i), nchar(barpairid) > 5), aes(x=x.coord, y=y.coord)) +
        geom_shape(aes(group=barpairid, size=as.factor(n_reps_pair)), color='black', alpha=0.8, fill=NA)
    }
    
    p <- p + 
      geom_point(color='white', alpha=1, size=4) + theme_void() +
      theme(legend.position='none', aspect.ratio = 1)  + 
      geom_text(data= net  %>% filter(replicate %in% c(i), nchar(barpairid) > 5) %>% dplyr::select(barnum, barcolor, x.coord, y.coord) %>% distinct(), 
                aes(x=x.coord, y=y.coord, label = barnum, color=barcolor), size = 2.5, fontface = "bold")  +
      scale_color_manual(values=c('fusion'='goldenrod','gfp'='forestgreen','mcherry'='red3')) + 
      ggtitle(paste0(cellline), subtitle=paste0('Replicate ', i)) + 
      ylim(axmin,axmax) + xlim(axmin,axmax)  + 
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
    
    plots[[1+i]] <- p
  }
  
  if(scale=='none'){
    p2 <- ggplot(net %>% filter(replicate %in% c(0,1,2,3,4)) %>% filter(nchar(barpairid) > 5), aes(x=x.coord, y=y.coord)) +
      geom_shape(aes(group=barpairid), color='black', fill=NA, alpha=0.8)  +
      geom_point(color='white', alpha=1, size=4) + theme_void() +
      theme(legend.position='none', aspect.ratio = 1)  + 
      geom_text(data= net  %>% filter(nchar(barpairid) > 5) %>% dplyr::select(barnum, barcolor, x.coord, y.coord) %>% distinct(), 
                aes(x=x.coord, y=y.coord, label = barnum, color=barcolor), size = 2.5, fontface = "bold") +
      scale_color_manual(values=c('fusion'='goldenrod','gfp'='forestgreen','mcherry'='red3')) 
    #+ scale_size_manual(values=c('1'=0.5, '2'=2, '3'=3,'4'=4)) + ggtitle(paste0(cellline), subtitle='all replicates') 
  }
  
  if(scale=='abundance'){
    p2 <- ggplot(net %>% filter(replicate %in% c(0,1,2,3,4)) %>% filter(nchar(barpairid) > 5), aes(x=x.coord, y=y.coord)) +
      geom_shape(aes(group=barpairid, size=prop_bar), color='black', alpha=0.5, fill=NA) +
      geom_point(color='white', alpha=1, size=7)
  }
  
  if(scale=='reps'){
    p2 <- ggplot(net %>% filter(replicate %in% c(0,1,2,3,4)) %>% filter(nchar(barpairid) > 5), aes(x=x.coord, y=y.coord)) +
      geom_shape(aes(group=barpairid, size=as.factor(n_reps_pair)), color='black', alpha=0.5, fill=NA) +
      geom_point(color='white', alpha=1, size=7) #+ scale_size_manual(values=c('1'=0.5, '2'=2, '3'=3,'4'=4))
  }
  
  p2 <- p2 +
      theme_void() +
      theme(legend.position='none', aspect.ratio = 1)  + 
      scale_color_manual(values=c('fusion'='goldenrod','gfp'='forestgreen','mcherry'='red3'))  + 
    ylim(axmin,axmax) + xlim(axmin,axmax)  +
    ggtitle(paste0(cellline), subtitle='')  + 
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  
  
  plots[[1+i+1]] <- p2
  
  return(plots)
}

############################
############################
############################
############################



##### CELL CLASSIFICATION
classify_cells_fusion_and_parent <- function(obj, min_fusion_clone_cells_per_rep){
  ###### IDENTIFY FUSION CLONES
  ### will only use if at least `min_fusion_clone_cells_per_rep` cells within the fusion sample contain the same barcode set
  confident_fusions <- obj@meta.data %>% filter(condition == 'fusion', clone_color == 'fusion') %>% 
    rownames_to_column('cellid') %>% 
    group_by(gbarid, rbarid, barpairid, replicate, n_gfp_bars, n_mch_bars) %>% 
    summarize(rep_count = n()) %>% ungroup() %>% 
    filter(rep_count >= min_fusion_clone_cells_per_rep) %>% 
    mutate(isConfidentFusion = TRUE) %>% 
    mutate(condition='fusion') %>% 
    ungroup()
  
  ###### IDENTIFY PARENTAL CLONES AND COUNT HOW MANY PARTNERS THEY HAVE
  parent_fusion_events <- confident_fusions %>% 
    dplyr::select(gbarid, rbarid, barpairid) %>% 
    distinct() %>% 
    pivot_longer(c(gbarid, rbarid), values_to = 'parent_id') %>% 
    group_by(parent_id) %>% 
    summarize(n_isFusionParent = n()) %>% 
    arrange(desc(n_isFusionParent)) %>% 
    ungroup()
  
  ### if a parental clone has two barcodes, then keep its tally and also add to its matched single barcode's tallies
  ## find multi-barcoded parents
  parent_fusion_events <- parent_fusion_events %>% mutate(n_parent_bars = 1+str_count(parent_id, "_"))
  
  ## break into single barcodes
  mod_counts <- parent_fusion_events[parent_fusion_events$n_parent_bars > 1,] %>% 
    separate(parent_id, into=c('b1','b2'), sep='_') %>% 
    pivot_longer(c(b1,b2), values_to = 'parent_id') %>% drop_na() %>% 
    dplyr::select(parent_id, n_isFusionParent)
  
  # add to data frame
  for(i in seq(1, nrow(mod_counts))){
    if(mod_counts$parent_id[i] %in% parent_fusion_events$parent_id){
      parent_fusion_events$n_isFusionParent[parent_fusion_events$parent_id==mod_counts$parent_id[i]] <- mod_counts$n_isFusionParent[i] + parent_fusion_events$n_isFusionParent[parent_fusion_events$parent_id==mod_counts$parent_id[i]]
    }
    if(mod_counts$parent_id[i] %notin% parent_fusion_events$parent_id){
      add <- data.frame(parent_id=mod_counts$parent_id[i], n_isFusionParent=mod_counts$n_isFusionParent[i], n_parent_bars = NA)
      parent_fusion_events <- rbind(parent_fusion_events, add)
    }
  }
  
  parent_fusion_events <- parent_fusion_events %>% 
    dplyr::select(parent_id, n_isFusionParent) %>% 
    dplyr::rename('barpairid'='parent_id') %>% 
    rowwise() %>% mutate(clone_color = case_when(str_count(barpairid, "g") >= 1~ 'green', 
                                                 str_count(barpairid, "r") >= 1~ 'red'))
  
  obj@meta.data <- obj@meta.data %>% rownames_to_column('cellid') %>% left_join(confident_fusions) %>% left_join(parent_fusion_events) %>% column_to_rownames('cellid')
  
  return(obj)
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### save figure function
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {png(paste0(fig_path, name, ".", type),width = width, height = height, units = "in", res = res)} 
  else {pdf(paste0(fig_path, name, ".", type),width = width, height = height)}
  print(plots)
  dev.off()
}