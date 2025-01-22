library(tidyverse)
library(ggrepel)

########################################################
# plot summary of n fusion partners in population
n_fusion_partners_plot <- function(h1806_classified, mb231_classified, x.text.size = 12, y.text.size=8, label.font.size=3){
  h1806_stats <- h1806_classified@meta.data %>% 
    filter(n_bar_cell >=1, condition %in% c('initial')) %>% 
    dplyr::select(barpairid, n_isFusionParent) %>% distinct() %>% 
    mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
    mutate(n_clones_total = n()) %>% group_by(n_isFusionParent, n_clones_total) %>% 
    summarize(n_with_n_edges = n()) %>% 
    mutate(perc_n_with_n_edges = n_with_n_edges/n_clones_total) %>% mutate(cellline='HCC1806')
  
  mb231_stats <- mb231_classified@meta.data %>% 
    filter(n_bar_cell >=1, condition %in% c('initial')) %>% 
    dplyr::select(barpairid, n_isFusionParent) %>% distinct() %>% 
    mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
    mutate(n_clones_total = n()) %>% group_by(n_isFusionParent, n_clones_total) %>% 
    summarize(n_with_n_edges = n()) %>% 
    mutate(perc_n_with_n_edges = n_with_n_edges/n_clones_total) %>% mutate(cellline='MDA-MB-231')
  
  stats <- rbind(mb231_stats, h1806_stats)
  
  stats <- stats %>% 
    mutate(group = case_when(n_isFusionParent <=2 ~ as.character(n_isFusionParent), n_isFusionParent >2 ~ '3+')) %>% 
    group_by(group, cellline) %>% 
    mutate(perc_n_with_n_edges_group = sum(perc_n_with_n_edges), n_with_n_edges_group=sum(n_with_n_edges))
  
  n_fusion_partners_plot <- ggplot() +
    theme_classic() + theme(strip.background = element_blank()) +
    scale_fill_manual(values=cellline_colors) +
    geom_bar(data= stats, aes(x=group, y=perc_n_with_n_edges_group, group=cellline, fill=cellline),
             stat='identity', position="dodge", color='black', size=0.1) + 
    geom_text(
      data = stats %>% dplyr::select(cellline, group,n_clones_total, n_with_n_edges_group, perc_n_with_n_edges_group) %>% distinct(),
      position = position_dodge(width=0.9),
      vjust = 0.5,
      size=label.font.size,
      aes(x=group, y=perc_n_with_n_edges_group+0.07, group=cellline, 
          label=paste0(100*round(perc_n_with_n_edges_group,3),"%","\n","(",n_with_n_edges_group,'/',n_clones_total,")"))) +
    ylab('Proportion of population with n fusion partners') +
    xlab('Number of fusion partners') + labs(fill='Population') + 
    theme(axis.text.x = element_text(size = x.text.size), 
          axis.text.y = element_text(size = y.text.size),
          axis.title.x = element_text(size = x.text.size+2),
          axis.title.y = element_text(size = y.text.size+2)) 
  
  return(n_fusion_partners_plot)
}




########################################################
#### plot of n partners colored by random chance
binom_test_n_partners <- function(obj, pcut = 0.05, label.font.size=3, x.text.size=10, y.text.size=10, vertical=TRUE){
  
  # find abundance of each clone within the red and green populations
  abund_int <- obj@meta.data %>% 
    filter(n_bar_cell >=1, condition %in% c('initial'), clone_color %in% c('red','green')) %>% 
    group_by(clone_color) %>% 
    mutate(total_int_cells = n()) %>% 
    group_by(barpairid, total_int_cells, clone_color) %>% 
    mutate(n_clone=n()) %>% rowwise() %>% 
    mutate(abundance = n_clone/total_int_cells) %>% ungroup() %>% 
    mutate(group = case_when(n_isFusionParent <=2 ~ as.character(n_isFusionParent), n_isFusionParent >2 ~ '3+', .default = '0')) %>%
    dplyr::select(barpairid, n_isFusionParent, group, abundance, clone_color) %>% 
    distinct() %>% ungroup()
  
  abund_int$n_isFusionParent[is.na(abund_int$n_isFusionParent)==TRUE] <- 0
  
  n_fusion_clones <- obj@meta.data %>% filter(condition %in% c('fusion')) %>% 
    filter(isConfidentFusion==TRUE) %>% pull(barpairid) %>% unique() %>% length()
  
  # Asks if the true probability of success is greater than expected from proportion
  # low p-value = more likely
  binom_test <- abund_int %>% 
    filter(n_isFusionParent > 0) %>% 
    rowwise() %>% 
    mutate(binom_test = binom.test(x = n_isFusionParent, n = n_fusion_clones, p = abundance, alternative = "greater")$p.value) %>%
    arrange((binom_test)) %>% 
    mutate(binom_test = round(binom_test, 5)) %>% 
    ungroup() 
  
  abund_int <- abund_int %>% 
    left_join(binom_test) %>% 
    mutate(sig_events = case_when(binom_test <= pcut ~ 'p <= 0.05', .default='not sig')) %>% 
    mutate(label = case_when(sig_events=='p <= 0.05' ~ barpairid)) 
  
  if(vertical==FALSE){
    p <- 
      ggplot(abund_int, aes(x=as.factor(group), y=abundance)) +
      theme_classic() +
      geom_violin(fill='grey80', alpha=1, color=NA, aes(group=as.factor(group)) ) + 
      geom_point(position=position_jitter(seed = 1, width=0.1, height=0), size=1,  aes(color=sig_events)) + 
      geom_text_repel(aes(label=label),
                      bg.color = "white", bg.r = 0.05, box.padding = 0.25,
                      size=label.font.size, seed=123, min.segment.length = 0,
                      max.overlaps = Inf,
                      position=position_jitter(seed=1, width=0.1, height=0)) +
      scale_color_manual(values=c('p <= 0.05'='magenta','not sig'='grey20'),
                         labels=c('p <= 0.05'=expression("p"<=0.05), 'not sig'='not sig')) + 
      theme(axis.text.x = element_text(size = x.text.size), 
            axis.text.y = element_text(size = y.text.size),
            axis.title.x = element_text(size = x.text.size+2),
            axis.title.y = element_text(size = y.text.size+2)) +
      xlab('Number of fusion partners') +
      ylab('Clonal abundance in initial population') + labs(color='p(random)')
  }
  if(vertical==TRUE){
    p <- ggplot(abund_int, aes(y=as.factor(n_isFusionParent), x=abundance)) +
      theme_classic() +
      geom_violin(fill='grey80', alpha=1, color=NA, aes(group=as.factor(n_isFusionParent)) ) + 
      geom_point(position=position_jitter(seed = 1, width=0, height=0.1), size=1,  aes(color=sig_events)) + 
      geom_text_repel(aes(label=label),
                      bg.color = "white", bg.r = 0.05, box.padding = 0.25,
                      size=label.font.size, seed=123, min.segment.length = 0,
                      max.overlaps = Inf,
                      position=position_jitter(seed=1, width=0, height=0.1)) +
      scale_color_manual(values=c('p <= 0.05'='magenta','not sig'='grey20'),
                         labels=c('p <= 0.05'=expression("p"<=0.05), 'not sig'='not sig'))
      theme(axis.text.x = element_text(size = x.text.size), 
            axis.text.y = element_text(size = y.text.size),
            axis.title.x = element_text(size = x.text.size+2),
            axis.title.y = element_text(size = y.text.size+2)) +
      ylab('Number of fusion partners') +
      xlab('Clonal abundance in initial population') + labs(color='p(random)')
  }
  return(list(p, abund_int))
}



