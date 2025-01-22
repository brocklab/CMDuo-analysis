library(ggrepel)
DE_volcano_pick_geneset <- function(res = fuse_markers, 
                                    genes = c('GAPDH','ESAM'),
                                    n_genes = Inf,
                                    title=paste0(''), 
                                    xmax=NULL, 
                                    ymax=NULL,
                                    colrs=c('blue','red','black'),
                                    set_max = 300,
                                    plot_type = 'p_val_adj')
{
  
  plot_volc <- res %>% 
    select(gene, log2FC, p_val_adj, p_val) %>% 
    mutate(of_interest = case_when(gene %in% genes ~ 'test', .default = "other" )) %>% 
    mutate(neg.log10pval = case_when(plot_type == 'p_val' ~ -log10(p_val), plot_type == 'p_val_adj' ~ -log10(p_val_adj))) %>% 
    mutate(neg.log10pval_adj = case_when(neg.log10pval == Inf ~ set_max,
                                         neg.log10pval != Inf ~ neg.log10pval))
  
  if(is.null(xmax)==TRUE){
    xmax = round(max(abs(res$log2FC)),1) + 0.1
  }
  xmin = -xmax
  
  if(is.null(ymax)==TRUE){
    ymax = round(-log10(min(res$p_val_adj))) + 1
  }
  
  volc <- ggplot(plot_volc %>% filter(of_interest!='test'), aes(x=log2FC, y=neg.log10pval_adj)) + 
    geom_point(alpha=0.5, size=1, color='grey', data=plot_volc) +
    geom_hline(yintercept = -log10(0.05), size=0.1) +
    geom_vline(xintercept = c(-0.2, 0.2), size=0.1) +
    geom_point(size=1, color=colrs[1], alpha=1, data=plot_volc %>% filter(of_interest=='test', log2FC < 0)) +
    geom_point(size=1, color=colrs[2], alpha=1, data=plot_volc %>% filter(of_interest=='test', log2FC > 0)) +
    theme_classic() + 
    ggtitle(title) + xlab('log2FC') + ylab('-log10(padj)') +
    xlim(xmin, xmax) +
    ylim(0, ymax)
  
  if(length(plot_volc$of_interest[plot_volc$of_interest=='test'] <= n_genes)){
    markgenes_hi <- plot_volc %>% 
      filter(of_interest=='test', log2FC > 0) %>% 
      mutate(score = log2FC*-log10(p_val_adj)) %>% 
      arrange(desc(score)) %>% 
      filter(score > 0)  %>% 
      head(n_genes) %>% 
      pull(gene)
    
    markgenes_lo <- plot_volc %>% 
      filter(of_interest=='test', log2FC < 0) %>% 
      mutate(score = abs(log2FC)*-log10(p_val_adj)) %>% 
      arrange(desc(score)) %>% 
      filter(score > 0)  %>% 
      head(n_genes) %>% 
      pull(gene)
    
    volc <- volc + geom_text_repel(data=plot_volc %>% filter(gene %in% unique(c(markgenes_hi, markgenes_lo))),
                                   bg.color = "white", bg.r = 0.2, box.padding = 0.25,
                                   max.overlaps = Inf, aes(label=gene), size=2, color='black', segment.size = 0.2) 
  }
  return(volc)
}


################################################################################

DE_volcano_top_and_pick_geneset <- function(res = fuse_markers, 
                                            shared_genes = c('GAPDH','ESAM'),
                                            highlight_shared_genes = c('GAPDH','ESAM'),
                                            n_genes = Inf,
                                            title=paste0(''), 
                                            p_cut = 0.05, 
                                            l2fc_cut = 0.2, 
                                            n_genes_shared = Inf, 
                                            n_genes_top = 10, 
                                            gene_label_size = 2.5,
                                            label.padding = 0.25,
                                            axis_title_font_size = 10,
                                            plot_title_font_size = 14,
                                            pick_seed = 1212,
                                            xmin=NULL,
                                            xmax=NULL, 
                                            ymin = 0,
                                            ymax=NULL,
                                            colrs=c('blue','red','black'),
                                            highlight_colrs = c('red','red'),
                                            highlight_shape = 3, 
                                            colrs_main=c('dodgerblue','goldenrod4'),
                                            top_color = 'goldenrod4',
                                            background_point_color = 'grey80', 
                                            plot_type = 'p_val_adj')
{
  
  top_hi <- res %>% 
    filter(log2FC > 0, gene %notin% highlight_shared_genes) %>% 
    mutate(score = log2FC*-log10(p_val_adj)) %>% 
    arrange(desc(score)) %>% 
    filter(score > 0) %>% 
    head(n_genes_top) %>% 
    pull(gene) 
  
  top_lo <- res %>% 
    filter(log2FC < 0, gene %notin% highlight_shared_genes) %>% 
    mutate(score = abs(log2FC)*-log10(p_val_adj)) %>% 
    arrange(desc(score)) %>% 
    filter(score > 0) %>% 
    head(n_genes_top) %>% 
    pull(gene)
  
  
  plot_volc <- res %>% 
    select(gene, log2FC, p_val_adj, p_val) %>% 
    mutate(is_shared = case_when(gene %in% shared_genes ~ 'shared', .default = 'not')) %>% 
    mutate(of_interest = case_when(gene %in% highlight_shared_genes ~ 'shared', 
                                   gene %in% c(top_hi, top_lo) ~ 'top',
                                   .default = "other" )) %>% 
    mutate(neg.log10pval = case_when(plot_type == 'p_val' ~ -log10(p_val), plot_type == 'p_val_adj' ~ -log10(p_val_adj))) %>% 
    mutate(neg.log10pval_adj = case_when(neg.log10pval > ymax ~ ymax, 
                                         .default = neg.log10pval)) %>% 
    mutate(log2FC = case_when(log2FC > xmax ~ xmax, 
                              log2FC < -xmax ~ -xmax,
                              .default = log2FC)) 
  
  if(is.null(xmax)==TRUE){
    xmax = round(max(abs(res$log2FC)),1) + 0.1
  }
  
  if(is.null(xmin)==TRUE){
    xmin = -xmax
  }else{
    xmin = xmin
  }
  
  
  if(is.null(ymax)==TRUE){
    ymax = round(-log10(min(res$p_val_adj))) + 1
  }
  
  volc <- ggplot(plot_volc, aes(x=log2FC, y=neg.log10pval_adj) )  +
    geom_point(data=plot_volc %>% filter(of_interest %notin% c('shared','top')), 
               alpha=0.5, size=1, color=background_point_color) +
    geom_point(data=plot_volc %>% filter(of_interest %notin% c('shared','top'), log2FC <= -l2fc_cut, p_val_adj <= p_cut), 
               alpha=0.5, size=1, color=colrs_main[1]) +
    geom_point(data=plot_volc %>% filter(of_interest %notin% c('shared','top'), log2FC >=  l2fc_cut, p_val_adj <= p_cut), alpha=0.5, size=1, color=colrs_main[2]) +
    geom_point(size=1, shape=16, color=colrs[1], alpha=1, data=plot_volc %>% filter(is_shared=='shared', log2FC < 0)) +
    geom_point(size=1, shape=16, color=colrs[2], alpha=1, data=plot_volc %>% filter(is_shared=='shared', log2FC > 0)) +
    geom_point(size=1, color=colrs_main[1], alpha=1, data=plot_volc %>% filter(of_interest=='top', log2FC < 0)) +
    geom_point(size=1, color=colrs_main[2],   alpha=1, data=plot_volc %>% filter(of_interest=='top', log2FC > 0)) +
    geom_point(size=1, shape=highlight_shape, color=highlight_colrs[1], alpha=1, data=plot_volc %>% filter(of_interest=='shared', log2FC < 0)) +
    geom_point(size=1, shape=highlight_shape, color=highlight_colrs[2], alpha=1, data=plot_volc %>% filter(of_interest=='shared', log2FC > 0)) +
    geom_hline(yintercept = -log10(p_cut), size=0.1) +
    geom_vline(xintercept = c(-l2fc_cut, l2fc_cut), size=0.1) +
    theme_classic() + 
    ggtitle(title) + xlab('log2FC') + ylab('-log10(padj)') +
    xlim(xmin, xmax) +
    ylim(ymin, ymax)
  
  
  markgenes_hi <- plot_volc %>% 
    filter(of_interest=='shared', log2FC > 0) %>% 
    mutate(score = log2FC*-log10(p_val_adj)) %>% 
    arrange(desc(score)) %>% 
    filter(score > 0)  %>% 
    head(n_genes_shared) %>% 
    pull(gene)
  
  markgenes_lo <- plot_volc %>% 
    filter(of_interest=='shared', log2FC < 0) %>% 
    mutate(score = abs(log2FC)*-log10(p_val_adj)) %>% 
    arrange(desc(score)) %>% 
    filter(score > 0)  %>% 
    head(n_genes_shared) %>% 
    pull(gene)
  
  # if(length(plot_volc$of_interest[plot_volc$of_interest=='test'] <= n_genes_shared)){
  
  set.seed(pick_seed)
  volc <- volc + geom_text_repel(data=plot_volc %>% filter(gene %in% unique(c(top_hi, top_lo, markgenes_hi, markgenes_lo))),
                                 bg.color = "white", bg.r = 0.2, box.padding = label.padding,
                                 max.overlaps = Inf, aes(label=gene, color=of_interest), size=gene_label_size, 
                                 min.segment.length = 0, segment.size = 0.15) + 
    scale_color_manual(values=c('shared'='black','top'=top_color)) +
    theme_bw() + 
    theme(legend.position = 'none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    labs(y=expression(paste("-",log[10],"(", p[adj],")",sep="")), x=expression(paste(log[2],"FC",sep=""))) + 
    theme(axis.title = element_text(size = axis_title_font_size),
          plot.title = element_text(size = plot_title_font_size) )
  
  
  return(volc)
}