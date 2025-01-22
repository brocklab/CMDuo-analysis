library(tidyverse)
library(fgsea)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### save figure
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {png(paste0(fig_path, name, ".", type),width = width, height = height, units = "in", res = res)} 
  else {pdf(paste0(fig_path, name, ".", type),width = width, height = height)}
  print(plots)
  dev.off()
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### import data
import_de_results <- function(de_tsv_name){
  # load in outputs from de test
  res <- list()
  res[['mb231']] <- read.csv(paste0('outdata/de_results/nebula/nebula_0.1pct_mb231_',de_tsv_name,'.tsv'), sep='\t') %>% 
    remove_rownames() %>% arrange(desc(log2FC))
  res[['h1806']] <- read.csv(paste0('outdata/de_results/nebula/nebula_0.1pct_h1806_',de_tsv_name,'.tsv'), sep='\t') %>% 
    remove_rownames() %>% arrange(desc(log2FC))
  return(res)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### save outputs of fgsea

save_gsea_outs <- function(outs, filename='filename.csv'){
  df <- do.call(rbind, lapply(outs, as.data.frame))
  df <- df %>% rowwise() %>% mutate(leadingEdge = paste(unlist(leadingEdge), collapse=', '))
  write.table(df, file=paste0('outdata/gsea_results/gsea_',filename), quote=FALSE, sep='\t', col.names = TRUE)
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### wrapper to call functions to import DE results, run fgsea, then save results
run_and_save_gsea <- function(de_tsv_name = 'binom_fusers_p05', 
                              min_geneset_size = 20,
                              gs_list = list(h.hallmarks),
                              save_file_name = 'match'){
  
  # import mixed model DE results from nebula
  res <- import_de_results(de_tsv_name)
  
  # run GSEA
  out.gsea <- run_gsea_on(res['mb231'][[1]], res['h1806'][[1]], 
                          min_geneset_size = 20, 
                          gs_list = gs_list)
  if(save_file_name == 'match'){
    save_file_name <- de_tsv_name
  }
  
  # save outputs
  save_gsea_outs(out.gsea, filename=paste0('gsea_',save_file_name,'_new.tsv'))
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### load DE and GSEA outputs for analysis and plotting
load_de_gsea_for_analysis <- function(de_tsv_name){
  outs <- list()
  
  res <- import_de_results(paste0(de_tsv_name))
  
  outs[['mb231.res']] <- res['mb231'][[1]]
  outs[['h1806.res']] <- res['h1806'][[1]]
  
  outs[['gsea']] <- read.csv(paste0('outdata/gsea_results/gsea_',de_tsv_name,'_new.tsv'), sep='\t') %>%
    rowwise() %>% 
    mutate(gs = paste((str_split(pathway, pattern='_')[[1]][1:1]),  collapse = " ")) %>%
    mutate(data=paste0(de_tsv_name)) %>% 
    ungroup()
  
  return(outs)
}



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### load gene sets
h.hallmarks <- gmtPathways('data/gsea_gmt/h.all.v2023.2.Hs.symbols.gmt')
c2.KEGG <- gmtPathways('data/gsea_gmt/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt')
c3.transcription_factors <- gmtPathways('data/gsea_gmt/c3.tft.v2023.2.Hs.symbols.gmt')
c5.go.biological_process <- gmtPathways('data/gsea_gmt/c5.go.bp.v2023.2.Hs.symbols.gmt')
c5.go.molecular_function <- gmtPathways('data/gsea_gmt/c5.go.mf.v2023.2.Hs.symbols.gmt')

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### function to call fgsea
run_fgsea <- function(de_data, pathways=gs, minSetSize=20){
  # first need to add scores to rank DE analysis outputs in order to rank genes. 
  # For this we will use avg log fold change
  res <- de_data %>% mutate(score=log2FC)%>% arrange(desc(score))
  ranks <- deframe(res %>% select(gene, score))
  # run fgsea
  fgseaRes <- fgsea(pathways=pathways, 
                    stats=ranks,
                    scoreType='std',
                    minSize  = minSetSize,
                    nPermSimple = 10000)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### run gsea with gene set and de results
library(fgsea)
fgsea_geneset <- function(res_list, gs, min_geneset_size){
  fgsea.outs <- data.frame()
  for(i in 1:length(res_list)){
    new.outs <- run_fgsea(de_data    = res_list[[i]], 
                          pathways   = gs, 
                          minSetSize = min_geneset_size) %>% 
      mutate(cellline=names(res_list)[i]) %>% arrange(desc(NES))
    fgsea.outs <- rbind(fgsea.outs, new.outs)
  }
  return(fgsea.outs)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### wrapper to run gsea tests on multiple gene sets
run_gsea_on <- function(mb231.res, h1806.res, 
                        min_geneset_size = 20,  
                        gs_list = list(h.hallmarks, c5.go.biological_process, c2.KEGG, c5.go.molecular_function)){
  gsea_outs <- list()
  for(i in seq(1,length(gs_list))){
    gsea_outs[[i]] <- fgsea_geneset(res_list = list('MDA-MB-231' = mb231.res, 
                                                    'HCC1806' = h1806.res),
                                    gs = gs_list[[i]],
                                    min_geneset_size = min_geneset_size)
  }
  return(gsea_outs)
}




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### save outputs of GSEA
save_gsea_outs <- function(outs, filename='filename.csv'){
  df <- do.call(rbind, lapply(outs, as.data.frame))
  df <- df %>% rowwise() %>% mutate(leadingEdge = paste(unlist(leadingEdge), collapse=', '))
  write.table(df, file=paste0('outdata/gsea_results/',filename), quote=FALSE, sep='\t', col.names = TRUE)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### load outputs of GSEA
load_gsea_outs <- function(filename, data_col_id){
  outs <- read.csv(paste0('outdata/gsea_results/',filename,'.tsv'), sep='\t') %>%
    rowwise() %>% 
    mutate(gs = paste((str_split(pathway, pattern='_')[[1]][1:1]),  collapse = " ")) %>%
    mutate(data=data_col_id) %>% 
    ungroup()
  return(outs)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### plot results
dot_plot_gsea <- function(gsea.outs, pick_gs='HALLMARK', pval_threshold  = 0.25, padj_threshold  = 0.25, NES_min = 1, 
                          n_pathway_plot = 20, min_celllines_in=2, gs_name_length=1,
                          data_title = '', only_pos = FALSE){
  # relabel pathway to not include underscores and gene set name, create new column with gene set name
  # only plot pathways which are significant in at least one cell line
  filterpaths <- gsea.outs %>% 
    filter(gs == pick_gs) %>% 
    rowwise() %>% 
    mutate(gs = paste((str_split(pathway, pattern='_')[[1]][1:gs_name_length]),  collapse = " ")) %>% 
    mutate(pathway = paste((str_split(pathway, pattern='_')[[1]][-1:-gs_name_length]),  collapse = " ")) %>% 
    group_by(pathway) %>% 
    mutate(n_samples = n()) %>% ungroup() %>% 
    filter(n_samples >= min_celllines_in) %>% 
    mutate(keep_path = case_when(padj <= pval_threshold & padj <= padj_threshold & abs(NES) >= NES_min ~ 1, 
                                 .default=0)) %>% 
    group_by(pathway) %>% 
    mutate(keep_score = sum(keep_path), mean_sign = mean(sign(NES))) %>% 
    ungroup() %>% 
    filter(keep_score >= 1) %>% 
    group_by(pathway) %>% 
    mutate(sumNES = sum(NES, na.rm=TRUE)) %>% 
    arrange(desc(sumNES))
  
  # extract pathways that are high in both and significant in at least one
  high <- filterpaths %>% 
    dplyr::select(pathway, mean_sign, sumNES) %>% 
    distinct() %>% 
    filter(mean_sign > 0) %>% 
    head(n=n_pathway_plot) %>% 
    pull(pathway)
  # extract pathways that are low in both and significant in at least one
  lows <- filterpaths %>% 
    dplyr::select(pathway, mean_sign, sumNES) %>% 
    distinct() %>% filter(mean_sign < 0) %>% 
    tail(n=n_pathway_plot) %>% 
    pull(pathway)
  
  if(only_pos == TRUE){
    lows <- c()
  }
  
  # filter data for plotting
  fgsea.outs.plot <- filterpaths %>% 
    filter(pathway %in% c(high, lows)) %>% 
    arrange(sumNES)
  
  # set plot order
 # fgsea.outs.plot$pathway <- factor(x = fgsea.outs.plot$pathway, levels = unique(fgsea.outs.plot$pathway))
  means <- fgsea.outs.plot %>% 
    group_by(pathway) %>% 
    summarize(, padj = mean(padj), NES = mean(NES)) %>% 
    mutate(cellline = 'mean', gs=pick_gs, sumNES=0)
  
  fgsea.outs.plot <- rbind(fgsea.outs.plot, means)
  
  fgsea.outs.plot$pathway <- factor(x = fgsea.outs.plot$pathway, 
                                    levels = fgsea.outs.plot %>% 
                                      filter(cellline=='mean') %>% 
                                      arrange((NES)) %>% 
                                      pull(pathway)
                                    )
  
  dot_plot <- ggplot(fgsea.outs.plot, aes(y=pathway, x=NES)) +
    geom_point(aes(color=as.factor(cellline), size=-log10(padj)), alpha=0.8) +
    theme_bw() + 
    scale_color_manual(values=c('MDA-MB-231'='#F8766D', 'HCC1806'='#00BFC4','mean' = 'black')) +
    ggtitle(paste0(pick_gs,' ',data_title)) + 
    ylab('') + 
    labs(color = "cell line", size=expression(paste("-",log[10],"(", p[adj],")",sep=""))) + 
    geom_vline(xintercept=0, linewidth=0.1)
  
  return(dot_plot)
}