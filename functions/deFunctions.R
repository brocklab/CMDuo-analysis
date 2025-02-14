### subset data to only have cells from the sorted fusion and sorted control samples
### further subset data to only contain clones for which there are at least 3 cells
subset_on_fusion_v_ctrl <- function(dat, n_cell_min=3){
  # Find high confidence fusion clones and related parental barcodes
  fusion_pairs <- dat@meta.data %>% dplyr::filter(isConfidentFusion==TRUE) %>% 
    dplyr::select(gbarid, rbarid, barpairid) %>% distinct()
  
  # subset data to only contain high confidence fusion cells or sorted control cells
  dat_fuse_ctrl <- subset(dat, subset = barpairid %in% fusion_pairs$barpairid & condition=='fusion' & isConfidentFusion==TRUE |
                            condition=='ctrl'   & fusion!='fusion')
  
  # classify cells as fusion parents or fusion progeny
  dat_fuse_ctrl@meta.data <- dat_fuse_ctrl@meta.data %>% 
    mutate(cell_class = case_when(condition=='fusion' & clone_color=='fusion' ~ 'fusion', 
                                  condition=='ctrl' & clone_color %in% c('green','red') & exp_color!='gfp_mch' ~ 'control',
                                  .default='not_classified'))
  
  
  # find clones that have at least n_cells
  min_cell_pass_clones <- dat_fuse_ctrl@meta.data %>% 
    dplyr::select(barpairid) %>% 
    pivot_longer(c(barpairid), values_to = 'barid') %>% 
    drop_na() %>% 
    group_by(barid) %>% 
    summarize(cells_per_clone = n()) %>% 
    dplyr::filter(cells_per_clone >= n_cell_min) %>% 
    pull(barid) %>% 
    unique()
  
  # subset to only contain high confidence fusion cells with >= n cels per clone
  dat_fuse_ctrl_mincells <- subset(dat_fuse_ctrl, subset = (barpairid %in% min_cell_pass_clones) & (cell_class %in% c('fusion','control')) )
  
  # plot subsetted data
  p1 <- ggumap(get_embedding(dat_fuse_ctrl)) + 
    geom_point(aes(color=as.factor(hclust)),alpha=0.5,size=0.2) + 
    labs(color='cluster') + 
    ggtitle(paste0('fusion sets, n cells = ',nrow(get_embedding(dat_fuse_ctrl))), 
            subtitle = paste0('n clones = ',length(unique(get_embedding(dat_fuse_ctrl)$barpairid)))) + 
    theme(aspect.ratio = 0.8)
  
  p2 <- ggumap(get_embedding(dat_fuse_ctrl_mincells)) + 
    geom_point(aes(color=as.factor(hclust)),alpha=0.5,size=0.2) + 
    labs(color='cluster') + 
    ggtitle(paste0('fusion sets ',n_cell_min,' cell min, n cells = ',nrow(get_embedding(dat_fuse_ctrl_mincells))), 
            subtitle = paste0('n clones = ',length(unique(get_embedding(dat_fuse_ctrl_mincells)$barpairid)))) +
    theme(aspect.ratio = 0.8)
  
  p3 <- ggumap(get_embedding(dat_fuse_ctrl_mincells)) + 
    geom_point(aes(color=as.factor(color)),alpha=0.5,size=0.2)  + 
    labs(color='color') + 
    ggtitle(paste0('fusion sets ',n_cell_min,' cell min, n cells = ',nrow(get_embedding(dat_fuse_ctrl_mincells))), 
            subtitle = paste0('n clones = ',length(unique(get_embedding(dat_fuse_ctrl_mincells)$barpairid)))) +
    theme(aspect.ratio = 0.8) + 
    scale_color_manual(values=c('gfp'='forestgreen','mcherry'='red3','fusion'='goldenrod'))
  
  p4 <- ggumap(get_embedding(dat_fuse_ctrl_mincells)) + 
    geom_point(aes(color=as.factor(cell_class)),alpha=0.5,size=0.2) + 
    labs(color='cell_class') + 
    ggtitle(paste0('fusion sets ',n_cell_min,' cell min, n cells = ',nrow(get_embedding(dat_fuse_ctrl_mincells))), 
            subtitle = paste0('n clones = ',length(unique(get_embedding(dat_fuse_ctrl_mincells)$barpairid)))) +
    theme(aspect.ratio = 0.8) + 
    scale_color_manual(values=c('control'='dodgerblue','fusion'='goldenrod'))
  
  return(list(dat_fuse_ctrl_mincells, p1, p2, p3, p4))
}


###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### subset data on fusogenic versus never fuser
subset_on_fusogenic <- function(dat, cellline='h1806', n_cell_min = 5, hc_thresh = 2, multifuser_confidence = 'high'){
  
  dat@meta.data$n_isFusionParent[is.na(dat@meta.data$n_isFusionParent)==TRUE] <- 0
  
  # remove ESAM-neg cluster from mb231 analysis since cell sorting removed these from control
  if(cellline=='mb231'){
    # Calculate average expression per cluster
    esam_expression <- t(AverageExpression(dat, features = c('ESAM'), group.by = 'seurat_clusters')$RNA)
    esam_pos_clust <- data.frame(esam_expression) %>% rownames_to_column('clust') %>% mutate(clust = substring(clust, 2)) %>% filter(ESAM > 0.4) %>% pull(clust)
    dat <- subset(dat, subset = seurat_clusters %in% esam_pos_clust)
  }
  
  #### classify cells
  high_confidence_parents <- dat@meta.data %>% 
    dplyr::filter(n_isFusionParent >= hc_thresh, condition %in% c('presort','initial')) %>%
    dplyr::select(rbarid, gbarid) %>% 
    pivot_longer(c(rbarid,gbarid)) %>% 
    dplyr::select(value) %>% 
    drop_na() %>% 
    group_by(value) %>% 
    summarize(Freq = n()) %>% 
    dplyr::filter(Freq >= n_cell_min) %>% 
    pull(value) %>% unique() %>% sort()
  
  low_confidence_parents <- dat@meta.data  %>% 
    dplyr::filter(n_isFusionParent < hc_thresh, n_isFusionParent > 0, condition %in% c('presort','initial')) %>% 
    dplyr::select(rbarid, gbarid) %>% 
    pivot_longer(c(rbarid,gbarid)) %>% 
    dplyr::filter(value %notin% high_confidence_parents) %>% 
    drop_na() %>% 
    group_by(value) %>% summarize(Freq = n()) %>% 
    dplyr::filter(Freq >= n_cell_min) %>% 
    pull(value) %>% unique() %>% sort()
  
  not_parents <- dat@meta.data  %>% 
    dplyr::filter(n_isFusionParent == 0, condition %in% c('presort','initial')) %>% 
    dplyr::select(rbarid, gbarid) %>% 
    pivot_longer(c(rbarid,gbarid)) %>% 
    dplyr::filter(value %notin% high_confidence_parents) %>% 
    dplyr::select(value) %>% 
    drop_na() %>% 
    group_by(value) %>% 
    summarize(Freq = n()) %>% 
    dplyr::filter(Freq >= n_cell_min) %>% 
    pull(value) %>% unique() %>% sort()
  ####
  
  # Add classifications to meta data
  # classify cells as fusion parents or progeny or not_parent
  dat@meta.data <- dat@meta.data %>% 
    mutate(cell_class = case_when((gbarid %in% high_confidence_parents & condition %in% c('presort','initial')) ~ 'parent',
                                  (rbarid %in% high_confidence_parents & condition %in% c('presort','initial')) ~ 'parent', 
                                  (gbarid %in%  low_confidence_parents & condition %in% c('presort','initial')) ~ 'parent_lc',
                                  (rbarid %in%  low_confidence_parents & condition %in% c('presort','initial')) ~ 'parent_lc', 
                                  (gbarid %in% not_parents             & condition %in% c('presort','initial')) ~ 'not_parent',
                                  (rbarid %in% not_parents             & condition %in% c('presort','initial')) ~ 'not_parent',
                                  .default = 'other'))
  
  # subset data
  if(multifuser_confidence == 'high'){
    dat_fusogenic <- subset(dat, subset = cell_class %in% c('parent', 'not_parent') )
  }
  if(multifuser_confidence=='any'){
    dat_fusogenic <- subset(dat, subset = cell_class %in% c('parent','parent_lc', 'not_parent') )
    dat_fusogenic$cell_class[dat_fusogenic$cell_class=='parent_lc'] <- 'parent'
  }
  
  # plot subsetted data
  p1 <- ggumap(get_embedding(dat_fusogenic)) + 
    geom_point(aes(color=as.factor(hclust)),alpha=0.5,size=0.2) + 
    labs(color='cluster') + 
    ggtitle(paste0('fusogenic vs not, ',n_cell_min,' cell min'), 
            subtitle = paste0('n cells = ',nrow(get_embedding(dat_fusogenic)), '\nn clones = ',length(unique(get_embedding(dat_fusogenic)$barpairid)))) +
    theme(aspect.ratio = 0.8)
  
  p2 <- ggumap(get_embedding(dat_fusogenic)) + 
    geom_point(aes(color=as.factor(clone_color)),alpha=0.5,size=0.2)  + 
    labs(color='color') + 
    ggtitle(paste0('fusogenic vs not, ',n_cell_min,' cell min'), 
            subtitle = paste0('n cells = ',nrow(get_embedding(dat_fusogenic)), '\nn clones = ',length(unique(get_embedding(dat_fusogenic)$barpairid))))+
    theme(aspect.ratio = 0.8) + 
    scale_color_manual(values=c('green'='forestgreen','red'='red3','fusion'='goldenrod'))
  
  p3 <- ggumap(get_embedding(dat_fusogenic)) + 
    geom_point(aes(color=as.factor(cell_class)),alpha=0.8,size=0.2) + 
    labs(color='cell_class') + 
    ggtitle(paste0('fusogenic vs not, ',n_cell_min,' cell min'), 
            subtitle = paste0('n cells = ',nrow(get_embedding(dat_fusogenic)), '\nn clones = ',length(unique(get_embedding(dat_fusogenic)$barpairid)))) +
    theme(aspect.ratio = 0.8) + 
    scale_color_manual(values=c('parent'='purple','not_parent'='black')) 
  
  return(list(dat_fusogenic, p1, p2, p3))
}


##########################
### plot clonal composition by class
plot_clonal_composition_byclass <- function(dat, text_position = -150, class_levels = NA){
  if(length(class_levels) == 1){
    class_levels = unique(dat@meta.data$cell_class) %>% sort()
  }
  
  clones_all <- dat@meta.data %>% 
    dplyr::select(barpairid, cell_class) %>% 
    group_by(barpairid) %>% 
    mutate(clone_count_cells = n()) %>% 
    distinct() %>% 
    group_by(cell_class) %>% 
    mutate(class_count_clones = n()) 
  
  clones_all$cell_class <- factor(clones_all$cell_class, levels = class_levels)
  
  bar_order <- clones_all %>% arrange(desc(clone_count_cells)) %>% pull(barpairid) %>% unique()
  
  clones_all$barpairid <- factor(clones_all$barpairid, levels = bar_order)
  
  ggplot() +
    geom_bar(data = clones_all, 
             aes(x=cell_class, y=clone_count_cells, fill=barpairid),
             position='stack', 
             stat='identity', 
             color='black', 
             size=0.1) + 
    theme_classic() + 
    xlab('cell class') + ylab('cells per clone') +
    theme(legend.position = 'none', 
          strip.background = element_blank(),
          aspect.ratio = 1) +
    scale_fill_manual(values=rep(c(wes_palette("Moonrise3"), wes_palette("GrandBudapest2"), 
                                   wes_palette("Royal2"), wes_palette("Royal1"), wes_palette("FrenchDispatch")), 1000)) +
    geom_text(data = clones_all %>% dplyr::select(cell_class, class_count_clones) %>% distinct(),
              aes(x=cell_class, 
                  y=text_position, 
                  label=paste0('clones = ',class_count_clones))) +
    geom_text(data = clones_all %>% dplyr::select(cell_class, clone_count_cells) %>% group_by(cell_class) %>% summarize(total_cells = sum(clone_count_cells)),
              aes(x=cell_class, 
                  y=total_cells+150,
                  label=paste0('cells = ',total_cells))) 
  
}
###########



### get data embedding from seurat
get_embedding <- function(plot_dat){
  umapcell <- plot_dat@reductions$umap@cell.embeddings %>% data.frame()
  if('cellid' %in% colnames(plot_dat@meta.data)==FALSE){
    embedding <- left_join(plot_dat@meta.data %>% rownames_to_column('cellid'), umapcell %>% rownames_to_column('cellid'), by='cellid')
  }
  if('cellid' %in% colnames(plot_dat@meta.data)==TRUE){
    embedding <- left_join(plot_dat@meta.data, umapcell %>% rownames_to_column('cellid'))
  }
  embedding <- embedding %>% dplyr::rename('UMAP_1'='umap_1', 'UMAP_2'='umap_2')
  return(embedding)
}

### make minimal ggplot object from data embedding
ggumap <- function(embedding){
  ggplot(data=embedding, aes(x=UMAP_1, y=UMAP_2)) +
    theme_classic() + xlab('UMAP 1') + ylab('UMAP 2') + 
    theme(text=element_text(size=14), 
          axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
    guides(color = guide_legend(override.aes = list(size = 4)))
}


### make standard umaps 
plot_sc_umap_condition <- function(dat, facet_on_rep = FALSE){
  if(facet_on_rep==TRUE){
    return(ggumap(get_embedding(mb231) %>% filter(replicate !=0)) + 
             geom_point(aes(color=as.factor(condition)), alpha=0.8, size=0.2) + 
             theme_void() + labs(color='cluster') + theme(aspect.ratio = 0.8) + 
             scale_color_manual(values=c('fusion'='goldenrod','ctrl'='dodgerblue','presort'='grey80', 'initial'='grey40')) + 
             facet_wrap(~replicate))
  }
  
  if(facet_on_rep==FALSE){
    return(ggumap(get_embedding(mb231)) + 
             geom_point(aes(color=as.factor(condition)), alpha=0.8,size=0.2) + 
             labs(color='cluster') + theme(aspect.ratio = 0.8) + 
             scale_color_manual(values=c('fusion'='goldenrod','ctrl'='dodgerblue','presort'='grey80', 'initial'='grey40')) )
  }
}


### run NEBULA for DE
## function to run NEBULA mixed effects model for DE analysis
library(nebula)
run_nebula_de <- function(dat, 
                          method = 'LN', 
                         min.pct = 0.1, 
                         n_cores = 12, 
                         min.cpc = 0.005,
                             sid = 'barpairid',
                           group = 'cell_class'
                          ){
  
  neb <- scToNeb(obj = dat, 
               assay = "RNA", 
                  id = sid, 
                pred = group, 
              offset = "nCount_RNA"
              )
  
  neb$id[is.na(neb$id)==TRUE] <- 'other'
  
  model.mat = model.matrix(~cell_class, data=neb$pred)
  
  grouped.neb = group_cell(count = neb$count,
                              id = neb$id,
                            pred = model.mat,
                          offset = neb$offset
                          )
  
  neb_out = nebula(count = grouped.neb$count,
                      id = grouped.neb$id,
                    pred = grouped.neb$pred,
                  offset = grouped.neb$offset,
                  method = method,
                   ncore = n_cores,
                   mincp = round(min.pct*nrow(model.mat)),
                     cpc = min.cpc
                   )
  return(neb_out)
}


#### extract Pearson residuals from Nebula out
get_pearson_nebula <- function(neb_res,
                               dat, 
                          method = 'LN', 
                          min.pct = 0.1, 
                          n_cores = 12, 
                          min.cpc = 0.005,
                          sid = 'barpairid',
                          group = 'cell_class'
){
  
  neb <- scToNeb(obj = dat, 
                 assay = "RNA", 
                 id = sid, 
                 pred = group, 
                 offset = "nCount_RNA"
  )
  
  neb$id[is.na(neb$id)==TRUE] <- 'other'
  
  model.mat = model.matrix(~cell_class, data=neb$pred)
  
  grouped.neb = group_cell(count = neb$count,
                           id = neb$id,
                           pred = model.mat,
                           offset = neb$offset
  )
  
  pres = nbresidual(neb_res,
                   count = grouped.neb$count,
                   id = grouped.neb$id,
                   pred = grouped.neb$pred,
                   offset = grouped.neb$offset)
  return(pres)
}

##### Format NEBULA outputs for downstream analysis
format_neb_outs <- function(neb_out, testgroup = 'fusion'){
  neb_reformat <- neb_out[[1]] %>% 
    dplyr::rename('logFC'=paste0('logFC_cell_class',testgroup), 'p_val'=paste0('p_cell_class',testgroup)) %>% 
    mutate(log2FC = logFC/log(2)) %>% # convert from log(e) to log(2)
    dplyr::select(gene, log2FC, p_val) %>% 
    filter(gene %notin% c('EGFP','MCHERRY')) %>% # not interested in DE of reporter genes
    mutate(p_val_adj = p.adjust(p_val, method = "BH")) # adjust p values for multiple hypothesis testing
  return(neb_reformat)
}






###### create pseudobulked data
create_pseduo_bulk<- function(obj, min.pct=0.1, min_cells_per_clone = 3){
  ### subset data to ensure minimum clonal representation prior to bulking
  # find clones
  min_cell_clones <- obj@meta.data %>% group_by(cell_class, barpairid) %>% 
    summarize(clone_count = n()) %>% filter(clone_count >= min_cells_per_clone) %>% pull(barpairid)
  # subset data
  obj <- subset(obj, subset = barpairid %in% min_cell_clones)
  
  ## filter out genes not expressed in at least 10% of cells
  allcounts <- obj@assays$RNA$counts
  n_cells <- ncol(allcounts)
  
  ## count the number of cells each gene is found in with at least 1 count
  geneTotals <- rowSums(allcounts > 0)
  
  ## identify genes found in at least 10% of cells
  keep_genes <- names(geneTotals[geneTotals > round(min.pct*n_cells)])
  
  filtered_allcounts <- allcounts[rownames(allcounts) %in% keep_genes, ]
  
  data_collapsed <- presto::collapse_counts(filtered_allcounts, 
                                            obj@meta.data, 
                                            c('barpairid','cell_class'))
  meta_data<- data_collapsed$meta_data 
  
  mat<- data_collapsed$counts_mat
  
  colnames(mat)<- paste0(meta_data$barpairid)#,'-',colnames(mat))
  return(list(mat = mat, meta_data = meta_data))
}


########## Run DESeq2 on pseudobulk
library(DESeq2)
run_deseq2_ps <- function(ps, min_total_gene_count = 500, ref_group = "control"){
  counts <- ps[[1]]
  meta <- ps[[2]]
  
  # Remove genes with low total counts across all pseudobulked cells
  keep_genes <- names(rowSums(counts)[rowSums(counts) > min_total_gene_count])
  filtered_counts <- counts[rownames(counts) %in% keep_genes, ]
  
  dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                                colData = meta,
                                design = ~cell_class)
  
  dds$cell_class <- relevel(dds$cell_class, ref = ref_group)
  
  dds <- DESeq(dds)
  res <- results(dds)
  
  test_obj <- CreateSeuratObject(counts(dds, normalized=TRUE))
  test_obj@meta.data <- left_join(test_obj@meta.data %>% rownames_to_column('barpairid'), meta) %>% column_to_rownames('barpairid')
  Idents(test_obj) <- 'cell_class'
  #Idents(test_obj) <- factor(Idents(test_obj), levels = c('not_parent','parent'))
  return(list(test_obj, dds, res))
}


#######################
# Volcano Plot
library(ggrepel)
DE_volcano_pick_geneset <- function(res = fuse_markers, 
                                    test = 'nebula',
                                    genes = c('GAPDH','ESAM'),
                                    n_genes = Inf,
                                    title=paste0(''), 
                                    colrs=c('blue','red','black'),
                                    set_max = 300)
{
  
  if(test == 'nebula'){
    plot_volc <- res %>% dplyr::select(gene, avg_logFC, p_val_adj) %>% 
      mutate(of_interest = case_when(gene %in% genes ~ 'test', 
                                     .default = "other" )) %>% 
      mutate(neg.log10pval = -log10(p_val_adj)) %>% 
      mutate(neg.log10pval_adj = case_when(neg.log10pval == Inf ~ set_max,
                                           neg.log10pval != Inf ~ neg.log10pval))
    
    volc <- ggplot(plot_volc %>% filter(of_interest!='test'), aes(x=avg_logFC, y=neg.log10pval_adj)) + 
      geom_point(alpha=0.5, size=1, color='grey',   data=plot_volc) +
      geom_point(size=1, color=colrs[1], alpha=1, data=plot_volc %>% filter(of_interest=='test', avg_logFC < 0)) +
      geom_point(size=1, color=colrs[2], alpha=1, data=plot_volc %>% filter(of_interest=='test', avg_logFC > 0)) +
      theme_bw() + 
      ggtitle(title) + xlab('logFC') + ylab('-log10(padj)')
    
    if(length(plot_volc$of_interest[plot_volc$of_interest=='test'] <= n_genes)){
      markgenes_hi <- plot_volc %>% filter(of_interest=='test', avg_logFC > 0) %>% 
        mutate(score = avg_logFC*-log10(p_val_adj)) %>% arrange(desc(score)) %>% 
        filter(score > 0)  %>% head(n_genes) %>% pull(gene)
      
      markgenes_lo <- plot_volc %>% filter(of_interest=='test', avg_logFC < 0) %>% 
        mutate(score = avg_logFC*-log10(p_val_adj)) %>% arrange((score)) %>% 
        filter(score < 0)  %>% head(n_genes) %>% pull(gene)
      
      volc <- volc + geom_text_repel(data=plot_volc %>% filter(gene %in% unique(c(markgenes_hi, markgenes_lo))),
                                     bg.color = "white", bg.r = 0.05, box.padding = 0.2,
                                     max.overlaps = Inf, aes(label=gene), size=2.5) 
      return(volc)
    }
  }
  
  if(test == 'pseudo'){
    plot_volc <- res %>% dplyr::select(log2FoldChange, stat, pvalue, padj) %>% 
      rownames_to_column('gene') %>% 
      mutate(of_interest = case_when(gene %in% genes ~ 'test', 
                                     .default = "other" )) %>% 
      mutate(neg.log10pval = -log10(padj)) %>% 
      mutate(neg.log10pval_adj = case_when(neg.log10pval == Inf ~ set_max,
                                           neg.log10pval != Inf ~ neg.log10pval))
    
    volc <- ggplot(plot_volc %>% filter(of_interest!='test'), aes(x=log2FoldChange, y=neg.log10pval_adj)) + 
      geom_point(alpha=0.5, size=1, color='grey',   data=plot_volc) +
      geom_point(size=1, color=colrs[1], alpha=1, data=plot_volc %>% filter(of_interest=='test', log2FoldChange < 0)) +
      geom_point(size=1, color=colrs[2], alpha=1, data=plot_volc %>% filter(of_interest=='test', log2FoldChange > 0)) +
      theme_bw() + 
      ggtitle(title) + xlab('log2FoldChange') + ylab('-log10(padj)')
    
    if(length(plot_volc$of_interest[plot_volc$of_interest=='test'] <= n_genes)){
      markgenes_hi <- plot_volc %>% filter(of_interest=='test', log2FoldChange > 0) %>% arrange(desc(stat)) %>% 
        filter(stat > 0)  %>% head(n_genes) %>% pull(gene)
      
      markgenes_lo <- plot_volc %>% filter(of_interest=='test', log2FoldChange < 0) %>% arrange((stat)) %>% 
        filter(stat < 0)  %>% head(n_genes) %>% pull(gene)
      
      volc <- volc + geom_text_repel(data=plot_volc %>% filter(gene %in% unique(c(markgenes_hi, markgenes_lo))),
                                     bg.color = "white", bg.r = 0.05, box.padding = 0.2,
                                     max.overlaps = Inf, aes(label=gene), size=2.5) 
      return(volc)
    }
  }
  
  if(test == 'standard'){
    plot_volc <- res %>% dplyr::select(avg_log2FC, p_val_adj, p_val) %>% rownames_to_column('gene') %>% 
      mutate(of_interest = case_when(gene %in% genes ~ 'test', 
                                     .default = "other" )) %>% 
      mutate(neg.log10pval = -log10(p_val_adj)) %>% 
      mutate(neg.log10pval_adj = case_when(neg.log10pval == Inf ~ set_max,
                                           neg.log10pval != Inf ~ neg.log10pval))
    
    volc <- ggplot(plot_volc %>% filter(of_interest!='test'), aes(x=avg_log2FC, y=neg.log10pval_adj)) + 
      geom_point(alpha=0.5, size=1, color='grey',   data=plot_volc) +
      geom_point(size=1, color=colrs[1], alpha=1, data=plot_volc %>% filter(of_interest=='test', avg_log2FC < 0)) +
      geom_point(size=1, color=colrs[2], alpha=1, data=plot_volc %>% filter(of_interest=='test', avg_log2FC > 0)) +
      theme_bw() + 
      ggtitle(paste(as.character(title[1]))) + xlab('avg(log2FC)') + ylab('-log10(padj)')
    
    if(length(plot_volc$of_interest[plot_volc$of_interest=='test'] <= n_genes)){
      markgenes_hi <- plot_volc %>% filter(of_interest=='test', avg_log2FC > 0) %>% mutate(score = avg_log2FC*-log10(p_val)) %>% arrange(desc(score)) %>% 
        filter(score > 0)  %>% head(n_genes) %>% pull(gene)
      markgenes_lo <- plot_volc %>% filter(of_interest=='test', avg_log2FC < 0) %>% mutate(score = avg_log2FC*-log10(p_val)) %>% arrange((score)) %>% 
        filter(score < 0)  %>% head(n_genes) %>% pull(gene)
      volc <- volc + geom_text_repel(data=plot_volc %>% filter(gene %in% unique(c(markgenes_hi, markgenes_lo))),
                                     bg.color = "white", bg.r = 0.05, box.padding = 0.2,
                                     max.overlaps = Inf, aes(label=gene), size=2.5) 
      return(volc)
    }
  }
}











