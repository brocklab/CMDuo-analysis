find_marker_genes_ppf_lowhi <- function(obj_trio_min, fusepairs, n_genes=10, min.pct=0.8){
  
  stor <- data.frame()
  for(i in seq(1, nrow(fusepairs))){ 
    sub <- subset_parent_parent_fusion(obj_trio_min, pick_clones = fusepairs[i,], parent_pairs = fusepairs)
    
    sub@assays$RNA$scale.data <- NULL
    sub <- ScaleData(sub, vars.to.regress = c("S.Score", "G2M.Score"))
    
    Idents(sub) <- 'barpairid'
    
    subset_markers_p <- FindMarkers(sub,
                                    ident.1 = c(fusepairs[i,]$gbarid),
                                    ident.2 = c(fusepairs[i,]$rbarid),
                                    min.pct=min.pct,
                                    only.pos = FALSE) %>% arrange(desc(avg_log2FC))
    
    subset_markers_f <- FindMarkers(sub,
                                    ident.1 = c(fusepairs[i,]$barpairid),
                                    ident.2 = c(fusepairs[i,]$gbarid, fusepairs[i,]$rbarid),
                                    min.pct=min.pct,
                                    only.pos = FALSE) %>% arrange(desc(avg_log2FC))
    
    subset_markers_g <- subset_markers_p %>% filter(avg_log2FC > 0)
    subset_markers_r <- subset_markers_p %>% filter(avg_log2FC < 0) %>% mutate(avg_log2FC=-avg_log2FC)
    subset_markers_notf <- subset_markers_f %>% filter(avg_log2FC < 0) %>% mutate(avg_log2FC=-avg_log2FC)
    
    # exclude genes which get turned on in control cells in from cell sorting
    sort_signature <- read.csv('outdata/de_results/archive/sort_stress_signature_genes.csv') %>% pull(x)
    
    rank_genes <- function(subset_markers, clone_color = 'na'){
      ranked <- subset_markers %>% 
        rownames_to_column('gene') %>% 
        filter(gene %notin% sort_signature) %>% 
        mutate(score = avg_log2FC*-log10(p_val_adj)) %>% 
        filter(gene %notin% c('EGFP','MCHERRY')) %>% 
        arrange(desc(score)) %>% 
        mutate(gene_rank = 1:n()) %>% 
        dplyr::select(gene, gene_rank) %>% 
        mutate(clone_color = clone_color) 
      
      return(ranked)
    }
    
    ranked_g <- rank_genes(subset_markers_g, clone_color = 'green')
    ranked_r <- rank_genes(subset_markers_r, clone_color = 'red')
    ranked_f <- rank_genes(subset_markers_f, clone_color = 'fusion')
    ranked_nf <- rank_genes(subset_markers_notf, clone_color = 'not_fusion')
    
    clone_set_markers <- rbind(ranked_g, ranked_r, ranked_f, ranked_nf) %>% filter(gene_rank <= n_genes) %>% mutate(cluster = fusepairs[i,]$barpairid)
    stor <- rbind(stor, clone_set_markers)
  }
  
  return(stor)
}
