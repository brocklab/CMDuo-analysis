### get data embedding from seurat
get_embedding <- function(plot_dat){
  umapcell <- plot_dat@reductions$umap@cell.embeddings %>% data.frame()
  if('cellid' %in% colnames(plot_dat@meta.data)==FALSE){
    embedding <- left_join(plot_dat@meta.data %>% rownames_to_column('cellid'), umapcell %>% rownames_to_column('cellid'), by='cellid')
  }
  if('cellid' %in% colnames(plot_dat@meta.data)==TRUE){
    embedding <- left_join(plot_dat@meta.data, umapcell %>% rownames_to_column('cellid'))
  }
  if('umap_1' %in% colnames(embedding)==TRUE){
    embedding <- embedding %>% dplyr::rename('UMAP_1'='umap_1', 'UMAP_2'='umap_2')
  }
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



### save figure function
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {png(paste0(fig_path, name, ".", type),width = width, height = height, units = "in", res = res)} 
  else {pdf(paste0(fig_path, name, ".", type),width = width, height = height)}
  print(plots)
  dev.off()
}