## cluster colors
clust_colors <- c('skyblue4','plum3','navyblue','slateblue3','maroon1','#00c6f2','darkorchid2','darkcyan','darkorchid','tan1')
#clust_colors <- c('skyblue4','plum3','#00c6f2','slateblue3','maroon1','navyblue','darkorchid2','darkcyan','darkorchid','tan1')
names(clust_colors) <- LETTERS[1:10]

## condition colors
#condition_colors <- c('initial'='grey30','presort'='grey80','fusion'='goldenrod3','ctrl'='dodgerblue3')
condition_colors <- c('initial'='grey30','presort'='grey80','pre-sort'='grey80',
                     'fusion-enriched'='goldenrod3','fusion'='goldenrod3','fusion-\nenriched'='goldenrod3',
                     'ctrl'='dodgerblue3','control'='dodgerblue3')

# parent-progeny class colors
pp_colors  <- c("fusion"="#c27100","parent"="grey30", 'other'='grey90')
pp_colors2 <- c("progeny"="#c27100","parent"="grey30", 'other'='grey90')

## index colors
index_colors <- c('green'='forestgreen','red'='red4','fusion'='goldenrod',
                  'green \nmarker genes'='forestgreen','red \nmarker genes'='red3','fusion \nmarker genes'='goldenrod',
                  'green markers'='forestgreen','red markers'='red3','fusion markers'='goldenrod')

## clone colors
# stored in seurat object
get_clone_colors <- function(obj){
  colors_df <- obj@meta.data %>% 
    filter(n_bar_cell >= 1) %>% 
    dplyr::select(barpairid, barpairid_color) %>% 
    distinct()
  
  color_map <- colors_df$barpairid_color
  names(color_map) <- colors_df$barpairid
  return(color_map)
}


## cell line colors
cellline_colors <- c('MDA-MB-231'='#F8766D', 'HCC1806'='#00BFC4','mean' = 'black')



### save figure function
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {png(paste0(fig_path, name, ".", type),width = width, height = height, units = "in", res = res)} 
  else {pdf(paste0(fig_path, name, ".", type),width = width, height = height)}
  print(plots)
  dev.off()
}
