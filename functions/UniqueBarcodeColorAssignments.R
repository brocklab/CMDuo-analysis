###########################################################
### Function to assign unique colors to each barcode
# green barcodes will be green-ish
# red barcodes will be red-ish
# fusion clones will be the average of their parental clone colors mixed with yellow
give_color_to_each_barcode <- function(obj){
  green_clones <- unique(obj@meta.data$gbarid[obj@meta.data$n_bar_cell >=1 & obj@meta.data$clone_color =='green']) %>% sort()
  red_clones <- unique(obj@meta.data$rbarid[obj@meta.data$n_bar_cell >=1 & obj@meta.data$clone_color =='red']) %>% sort()
  nd_clones <- unique(obj@meta.data$ubarid[obj@meta.data$n_bar_cell >=1 & is.na(obj@meta.data$clone_color) ==TRUE]) %>% sort()
  
  # Generate colors
  set.seed(8786)
  green_colors <- colorRampPalette(c("#d7fada","#90c447","#9ab39b","forestgreen","#0a4543"))(length(green_clones))
  pick_greens <- sample(green_colors, length(green_clones), replace=FALSE)
  
  set.seed(1924)
  red_colors <- colorRampPalette(c("#ffccd4","#b05180","#d96f68","red3","#5c0600","#bf5f3f","#4d1036"))(length(red_clones))
  red_colors <- colorRampPalette(c("#ffddd4","#c7674a","#d96f68","red3","#47040d"))(length(red_clones))
  pick_reds <- sample(red_colors, length(red_clones), replace=FALSE)
  
  nd_colors <- colorRampPalette(c("grey80","grey20"))(length(nd_clones))
  pick_nds <- sample(nd_colors, length(nd_clones), replace=FALSE)
  
  gmap <- data.frame(gbarid = green_clones, gbarid_color = pick_greens)
  rmap <- data.frame(rbarid = red_clones,   rbarid_color = pick_reds)
  umap <- data.frame(ubarid = nd_clones, ubarid_color = pick_nds)

  
  if(nrow(umap)>0){
    assign_colors <- obj@meta.data %>% filter(n_bar_cell >=1) %>% 
      dplyr::select(gbarid, rbarid, ubarid, barpairid) %>% 
      left_join(gmap, by='gbarid') %>% 
      left_join(rmap, by='rbarid') %>% 
      left_join(umap, by='ubarid') %>% 
      rowwise() %>% mutate(barpairid_color = average_colors(gbarid_color, rbarid_color)) %>% 
      mutate(barpairid_color = average_colors(barpairid_color, ubarid_color)) %>% ungroup()
  }
  
  if(nrow(umap)==0){
    assign_colors <- obj@meta.data %>% filter(n_bar_cell >=1) %>% 
      dplyr::select(gbarid, rbarid, barpairid) %>% 
      left_join(gmap, by='gbarid') %>% 
      left_join(rmap, by='rbarid') %>% 
      rowwise() %>% mutate(barpairid_color = average_colors(gbarid_color, rbarid_color)) %>% ungroup()
  }
  
  color_df <- assign_colors %>% select(barpairid, barpairid_color) %>% distinct()
  color_map <- color_df$barpairid_color
  names(color_map) <- color_df$barpairid
  
  obj@meta.data <- obj@meta.data %>% rownames_to_column('cellid') %>% left_join(color_df, by='barpairid') %>% column_to_rownames('cellid')
  
  return(list(obj, color_df, color_map))
}


###########################################################
### Function to average two hex colors and add in yellow
average_colors <- function(color1, color2) {
  if(is.na(color1)==FALSE & is.na(color2)==FALSE){
    avg_hex <- NA
  }
  if(is.na(color1)==TRUE){
    avg_hex <- color2
  }
  if(is.na(color2)==TRUE){
    avg_hex <- color1
  }
  if(is.na(color1)==FALSE & is.na(color2)==FALSE){
    # Ensure the colors have the # prefix
    if (!startsWith(color1, "#")) color1 <- paste0("#", color1)
    if (!startsWith(color2, "#")) color2 <- paste0("#", color2)
    
    # Convert hex to RGB
    rgb1 <- col2rgb(color1)
    rgb2 <- col2rgb(color2)
    rgb3 <- col2rgb("#f7d84a")
    
    # Average the RGB components
    avg_rgb <- (rgb1 + rgb2 + rgb3) / 3 # make more yellow with rgb3
    
    # Convert averaged RGB back to hex
    avg_hex <- rgb(avg_rgb[1] / 255, avg_rgb[2] / 255, avg_rgb[3] / 255, maxColorValue = 1)
  }
  return(avg_hex)
}