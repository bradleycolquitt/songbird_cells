get_aba_panel_ids <- function() {
  get_ids <- "http://api.brain-map.org/api/v2/data/query.csv?criteria=
     model::AtlasImage,
     rma::criteria,atlas_data_set(atlases[id$eq602630314]),graphic_objects(graphic_group_label[id$eq28]),
     rma::options[tabular$eq'sub_images.id'][order$eq'sub_images.id']
     &num_rows=all&start_row=0"
  
  get_ids <- gsub("[ \n]+","",get_ids)
  
  read.csv(url(get_ids))$id
}

get_aba_adult_mouse_ids = function() {
  get_ids = "http://api.brain-map.org/api/v2/data/query.xml?criteria=model::AtlasImage,
  rma::criteria,
  [annotated$eqtrue],
  atlas_data_set(atlases[id$eq1]),
  alternate_images[image_type$eq'Atlas+-+Adult+Mouse'],
  rma::options[order$eq'sub_images.section_number'][num_rows$eqall]"
  
  get_ids = "http://api.brain-map.org/api/v2/data/query.csv?criteria=
model::AtlasImage,
rma::criteria,
atlas_data_set(atlases[id$eq1]),
graphic_objects(graphic_group_label[id$eq28]),
rma::options[tabular$eq'sub_images.id'][order$eq'sub_images.id']
&num_rows=all&start_row=0"
  
  get_ids = gsub("[ \n]+", "", get_ids)
  read.csv(url(get_ids))$id
}

get_aba_panels <- function() {
  get_ids <- "http://api.brain-map.org/api/v2/data/query.csv?criteria=
     model::AtlasImage,
     rma::criteria,atlas_data_set(atlases[id$eq602630314]),graphic_objects(graphic_group_label[id$eq28]),
     rma::options[tabular$eq'sub_images.id'][order$eq'sub_images.id']
     &num_rows=all&start_row=0"
  
  get_ids <- gsub("[ \n]+","",get_ids)
  
  read.csv(url(get_ids))
}

get_aba_svg <- function(id,
                        downsample) {
  downsample <- as.character(downsample)
  
  id_urls <- paste0("http://api.brain-map.org/api/v2/svg/",
                   id,
                   "?downsample=",downsample)
  
  # id_urls = paste0("http://api.brain-map.org/api/v2/svg_download",
  #                  id,
  #                  "?downsample=",downsample)
  
  in_con <- curl(id_urls, open = "r")
  svg_lines <- suppressWarnings(readLines(in_con))
  close(in_con)
  
  return(svg_lines)
}

save_aba_svgs <- function(ids = NULL,
                          out_dir,
                          downsample = 4,
                          remove_colors = FALSE) {
  
  if(is.null(ids)) {
    ids <- get_aba_panel_ids()
  }
  
  for(i in 1:length(ids)) {
    print(ids[i])
    svg_lines <- get_aba_svg(ids[i], downsample=downsample)
    
    if(remove_colors) {
      svg_lines <- gsub("fill:#.{6}","fill:#ffffff",svg_lines)
    }
    
    out_con <- file(file.path(out_svg_dir, paste0(ids[i],".svg")), open = "w")
    writeLines(svg_lines, out_con)
    close(out_con)
  }
  
}

svg_to_tags <- function(x) {
  unlist(strsplit(x, "><"))
}

aba_svg_tags_to_list <- function(x) {
  split_on_space <- unlist(strsplit(x,"\" "))
  no_quotes <- gsub("\"","",split_on_space)
  no_brackets <- sub("<g |path ","",no_quotes)
  split_on_eq <- strsplit(no_brackets, "=")
  
  out_list <- map(split_on_eq,
                  function(z) {
                    if(length(z) == 2) {
                      z[2]
                    } else {
                      out <- NA
                    }
                  })
  names(out_list) <- map_chr(split_on_eq, 1)
  
  out_list
}

aba_svg_list_to_coords <- function(x) {
  if("d" %in% names(x)) {
    d <- sub("M ","",x$d)
    p <- unlist(strsplit(d, " L "))
    df <- map_dfr(p,
                  function(point) {
                    xy <- unlist(strsplit(point, ","))
                    data.frame(x = as.numeric(xy[1]),
                               y = as.numeric(xy[2]))
                  })
    
    # remove redundant points
    df <- df %>%
      filter(!(x == lag(x) & x == lead(x))) %>%
      filter(!(y == lag(y) & y == lead(y)))
    
    df
  } else {
    NULL
  }
}

aba_svg_list_to_attr <- function(x) {
  keep <- names(x) != "d"
  out_list <- x[keep]
  if("style" %in% names(out_list)) {
    style <- out_list$style
    style <- gsub("/","",style)
    styles <- unlist(strsplit(style, ";"))
    color <- sub("stroke:","",styles[grepl("stroke:",styles)])
    fill <- sub("fill:","",styles[grepl("fill:",styles)])
    out_list$color <- color
    out_list$fill <- fill
  }
  out_list
}

aba_svg_coords_to_segs <- function(df) {
  if(!is.null(df)) {
    
    points <- df
    
    segs <- data.frame(x = points$x,
                       y = points$y,
                       xend = lead(points$x),
                       yend = lead(points$y))
    
    segs <- segs[-nrow(segs),]
    
    last_seg <- data.frame(x = segs$xend[nrow(segs)],
                           y = segs$yend[nrow(segs)],
                           xend = segs$x[1],
                           yend = segs$y[1])
    
    segs <- rbind(segs, last_seg)
    segs
  }
}

plot_aba_svg_coords <- function(svg_coords,
                                svg_attr,
                                min_pts = 20) {
  # remove nulls
  keep_coords <- !map_lgl(svg_coords, is.null)
  svg_coords <- svg_coords[keep_coords]
  svg_attr <- svg_attr[keep_coords]
  
  keep_coords <- map_int(svg_coords, nrow) >= min_pts
  svg_coords <- svg_coords[keep_coords]
  svg_attr <- svg_attr[keep_coords]
  
  plot_data <- map2_dfr(svg_coords,
                        svg_attr,
                        function(x, y) {
                          out <- x
                          out$fill <- y$fill
                          out$color <- y$color
                          out$id <- y$id
                          out$order <- y$order
                          out
                        })
  
  plot_list <- split(plot_data, plot_data$order)
  
  p <- ggplot() +
    scale_fill_identity() +
    scale_color_identity() +
    scale_y_reverse() +
    theme_void()
  
  for(i in 1:length(plot_list)) {
    p <- p +
      geom_polygon(data = plot_list[[i]],
                   aes(x = x,
                       y = y,
                       group = id,
                       fill = fill,
                       color = color))
  }
  
  p
}

ish_slice_heatmap_flat <- function(mat,
                              anno = NULL,
                              taxon = NULL,
                              slice_num,
                              plane = "coronal",
                              normalize = "slice",
                              colorset = c("darkblue","gray90","red")) {
  
  library(rbokeh)
  library(dplyr)
  library(reshape2)
  library(scrattch.vis)
  
  slice_mat <- slice_ccf_arr(mat, slice_num, plane)
  
  if(normalize == "slice") {
    max_val <- max(slice_mat)
  } else if(normalize == "all") {
    max_val <- max(c(mat), na.rm=T)
  }
  max_val = round(max_val, 1)
  slice_flat <- reshape2::melt(slice_mat)
  names(slice_flat) <- c("y","x","value")
  
  slice_flat$value[slice_flat$value < 0] <- 0
  slice_flat$color <- scrattch.vis::values_to_colors(slice_flat$value,
                                                     colorset = colorset,
                                                     min_val = 0,
                                                     max_val = max_val)
  
  if(is.null(anno)) {
    hover_list <- list("Value" = "value")
  } else {
    anno_flat <- reshape2::melt(slice_ccf_arr(anno, slice_num, plane)) 
    names(anno_flat) <- c("y","x","id")
    anno_flat = anno_flat %>% left_join(taxon %>% select(id, acronym), by="id")
    
    slice_flat <- dplyr::left_join(slice_flat, anno_flat, by = c("x","y"))
    hover_list <- list("Value" = "value",
                       "Annotation" = "annotation")
  }
  
  gg = ggplot(slice_flat) + 
    theme_void() + 
    #scale_fill_identity() 
    scale_fill_gradient(low="grey90",
                        high=colorset[length(colorset)],
                        na.value="white",
                        breaks=c(0, max_val),
                        limits=c(0, max_val)) + 
    theme(legend.position = c(.95, .5))

  if(plane == "coronal") {

    gg = gg +
      geom_raster(aes(x=x, y=-y, fill=value)) 
      # rbokeh::ly_crect(data = slice_flat,
      #                  x = x,
      #                  y = -y,
      #                  fill_color = color,
      #                  line_color = NA,
      #                  fill_alpha = 1,
      #                  hover = hover_list)
  } else if(plane == "horizontal") {
    gg = gg +
      geom_raster(aes(x=y, y=x, fill=color, line=NA)) 
    # f <- rbokeh::figure(width = dim(slice_mat)[1]*10,
    #                     height = dim(slice_mat)[2]*10) %>%
    #   rbokeh::ly_crect(data = slice_flat,
    #                    x = y,
    #                    y = x,
    #                    fill_color = color,
    #                    line_color = NA,
    #                    fill_alpha = 1,
    #                    hover = hover_list)
  } else if(plane == "saggital") {
    gg = gg +
      geom_raster(aes(x=y, y=-x, fill=color, line=NA)) 
    # f <- rbokeh::figure(width = dim(slice_mat)[1]*10,
    #                     height = dim(slice_mat)[2]*10) %>%
    #   rbokeh::ly_crect(data = slice_flat,
    #                    x = y,
    #                    y = -x,
    #                    fill_color = color,
    #                    line_color = NA,
    #                    fill_alpha = 1,
    #                    hover = hover_list)
  }
  
  return(gg)
  
  
}