generate_plasmid <- function(plasmid_size,
                             plasmid_name = "",
                             name_location = 0,
                             size_factor = 100000,
                             KPC_start = NULL,
                             Tn_IS_start = NULL,
                             Tn_IS_end = NULL,
                             backbone_color = "black",
                             Tn_IS_color = "red",
                             KPC_color = "blue",
                             plasmid_width_size = 4,
                             text_size = 4,
                             add_plasmid_size = NULL,
                             add_gene_locations = NULL,
                             add_gene_colors = "green"
                             ) {


  #Set size ratio for plasmid (eg. 100kb plasmid will have radius of 1)
  plasmid_radius <- plasmid_size / size_factor
  
  #KPC location
  denominator_KPC <- (KPC_start / plasmid_size) + 0.75 #0.75 so that origin starts at top of plasmid
  x_KPC_location <- plasmid_radius * cos((-2 *pi) * denominator_KPC) 
  y_KPC_location <- plasmid_radius * sin((-2 *pi) * denominator_KPC)

  #Tn/IS location
  denominator_TN_IS_start <- 1 / (Tn_IS_start / plasmid_size)
  denominator_TN_IS_end <- 1 / (Tn_IS_end / plasmid_size)
  Tn_IS_arc_start <- (2 * pi) / denominator_TN_IS_start
  Tn_IS_arc_end <- (2 * pi) / denominator_TN_IS_end
  
  #Create backbone
  backbone <- ggplot() + geom_circle(aes(x0=0, y0=0, r=plasmid_radius), 
                                     size = plasmid_width_size, 
                                     color=backbone_color) + coord_fixed()
  
  #Set axis lengths
  backbone <- backbone + xlim(-4, 4) + ylim(-4,4)
  
  #Draw arcs for Tn/IS elements
  if (!is.null(Tn_IS_start) && !is.null(Tn_IS_end)) {  #Check if argument is given
    backbone <- backbone + 
    geom_arc(aes(x0=0, y0=0, r=plasmid_radius, start=Tn_IS_arc_start, end=Tn_IS_arc_end), 
    lineend = "round", 
    color = Tn_IS_color,
    size=5)
  }
  
  #Draw KPC gene circle
  if (!is.null(KPC_start)) {  #Check if argument is given
    backbone <- backbone + 
    geom_point(aes(x=x_KPC_location, y=y_KPC_location),
               shape=21,
               size=6, 
               color="black", #color of outer ring
               fill=KPC_color,
               stroke=1) #thickness of outer ring
  }
  
  #Add additional gene circles
  if (!is.null(add_gene_locations)){
    gene_circles_x <- plasmid_radius * cos((-2 * pi) * ((add_gene_locations / plasmid_size) + 0.75))
    gene_circles_y <- plasmid_radius * sin((-2 * pi) * ((add_gene_locations / plasmid_size) + 0.75))
    backbone <- backbone + geom_point(aes(x=c(gene_circles_x), y=c(gene_circles_y)), shape = 21, size = 6, color = "black", fill = add_gene_colors, stroke = 1)
  }
  
  
  
  #Add plasmid name
  backbone <- backbone + 
                      annotate("text",
                      x=0, y=(-0.3 + name_location - plasmid_radius), #location centered below plasmid
                      label=plasmid_name,
                      size = text_size)
  
  #Add plasmid size in center (in kb)
  if (!is.null(add_plasmid_size)) {
    size_in_kb <- format(round((plasmid_size / 1000), 1), nsmall =1)
    backbone <- backbone + 
    annotate("text", x=0, y=0,
    label = paste(toString(size_in_kb), "kb", sep = " "),
    size = text_size)
      
  }
  
  #remove background
  backbone <-  backbone + theme_void()

  #Draw final image
  backbone
}
