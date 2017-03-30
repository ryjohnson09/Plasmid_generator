library(ggforce)

generate_plasmid <- function(plasmid_size,
                             plasmid_size_location = 0,
                             plasmid_name = "",
                             name_location = 0,
                             size_factor = 100000,
                             KPC_start = NULL,
                             Tn_IS_start = NULL,
                             Tn_IS_end = NULL,
                             backbone_color = "black",
                             Tn_IS_color = "red",
                             Tn_IS_size = 5,
                             KPC_color = "blue",
                             plasmid_width_size = 4,
                             gene_circle_size = 6,
                             name_text_size = 4,
                             size_text_size = 4,
                             add_plasmid_size = NULL,
                             add_gene_locations = NULL,
                             add_gene_colors = "green",
                             add_arc_start = NULL,
                             add_arc_end = NULL,
                             add_arc_color = "purple",
                             add_arc_size = 5
                             ) {


  #Set size ratio for plasmid (eg. 100kb plasmid will have radius of 1)
  plasmid_radius <- plasmid_size / size_factor
  
  #Create backbone
  backbone <- ggplot() + 
    geom_circle(aes(x0=0, y0=0, r=plasmid_radius), 
    size = plasmid_width_size, 
    color=backbone_color) + coord_fixed()
  
  #Set axis lengths
  backbone <- backbone + xlim(-4, 4) + ylim(-4,4)
  
  #Draw arcs for Tn/IS elements
  if (!is.null(Tn_IS_start) && !is.null(Tn_IS_end)) {
    backbone <- backbone + 
    geom_arc(aes(x0=0, y0=0, r=plasmid_radius, 
    start=((2 * pi) * (Tn_IS_start/plasmid_size)), 
    end=((2 * pi) * (Tn_IS_end/plasmid_size))), 
    lineend = "round", 
    color = Tn_IS_color,
    size=Tn_IS_size)
  }
  
  #Draw any additional arcs
  if (!is.null(add_arc_start) && !is.null(add_arc_end)) {
    backbone <- backbone +
    geom_arc(aes(x0=0, y0=0, r=plasmid_radius, 
    start = ((2 * pi) * (add_arc_start/plasmid_size)),
    end = ((2 * pi) * (add_arc_end/plasmid_size))),
    lineend = "round",
    color = add_arc_color,
    size = add_arc_size)
  }
  
  #Add additional gene circles
  if (!is.null(add_gene_locations)){
    gene_circles_x <- plasmid_radius * cos((-2 * pi) * ((add_gene_locations / plasmid_size) + 0.75))
    gene_circles_y <- plasmid_radius * sin((-2 * pi) * ((add_gene_locations / plasmid_size) + 0.75))
    backbone <- backbone + 
    geom_point(aes(x=c(gene_circles_x), y=c(gene_circles_y)), 
    shape = 21, 
    size = gene_circle_size, 
    color = "black", 
    fill = add_gene_colors, 
    stroke = 1)
  }
  
  #Draw KPC gene circle
  if (!is.null(KPC_start)) {
    backbone <- backbone + 
    geom_point(aes(x=plasmid_radius * cos((-2 * pi) * ((KPC_start / plasmid_size) + 0.75)),
    y=plasmid_radius * sin((-2 * pi) * ((KPC_start / plasmid_size) + 0.75))),
    shape=21,
    size=gene_circle_size, 
    color="black",
    fill=KPC_color,
    stroke=1)
  }
  
  #Add plasmid name
  backbone <- backbone + 
    annotate("text",
    x=0, y=(-0.4 + name_location - plasmid_radius),
    label=plasmid_name,
    size = name_text_size)
  
  #Add plasmid size in center (in kb)
  if (!is.null(add_plasmid_size)) {
    size_in_kb <- format(round((plasmid_size / 1000), 1), nsmall =1)
    backbone <- backbone + 
    annotate("text", x=0, y=(0 + plasmid_size_location),
    label = paste(toString(size_in_kb), "kb", sep = " "),
    size = size_text_size)
  }
  
  #remove background
  backbone <-  backbone + theme_void()

  #Draw final image
  backbone
}
