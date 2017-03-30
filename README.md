# Plasmid Generator  ![pKPC272](https://ibin.co/3HMQZ4JVATnC.png)

## Introduction  
Generate basic plasmid maps with user-defined aesthetics. This was written in the R programming language. A step-by-step tutorial is provided and encompasses all attributes of the function.

### Install  
Install necessary packages (ggforce) and the function itself  

```{r}
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
```


### Generate basic plasmid backbone  
To generate a basic ring shape, simply provide the function with the plasmid size (in basepairs). In this example, the plasmid is 100,000 basepairs in length. You can also change the color of the backbone by modifying the `backbone_color` argument:   
```{r}
generate_plasmid(plasmid_size = 100000)
generate_plasmid(plasmid_size = 100000, backbone_color = "red")
```

This is not very interesting, but it will provide the scaffold for the rest of the program. As a note, the program is scalled so that a plasmid of size **100,000**bp will have a radius of **1**. In other words, a plasmid of size 50,000bp will have a radius of 0.5 and a plasmid of 200,000bp will have a radius of 2. To show the size of the plasmid, set `add_plasmid_size` to `TRUE`:  

```{r}
generate_plasmid(plasmid_size = 100000, add_plasmid_size = TRUE)
```


Plasmids larger than 350,000bp may not render correctly. Similaryly, plasmids smaller than 50,000bp will be a bit cramped. To fix this, find a good middle ground for your plasmids of interest, and set that to `size_factor`. Example, if you want to plot 3 seprate plasmids of sizes 100kb, 60kb, and 40kb, try setting the `size_factor` to `40000`. Example:  
```{r}
generate_plasmid(plasmid_size = 100000, add_plasmid_size = TRUE)
generate_plasmid(plasmid_size = 60000, add_plasmid_size = TRUE)
generate_plasmid(plasmid_size = 40000, add_plasmid_size = TRUE)
```


Now let's change the `size_factor` to `40000`:

```{r}
generate_plasmid(plasmid_size = 100000, add_plasmid_size = TRUE, size_factor = 40000)
generate_plasmid(plasmid_size = 60000, add_plasmid_size = TRUE, size_factor = 40000)
generate_plasmid(plasmid_size = 40000, add_plasmid_size = TRUE, size_factor = 40000)
```

To change the size of the plasmid-size text, modify the `size_text_size` argument (default = 4):  
```{r}
generate_plasmid(plasmid_size = 100000, add_plasmid_size = TRUE, size_text_size = 6)
```


### Add label to plasmid  

This will add a label directly underneath the plasmid. To change the size of the label, modify the `name_text_size` argument (default = 4):    
```{r}
generate_plasmid(plasmid_size = 100000, plasmid_name = "pTEST-123")
generate_plasmid(plasmid_size = 100000, plasmid_name = "pTEST-123", name_text_size = 6)
```

If you need to move the plasmid name up or down, simply modify the `name_location` parameter. Positive values will move it up, negative values will move it down:  
```{r}
generate_plasmid(plasmid_size = 100000, plasmid_name = "pTEST-123", name_location = 1)
generate_plasmid(plasmid_size = 100000, plasmid_name = "pTEST-123", name_location = -1)
```


### Add KPC gene  

For this program, we are interested in plotting the location of the **KPC** gene. To do this, you must know the start site of the KPC gene on the plasmid sequence. For example, if the KPC gene starts at nucleotide position 25,000:  
```{r}
generate_plasmid(plasmid_size = 100000, KPC_start = 25000)
```

You may have already deteremined this, but nucleotide position 1 is the very top of the plasmid, and increases in a clock-wise manner. So, the KPC gene at position 25,000 (exactly 1/4 of the way around the plasmid), shows up on the right side of the plasmid.  

To change the color of the KPC gene, modify the `KPC_color` argument:  
```{r}
generate_plasmid(plasmid_size = 100000, KPC_start = 25000, KPC_color = "yellow")
```

### Add Tn or IS sequence  
In most cases, the KPC gene is within the context of a Transposon (Tn) or Insertion Sequence (IS). To add an "arc" that denotes this sequence, you must know the start and end site. Then plot in the following manner:  
```{r}
generate_plasmid(plasmid_size = 100000, KPC_start = 25000, Tn_IS_start = 20000, Tn_IS_end = 30000)
```

As with the KPC gene, you can change the color of this arc:  
```{r}
generate_plasmid(plasmid_size = 100000, KPC_start = 25000, Tn_IS_start = 20000, Tn_IS_end = 30000, Tn_IS_color = "green")
```

### Add additional genes of interest  
If there are other genes that you would like to plot (eg. other antibiotic resistance genes), you can pass a vector of start sites to the `add_gene_locations`. For example, you want to plot three additional genes at locations 30, 60000, and 75000:  
```{r}
generate_plasmid(plasmid_size = 100000, add_gene_locations = c(30, 60000, 75000))
```

Change the colors of the genes by passing an equally sized vector of colors to `add_gene_colors`:  
```{r}
generate_plasmid(plasmid_size = 100000, add_gene_locations = c(30, 60000, 75000), add_gene_colors = c("green", "red", "purple"))
```

### Bringing it all together  

As an example, let's graph the plasmid [pKPC-272](https://www.ncbi.nlm.nih.gov/nuccore/CP008825). Here's what we know about this plasmid:  

Element | size/location
--------|--------------
Plasmid Size | 282,439 bp
KPC location | 126,219 bp
Tn4401b start | 119,018 bp
Tn4401b end | 129,023 bp
Other abx genes | 116,356; 138,347 bp  

Let's make it!  
```{r}
generate_plasmid(plasmid_size = 282439, plasmid_name = "pKPC-272", add_plasmid_size = TRUE, KPC_start = 126219, Tn_IS_start = 119018, Tn_IS_end = 129023, add_gene_locations = c(116356, 138347))
```

That looks good, but let's increase the size of the text a bit to make it a bit more pleasing to the eye. We can also lower the label a bit so it's not so close to the plasmid itself:  

```{r}
generate_plasmid(plasmid_size = 282439, plasmid_name = "pKPC-272", add_plasmid_size = TRUE, KPC_start = 126219, Tn_IS_start = 119018, Tn_IS_end = 129023, add_gene_locations = c(116356, 138347), size_text_size = 8, name_text_size = 8, name_location = -0.2)
```

### Special cases  
#### Small Plasmids  
If you need to plot a very small plasmid (< 15 kb), here is work around.
```{r}
generate_plasmid(plasmid_size = 12000,
                 add_plasmid_size = TRUE,
                 plasmid_name = "pSMALL", #place size of plasmid under the name
                 plasmid_size_location = 0.4, #can raise/lower (positive or negative values, respectively) the plasmid size text location. This may require some playing around with to get it where you want.
                 plasmid_width_size = 0.5, #decrease width of the backbone (default = 4)
                 KPC_start = 10000,
                 Tn_IS_start = 9000,
                 Tn_IS_end = 11000,
                 Tn_IS_size = 0.6, #decrease width of Tn_IS arc, but make slightly larger than plasmid width
                 add_gene_locations = c(4000, 5000),
                 gene_circle_size = 0.7) #decrease the size of all gene circles
```

#### Need more than one arc 

If in addition to the Tn_IS element, you would like to add another "arc" to your plasmid, use the __add_arc__ arguments. It has all the same attributes as the __Tn_IS__ arguments (i.e. color, size).


```{r}
generate_plasmid(plasmid_size = 150000, Tn_IS_start = 20000, Tn_IS_end = 50000, add_arc_start = 60000, add_arc_end = 70000)
```


