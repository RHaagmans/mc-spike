library(ggplot2)
library(ggpubr)

# Default color palettes
pal_cat <- c("#b10609","#157cba","#e4c458","#009a5b","#d6389f","#474747")
pal_cat_light <- c("#e32527","#1a9feb","#f9db76","#0fce81","#ef5bb7","#535353")
pal_cat_dark <- c("#850607","#095785","#b99e47","#027344","#892a6a","#323232") 

# Colour sequence generator
pal_seq <- function(...){return(colorRampPalette(c(...)))}
pal_seq1 <- colorRampPalette(c(pal_cat[1],pal_cat_light[3]))
pal_seq2 <- colorRampPalette(c(pal_cat[2],pal_cat_light[4]))

nice_uniform_labels <- function(n=8){
  return(
    theme(
      plot.title = element_text(size=n),
      axis.title = element_text(size=n),
      axis.text = element_text(size=n),
      legend.title = element_text(size=n),
      legend.text = element_text(size=n)
    )
  )
}

nice_legend_tight_y <- theme(
    legend.margin = margin(0,0,0,0),
    legend.box.spacing = unit(0.4, "line"),
    legend.box.margin = margin(0,0,0,0),
    legend.spacing.y = unit(0, "line")
  )


nice_theme <- theme(
  plot.title = element_blank(),
  legend.background = element_rect(fill = "white", linewidth = 0, colour = "white"),
  panel.grid = element_blank(),
  strip.placement = "outside",
  strip.background = element_rect(fill = "white", colour = NA),
)

nice_plot <- function(labsize=8){
  return(
    theme_bw() +
    nice_theme +
    labs_pubr(labsize)+
    nice_uniform_labels(labsize)
  )
}

# Scales
nice_scale_expansion <- expansion(add=0, mult=c(0,0.1))
nice_scale_expansion_zero <- expansion(add=0, mult=0)

nice_scale_x <- scale_x_continuous(expand = nice_scale_expansion)
nice_scale_y <- scale_y_continuous(expand = nice_scale_expansion)

# Factor ordering
base_sample_order <- c("S06", "S07", "S08", "SMC", "SBL")
spike_order <- c("NO", "LO", "HI")
method_order <- c("WTA2", "SISPA")

sample_order <- c(
  "S06-NO",
  "S06-LO",
  "S06-HI",
  "S07-NO",
  "S07-LO", 
  "S07-HI",
  "S08-NO",
  "S08-LO",
  "S08-HI",
  "SMC-LO",
  "SMC-HI", 
  "SBL-NO"
)

full_sample_order <- expand.grid(sample_order, method_order)
full_sample_order <-
  sprintf('%s-%s', full_sample_order[, 1], full_sample_order[, 2])

virus_order <- c(
  "T5", 
  "M13", 
  "P22", 
  "Det7", 
  "MHV", 
  "RV-A", 
  "BVDV"
)

segment_order <- c(
  "T5", 
  "M13", 
  "P22", 
  "Det7", 
  "MHV", 
  "RV-A_s1", 
  "RV-A_s2", 
  "RV-A_s3", 
  "RV-A_s4", 
  "RV-A_s5", 
  "RV-A_s6", 
  "RV-A_s7", 
  "RV-A_s8", 
  "RV-A_s9", 
  "RV-A_s10", 
  "RV-A_s11",
  "BVDV"
)

segment_to_virus <- c(
  "T5"="T5", 
  "M13"="M13", 
  "P22"="P22", 
  "Det7"="Det7", 
  "MHV"="MHV", 
  "RV-A_s1"="RV-A", 
  "RV-A_s2"="RV-A", 
  "RV-A_s3"="RV-A", 
  "RV-A_s4"="RV-A", 
  "RV-A_s5"="RV-A", 
  "RV-A_s6"="RV-A", 
  "RV-A_s7"="RV-A", 
  "RV-A_s8"="RV-A", 
  "RV-A_s9"="RV-A", 
  "RV-A_s10"="RV-A", 
  "RV-A_s11"="RV-A",
  "BVDV"="BVDV"
) %>% 
  enframe(name = "ref", value = "virus") %>% 
  mutate(
    ref = factor(ref, levels = segment_order), 
    virus = factor(virus, levels = virus_order)
  )

virus_host <- c(
  "T5" = "Prokaryotic", 
  "M13" = "Prokaryotic", 
  "P22" = "Prokaryotic", 
  "Det7" = "Prokaryotic", 
  "MHV" = "Eukaryotic", 
  "RV-A" = "Eukaryotic", 
  "BVDV" = "Eukaryotic"
)


base_sample_faecal <- c("S06", "S07", "S08")
base_sample_control <- c("SMC", "SBL")

# Factor Aesthetics
base_sample_palette <- c('#1b9e77',  '#d95f02',  '#7570b3',  '#e7298a', "grey23")
names(base_sample_palette) <- c(base_sample_faecal, base_sample_control)

sample_palette <- c(
  '#1b9e77',  '#d95E01',  '#7570b3',  '#ec53a1', 
  '#54B699','#E38641','#9894c6', '#f17eb9', 
  '#8DCEBB','#ECAE80','#bab7d9', 
  'grey23')

names(sample_palette) <- c(
  "S06-HI", "S07-HI","S08-HI","SMC-HI",
  "S06-LO", "S07-LO", "S08-LO","SMC-LO",
  "S06-NO", "S07-NO", "S08-NO",
  "SBL-NO"
)


virus_palette <- c(pal_cat, c("orange"))
names(virus_palette) <- virus_order

method_palette <- c("SISPA"="#BC3C29FF", "WTA2" = "#0072B5FF")

checkv_quality_palette <- c(pal_cat[5], pal_cat[2], pal_cat[4], pal_cat[3], pal_cat[6])
names(checkv_quality_palette) <- c(
  "Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"
)

base_sample_shapes <- c("SMC" = 16, "S06"=0, "S07"=1, "S08"=2)
spike_shapes <- c("NO"=1, "LO"=16,"HI"=17)

assembly_qc_stage_colors <- c(
  "duplicate" = pal_cat[6],
  "unique_nonviral" = pal_cat[1],
  "unique_viral" = pal_cat[4]
)



contig_classifications = c(
  "unsure" = "Unclassified",
  "endo" = "Faecal viruses",
  "mc" = "Mock community viruses"
  )

virclass_palette <- c(pal_cat[4],pal_cat[1],pal_cat[6])
names(virclass_palette) <- c(
  contig_classifications[["endo"]],
  contig_classifications[["mc"]],
  contig_classifications[["unsure"]]
)

# Color schemes
colorscheme_base_sample_color <- scale_color_manual(values = base_sample_palette)
colorscheme_base_sample_fill <- scale_fill_manual(values = base_sample_palette)

colorscheme_sample_color <- scale_color_manual(values = sample_palette)
colorscheme_sample_fill <- scale_fill_manual(values = sample_palette)

colorscheme_virus_color <- scale_color_manual(values = virus_palette)
colorscheme_virus_fill <- scale_fill_manual(values = virus_palette)

serial_comma <- function(vec){
  if(length(vec) <= 1){
    return(vec)
  }else if(length(vec) == 2){
    return(
      paste(vec[1], vec[2], sep = " and ")
      )
  }else{
    return(
      
      paste(paste(vec[1:length(vec)-1], collapse=", "), vec[length(vec)], sep = ", and ")
      
      )
  }
}

