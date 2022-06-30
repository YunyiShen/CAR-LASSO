library("rasterVis")
library("raster")
library(ggsn)
library(ggplot2)
library(maptools)
library(tmap)
library(extrafont)

generate_raster <- function(f){
  temp <- read.csv(f,row.names = 1)/50
  raster(as.matrix(temp))
}

generate_plot <- function(thestack,name){
  names(thestack) <- paste0("model", 1:5)
  gplot(thestack) +
    geom_tile(aes(fill = value), color="grey65") +
    facet_grid(~ variable) +
    scale_fill_gradient2(#low = "#0339ce",
      limits = c(-1,1),
      #high =  "#fdcf18",
      #low = "blue",
      low = "#0072B2",
      mid = "white",
      #mid = "grey95",
      #high = "red",
      high = "#990000",
      midpoint = 0,
      na.value="transparent") +
    labs(fill = "value", x = "", y=name) +
    theme(text = element_text(size=20), 
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(.15, .15, .15, .15, "cm"),
          #legend.position = "bottom",
          legend.title = element_blank())+
    theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      legend.background = element_rect(fill = "transparent")#, # get rid of legend bg
      #legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )
  
}


all_resCAR_ACAR <- list.files("./tests/Formal/Accurancy/k10/results/graph_visual/Omega",
                      pattern = "CARO_beta_s0.5_design_p10_n50",
                      full.names = T)
all_resCAR <- all_resCAR_ACAR[1:5+6]
all_resACAR <- all_resCAR_ACAR[1:5]

mat_CAR <- lapply(all_resCAR, generate_raster)
stack_CAR <- Reduce(stack, mat_CAR)
mat_ACAR <- lapply(all_resACAR, generate_raster)
stack_ACAR <- Reduce(stack, mat_ACAR)

all_resSRG <- list.files("./tests/Formal/Accurancy/k10/results/graph_visual/Omega",
                              pattern = "SRGO1_beta_s0.5_design_p10_n50",
                              full.names = T)[1:5]
mat_SRG <- lapply(all_resSRG, generate_raster)
stack_SRG <- Reduce(stack, mat_SRG)


all_resGlasso <- list.files("./tests/Formal/Accurancy/k10/results/graph_visual/Omega",
                         pattern = "GlassoaugO_beta_s0.5_design_p10_n50",
                         full.names = T)[1:5+6]
mat_Glasso <- lapply(all_resGlasso, generate_raster)
stack_Glasso <- Reduce(stack, mat_Glasso)

all_resAGlasso <- list.files("./tests/Formal/Accurancy/k10/results/graph_visual/Omega",
                            pattern = "GlassoO_beta_s0.5_design_p10_n50",
                            full.names = T)[1:5]
mat_AGlasso <- lapply(all_resAGlasso, generate_raster)
stack_AGlasso <- Reduce(stack, mat_AGlasso)

all_resmultireg <- list.files("./tests/Formal/Accurancy/k10/results/graph_visual/Omega",
                            pattern = "multiregO_beta_s0.5_design_p10_n50",
                            full.names = T)[1:5]
mat_multireg <- lapply(all_resmultireg, generate_raster)
stack_multireg <- Reduce(stack, mat_multireg)


ACAR_res <- generate_plot(stack_ACAR, "CG-A")
CAR_res <- generate_plot(stack_CAR, "CG")
SRG_res <- generate_plot(stack_SRG, "SRG")
multireg_res <- generate_plot(stack_multireg,"multireg")
Glasso_res <- generate_plot(stack_Glasso,"Glasso-aug")
AGlasso_res <- generate_plot(stack_AGlasso,"AGlasso")


ggpubr::ggarrange(ACAR_res, CAR_res, SRG_res, Glasso_res, multireg_res,
                  label.x = 0,
                  ncol = 1,
                  legend = "right",align = "hv",
                  common.legend = T)

ggsave("./tests/Formal/Accurancy/Figs/graph_reconstruct.5_10.pdf",width = 10, height = 11)
ggsave("./tests/Formal/Accurancy/Figs/graph_reconstruct.5_10.jpg",width = 10, height = 11)
