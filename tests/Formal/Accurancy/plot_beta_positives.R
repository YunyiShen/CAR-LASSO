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
  names(thestack) <- c("AR1","AR2","Block","Star","Circle")
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

s <- 0.2
p <- 5

all_resCAR_ACAR <- list.files("./tests/Formal/Accurancy/k10/results/graph_visual/B",
                      pattern = paste0( "CARB_beta_s",s,"_design_p",p,"_n50"),
                      full.names = T)
all_resCAR <- all_resCAR_ACAR[1:6+6]
all_resACAR <- all_resCAR_ACAR[1:6]

mat_CAR <- lapply(all_resCAR, generate_raster)
stack_CAR <- Reduce(stack, mat_CAR)
mat_ACAR <- lapply(all_resACAR, generate_raster)
stack_ACAR <- Reduce(stack, mat_ACAR)

all_resSRG <- list.files("./tests/Formal/Accurancy/k10/results/graph_visual/B",
                              pattern = paste0("SRGB_beta_s",s,"_design_p",p,"_n50"),
                              full.names = T)[1:6]
mat_SRG <- lapply(all_resSRG, generate_raster)
stack_SRG <- Reduce(stack, mat_SRG)



all_resmultireg <- list.files("./tests/Formal/Accurancy/k10/results/graph_visual/B",
                            pattern = paste0("multiregB_beta_s",s,"_design_p",p,"_n50"),
                            full.names = T)[1:6]
mat_multireg <- lapply(all_resmultireg, generate_raster)
stack_multireg <- Reduce(stack, mat_multireg)


ACAR_res <- generate_plot(stack_ACAR, "CG-A")
CAR_res <- generate_plot(stack_CAR, "CG")
SRG_res <- generate_plot(stack_SRG, "SRG")
multireg_res <- generate_plot(stack_multireg,"multireg")



ggpubr::ggarrange(ACAR_res, CAR_res, SRG_res, multireg_res,
                  label.x = 0,
                  ncol = 1,
                  legend = "right",align = "hv",
                  common.legend = T)

ggsave(paste0("./tests/Formal/Accurancy/Figs/beta_reconstruct",sub("0","",s),"_",p,".pdf"),width = 10, height = 9)
#ggsave("./tests/Formal/Accurancy/Figs/beta_reconstruct.5_5.jpg",width = 10, height = 9)
