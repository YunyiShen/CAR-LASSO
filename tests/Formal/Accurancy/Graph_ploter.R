library("rasterVis")
library("raster")
library(ggsn)
library(ggplot2)
library(maptools)
library(tmap)
library(extrafont)

source("./tests/Formal/Accurancy/Graph_generator.R")

graphss <- lapply(1:6, function(i){
  temp <- do.call(paste0("g_model",i),list(k=10))
  raster(temp$Omega)
})

graphss <- Reduce(stack, graphss)
names(graphss) <- paste0("model",1:6)

gplot(graphss) +
  geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_gradient2(#low = "#0339ce",
                      #high =  "#fdcf18",
                      low = "red",
                      #mid = "white",
                      mid = "grey95",
                      high = "blue",
                      midpoint = 0,
                      na.value="transparent") +
  labs(fill = "value", x = "", y="") +
  theme(text = element_text(size=14), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(.15, .15, .15, .15, "cm"))+
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent")#, # get rid of legend bg
    #legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

ggsave("./simu_graphs.pdf",width = 10.5,height=7, scale = 0.8)
