library(CARlasso)
library(ggplot2)
library(hexSticker)

example(CARlasso)

car_res$nodes$predictors <- rep("", 5)
car_res$nodes$response <- rep("", 5)
p <- plot(car_res)
pp <- p + theme(legend.position = "none") +
  theme(
  panel.background = element_rect(fill = "transparent",color = NA), # bg of the panel
  plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
  panel.grid.major = element_blank(), # get rid of major grid
  panel.grid.minor = element_blank(), # get rid of minor grid
  legend.background = element_rect(fill = "transparent"), # get rid of legend bg
  legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
  legend.key = element_rect(fill = "transparent", colour = NA), # get rid of key legend fill, and of the surrounding
  axis.line = element_line(colour = NA) # adding a black line for x and y axis
)


ww <- sticker(pp, package = "CARlasso", s_width = 1.83, s_height=1.83, s_x=0.96, s_y=0.7, h_fill = "#4e8ff0", p_y = 1.45)
ww

ggsave("CARlasso.png",width = 1,height = 1, unit = "in", dpi = 500, scale = 1)

sticker(pp, package = "CARlasso", s_width = 2.3, s_height=2.3, s_x=0.96, s_y=0.65, h_fill = "#4e8ff0", p_y = 1.5, p_size = 15.5) 


## New sticker:
# CAR-LASSO
library(CARlasso)
library(ggplot2)
library(hexSticker)

#example(CARlasso)
set.seed(1234)
dt <- simu_AR1(n=100,k=4, rho=0.7)
car_res <- CARlasso(y1+y2+y3+y4~x1+x2+x3, data = dt, adaptive = TRUE)

car_res$nodes$predictors <- rep("", 3)
car_res$nodes$response <- rep("", 4)
p <- plot(car_res)
pp <- p + theme(legend.position = "none") +
  theme(
    panel.background = element_rect(fill = "transparent",color = NA), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.key = element_rect(fill = "transparent", colour = NA), # get rid of key legend fill, and of the surrounding
    axis.line = element_line(colour = NA) # adding a black line for x and y axis
  )


ww <- sticker(pp, package = "CARlasso", s_width = 1.83, s_height=1.83, s_x=1.05, s_y=0.8, h_fill = "#edeff2", p_y = 1.6, p_color="black", p_size=15, h_color="#9d9d9e")
ww