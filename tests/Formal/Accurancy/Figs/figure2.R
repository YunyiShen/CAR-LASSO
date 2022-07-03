library(igraph)
library(RColorBrewer)

col_list <- rev(brewer.pal(n = 8, name = "RdYlBu")) # color scale

plot_Omega <- function(Omega, plot_title, color = TRUE){
  rescaled_Omega <- Omega <=0 #scales::rescale(Omega, 
                    #                from = c(-1,1)* max(abs(Omega)),
                    #                to = c(0,1))
  plot(1, type = "n", xlim = c(0.5,q+ 0.5), ylim = c(0.5,q+0.5), 
       main = plot_title, xlab = "", ylab = "",
       xaxt = "n", yaxt = "n", bty = "n")
  for(i in 1:q){
    for(j in 1:q){
      if(abs(Omega[i,j]) < 1e-12){
        rect(i-0.5, (q-j+1) - 0.5, i + 0.5, (q - j + 1) + 0.5,
             col = "white", border = "lightgray")
      } else{
        
        if(color){
          if(i!=j)
          rect(i-0.5, (q-j+1) - 0.5, i + 0.5, (q-j+1) + 0.5, 
               col = ifelse(rescaled_Omega[i,j], "black","gray"), #rgb(colorRamp(col_list, bias = 1)(rescaled_Omega[i,j])/255),
               border = "lightgray")
          else{
            rect(i-0.5, (q-j+1) - 0.5, i + 0.5, (q-j+1) + 0.5, 
                 col = "lightgray", #rgb(colorRamp(col_list, bias = 1)(rescaled_Omega[i,j])/255),
                 border = "lightgray")
          }
        } else{
          rect(i-0.5, (q-j+1) - 0.5, i + 0.5, (q-j+1) + 0.5, 
               col = "gray",
               border = "lightgray")
        }
      }
    }
  }
  rect(0.5, 0.5, q+0.5, q + 0.5)
}

plot_graph <- function(Omega, star_layout = FALSE){
  if(nrow(Omega) != 10) stop("hardcoded q = 10")
  A <- 1*(Omega > 0) - 1 * (Omega<0)
  diag(A) <- 0
  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", weighted = T)
  V(g)$name <- c(expression(Y[1]), expression(Y[2]), expression(Y[3]), 
                 expression(Y[4]), expression(Y[5]), expression(Y[6]),
                 expression(Y[7]), expression(Y[8]), expression(Y[9]),
                 expression(Y[10]))
  V(g)$color <- rep("lightgray", times = nrow(Omega))
  E(g)$color <- ifelse(E(g)$weight<0, "black","gray")
  
  #par(mar = c(1,1,1,1), mgp = c(1.8, 0.5, 0))
  if(star_layout){
    plot(g, vertex.size = 40, vertex.label.color = "black", vertex.label.cex = 2,
         layout = layout_as_star)
  } else{
    plot(g, vertex.size = 40, vertex.label.color = "black", vertex.label.cex = 2,
         layout = layout_in_circle)
  }
  
}



q <- 10
# AR1 
Sigma1 <- 0.7^(abs(outer(1:q, 1:q, FUN = "-")))
Omega1 <- solve(Sigma1)
Omega1[abs(Omega1) < 1e-12] <- 0

png("omega_AR1.png", width = 4.5, height = 4.5, units = "in", res = 300)
par(mar = c(3,1,2,1), mgp = c(1.8, 0.5, 0), cex.main = 3)
plot_Omega(Omega1, "AR1", color = T)
dev.off()

png("graph_AR1.png", width = 4.5, height = 4.5, units = "in", res = 300)
par(mar = c(1,1,1,1), mgp = c(1.8, 0.5, 0))
plot_graph(Omega1)
dev.off()



#########
Omega2 <- matrix(NA, nrow = q, ncol = q)
for(i in 1:q){
  for(j in 1:q){
    if(abs(i-j) > 2) Omega2[i,j] <- 0
    else if(abs(i-j) == 2) Omega2[i,j] <- 0.25
    else if(abs(i-j) ==1) Omega2[i,j] <- 0.5
    else Omega2[i,j] <- 1
  }
}


png("omega_AR2.png", width = 4.5, height = 4.5, units = "in", res = 300)
par(mar = c(3,1,2,1), mgp = c(1.8, 0.5, 0), cex.main = 3)
plot_Omega(Omega2, "AR2")
dev.off()

png("graph_AR2.png", width = 4.5, height = 4.5, units = "in", res = 300)
par(mar = c(1,1,1,1), mgp = c(1.8, 0.5, 0))
plot_graph(Omega2)
dev.off()

######################
# Block model
Omega3 <- matrix(0, nrow = q, ncol = q)
Omega3[1:(q/2),1:(q/2)] <- 0.5
Omega3[q/2 + 1:(q/2),1:(q/2) + q/2] <- 0.5
diag(Omega3) <- 1

png("omega_block.png", width = 4.5, height = 4.5, units = "in", res = 300)
par(mar = c(3,1,2,1), mgp = c(1.8, 0.5, 0), cex.main = 3)
plot_Omega(Omega3, "Block")
dev.off()

png("graph_block.png", width = 4.5, height = 4.5, units = "in", res = 300)
par(mar = c(1,1,1,1), mgp = c(1.8, 0.5, 0))
plot_graph(Omega3)
dev.off()
#######################

Omega4 <- matrix(0, nrow = q, ncol = q)
Omega4[1,] <- 0.75
Omega4[,1] <- 0.75
diag(Omega4) <- 1


png("omega_star.png", width = 4.5, height = 4.5, units = "in", res = 300)
par(mar = c(3,1,2,1), mgp = c(1.8, 0.5, 0), cex.main = 3)
plot_Omega(Omega4, "Star")
dev.off()

png("graph_star.png", width = 4.5, height = 4.5, units = "in", res = 300)
par(mar = c(1,1,1,1), mgp = c(1.8, 0.5, 0))
plot_graph(Omega4, star_layout = TRUE)
dev.off()


#####################
Omega5 <- Omega1
Omega5[1,10] <- Omega5[10,1] <- -1.37
diag(Omega5) <- 1

png("omega_circle.png", width = 4.5, height = 4.5, units = "in", res = 300)
par(mar = c(3,1,2,1), mgp = c(1.8, 0.5, 0), cex.main = 3)
plot_Omega(Omega5, "Circle")
dev.off()

png("graph_circle.png", width = 4.5, height = 4.5, units = "in", res = 300)
par(mar = c(1,1,1,1), mgp = c(1.8, 0.5, 0))
plot_graph(Omega5)
dev.off()


#####################
Omega6 <- matrix(0.75, nrow = q, ncol = q)
diag(Omega6) <- 1

png("omega_full.png", width = 4.5, height = 4.5, units = "in", res = 300)
par(mar = c(3,1,2,1), mgp = c(3, 0.5, 0), cex.main = 3)
plot_Omega(Omega6, "Dense")
dev.off()

png("graph_full.png", width = 4.5, height = 4.5, units = "in", res = 300)
par(mar = c(1,1,1,1), mgp = c(1.8, 0.5, 0))
plot_graph(Omega6)
dev.off()








