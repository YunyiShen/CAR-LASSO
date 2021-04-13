library(igraph)
library(GGally)
source("./R/misc.R")
# Some resources on network plots: http://mr.schochastics.net/netVizR.html

# for human
load("./real_data/Human/res/CAR_full_design_genus_.005_50_without_unclass_long_chain_lambdaprior_1e-2_1e-6.RData")

# for soil
load("./real_data/Soil/res/CAR_full_design_genus_.01_50_without_unclass_lambda_prior_1e-2_1e-6_long_chain.RData")

B_marginal <- multireg_beta %*% solve(multireg_Graph)

B_binary <- abs(B_marginal) > 1e-3
Graph_binary <- abs(multireg_Graph) > 1e-3
diag(Graph_binary) <- 1

CAR <- get_CAR_MB(B_marginal*B_binary,Graph_binary*multireg_Graph)

among_spp <- graph_from_adjacency_matrix(CAR$C,mode = "directed",weighted = T,diag = F)
abs_graph <- graph_from_adjacency_matrix(abs(Graph_binary*A_Graph),mode = "undirected",weighted = T,diag = F)
linear_reg_graph <- graph_from_adjacency_matrix(t(CAR$C),weighted = T, mode = "directed")

col_pn <- c("lightblue","pink")
l <-layout_with_fr(among_spp)#, repulserad=vcount(among_spp)^3,area=vcount(among_spp)^2.4)
plot(among_spp,edge.arrow.size=.1,
     vertex.label=colnames(comp_mat)[-ncol(comp_mat)],
     vertex.size = (3* alpha_centrality(among_spp)),
     #vertex.label=1:35,
     layout = layout_with_gem,
     edge.color = col_pn[1+(sign(E(among_spp)$weight)+1)/2],
     edge.curved=0)


vertices_df <- data.frame(id=c(paste0("M",1:(ncol(comp_mat)-1)),paste0("E",1:ncol(Design_dummy))),group = c(rep("microbe",(ncol(comp_mat)-1)),rep("env",ncol(Design_dummy))))

ind_mat_micro <- expand.grid(from = 1:(ncol(comp_mat)-1),to = 1:(ncol(comp_mat)-1))
ind_mat_micro <- ind_mat_micro[ind_mat_micro$from!=ind_mat_micro$to,]
ind_mat_micro$weight <- sapply(1:nrow(ind_mat_micro),function(i,indmat,mat){
  mat[indmat$from[i],indmat$to[i]]
},ind_mat_micro,CAR$C)

ind_mat_env <- expand.grid(from = 1:ncol(Design_dummy),to = 1:(ncol(comp_mat)-1))
ind_mat_env$weight <- sapply(1:nrow(ind_mat_env),function(i,indmat,mat){
  mat[indmat$from[i],indmat$to[i]]
},ind_mat_env,B_marginal * B_binary)

ind_mat_env$from <- paste0("E",ind_mat_env$from)
ind_mat_env$to <- paste0("M",ind_mat_env$to)
ind_mat_micro$from <- paste0("M",ind_mat_micro$from)
ind_mat_micro$to <- paste0("M",ind_mat_micro$to)

edge_df <- rbind(ind_mat_micro,ind_mat_env)
edge_df <- edge_df[edge_df$weight!=0,]

edge_abs_df <- edge_df
edge_abs_df$weight <- abs(edge_abs_df$weight)

full_graph <- graph.data.frame(edge_df,vertices_df,directed=T)
full_graph_abs <- graph.data.frame(edge_abs_df,vertices_df,directed=T)
rev_edge_df <- edge_df
rev_edge_df$from <- edge_df$to
rev_edge_df$to <- edge_df$from
rev_graph <- graph.data.frame(rev_edge_df,vertices_df,directed=T)


col_ER <- c("orange","darkgreen")
shape_ER <- c("square","circle")
type <- c("predictors","microbe")
direction <- c("negative","positive")

E(full_graph)$edge.color <- col_pn[(sign(E(full_graph)$weight)+1)/2+1]
E(full_graph)$direction. <- direction[(sign(E(full_graph)$weight)+1)/2+1]
E(full_graph)$abs_weight <- abs( E(full_graph)$weight)

# IF
# for human:
V(full_graph)$name <- c(colnames(comp_mat)[-ncol(comp_mat)],colnames(Design_dummy))
# for soil
V(full_graph)$name <- c(colnames(comp_mat)[-ncol(comp_mat)],"fertilizer","crop agriculture","poorly drained","total N","October measure")
# ENDIF

V(full_graph)$alpha_centrality <- alpha_centrality(full_graph)
V(full_graph)$type <- type[c(rep(2,(ncol(comp_mat)-1)),rep(1,ncol(Design_dummy)))]

#cbPalette_edge <- c( "#0072B2", "#D55E00")
cbPalette_edge <- c( "#0072B2", "#990000")
cbPalette_node <- c( "#0815d3", "#682d01")

# IF
## FOR HUMAN
V(full_graph)$name[17] <- "StratumDayHospital"
V(full_graph)$name[18] <- "StratumLong-term"
V(full_graph)$name[19] <- "StratumRehab"
V(full_graph)$name[20] <- "Diet1"
V(full_graph)$name[21] <- "Diet2"
V(full_graph)$name[22] <- "Diet3"
V(full_graph)$name[23] <- "Diet4"
V(full_graph)$name[24] <- "DietPEG"
## FOR SOIL
V(full_graph)$name[4] <- "C.Koribacter"
V(full_graph)$name[5] <- "C.Nitrososphaera"
V(full_graph)$name[6] <- "C.Solibacter"
# END IF

library(ggplot2)
library(ggraph)
set_graph_style(plot_margin = margin(10,10,10,10))
p <- ggraph(full_graph,layout = "circle")+
  geom_edge_link(aes(color = direction.,width = abs_weight,alpha = abs_weight)) + 
  scale_edge_color_manual(values = (  cbPalette_edge))+
  #geom_node_point(mapping = aes(shape = type,size = alpha_centrality,stroke = 6),col = "#696969",alpha = 1) +
  geom_node_point(mapping = aes(shape = type,size = alpha_centrality,stroke = 1.5),col = "#000000",fill= "white", alpha = 1) +
  scale_shape_manual(values = c(21,24))+
  #scale_color_manual(values = cbPalette)+
  coord_fixed(clip = 'off')+
  guides(width = guide_legend(order = 1),
         #color = guide_legend(order=3),
         size = guide_legend(order=2),
         shape = FALSE, #guide_legend(order=1),
         edge_color = FALSE #guide_legend(order = 3)
  )

## FOR HUMAN: 
dd = c(rep(0,6),0.1,rep(0,12),-0.1,rep(0,5))
## FOR SOIL:
dd = c(rep(0,5),0.1,rep(0,16))

## to put labels names outside the circle
p2 <- p  + geom_node_text(aes(label = name),nudge_x = p$data$x * .38, nudge_y = p$data$y * .2+dd, family = "")+  #repel = T,check_overlap = T)+
  theme_graph(base_family = 'Helvetica')+
  theme(legend.text=element_text(size=9),
        legend.position = "bottom")


p2
#ggsave("./real_data/humangut.pdf",width = 8,height = 6)
ggsave("./real_data/humangut2_margin.pdf",p2,width = 8,height = 7)

ggsave("./real_data/soil2_margin.pdf",p2,width = 8,height = 7)

