library(igraph)
library(GGally)
source("./R/misc.R")
B_binary <- abs(A_beta/multireg_beta) > .5
Graph_binary <- abs(A_Graph/multireg_Graph) > .5
diag(Graph_binary) <- 1

diag(Graph_binary) <- 1
CAR <- get_CAR_MB(A_beta*B_binary,Graph_binary*A_Graph)

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
},ind_mat_env,CAR$B)

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

cbPalette_edge <- c( "#0072B2", "#D55E00")
cbPalette_node <- c( "#0815d3", "#682d01")


library(ggplot2)
library(ggraph)
set_graph_style(plot_margin = margin(10,10,10,10))
ggraph(full_graph,layout = "circle")+
  geom_edge_link(aes(color = direction.,width = abs_weight,alpha = abs_weight)) + 
  scale_edge_color_manual(values = (  cbPalette_edge))+
  geom_node_point(mapping = aes(shape = type,size = alpha_centrality,stroke = 6),col = "#696969",alpha = 1) +
  #scale_color_manual(values = cbPalette) + 
  geom_node_text(aes(label = name),nudge_y = 0.0,family = "",repel = T,check_overlap = T)+
  coord_fixed(clip = 'off')+
  theme(legend.text=element_text(size=10))+
  guides(width = guide_legend(order = 4),
         #color = guide_legend(order=3),
         size = guide_legend(order=2),
         shape = guide_legend(order=1),
         edge_color = guide_legend(order = 3))+
  theme_void()

ggsave("./real_data/humangut.pdf",width = 8,height = 6)

