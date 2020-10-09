library(igraph)
library(GGally)
source("./R/misc.R")
B_binary <- multireg_beta/A_beta<.5
diag(Graph_binary) <- 1
CAR <- get_CAR_MB(A_beta*B_binary,Graph_binary*A_Graph)

among_spp <- graph_from_adjacency_matrix(CAR$C,mode = "directed",weighted = T,diag = F)
abs_graph <- graph_from_adjacency_matrix(abs(Graph_binary*A_Graph),mode = "undirected",weighted = T,diag = F)
linear_reg_graph <- graph_from_adjacency_matrix(t(CAR$C),weighted = T, mode = "directed")

col_pn <- c("pink","lightblue")
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
rev_edge_df <- edge_abs_df
rev_edge_df$from <- edge_df$to
rev_edge_df$to <- edge_df$from
rev_graph <- graph.data.frame(rev_edge_df,vertices_df,directed=T)


col_ER <- c("orange","darkgreen")
shape_ER <- c("square","circle")

E(full_graph)$edge.color <- col_pn[(sign(E(full_graph)$weight)+1)/2+1]

l <-layout_with_fr(full_graph)
par(mar=c(0,0,0,0)+.1)
set.seed(428)
plot(full_graph,edge.arrow.size=.2,
     #vertex.label=colnames(comp_mat)[1:35],
     vertex.label=c(colnames(comp_mat)[-ncol(comp_mat)],colnames(Design_dummy)),
     vertex.shape = shape_ER[c(rep(2,(ncol(comp_mat)-1)),rep(1,ncol(Design_dummy)))],
     vertex.size = 2*(alpha_centrality(full_graph)),
     vertex.color = col_ER[c(rep(1,(ncol(comp_mat)-1)),rep(2,ncol(Design_dummy)))],
     layout = layout_with_gem,
     #layout = layout_in_circle,
     edge.color = col_pn[(sign(E(full_graph)$weight)+1)/2+1],
     edge.curved=0)
    