#' plot the chain graph estimated by CAR-LASSO with threshold method using ggraph
#'
#' @param obj The carlasso_out object
#' @param tol threshold default 0.01
#' @return A `ggplot` object
#' @export



plot.carlasso_out <- function(obj, tol = 0.01) {
    col_pn <- c("lightblue","pink")
    # graph structure using threshold:
    B_binary <- abs(obj$point_est$beta) > tol
    Graph_binary <- abs(obj$point_est$Omega) > tol

    diag(Graph_binary) <- 1
    CAR <- get_CAR_MB(obj$point_est$beta*B_binary,
        Graph_binary*obj$point_est$Omega)
    
    n_resp <- length(obj$nodes$response)
    n_pred <- length(obj$nodes$predictors)

    vertices_df <- data.frame(id = c(paste0("resp", 1:n_resp), 
        paste0("pred", 1:n_pred)),
        group = c(rep("resp", n_resp),
        rep("pred", n_pred)))

    ind_mat_resp <- expand.grid(from = 1:n_resp, to = 1:n_resp)
    ind_mat_resp <- ind_mat_resp[ind_mat_resp$from != ind_mat_resp$to, ]
    ind_mat_resp$weight <- sapply(1:nrow(ind_mat_resp),
            function(i, indmat, mat) {
                mat[indmat$from[i], indmat$to[i]]
            }
        ,ind_mat_resp, CAR$C)

    ind_mat_pred <- expand.grid(from = 1:n_pred, to = 1:n_resp)
    ind_mat_pred$weight <- sapply(1:nrow(ind_mat_pred), 
            function(i, indmat, mat) {
                mat[indmat$from[i], indmat$to[i]]
            }
        , ind_mat_pred, CAR$B)

    ind_mat_pred$from <- paste0("pred", ind_mat_pred$from)
    ind_mat_pred$to <- paste0("resp", ind_mat_pred$to)
    ind_mat_resp$from <- paste0("resp", ind_mat_resp$from)
    ind_mat_resp$to <- paste0("resp", ind_mat_resp$to)

    edge_df <- rbind(ind_mat_resp, ind_mat_pred)
    edge_df <- edge_df[edge_df$weight != 0, ]

    edge_abs_df <- edge_df
    edge_abs_df$weight <- abs(edge_abs_df$weight)

    full_graph <- graph.data.frame(edge_df, vertices_df, directed = T)


    col_ER <- c("orange", "darkgreen")
    shape_ER <- c("square", "circle")
    type <- c("predictors", "microbe")
    direction <- c("negative", "positive")

    E(full_graph)$edge.color <- col_pn[(sign(E(full_graph)$weight) + 1) / 2 + 1]
    E(full_graph)$direction. <- direction[(sign(E(full_graph)$weight) + 1) / 2 + 1]
    E(full_graph)$abs_weight <- abs(E(full_graph)$weight)


    V(full_graph)$name <- c(obj$nodes$responses, obj$nodes$predictors)

    V(full_graph)$alpha_centrality <- alpha_centrality(full_graph)
    V(full_graph)$type <- type[c(rep(2, n_resp), rep(1, n_pred))]

    cbPalette_edge <- c("#0072B2", "#990000")
    cbPalette_node <- c("#0815d3", "#682d01")

    set_graph_style(plot_margin = margin(10, 10, 10, 10))
    p <- ggraph(full_graph, layout = "circle") +
        geom_edge_link(aes(color = direction.,
            width = abs_weight, alpha = abs_weight)) +
        scale_edge_color_manual(values = (cbPalette_edge)) +
        geom_node_point(mapping = aes(shape = type,
                size = alpha_centrality, stroke = 1.5),
            col = "#000000", fill = "white", alpha = 1) +
        scale_shape_manual(values = c(21, 24)) +
        coord_fixed(clip = "off") +
        guides(
            width = guide_legend(order = 1),
            size = guide_legend(order = 2),
            shape = FALSE, 
            edge_color = FALSE 
        )


    dd <- rep(0, length(V(full_graph)$name))

    p <- p + geom_node_text(aes(label = name), nudge_x = p$data$x * .38, nudge_y = p$data$y * .2 + dd, family = "") + # repel = T,check_overlap = T)+
        theme_graph(base_family = "Helvetica") +
        theme(
            legend.text = element_text(size = 9),
            legend.position = "bottom"
        )
    p
}