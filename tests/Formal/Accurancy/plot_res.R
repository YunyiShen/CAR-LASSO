library(ggplot2)
all_res <- list.files("./tests/Formal/Accurancy/results",full.names = T)

all_loss_files <- all_res[grep("loss",all_res)]
all_loss <- lapply(all_loss_files,read.csv,row.names=1)
all_loss <- Reduce(rbind,all_loss)

all_loss_k30 <- all_loss[!is.na(all_loss$steinOmega),]
all_loss_k30 <- all_loss_k30[all_loss_k30$k==30,]

all_loss_k30$beta.sparsity <- as.factor( all_loss_k30$s)
all_loss_k30$p <- paste0("#.predictors:",all_loss_k30$p)
all_loss_k30$mod <- paste0("Cov.model:",all_loss_k30$mod)

Stein_k30 <- ggplot(data = all_loss_k30,aes(x=algo,y=log(steinOmega))) + 
    geom_boxplot(aes(fill = beta.sparsity)) + 
    facet_grid(mod~p) + 
    ylab("log Stein's loss of Omega") + 
    xlab("") +
    theme(legend.position="top") + 
    theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=45,hjust = 1,vjust=1),
        plot.margin = margin(.15, .15, .15, .15, "cm"))
Stein_k30
ggsave("./tests/Formal/Accurancy/Figs/Stein_k30.pdf",Stein_k30,width = 10,height = 10,unit = "in")
ggsave("./tests/Formal/Accurancy/Figs/Stein_k30.jpg",Stein_k30,width = 10,height = 10,unit = "in")


beta_k30 <- ggplot(data = all_loss_k30[!is.na(all_loss_k30$logL2beta),],
        aes(x=algo,y=.5*logL2beta)) + 
    geom_boxplot(aes(fill = beta.sparsity)) + 
    facet_grid(mod~p) + 
    ylab("log L2 loss of beta") + 
    xlab("") +
    theme(legend.position="top") + 
    theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=45,hjust = 1,vjust=1),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

beta_k30
ggsave("./tests/Formal/Accurancy/Figs/beta_k30.pdf",beta_k30,width = 8,height = 8,unit = "in")
ggsave("./tests/Formal/Accurancy/Figs/beta_k30.jpg",beta_k30,width = 8,height = 8,unit = "in")