library(ggplot2)
all_res <- list.files("./tests/Formal/Accurancy/k10/results",full.names = T)
##all_res <- list.files("results",full.names = T)

all_loss_files <- all_res[grep("loss",all_res)]
all_loss <- lapply(all_loss_files,read.csv,row.names=1)
all_loss <- Reduce(rbind,all_loss)

all_loss_k30 <- all_loss[!is.na(all_loss$steinOmega),]
all_loss_k30 <- all_loss_k30[all_loss_k30$k==10,]

all_loss_k30$beta.sparsity <- factor( 1-all_loss_k30$s,levels = c(0.8,0.5))
#all_loss_k30$p <- paste0("#.predictors:",all_loss_k30$p)
all_loss_k30$p <- paste0(all_loss_k30$p, " predictors")
all_loss_k30 <- within(all_loss_k30, p<-factor(p, levels=c("5 predictors", "10 predictors")))
all_loss_k30$mod <- paste0("model ",all_loss_k30$mod)

## Original plot:
Stein_k30 <- ggplot(data = all_loss_k30,aes(x=algo,y=log(steinOmega))) + 
    geom_boxplot(aes(fill = beta.sparsity)) + 
    facet_grid(mod~p) + 
    ylab("log Stein's loss of Omega") + 
    xlab("") +
    theme(legend.position="top") + 
    theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=45,hjust = 1,vjust=1),
        plot.margin = margin(.15, .15, .15, .15, "cm"))
ggsave("./tests/Formal/Accurancy/Figs/Stein_k30_Asmallprior.pdf",Stein_k30,width = 10,height = 6,unit = "in")
ggsave("./tests/Formal/Accurancy/Figs/Stein_k30_Asmallprior.jpg",Stein_k30,width = 10,height = 6,unit = "in")

## New plot:
all_loss_k30 <- within(all_loss_k30, algo<-factor(algo, levels=c("CAR-LASSO", "CAR-ALASSO", "SRG-LASSO","GLASSO-aug", "GALASSO-aug", "ad-hoc-aug" ,"GLASSO", "GALASSO","multireg_mu0-aug" ,"multireg", "multireg_mu0", "ad-hoc")))
Stein_k30 <- ggplot(data = all_loss_k30,aes(x=algo,y=log(steinOmega))) + 
  geom_point(aes(color = beta.sparsity), alpha=0.1, size=1)+
  geom_boxplot(aes(fill = beta.sparsity)) + 
  facet_grid(p~mod) + 
  ylab("log Stein's loss of Omega") + 
  xlab("") + theme_bw() +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=90,hjust = 1,vjust=.5),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

Stein_k30
ggsave("./tests/Formal/Accurancy/Figs/Stein_k10_Asmallprior2.pdf",Stein_k30,width = 12,height = 8,unit = "in")
ggsave("Figs/Stein_k10_Asmallprior2.jpg",Stein_k30,width = 10,height = 8,unit = "in")

Stein_k30_aggregate_mean <- aggregate(steinOmega~p+mod+algo+beta.sparsity,data = all_loss_k30,FUN = mean)
Stein_k30_aggregate_sd <- aggregate(steinOmega~p+mod+algo+beta.sparsity,data = all_loss_k30,FUN = sd)

Stein_k30_aggregate <- Stein_k30_aggregate_mean[,1:4]
Stein_k30_aggregate$mean <- Stein_k30_aggregate_mean$steinOmega
Stein_k30_aggregate$sd <- Stein_k30_aggregate_sd$steinOmega


## Original plot:
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
ggsave("./tests/Formal/Accurancy/Figs/beta_k30_Asmallprior.pdf",beta_k30,width = 10,height = 6,unit = "in")
ggsave("./tests/Formal/Accurancy/Figs/beta_k30_Asmallprior.jpg",beta_k30,width = 10,height = 6,unit = "in")


## New plot:
beta_k30 <- ggplot(data = all_loss_k30[!is.na(all_loss_k30$logL2beta),],
                   aes(x=algo,y=.5*logL2beta)) + 
  geom_point(aes(color = beta.sparsity), alpha=0.1, size=1)+
  geom_boxplot(aes(fill = beta.sparsity)) + 
  facet_grid(p~mod) + 
  ylab("log L2 loss of beta") + 
  xlab("") + theme_bw() +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=90,hjust = 1,vjust=1),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

beta_k30
ggsave("Figs/beta_k30_Asmallprior2.pdf",beta_k30,width = 10,height = 8,unit = "in")
ggsave("Figs/beta_k30_Asmallprior2.jpg",beta_k30,width = 10,height = 8,unit = "in")
    