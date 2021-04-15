library(ggplot2)
all_res <- list.files("./tests/Formal/Accurancy/k10/results",full.names = T)
#all_res <- list.files("results",full.names = T)

## Omegas
all_Omega_files <- all_res[grep("graph_Omega",all_res)]
all_Omega <- lapply(all_Omega_files,read.csv,row.names=1)
all_Omega <- Reduce(rbind,all_Omega)

all_Omega_k30 <- all_Omega[!is.na(all_Omega$k),]
all_Omega_k30 <- all_Omega_k30[all_Omega_k30$k==10,]

all_Omega_k30$beta.sparsity <- factor( 1-all_Omega_k30$s,levels = c(0.8,0.5))
all_Omega_k30$p <- paste0(all_Omega_k30$p, " predictors")
all_Omega_k30 <- within(all_Omega_k30, p<-factor(p, levels=c("5 predictors", "10 predictors")))
all_Omega_k30$mod <- paste0("model ",all_Omega_k30$mod)

all_Omega_bayes <- all_Omega_k30[!is.na(all_Omega_k30$TP_bayes) & all_Omega_k30$mod!="model 6",c(1:11,20)]
all_Omega_multireg <- all_Omega_k30[(all_Omega_k30$algo=="multireg"|all_Omega_k30$algo=="multireg_mu0"|all_Omega_k30$algo=="multireg_mu0-aug") & all_Omega_k30$mod!="model 6",c(1:7,16:19,20)]

all_Omega_bayes$Sensitivity <- all_Omega_bayes$TP_bayes/(all_Omega_bayes$TP_bayes+all_Omega_bayes$FN_bayes)
all_Omega_bayes$Specificity <- all_Omega_bayes$TN_bayes/(all_Omega_bayes$TN_bayes+all_Omega_bayes$FP_bayes)
all_Omega_bayes$MCC <- (all_Omega_bayes$TP_bayes*all_Omega_bayes$TN_bayes-all_Omega_bayes$FP_bayes*all_Omega_bayes$FN_bayes)/
  sqrt((all_Omega_bayes$TP_bayes+all_Omega_bayes$FP_bayes)*(all_Omega_bayes$TP_bayes+all_Omega_bayes$FN_bayes)*(all_Omega_bayes$TN_bayes+all_Omega_bayes$FP_bayes)*(all_Omega_bayes$TN_bayes+all_Omega_bayes$FN_bayes))

all_Omega_multireg$Sensitivity <- all_Omega_multireg$TP_thrl/(all_Omega_multireg$TP_thrl+all_Omega_multireg$FN_thrl)
all_Omega_multireg$Specificity <- all_Omega_multireg$TN_thrl/(all_Omega_multireg$TN_thrl+all_Omega_multireg$FP_thrl)
all_Omega_multireg$MCC <- (all_Omega_multireg$TP_thrl*all_Omega_multireg$TN_thrl-all_Omega_multireg$FP_thrl*all_Omega_multireg$FN_thrl)/
  sqrt((all_Omega_multireg$TP_thrl+all_Omega_multireg$FP_thrl)*(all_Omega_multireg$TP_thrl+all_Omega_multireg$FN_thrl)*(all_Omega_multireg$TN_thrl+all_Omega_multireg$FP_thrl)*(all_Omega_multireg$TN_thrl+all_Omega_multireg$FN_thrl))

colnames(all_Omega_bayes)[8:11] <- c("TP","TN","FP","FN")
colnames(all_Omega_multireg)[8:11] <- c("TP","TN","FP","FN")

Omega_learning_k30 <- rbind(all_Omega_bayes,all_Omega_multireg)
Omega_learning_k30[is.na(Omega_learning_k30)] <- 0
#Omega_learning_k30$beta.sparsity <- Omega_learning_k30$s

## Original plot:
Graph_MCC <- ggplot(data = Omega_learning_k30,aes(x=algo,y = MCC)) + 
  geom_boxplot(aes(fill = beta.sparsity)) + 
  facet_grid(mod~p) + 
  ylab("MCC on Omega") + 
  xlab("") +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=45,hjust = 1,vjust=1),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

Graph_MCC
ggsave("./tests/Formal/Accurancy/Figs/MCC_k30_Asmallprior_Omega.pdf",Graph_MCC,width = 9,height = 5.5,unit = "in")

## New plot
  Omega_learning_k30 <- within(Omega_learning_k30, algo<-factor(algo, levels=c("CAR-LASSO", "CAR-ALASSO", "SRG-LASSO","GLASSO-aug", "GALASSO-aug", "ad-hoc-aug" ,"GLASSO", "GALASSO","multireg_mu0-aug" ,"multireg", "multireg_mu0", "ad-hoc")))
Graph_MCC <- ggplot(data = Omega_learning_k30,aes(x=algo,y = MCC)) + 
  geom_point(aes(color = beta.sparsity), alpha=0.1, size=1)+
  geom_boxplot(aes(fill = beta.sparsity)) + 
  facet_grid(p~mod) + 
  ylab("MCC on Omega") + 
  xlab("") + theme_bw() +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

Graph_MCC
ggsave("./tests/Formal/Accurancy/Figs/MCC_k10_Asmallprior_Omega2.pdf",Graph_MCC,width = 10,height = 8,unit = "in")


Graph_Sensitivity <- ggplot(data = Omega_learning_k30,aes(x=algo,y = Sensitivity)) + 
  geom_boxplot(aes(fill = beta.sparsity)) + 
  facet_grid(mod~p) + 
  ylab("Sensitivity on Omega") + 
  xlab("") +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=45,hjust = 1,vjust=1),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

Graph_Sensitivity
ggsave("./tests/Formal/Accurancy/Figs/Sensitivity_k10_Asmallprior_Omega.pdf",Graph_Sensitivity,width = 9,height = 5.5,unit = "in")

Graph_Specificity <- ggplot(data = Omega_learning_k30,aes(x=algo,y = Specificity)) + 
  geom_boxplot(aes(fill = beta.sparsity)) + 
  facet_grid(mod~p) + 
  ylab("Specificity on Omega") + 
  xlab("") +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=45,hjust = 1,vjust=1),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

Graph_Specificity
ggsave("./tests/Formal/Accurancy/Figs/Specificity_k30_Asmallprior_Omega.pdf",Graph_Specificity,width = 9,height = 5.5,unit = "in")

