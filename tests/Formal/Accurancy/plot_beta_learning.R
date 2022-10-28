
all_res <- list.files("./tests/Formal/Accurancy/k10/results",full.names = T)
all_beta_files <- all_res[grep("graph_beta",all_res)]
all_beta <- lapply(all_beta_files,read.csv,row.names=1)
all_beta <- Reduce(rbind,all_beta)
all_beta$algo <- sub("CAR","CG",all_beta$algo)

all_beta_k30 <- all_beta[!is.na(all_beta$k),]
all_beta_k30 <- all_beta_k30[all_beta_k30$k==10,]

all_beta_k30$beta.sparsity <- factor( 1-all_beta_k30$s,levels = c(0.8,0.5))
all_beta_k30$p <- paste0(all_beta_k30$p, " predictors")
all_beta_k30 <- within(all_beta_k30, p<-factor(p, levels=c("5 predictors", "10 predictors")))
model_list <- c("AR1","AR2","Block","Star","Circle","Dense")
all_beta_k30$mod <- model_list[all_beta_k30$mod]
all_beta_k30$mod <- factor(all_beta_k30$mod, levels = model_list)

all_beta_bayes <- all_beta_k30[!is.na(all_beta_k30$TP_bayes) & all_beta_k30$mod!="model 6",c(1:11,20)]
all_beta_multireg <- all_beta_k30[all_beta_k30$algo=="multireg" & all_beta_k30$mod!="model 6",c(1:7,16:19,20)]

all_beta_bayes$Sensitivity <- all_beta_bayes$TP_bayes/(all_beta_bayes$TP_bayes+all_beta_bayes$FN_bayes)
all_beta_bayes$Specificity <- all_beta_bayes$TN_bayes/(all_beta_bayes$TN_bayes+all_beta_bayes$FP_bayes)
all_beta_bayes$MCC <- (all_beta_bayes$TP_bayes*all_beta_bayes$TN_bayes-all_beta_bayes$FP_bayes*all_beta_bayes$FN_bayes)/
  sqrt((all_beta_bayes$TP_bayes+all_beta_bayes$FP_bayes)*(all_beta_bayes$TP_bayes+all_beta_bayes$FN_bayes)*(all_beta_bayes$TN_bayes+all_beta_bayes$FP_bayes)*(all_beta_bayes$TN_bayes+all_beta_bayes$FN_bayes))

all_beta_multireg$Sensitivity <- all_beta_multireg$TP_thrl/(all_beta_multireg$TP_thrl+all_beta_multireg$FN_thrl)
all_beta_multireg$Specificity <- all_beta_multireg$TN_thrl/(all_beta_multireg$TN_thrl+all_beta_multireg$FP_thrl)
all_beta_multireg$MCC <- (all_beta_multireg$TP_thrl*all_beta_multireg$TN_thrl-all_beta_multireg$FP_thrl*all_beta_multireg$FN_thrl)/
  sqrt((all_beta_multireg$TP_thrl+all_beta_multireg$FP_thrl)*(all_beta_multireg$TP_thrl+all_beta_multireg$FN_thrl)*(all_beta_multireg$TN_thrl+all_beta_multireg$FP_thrl)*(all_beta_multireg$TN_thrl+all_beta_multireg$FN_thrl))

colnames(all_beta_bayes)[8:11] <- c("TP","TN","FP","FN")
colnames(all_beta_multireg)[8:11] <- c("TP","TN","FP","FN")

beta_learning_k30 <- rbind(all_beta_bayes,all_beta_multireg)


## New plot:
library(ggplot2)
beta_learning_k30 <- within(beta_learning_k30, algo<-factor(algo, levels=c("CG-ALASSO", "CG-LASSO",  "SRG-LASSO", "multireg")))
beta_MCC <- ggplot(data = beta_learning_k30,aes(x=algo,y = MCC)) + 
  geom_point(aes(color = beta.sparsity), alpha=0.1, size=1)+
  geom_boxplot(aes(fill = beta.sparsity)) + 
  facet_grid(p~mod) + 
  ylab("MCC on B") + 
  xlab("") + theme_bw() +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=90,hjust = 1,vjust=1),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

beta_MCC
ggsave("./tests/Formal/Accurancy/Figs/MCC_k10_Asmallprior_beta2.pdf",beta_MCC,width = 10,height = 8,unit = "in")


beta_Sensitivity <- ggplot(data = beta_learning_k30,aes(x=algo,y = Sensitivity)) + 
  geom_point(aes(color = beta.sparsity), alpha=0.1, size=1)+
  geom_boxplot(aes(fill = beta.sparsity)) + 
  facet_grid(p~mod) + 
  ylab("Sensitivity on B") + 
  xlab("") + theme_bw() +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=90,hjust = 1,vjust=1),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

beta_Sensitivity
ggsave("./tests/Formal/Accurancy/Figs/Sensitivity_k30_Asmallprior_beta.pdf",beta_Sensitivity,width = 8,height = 5,unit = "in")

beta_Specificity <- ggplot(data = beta_learning_k30,aes(x=algo,y = Specificity)) + 
  geom_point(aes(color = beta.sparsity), alpha=0.1, size=1)+
  geom_boxplot(aes(fill = beta.sparsity)) + 
  facet_grid(p~mod) + 
  ylab("Specificity on B") + 
  xlab("") + theme_bw() +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=90,hjust = 1,vjust=1),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

beta_Specificity
ggsave("./tests/Formal/Accurancy/Figs/Specificity_k30_Asmallprior_beta.pdf",beta_Specificity,width = 8,height = 5,unit = "in")
