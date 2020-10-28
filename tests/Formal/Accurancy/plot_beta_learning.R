all_beta_files <- all_res[grep("graph_beta",all_res)]
all_beta <- lapply(all_beta_files,read.csv,row.names=1)
all_beta <- Reduce(rbind,all_beta)

all_beta_k30 <- all_beta[!is.na(all_beta$k),]
all_beta_k30 <- all_beta_k30[all_beta_k30$k==30,]

all_beta_k30$beta.sparsity <- factor( 1-all_beta_k30$s,levels = c(0.8,0.5))
all_beta_k30$p <- paste0("#.predictors:",all_beta_k30$p)
all_beta_k30$mod <- paste0("model:",all_beta_k30$mod)

all_beta_bayes <- all_beta_k30[!is.na(all_beta_k30$TP_bayes) & all_beta_k30$mod!="model:6",c(1:11,20)]
all_beta_multireg <- all_beta_k30[all_beta_k30$algo=="multireg" & all_beta_k30$mod!="model:6",c(1:7,16:19,20)]

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


beta_MCC <- ggplot(data = beta_learning_k30,aes(x=algo,y = MCC)) + 
  geom_boxplot(aes(fill = beta.sparsity)) + 
  facet_grid(mod~p) + 
  ylab("MCC on B") + 
  xlab("") +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=45,hjust = 1,vjust=1),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

beta_MCC
ggsave("./tests/Formal/Accurancy/Figs/MCC_k30_Asmallprior_beta.pdf",beta_MCC,width = 8,height = 5,unit = "in")


beta_Sensitivity <- ggplot(data = beta_learning_k30,aes(x=algo,y = Sensitivity)) + 
  geom_boxplot(aes(fill = beta.sparsity)) + 
  facet_grid(mod~p) + 
  ylab("Sensitivity on B") + 
  xlab("") +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=45,hjust = 1,vjust=1),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

beta_Sensitivity
ggsave("./tests/Formal/Accurancy/Figs/Sensitivity_k30_Asmallprior_beta.pdf",beta_Sensitivity,width = 8,height = 5,unit = "in")

beta_Specificity <- ggplot(data = beta_learning_k30,aes(x=algo,y = Specificity)) + 
  geom_boxplot(aes(fill = beta.sparsity)) + 
  facet_grid(mod~p) + 
  ylab("Specificity on B") + 
  xlab("") +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=45,hjust = 1,vjust=1),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

beta_Specificity
ggsave("./tests/Formal/Accurancy/Figs/Specificity_k30_Asmallprior_beta.pdf",beta_Specificity,width = 8,height = 5,unit = "in")
