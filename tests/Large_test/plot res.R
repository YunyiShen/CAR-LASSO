big_test_res <- list.files("./tests/Large_test/Res",full.names = T)

all_res <- lapply(big_test_res,read.csv)
all_res <- Reduce(rbind, all_res)

agg_res <- aggregate(formula = loss~k+n+sparse+algo,FUN = "median", data = all_res)
low.025_res <- aggregate(formula = loss~k+n+sparse+algo,FUN = "quantile", data = all_res,probs = .05)
high.975_res <- aggregate(formula = loss~k+n+sparse+algo,FUN = "quantile", data = all_res,probs = .95)

agg_res$lw <- low.025_res$loss
agg_res$hg <- high.975_res$loss

agg_res$k <- paste0("nodes:",agg_res$k)
agg_res$sparse <- paste0("sparsity:",agg_res$sparse)

require(ggplot2)
ggplot(agg_res,
       aes(log10(n), log10(loss), col=algo)) + 
  geom_line()+
  geom_errorbar(aes(ymin = log10(lw),ymax = log10(hg), col = factor(algo)))+
  xlab("log10 sample size")+
  ylab("log10 Stein's loss")+
  labs(color = "Algorithm") +
  theme(legend.position="top") + 
  scale_fill_brewer()+
  theme(text = element_text(size=12), 
        axis.text.x = element_text(angle=0,size = 12),
        plot.margin = margin(.15, .15, .15, .15, "cm"))+
  facet_grid(sparse~k,labeller = label_parsed, scales = 'free')

ggsave("loss_samle_size_node_sp.pdf",width = 6, height = 5, unit = "in")
ggsave("loss_samle_size_node_sp.png",width = 6, height = 5, unit = "in", dpi = 500)


