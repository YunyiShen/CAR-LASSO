library(ggplot2)
soil <- read.csv("./real_data/soil_alpha_cent.csv")
human <- read.csv("./real_data/human_alpha_cent.csv")

soil <- soil[order(soil$mu),]
soil$genus <- factor(soil$X, levels = (as.character(soil$X)))
human <- human[order(human$mu),]
human$genus <- factor(human$X, levels = (as.character(human$X)))

gsoil <- ggplot(data = soil,aes(x=genus,y=alpha_cent)) + 
  geom_bar( stat = "identity", fill="chocolate") + 
  ylab("alpha centrality") +
  xlab("")+ ggtitle("Soil") +
  theme(
    plot.title = element_text(hjust=0.5, size=rel(1.5)),
    axis.title.x = element_text(size=rel(1.2)),
    axis.title.y = element_text(size=rel(1.2), angle=90, vjust=0.5, hjust=0.5),
    axis.text.x = element_text(colour="grey49", size=rel(1.0), angle=0, face="plain"),
    axis.text.y = element_text(colour="grey49", size=rel(1.2), angle=0, face="plain"),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    strip.text = element_text(size = 17)
  ) +
  coord_flip()

ghuman <- ggplot(data = human,aes(x=genus,y=alpha_cent)) + 
  geom_bar( stat = "identity", fill="gold3") + 
  ylab("alpha centrality") +
  xlab("") + ggtitle("Human gut") +
  theme(
    plot.title = element_text(hjust=0.5, size=rel(1.5)),
    axis.title.x = element_text(size=rel(1.2)),
    axis.title.y = element_text(size=rel(1.2), angle=90, vjust=0.5, hjust=0.5),
    axis.text.x = element_text(colour="grey49", size=rel(1.0), angle=0, face="plain"),
    axis.text.y = element_text(colour="grey49", size=rel(1.2), angle=0, face="plain"),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    strip.text = element_text(size = 17)
  ) +
  coord_flip() 


ggpubr::ggarrange(ghuman,gsoil,
          #labels = c("gut", "soil"),
          ncol = 2, nrow = 1)


#ggsave("./real_data/alpha_cent_vs_mu.pdf", width = 9, height = 4, scale = 0.8)
ggsave("alpha_cent_vs_mu2.pdf", width = 10, height = 5, scale = 0.8)  


cor.test(human$mu, human$alpha_cent,
         method = "spearman")
