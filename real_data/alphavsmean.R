library(ggplot2)
soil <- read.csv("./real_data/soil_alpha_cent.csv")
human <- read.csv("./real_data/human_alpha_cent.csv")

soil <- soil[order(soil$mu),]
soil$genus <- factor(soil$X, levels = (as.character(soil$X)))
human <- human[order(human$mu),]
human$genus <- factor(human$X, levels = (as.character(human$X)))
gsoil <- ggplot(data = soil,aes(x=genus,y=alpha_cent)) + 
  geom_bar( stat = "identity") + 
  ylab("alpha centrality") +
  xlab("soil")+
  coord_flip()

ghuman <- ggplot(data = human,aes(x=genus,y=alpha_cent)) + 
  geom_bar( stat = "identity") + 
  ylab("alpha centrality") +
  xlab("gut") +
  coord_flip() 


ggpubr::ggarrange(ghuman,gsoil,
          #labels = c("gut", "soil"),
          ncol = 2, nrow = 1)


ggsave("./real_data/alpha_cent_vs_mu.pdf", width = 9, height = 4, scale = 0.8)
  


cor.test(human$mu, human$alpha_cent,
         method = "spearman")
