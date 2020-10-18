raw_annotation <- read.csv("./real_data/Soil/raw_data/abundance_profiling_e5_len16_id60_minabu_10.csv",sep = "\t", header = F)

genus_raw <- as.matrix( raw_annotation[,-c(1:5)])
#genus_raw <- genus_raw[,c(1,grep("-July",genus_raw[1,]))]

genus_raw <- genus_raw[,-grep(" / *.+",genus_raw[1,])]
genus_raw <- genus_raw[-grep("unclass",genus_raw[,1]),]

all_genus_data <- lapply(1:(ncol(genus_raw)-1),function(i,rawd){
  temp <- rawd[,c(1,i+1)]
  colnames(temp) <- temp[1,]
  temp <- temp[-1,]
  return(temp)
},genus_raw)

names(all_genus_data) <- (genus_raw)[1,-1]

source("./real_data/Human/misc.R")
genus_name <- lapply(all_genus_data,function(w){w[,1]})
genus_list <- Reduce(union,genus_name)
genus_thr_count <- get_taxa_abund_count(all_genus_data,genus_list,0.01)
genus_count_thr <- sapply(1:59,function(w,taxa_count_above_thr){sum(taxa_count_above_thr[,2]>w)},genus_thr_count)
plot(1:59,genus_count_thr)

genus_list_clean <- genus_thr_count$taxa[genus_thr_count$count>50]
#genus_list_clean <- genus_list_clean[-c(15:21)]
genus_cout_mat <- get_counting_data(all_genus_data,genus_list_clean)
#genus_cout_mat <- genus_cout_mat[-grep("c",rownames(genus_cout_mat)),]
Design <- read.csv("./real_data/Soil/raw_data/Design.csv",row.names = 1)
Design$Octo <- 0+is.na(Design$total.nitrogen.content.of.the.soil.Units.of.g.N.kg.soil)
Design$total.nitrogen.content.of.the.soil.Units.of.g.N.kg.soil[is.na(Design$total.nitrogen.content.of.the.soil.Units.of.g.N.kg.soil)]<-
Design$total.nitrogen.content.of.the.soil.Units.of.g.N.kg.soil[which(is.na(Design$total.nitrogen.content.of.the.soil.Units.of.g.N.kg.soil))-1]

#genus_cout_mat <- genus_cout_mat[rownames(Design),]
Design_clean <- Design[rownames(genus_cout_mat),]

#genus_cout_mat <- genus_cout_mat[,c(1:8,10:15,9)]

write.csv(genus_cout_mat,"./real_data/Soil/clean_data/genus_mat_.01_50_without_unclass.csv",row.names = T)
write.csv(Design_clean,"./real_data/Soil/clean_data/Design_.01_50_without_unclass.csv",row.names = T)
