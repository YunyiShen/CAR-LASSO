raw_annotation <- read.csv("./real_data/Human/raw_data/abundance profiling_MG_mgp154_e5_ind60_len15_min.ab10.tsv.csv", header = F)

genus_raw <- as.matrix( raw_annotation[,-c(1:5)])

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
genus_thr_count <- get_taxa_abund_count(all_genus_data,genus_list,0.005)
genus_count_thr <- sapply(1:90,function(w,taxa_count_above_thr){sum(taxa_count_above_thr[,2]>w)},genus_thr_count)
plot(1:90,genus_count_thr)

genus_list_clean <- genus_thr_count$taxa[genus_thr_count$count>50]
genus_list_clean <- genus_list_clean[-c(12:14)]
genus_cout_mat <- get_counting_data(all_genus_data,genus_list_clean)
rownames(genus_cout_mat) <- gsub("-T.","",rownames(genus_cout_mat))
genus_cout_mat <- genus_cout_mat[-grep("c",rownames(genus_cout_mat)),]
Design <- read.csv("./real_data/Human/raw_data/Design.csv",row.names = 1)
#genus_cout_mat <- genus_cout_mat[rownames(Design),]
Design_clean <- Design[rownames(genus_cout_mat),]


write.csv(genus_cout_mat,"./real_data/Human/clean_data/genus_mat_.005_50.csv",row.names = T)
write.csv(Design_clean,"./real_data/Human/clean_data/Design_.005_50.csv",row.names = T)
