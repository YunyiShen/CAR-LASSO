raw_annotation <- read.csv("./real_data/Soil/raw_data/abundance_profiling_e5_len16_id60_minabu_10.csv",sep = "\t", header = F,stringsAsFactors = F)
raw_annotation <- raw_annotation[,-grep(" / *.+",as.matrix(raw_annotation[1,]))]
raw_annotation <- raw_annotation[,-c(1:4,6)]
raw_family_list <- as.character( unique(raw_annotation$V5[-1]))

family_annotation <- matrix(NA,length(raw_family_list),ncol(raw_annotation))
#family_annotation[,1] <- raw_family_list
for(i in 1:length( raw_family_list)){
  #for(j in 1:(ncol(raw_annotation)-1)){
    taxa <- raw_family_list[i]
    temp <- raw_annotation[raw_annotation$V5==taxa,]
    if(ncol(temp)>1){
      temp1 <- temp[1,]
      temp1[,-1] <- colSums(matrix(as.numeric(as.matrix( temp[,-1])),nrow = nrow(temp)))
    }
    else {
      temp1 <- temp
      temp1[,-1] <- as.numeric(temp[,-1])
    }
    family_annotation[i,] <- as.matrix(temp1)
    
  #}
}
family_annotation <- rbind(as.matrix(raw_annotation[1,]),family_annotation)
family_annotation <- family_annotation[-grep("unclass",family_annotation[,1]),]

all_family_data <- lapply(1:(ncol(family_annotation)-1),function(i,rawd){
  temp <- rawd[,c(1,i+1)]
  colnames(temp) <- temp[1,]
  temp <- temp[-1,]
  return(temp)
},family_annotation)

names(all_family_data) <- (family_annotation)[1,-1]

source("./real_data/Human/misc.R")
family_name <- lapply(all_family_data,function(w){w[,1]})
family_list <- Reduce(union,family_name)
family_thr_count <- get_taxa_abund_count(all_family_data,family_list,0.01)
family_count_thr <- sapply(1:59,function(w,taxa_count_above_thr){sum(taxa_count_above_thr[,2]>w)},family_thr_count)
plot(1:59,family_count_thr)

