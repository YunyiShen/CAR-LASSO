get_taxa_abund_count <- function(raw_data, taxa_list, thr_comp = .005){
  composition <- lapply(raw_data,function(w,thr_comp){w <- data.frame(w,stringsAsFactors = F);w[,2] <- as.numeric(w[,2])/sum(as.numeric(w[,2]))>thr_comp;return(w)},thr_comp)
  
  taxa_count_above_thr <- data.frame(taxa = taxa_list,count = NA,stringsAsFactors = F)
  
  for(taxa in taxa_list){
    n <- 0
    for(j in 1:length(composition)){
      temp <- composition[[j]]
      taxa_in_temp <- temp[,1]
      if(!taxa %in% taxa_in_temp) {next}
      else {n <- n + (temp[temp[,1]==taxa,2])}
    }
    taxa_count_above_thr[taxa_count_above_thr$taxa==taxa,2] <- n
  }
  return(taxa_count_above_thr)
}


get_counting_data <- function(raw_data, taxa_list){
  n <- length(raw_data)
  k <- length(taxa_list)
  res <- matrix(0,nrow = n, ncol = k+1)
  colnames(res) <- c(taxa_list,"all_others")
  
  for(i in 1:n){
    temp <- raw_data[[i]]
    temp <- data.frame(temp,stringsAsFactors = F)
    temp[,2] <- as.numeric(temp[,2])
    temp[,1] <- as.character(temp[,1])
    taxa_list_temp <- intersect(temp[,1],taxa_list)
    for (taxa in taxa_list_temp) {
      res[i,taxa] <- temp[temp[,1]==taxa,2]
    }
    res[i,k+1] <- sum(temp[!temp[,1]%in%taxa_list,2])
    
  }
  rownames(res) <- names(raw_data)
  return(res)
}
