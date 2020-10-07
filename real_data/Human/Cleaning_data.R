load(file = "./real_data/Human/raw_data/mgp154_all_stat_data.RData")
source("./real_data/Human/misc.R")

all_genus_data <- lapply(all_data,function(w){data.frame(w$statistics$taxonomy$genus,stringsAsFactors = F)})
genus_name <- lapply(all_genus_data,function(w){w[,1]})
genus_list <- Reduce(union,genus_name)
genus_thr_count <- get_taxa_abund_count(all_genus_data,genus_list,0.05)
genus_count_thr <- sapply(1:90,function(w,taxa_count_above_thr){sum(taxa_count_above_thr[,2]>w)},genus_thr_count)
plot(1:90,genus_count_thr)
genus_list_clean <- genus_thr_count$taxa[genus_thr_count$count>90]
genus_list_clean <- genus_list_clean[-c(16:18)]
genus_cout_mat <- get_counting_data(all_genus_data,genus_list_clean)
write.csv(genus_cout_mat,"./real_data/Human/clean_data/genus_mat_.01_60.csv",row.names = F)

all_spp_data <- lapply(all_data,function(w){data.frame(w$statistics$taxonomy$species,stringsAsFactors = F)})
spp_name <- lapply(all_spp_data,function(w){w[,1]})
spp_list <- Reduce(union,spp_name)
spp_thr_count <- get_taxa_abund_count(all_spp_data,spp_list,0.001)
spp_count_thr <- sapply(1:90,function(w,taxa_count_above_thr){sum(taxa_count_above_thr[,2]>w)},spp_thr_count)
spp_list_clean <- spp_thr_count$taxa[spp_thr_count$count>90]




all_family_data <- lapply(all_data,function(w){data.frame(w$statistics$taxonomy$family,stringsAsFactors = F)})
family_name <- lapply(all_family_data,function(w){w[,1]})
family_list <- Reduce(union,family_name)
family_thr_count <- get_taxa_abund_count(all_family_data,family_list,0.05)
family_count_thr <- sapply(1:90,function(w,taxa_count_above_thr){sum(taxa_count_above_thr[,2]>w)},family_thr_count)
plot(1:90,family_count_thr)
family_list_clean <- family_thr_count$taxa[family_thr_count$count>20]
family_list_clean <- family_list_clean[-c(16:18)]
family_cout_mat <- get_counting_data(all_family_data,family_list_clean)
write.csv(family_cout_mat,"./real_data/Human/clean_data/genus_mat_.01_60.csv",row.names = F)





spp_cout_mat <- get_counting_data(all_spp_data,spp_list_clean)
write.csv(spp_cout_mat,"./real_data/Human/clean_data/spp_mat_.001_90.csv",row.names = F)

sample_meta <- read.csv("./real_data/Human/raw_data/metadata_samples.csv")
rownames( sample_meta) <- sample_meta$Subject
seq_sample_desc <- read.csv("./real_data/Human/raw_data/mgp154_metadata-1_sample.csv")

Design_raw <- sample_meta[seq_sample_desc$sample_name,]
Design <- Design_raw[,c("Age","Gender","Stratum","Diet.Group","BMI")]
write.csv(Design,"./real_data/Human/clean_data/Design.csv")
