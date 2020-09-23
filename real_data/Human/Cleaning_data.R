load(file = "mgp154_all_stat_data.RData")

all_genus_data <- lapply(all_data,function(w){w$statistics$taxonomy$genus})
genus_name <- lapply(all_genus_data,function(w){w[,1]})
genus_list <- Reduce(union,genus_name)

all_family_data <- lapply(all_data,function(w){w$statistics$taxonomy$family})
family_name <- lapply(all_family_data,function(w){w[,1]})
family_list <- Reduce(union,family_name)