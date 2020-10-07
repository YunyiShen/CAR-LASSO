get_all_data <- function(sample_id,url){
  temp <- GET(gsub("<id>",sample_id,url))
  json_stat <- rawToChar(temp$content)
  fromJSON(json_stat)
}



#Claesson, Marcus J., et al. "Gut microbiota composition correlates with diet and health in the elderly." Nature 488.7410 (2012): 178-184.

api_url <- "https://api-ui.mg-rast.org/metagenome/<id>?verbosity=stats&detail=taxonomy"

# get all ids
all_repos <- read.csv("./real_data/Human/raw_data/mgp154_metadata_library.csv")
all_ids <- as.character(all_repos$metagenome_id)


library(httr)
library(jsonlite)

all_data <- lapply(all_ids,get_all_data,api_url )
save(file = "mgp154_all_stat_data.RData",list = c("all_data"))

all_genus_data <- lapply(all_data,function(w){w$statistics$taxonomy$genus})
genus_name <- lapply(all_genus_data,function(w){w[,1]})
genus_list <- Reduce(union,genus_name)

all_family_data <- lapply(all_data,function(w){w$statistics$taxonomy$family})
family_name <- lapply(all_family_data,function(w){w[,1]})
family_list <- Reduce(union,family_name)
