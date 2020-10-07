download_data <- function(sample_id,url,pathsave = "./raw_seq/"){
  download.file(gsub("<id>",sample_id,url),paste0(pathsave,sample_id,".fna"))
}

url_api <- "https://api.mg-rast.org/download/<id>?file=100.1"

all_repos <- read.csv("./real_data/Human/raw_data/mgp154_metadata_library.csv")
all_ids <- as.character(all_repos$metagenome_id)