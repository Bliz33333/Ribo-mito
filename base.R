library(rentrez)
library(readr)
library(BiocManager)
#library(dplyr)
library(biomaRt)
library(utils)
library(R.utils)
library(rjson)
library(GenomicDataCommons)
setwd("~/R Workspace/Ribo_mito project")


cyt_rib_l <- read_csv("cyt_rib_l.csv", skip = 1)
cyt_rib_s <- read_csv("cyt_rib_s.csv", skip = 1)
mit_rib_l <- read_csv("mit_rib_l.csv", skip = 1)
mit_rib_s <- read_csv("mit_rib_s.csv", skip = 1)
mito_gene_data <- read_csv("mito_gene_data.csv")

#cyt_rib_hgnc <- c(unlist(cyt_rib_l[2]),unlist(cyt_rib_s[2]))
#names(cyt_rib_hgnc) <- c()
#write.table(cyt_rib_hgnc, file="cyt_rib_hgnc.txt", row.names=FALSE, col.names=FALSE, sep = ",")

raw_cyt_rib_ensembl <- read_csv("raw_cyt_rib_ensembl.txt")
cyt_rib_ensembl <- unique(raw_cyt_rib_ensembl$`Gene stable ID`)

save(cyt_rib_ensembl,file = "cyt_rib_ensemble_file")

datafiles_zipped_raw <- list.files(path = "C:/Users/blaze/Documents/R Workspace/Ribo_mito project/datafiles", recursive = TRUE,full.names = TRUE)
datafiles_zipped <- datafiles_zipped_raw[grep(pattern = ".gz$", fixed = F, datafiles_zipped_raw)]

#untar(tarfile = datafiles_zipped[1],exdir = "C:/Users/blaze/Documents/R Workspace/Ribo_mito project/unzipped", tar = "internal")

for (i in 1:length(datafiles_zipped)) 
{
  my_filename <- strsplit(   strsplit(datafiles_zipped[i],"/")[[1]][length(strsplit(datafiles_zipped[i],"/")[[1]])]   , ".gz")
  output_path <- paste("C:/Users/blaze/Documents/R Workspace/Ribo_mito project/unzipped/", my_filename, sep = "")
  gunzip(filename = datafiles_zipped[i], destname= output_path, overwrite=FALSE, remove=F)
}

unzipped_files <- list.files(path = "C:/Users/blaze/Documents/R Workspace/Ribo_mito project/unzipped", full.names = TRUE)
index_data <- matrix(nrow = length(unzipped_files),ncol = 3)
index_data[,1] <- unzipped_files

index_data[grep("htseq", index_data[,1], fixed = T),2] <- "htseq"
sum(index_data[,2] == "htseq", na.rm = T)
index_data[grep("star_gene", index_data[,1], fixed = T),2] <- "star_gene"
sum(index_data[,2] == "star_gene", na.rm = T)
index_data[grep("FPKM.txt", index_data[,1], fixed = T),2] <- "FPKM"
sum(index_data[,2] == "FPKM", na.rm = T)
index_data[grep("FPKM-UQ", index_data[,1], fixed = T),2] <- "FPKM-UQ"
sum(index_data[,2] == "FPKM-UQ", na.rm = T)

sum(is.na(index_data[,2]))


temp_filenames <- unlist(strsplit(index_data[,1], "/"))
temp_filenames <- temp_filenames[1:(length(temp_filenames)/8)*8]

for(i in 1:nrow(index_data))
{
  temp_index <- grep(temp_filenames[i], datafiles_zipped_raw, fixed = T)
  temp_index
  if(length(temp_index) > 1)
  {
    temp_index <- temp_index[-grep("parcel",datafiles_zipped_raw[temp_index], fixed = T)]
    temp_index 
  }
  temp_path <- unlist(strsplit(datafiles_zipped_raw[temp_index],split = "/"))
  temp_path
  temp_length <- length(temp_path)
  temp_length
  temp_case <- temp_path[temp_length-1]
  temp_case
  index_data[i,3] <- temp_case
  print(i)
}

named_index_data <- cbind(index_data,temp_filenames)

save(named_index_data,file = "named_index_data_file")

uq_data <- named_index_data[named_index_data[,2]== "FPKM-UQ",]

manifest <- fromJSON(file="files.2021-01-26.json")
caseUUIDs <- c()

for (i in 1:nrow(uq_data))
{
  print(i)
  record <- manifest[[grep(uq_data[i,4], manifest, fixed=TRUE, ignore.case=FALSE)]]
  if (record$file_name != paste(uq_data[i,4],".gz",sep = ""))
  {
    print("FALSE")
  }
  caseUUIDs[i] <- record$cases[[1]]$case_id
}

uq_full_index <- cbind(uq_data, caseUUIDs)

#save(uq_full_index,file = "uq_full_index_file")
#load("uq_full_index_file", verbose = T)

disease_data <- gdc_clinical(uq_full_index[,5])$main[,c(1,2,5)]
gdc_clinical(uq_full_index[1:3,5])$main[,c(1,2,5)]

temp_sorted <- matrix(nrow = nrow(uq_full_index),ncol = 3)
for(i in 1:nrow(uq_full_index))
{
  temp_sorted[i,] <- unlist(disease_data[which(disease_data[,1]==uq_full_index[i,5]),1:3])
}



uq_loc_data <- cbind(uq_full_index,temp_sorted[,c(2,3)])

save(uq_loc_data,file = "uq_loc_data_file")

#sum(uq_full_index[,5] != temp_sorted[,1])

#my_filename <- strsplit(   strsplit(datafiles_zipped[1],"/")[[1]][length(strsplit(datafiles_zipped[1],"/")[[1]])]   , ".gz")

#gunzip(filename = filename2, destname= "C:/Users/blaze/Documents/R Workspace/Ribo_mito project/unzipped/1", overwrite=FALSE, remove=F)
