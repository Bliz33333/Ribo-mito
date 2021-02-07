library(rentrez)
library(readr)
library(BiocManager)
#library(dplyr)
library(biomaRt)
library(utils)
library(R.utils)
library(rjson)
library(GenomicDataCommons)
library(data.table)
setwd("~/R Workspace/Ribo_mito project")

cyt_rib_l <- read_csv("cyt_rib_l.csv", skip = 1)
cyt_rib_s <- read_csv("cyt_rib_s.csv", skip = 1)
mit_rib_l <- read_csv("mit_rib_l.csv", skip = 1)
mit_rib_s <- read_csv("mit_rib_s.csv", skip = 1)
mito_gene_data <- read_csv("mito_gene_data.csv")

load("named_index_data_file")
load("uq_loc_data_file")
load("cyt_rib_ensemble_file")

transcript_data <- matrix(nrow = nrow(uq_loc_data), ncol = (1+length(cyt_rib_ensembl)))
colnames(transcript_data) <- c("File Path",cyt_rib_ensembl)

print(nrow(uq_loc_data))
#
for(i in 1:nrow(uq_loc_data))
{
  temp_full_data <- read_table2(file = uq_loc_data[i,1],col_names = F)
  temp_names_zeros <- unlist(temp_full_data$X1)
  temp_vals_zeros <- unlist(temp_full_data$X2)
  
  temp_non_zeros <- (temp_vals_zeros != 0)
  
  temp_vals <- temp_vals_zeros[temp_non_zeros]
  temp_names <- temp_names_zeros[temp_non_zeros]
  
  temp_names_clean <- unlist(strsplit(temp_names,split = ".", fixed = T))[2*(1:length(temp_names))-1]
  
  transcript_data[i,1] <- uq_loc_data[i,1]
  
  for(j in 1:length(cyt_rib_ensembl))
  {
    transcript_data[i,j+1] <- sum(temp_vals[which(temp_names_clean == cyt_rib_ensembl[j])])
  }
  
  
  print(i)
}

#save(transcript_data,file = "transcript_data_file")
load("transcript_data_file", verbose = T)

sum(uq_loc_data[,1]==transcript_data[,1])
nrow(transcript_data)

final_data <- as.data.frame((transcript_data[,]))
class(final_data$ENSG00000149806)


final_temp1 <- as.data.frame(uq_loc_data)
final_temp2 <- as.data.frame(transcript_data[,])

colnames(final_temp1)[1] <- colnames(final_temp2)[1]

final_data <- merge(final_temp1,final_temp2, by = colnames(final_temp2)[1])

temp_names <- c()





