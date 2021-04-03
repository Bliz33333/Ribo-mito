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
library(ggplot2)
setwd("~/R Workspace/Ribo_mito project")

cyt_rib_l <- read_csv("cyt_rib_l.csv", skip = 1)
cyt_rib_s <- read_csv("cyt_rib_s.csv", skip = 1)
mit_rib_l <- read_csv("mit_rib_l.csv", skip = 1)
mit_rib_s <- read_csv("mit_rib_s.csv", skip = 1)
mito_gene_data <- read_csv("mito_gene_data.csv")

load("named_index_data_file")
load("uq_loc_data_file")
load("cyt_rib_ensemble_file")

#________________________________________________________________________

transcript_data <- matrix(nrow = nrow(uq_loc_data), ncol = (1+length(cyt_rib_ensembl)))
colnames(transcript_data) <- c("File Path",cyt_rib_ensembl)

print(nrow(uq_loc_data))
#
for(i in 12106:nrow(uq_loc_data))
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

#__________________________________________________

load("transcript_data_file", verbose = T)

sum(uq_loc_data[,1]==transcript_data[,1])
nrow(transcript_data)

final_data <- as.data.frame((transcript_data[,]))
class(final_data$ENSG00000149806)


final_temp1 <- as.data.frame(uq_loc_data)
final_temp2 <- as.data.frame(transcript_data[,])

colnames(final_temp1)[1] <- colnames(final_temp2)[1]

final_data <- merge(final_temp1,final_temp2, by = colnames(final_temp2)[1])

temp_names <- c("File Type","File UUID","Type","Location")

colnames(final_data)[c(2,3,6,7)] <- temp_names

#final_data_test <- sapply(sapply(final_data[,-(1:7)],as.character),as.numeric)

for(i in 8:ncol(final_data))
{
  final_data[,i] <- as.numeric(as.character(final_data[,i]))
}

final_data <- final_data[!(is.na(final_data$Type)),]

final_data$`File UUID` <- as.character(final_data$`File UUID`)
final_data$temp_filenames <- as.character(final_data$temp_filenames)
final_data$caseUUIDs <- as.character(final_data$caseUUIDs)

#save(final_data,file = "final_data_file")

#________________________________________

load("final_data_file", verbose = T)


cancer_types <- names(sort(table(final_data$Type), decreasing = T))
gene_tables_gene_type <- list()
for(gcol in 8:(ncol(final_data)))
{
  temp_dat <- final_data[,c(6,gcol)]
  
  temp_list <- list()
  
  for (t in 1:length(cancer_types)) 
  {
    
    temp_sel_data <- temp_dat[temp_dat[,1]==cancer_types[t],2]
    if(length(temp_sel_data) == 0)
    {
      temp_list[[t]] <- NA
    } else {
      temp_list[[t]] <- temp_sel_data
    }
    
  }
  names(temp_list) <- cancer_types
  
  gene_tables_gene_type[[gcol-7]] <- temp_list
  
}
names(gene_tables_gene_type) <- colnames(final_data)[-(1:7)]

#save(gene_tables_gene_type,file = "gene_tables_gene_type_file")
#save(cancer_types, file = "cancer_types_file")



#________________________________________

load("final_data_file", verbose = T)
load("loss_cases_file")


temp_genes <- character(0)
temp_types <- character(0)
temp_vals <- numeric(0)



for (i in 8:ncol(final_data)) 
{
  temp_gene_vec <- rep(colnames(final_data)[i],times = nrow(final_data))
  temp_genes <- c(temp_genes, temp_gene_vec)
  
  temp_types <- c(temp_types, as.character(final_data$Type))
  
  temp_vals <- c(temp_vals, final_data[,i])
  
  print(i)
}

temp_genes <- as.factor(temp_genes)
temp_types <- as.factor(temp_types)
temp_types <- reorder(temp_types,temp_types,length)

transposed_final_data <- data.table()

transposed_final_data$Type <- temp_types
transposed_final_data$Gene <- temp_genes
transposed_final_data$Expression <- temp_vals





#IDK ABOUT THIS
loss_cases$Gene[loss_cases[,2] == "ENSG00000215472"] <- rep("ENSG00000265681", times = sum(loss_cases[,2] == "ENSG00000215472"))

is_outlier <- matrix(FALSE, nrow = nrow(final_data), ncol = ncol(final_data))

for (i in 1:nrow(loss_cases)) 
{
  temp_col <- which(colnames(final_data) %in% loss_cases[i,2])
  
  is_outlier[final_data$caseUUIDs %in% loss_cases[i,1],temp_col] <- TRUE
  
  #print(i)
  if(sum(final_data$caseUUIDs %in% loss_cases[i,1]) == 0)
  {
    print("ID missing")
  }
}

is_outlier_vec <- as.vector(is_outlier[,-(1:7)])

transposed_final_data$CNV_outl <- is_outlier_vec

transposed_final_data$LogExp <- log2(temp_vals)

#save(transposed_final_data,file= "transposed_final_data_file")

#______________________________________________________________

load("transposed_final_data_file", verbose= T)



load("gene_tables_gene_type_file", verbose = T)

boxplot((gene_tables_gene_type[[1]]), use.cols = TRUE)


boxplot(lapply(gene_tables_gene_type[[1]],log), use.cols = TRUE)






commonLoc <- names(sort(table(final_data$Location),decreasing = T)[1])
commonType <- names(sort(table(final_data[final_data$Location == commonLoc,]$Type),decreasing = T)[1])


averages <- numeric(nrow(final_data))
for (i in 1:nrow(final_data)) 
{
  averages[i] <- mean(unlist(final_data[i,-(1:7)]))
}

selection <- final_data$Location == commonLoc & final_data$Type == commonType

boxplot(final_data[selection,][,-(1:7)], use.cols = TRUE)

boxplot(final_data[selection,][,-(1:7)]/averages[selection], use.cols = TRUE)

boxplot(log(final_data[selection,][,-(1:7)]/averages[selection]), use.cols = TRUE)
boxplot(log(final_data[selection,][,-(1:7)]), use.cols = TRUE)

#boxplot(final_data[final_data$Type == levels(final_data$Type)[11],][,-(1:7)], use.cols = TRUE)
hist(log(final_data[selection,][,11]))





