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
library(gridExtra)
setwd("~/R Workspace/Ribo_mito project")

load("final_data_file", verbose = T)
load("transposed_final_data_file", verbose= T)
load("cancer_types_file", verbose = T)


lowMean <- function(dat)
{
  return(mean(dat)-1)
}

highMean <- function(dat)
{
  return(mean(dat)+1)
}

plots <- list()

for(i in 8:length(colnames(final_data)))
{
  gene_name <- colnames(final_data)[i]
  this_data <- transposed_final_data[transposed_final_data$Gene==gene_name,]
  this_data_out <- this_data[this_data$CNV_outl,]

  if(!(sum(this_data$Expression) == 0))
  {
    plots[[i-7]] <- ggplot(this_data,aes(x=LogExp, y=Type)) + geom_jitter(alpha = .25, shape = 1) + ggtitle(gene_name) + facet_grid(Type~., scales = "free", space = "free") + stat_summary(color = "red", fun.min = lowMean, fun = mean, fun.max = highMean ) + geom_point(this_data_out, mapping =  aes(x=LogExp, y=Type), color = "blue")
  }
  

  
  
  #ggplot(transposed_final_data[transposed_final_data$Gene==gene_name,],aes(x=LogExp, y=Type)) + geom_point(alpha = .25, shape = 1)
  #ggplot(transposed_final_data[transposed_final_data$Gene==gene_name,],aes(x=LogExp, y=Type)) + geom_violin()
  print(i)
}

for(i in i:length(plots))
{
  if(!(is.null(plots[[i]])))
  {
    ggsave(filename = (paste("plot", i, ".pdf", sep = "")), width = 20, height = 9, plot = plots[[i]])
  }
  
  
  print(i)
}

#__________________________________________________________

#out_transposed_fd <- transposed_final_data[transposed_final_data$CNV_outl == T,]

mean_chart <- matrix(nrow = length(cancer_types), ncol = ncol(final_data)-7)
rownames(mean_chart) <- cancer_types
colnames(mean_chart) <- colnames(final_data)[-(1:7)]

log_mean_chart <- mean_chart

#out_mean_chart <- mean_chart

for(i in 1:nrow(mean_chart))
{
  
  for(j in 1:ncol(mean_chart))
  {
    temp_sel <- ((transposed_final_data$Gene == colnames(mean_chart)[j]) & (transposed_final_data$Type == rownames(mean_chart)[i])) 
    #temp_out_sel <- temp_sel & transposed_final_data$CNV_outl
    
    temp_sel_2 <- log2(transposed_final_data$Expression[temp_sel])
    temp_sel_2[temp_sel_2==-Inf] <- NA
    
    log_mean_chart[i,j] <- mean(temp_sel_2, na.rm = T)
    
    #mean_chart[i,j] <- mean(transposed_final_data$Expression[temp_sel])
    #out_mean_chart[i,j] <- mean(transposed_final_data$Expression[temp_out_sel])
    
    print(j)
    
  }
  
  print(i)
  
}


out_ratios <- out_mean_chart / mean_chart

#save(log_mean_chart, file =  "log_mean_chart_file")
#save(mean_chart, file =  "mean_chart_file")
#save(out_mean_chart, file = "out_mean_chart_file")
#save(out_ratios, file = "out_ratios_file")

out_ratios[is.nan(out_ratios)] <- NA



temp_remove <- logical(length = nrow(out_ratios))

for (i in 1:nrow(out_ratios)) 
{
  if(sum(out_ratios[i,], na.rm = T) == 0)
  {
    temp_remove[i] <- T
  }
}

out_ratios_clean <- out_ratios[!temp_remove,]


temp_remove_2 <- logical(length = ncol(out_ratios_clean))

for (i in 1:ncol(out_ratios_clean)) 
{
  if(sum(out_ratios_clean[,i], na.rm = T) == 0)
  {
    temp_remove_2[i] <- T
  }
}

out_ratios_clean_2 <- out_ratios_clean[,!temp_remove_2]


heatmap(out_ratios_clean_2, na.rm = T)


ggplot(out_ratios, aes(colnames(out_ratios), rownames(out_ratios), fill= out_ratios)) + 
  geom_tile()


x <- LETTERS[1:20]
y <- paste0("var", seq(1,20))
data <- expand.grid(X=x, Y=y)
data$Z <- runif(400, 0, 5)

#------------------------------------------------------------------

load("log_mean_chart_file")

final_data_norm <- final_data

final_data_norm[,-(1:7)] <- log2(final_data[,-(1:7)])
final_data_norm[final_data_norm[,] == -Inf] <- NA

log_sd_chart <- log_mean_chart

for(i in 1:nrow(log_mean_chart))
{
  
  for(j in 1:ncol(log_mean_chart))
  {
    temp_sel <- ((final_data$Gene == colnames(mean_chart)[j]) & (transposed_final_data$Type == rownames(mean_chart)[i])) 
    
    
    
    print(j)
    
  }
  
  print(i)
  
}