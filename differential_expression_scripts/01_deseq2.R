library(DESeq2)
library(tidyverse)

setwd("/home/victoe/Genome_analysis/Data/DifferentialExpression")

metadata <- read.csv("/home/victoe/Genome_analysis/Data/Metadata/Metadata.csv")

metadata <- metadata[order(metadata$Run),]  

count_files <- list.files(pattern="*.counts.txt$")
sample_names <- gsub(".counts.txt", "", count_files)


count_data_list <- lapply(count_files, function(file) {
  read.table(file, header=FALSE, row.names=1)
})
names(count_data_list) <- sample_names
count_matrix <- do.call(cbind, count_data_list)

conditions <- factor(metadata$tissue)
cultivar <- factor(metadata$cultivar)  
coldata <- DataFrame(condition=conditions, cultivar=cultivar)
rownames(coldata) <- sample_names


dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ cultivar + condition)

#run the DESeq analysis
dds <- DESeq(dds)

results <- results(dds, contrast=c("condition", "treatment", "control"))


write.csv(as.data.frame(results), "DESeq2_results.csv")


plotMA(results, main="MA-plot, treatment vs control")
