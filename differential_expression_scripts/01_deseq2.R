library(DESeq2)
library(tidyverse)

# Set the working directory (change to your directory)
setwd("/home/victoe/Genome_analysis/Data/DifferentialExpression")

# Load the metadata
metadata <- read.csv("/home/victoe/Genome_analysis/Data/Metadata/Metadata.csv")

# Ensure the order of count files matches the metadata
metadata <- metadata[order(metadata$Run),]  # Ensure this matches your count data order

# List all count files
count_files <- list.files(pattern="*.counts.txt$")
sample_names <- gsub(".counts.txt", "", count_files)

# Read count data and create a data frame with count data from all samples
count_data_list <- lapply(count_files, function(file) {
  read.table(file, header=FALSE, row.names=1)
})
names(count_data_list) <- sample_names
count_matrix <- do.call(cbind, count_data_list)

# Create a DataFrame for the sample conditions using metadata
conditions <- factor(metadata$tissue)
cultivar <- factor(metadata$cultivar)  # Assuming 'cultivar' is the column name
coldata <- DataFrame(condition=conditions, cultivar=cultivar)
rownames(coldata) <- sample_names

# Construct a DESeqDataSet object with an interaction term if needed
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ cultivar + condition)

# Run the DESeq analysis
dds <- DESeq(dds)

# Getting results, possibly for each level of interaction
# For example, results for condition across each cultivar
results <- results(dds, contrast=c("condition", "treatment", "control"))

# Save the results to a CSV file
write.csv(as.data.frame(results), "DESeq2_results.csv")

# Creating plots
plotMA(results, main="MA-plot, treatment vs control")
