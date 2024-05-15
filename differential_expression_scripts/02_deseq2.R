library(DESeq2)
library(tidyverse)

setwd("/home/victoe/Genome_analysis/Data/DifferentialExpression")

#metadata
metadata_path <- "/home/victoe/Genome_analysis/Data/Metadata/Metadata.csv"
metadata <- read.csv(metadata_path)

#lists all count files
count_files <- list.files(pattern="*.txt$")
if (length(count_files) == 0) {
  stop("No count files found. Check the file path and pattern.")
}

# Read count data
count_data_list <- lapply(count_files, function(file) {
  filepath <- file.path(getwd(), file)
  if (file.info(filepath)$size == 0) {
    stop(paste("File is empty:", filepath))
  }
  counts <- read.table(filepath, header=FALSE, col.names=c("gene", "counts"))

  #remove the rows in the count files that are not necessary
  valid_genes <- !grepl("^__", counts$gene)  # genes not starting with "__"
  counts <- counts[valid_genes, ]

  # Ensure the counts are returned as a vector with gene names as names
  setNames(counts$counts, counts$gene)
})
names(count_data_list) <- gsub(".txt$", "", count_files)


metadata <- metadata[match(names(count_data_list), metadata$Run),]
if (anyNA(metadata)) {
  stop("Some samples in metadata do not have matching count files.")
}

#count matrix
count_matrix <- do.call(cbind, count_data_list)


conditions <- factor(metadata$tissue)
coldata <- DataFrame(condition = conditions)
rownames(coldata) <- names(count_data_list)


dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ condition)

# Run the DESeq analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds)  # This calls the results function and stores the output in 'res'


write.csv(as.data.frame(res), file.path(getwd(), "DESeq2_results.csv"))  # Convert 'res' to a dataframe and save it


res_df <- as.data.frame(res)  

#adding a 'diffexpressed' column based on log2 fold change and adjusted p-values
res_df$diffexpressed <- ifelse(res_df$log2FoldChange >= 0.6 & res_df$padj < 0.05, "Upregulated",
                               ifelse(res_df$log2FoldChange <= -0.6 & res_df$padj < 0.05, "Downregulated", "Not significant"))

res_df$delabel <- ifelse(res_df$diffexpressed == "Upregulated", as.character(res_df$gene), "")




# Generate the volcano plot using ggplot2
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

#plots the volcano plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("grey", "#00AFBB", "#bb0c00"),
                     labels = c("Not significant", "Downregulated", "Upregulated")) +
  coord_cartesian(ylim = c(0, 3), xlim = c(-10, 10)) +
  labs(color = 'Expression Status',
       x = expression("log"[2]*" fold change"), y = expression("-log"[10]*" adjusted p-value")) +
  scale_x_continuous(breaks = seq(-10, 10, 2)) +
  ggtitle('Durian fruit volcano plot') +
  geom_text_repel(max.overlaps = Inf)


ggsave("volcano_plot_best.png", plot = volcano_plot, width = 10, height = 8)

