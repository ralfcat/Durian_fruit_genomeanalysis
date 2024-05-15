library(DESeq2)
library(tidyverse)


original_metadata_path <- "/Users/victorenglof/Downloads/deseq2/Metadata.csv"
count_files_path <- "/Users/victorenglof/Downloads/deseq2/Musan King"

#sets the metadata
metadata <- read.csv(original_metadata_path, header = TRUE, row.names = 1)



#count files
count_files <- list.files(
  path = count_files_path,
  pattern = "*.txt",
  full.names = TRUE
)

# Check count files and their names
if (length(count_files) == 0) {
  stop("No count files found in the specified directory.")
}
print("Loaded count files:")
print(count_files)

# Read count data
count_data_list <- lapply(count_files, function(file) {
  counts <- read.table(file, header = FALSE, col.names = c("gene", "counts"))
  
  # Filter out rows that do not represent typical gene identifiers
  valid_genes <- !grepl("^__", counts$gene)
  counts <- counts[valid_genes, ]
  
  # Convert counts to a named vector (gene counts)
  setNames(counts$counts, counts$gene)
})


count_file_names <- gsub(".txt$", "", basename(count_files))
names(count_data_list) <- count_file_names

#check the names of count files
print("Count file names (sample identifiers):")
print(count_file_names)
print("Metadata sample names:")
print(row.names(metadata))

# Check if there is any mismatch between count file names and metadata sample names
mismatch <- setdiff(row.names(metadata), count_file_names)
if (length(mismatch) > 0) {
  print("Mismatch between metadata and count files:")
  print(mismatch)
  metadata <- metadata[row.names(metadata) %in% count_file_names,]
}

#count matrix
count_matrix <- do.call(cbind, count_data_list)


count_matrix <- count_matrix[, row.names(metadata), drop = FALSE]


metadata$tissue <- as.factor(metadata$tissue)

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata,
                              design = ~ tissue)

keep <- rowSums(counts(dds) >= 10) >= 2  
dds <- dds[keep,]
dds <- DESeq(dds)

res <- results(dds)

# Preparing data for the volcano plot
res_df <- as.data.frame(res) 


library(ggplot2)
max_padj <- max(-log10(res_df$padj), na.rm = TRUE)
min_padj <- min(-log10(res_df$padj), na.rm = TRUE)
upper_limit <- ceiling(max_padj + 1)
library(ggrepel)


res_df$gene <- rownames(res)

res_df$diffexpressed <- ifelse(res_df$padj < 0.05 & res_df$log2FoldChange > 1, "Upregulated",
                               ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < -1, "Downregulated", "Not significant"))
library(DESeq2)
library(pheatmap)

# Identify genes with significant differential expression
significant_genes <- res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]

# Filter out rows with NA values from significant_genes
significant_genes_filtered <- significant_genes[complete.cases(significant_genes), ]

# Sort genes by adjusted p-value and select the top 10
top_significant_genes <- significant_genes_filtered[order(significant_genes_filtered$padj), ][1:10, ]

# Get row names of the top 10 filtered significant genes
top_significant_genes_list <- rownames(top_significant_genes)

# Extract normalized counts from the DESeqDataSet
normalized_counts <- counts(dds, normalized = TRUE)

# Subset the normalized counts matrix to include only the top 10 filtered significant genes
filtered_normalized_counts <- normalized_counts[top_significant_genes_list, , drop = FALSE]

# Printing names
print('Top 10 significant gene names:')
print(head(rownames(top_significant_genes)))
print(top_significant_genes)

# Generate the heatmap
heatmap_plot <- pheatmap(filtered_normalized_counts,
                         scale = "row",  # Scale rows (genes) for better visualization
                         clustering_method = "average",  # Clustering method (e.g., "average", "ward.D2")
                         cluster_cols = TRUE,  # Cluster samples (columns)
                         cluster_rows = TRUE,  # Cluster genes (rows)
                         show_rownames = TRUE,  # Show gene names (row names)
                         show_colnames = TRUE,  # Show sample names (column names)
                         fontsize_row = 5,  # Font size for gene names
                         fontsize_col = 8,  # Font size for sample names
                         main = "Heatmap of Top 10 Normalized Counts")  # Title of the heatmap

output_dir <- "/path/to/your/output/directory"
output_heatmap_path <- file.path(output_dir, "01_heatmap_top_10_significant_genes.png")
png(output_heatmap_path, width = 10, height = 8, units = "in", res = 300)
print(heatmap_plot)


