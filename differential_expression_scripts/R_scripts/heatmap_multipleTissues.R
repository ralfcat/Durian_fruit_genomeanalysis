library(DESeq2)
library(tidyverse)


original_metadata_path <- "/Users/victorenglof/Downloads/deseq2/Metadata.csv"
count_files_path <- "/Users/victorenglof/Downloads/deseq2/Musan King"

# Load metadata from the CSV file
metadata <- read.csv(original_metadata_path, header = TRUE, row.names = 1)



# Load count files
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
  # Read counts data from each file
  counts <- read.table(file, header = FALSE, col.names = c("gene", "counts"))
  
  # Filter out rows that do not represent typical gene identifiers
  valid_genes <- !grepl("^__", counts$gene)
  counts <- counts[valid_genes, ]
  
  # Convert counts to a named vector (gene counts)
  setNames(counts$counts, counts$gene)
})

# Extract base file names (sample identifiers) and assign to count data list
count_file_names <- gsub(".txt$", "", basename(count_files))
names(count_data_list) <- count_file_names

# Check the names of count files and their correspondence to metadata samples
print("Count file names (sample identifiers):")
print(count_file_names)
print("Metadata sample names:")
print(row.names(metadata))

# Check if there is any mismatch between count file names and metadata sample names
mismatch <- setdiff(row.names(metadata), count_file_names)
if (length(mismatch) > 0) {
  print("Mismatch between metadata and count files:")
  print(mismatch)
  # Optional: Handle mismatched samples in metadata, or count data as per your needs
  metadata <- metadata[row.names(metadata) %in% count_file_names,]
}

# Create count matrix
count_matrix <- do.call(cbind, count_data_list)

# Filter count matrix to include only samples in the filtered metadata
count_matrix <- count_matrix[, row.names(metadata), drop = FALSE]

# Convert tissue column to a factor (if necessary)
metadata$tissue <- as.factor(metadata$tissue)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata,
                              design = ~ tissue)

keep <- rowSums(counts(dds) >= 10) >= 2  # Adjust parameters as needed
dds <- dds[keep,]
# Run DESeq analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# Optional: Plotting or additional analyses can be added here
# For example, creating an MA plot
plotMA(dds, main = "MA-plot")

# Optional: Save the filtered metadata if needed
# write.csv(metadata, "/path/to/filtered_metadata.csv")

# Creating plots, such as MA plot
plotMA(dds, main="MA-plot, treatment vs control")  # 'plotMA' function should take 'dds', not 'res'

# Preparing data for the volcano plot
res_df <- as.data.frame(res)  # Convert 'res' to a data frame for plotting

# Generate the volcano plot using ggplot2
library(ggplot2)
# Dynamically set the upper limit based on the maximum value found, add a small buffer
max_padj <- max(-log10(res_df$padj), na.rm = TRUE)
min_padj <- min(-log10(res_df$padj), na.rm = TRUE)
upper_limit <- ceiling(max_padj + 1)
library(ggrepel)
# Convert 'res' to a data frame for plotting


# Convert 'res' to a data frame for plotting


# Ensure that gene names are included as a column for plotting
res_df$gene <- rownames(res)

# Add a new column 'diffexpressed' to classify the genes
# Adjust the thresholds as per your specific criteria
res_df$diffexpressed <- ifelse(res_df$padj < 0.05 & res_df$log2FoldChange > 1, "Upregulated",
                               ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < -1, "Downregulated", "Not significant"))
library(DESeq2)
# heatmap 
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

# Printing names and details for verification
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

# Set the path for output directory, replace 'output_dir' with your actual directory path
output_dir <- "/path/to/your/output/directory"
output_heatmap_path <- file.path(output_dir, "01_heatmap_top_10_significant_genes.png")
png(output_heatmap_path, width = 10, height = 8, units = "in", res = 300)
print(heatmap_plot)
dev.off()  # Don't forget to close the PNG device

