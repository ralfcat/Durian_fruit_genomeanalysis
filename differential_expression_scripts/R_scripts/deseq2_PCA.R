library(DESeq2)
library(tidyverse)


original_metadata_path <- "/Users/victorenglof/Downloads/deseq2/Metadata.csv"
count_files_path <- "/Users/victorenglof/Downloads/deseq2"

#metadata
metadata <- read.csv(original_metadata_path, header = TRUE, row.names = 1)



#count files
count_files <- list.files(
  path = count_files_path,
  pattern = "*.txt",
  full.names = TRUE
)

#make sure the count files are correct
if (length(count_files) == 0) {
  stop("No count files found in the specified directory.")
}
print("Loaded count files:")
print(count_files)


count_data_list <- lapply(count_files, function(file) {
  
  counts <- read.table(file, header = FALSE, col.names = c("gene", "counts"))
  
  # Filter out rows that do not represent typical gene identifiers
  valid_genes <- !grepl("^__", counts$gene)
  counts <- counts[valid_genes, ]
  
  # Convert counts to a named vector (gene counts)
  setNames(counts$counts, counts$gene)

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
  metadata <- metadata[row.names(metadata) %in% count_file_names,]
}

#Count matrix
count_matrix <- do.call(cbind, count_data_list)


count_matrix <- count_matrix[, row.names(metadata), drop = FALSE]


metadata$Cultivar <- as.factor(metadata$Cultivar)


dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata,
                              design = ~ Cultivar)

keep <- rowSums(counts(dds) >= 10) >= 2  # Adjust parameters as needed
dds <- dds[keep,]
# Run DESeq analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# Save the results to a CSV file
write.csv(as.data.frame(res), file.path("01_DESeq2_results_mk.csv"))

# Preparing data for the volcano plot
res_df <- as.data.frame(res) 

# Generate the volcano plot using ggplot2
library(ggplot2)

max_padj <- max(-log10(res_df$padj), na.rm = TRUE)
min_padj <- min(-log10(res_df$padj), na.rm = TRUE)
upper_limit <- ceiling(max_padj + 1)
library(ggrepel)

res_df$gene <- rownames(res)
res_df$diffexpressed <- ifelse(res_df$padj < 0.05 & res_df$log2FoldChange > 1, "Upregulated",
                               ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < -1, "Downregulated", "Not significant"))


volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = diffexpressed)) +
  geom_point(size = 2) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  scale_color_manual(values = c("Not significant" = "grey", "Upregulated" = "#00AFBB", "Downregulated" = "#bb0c00")) +
  coord_cartesian(ylim = c(0, upper_limit), xlim = c(-12, 12)) +
  labs(color = 'Expression Status', x = expression(log[2]*" fold change"), y = expression("-log"[10]*" adjusted p-value")) +
  ggtitle('Aril comparison in different cultivars')


ggsave("aril_comparison_multipleCultivars_Volcano.png", plot = volcano_plot, width = 10, height = 8,path="/Users/victorenglof/Downloads/deseq2")
print(volcano_plot)

#count the number of "Upregulated" entries, ignoring NA values
upregulated_count <- sum(res_df$diffexpressed == "Upregulated", na.rm = TRUE)


print(upregulated_count)



# Load necessary libraries
library(DESeq2)
library(ggplot2)

# Ensure your DESeq2 object 'dds' has been setup and run
dds <- DESeq(dds)
# Use Variance Stabilizing Transformation instead of rlog
vst_data <- varianceStabilizingTransformation(dds, blind = FALSE)

# Check column names in metadata
print("Column names in metadata:")
print(colnames(colData(dds)))

# Perform PCA on the VST data, make sure to replace 'condition' and 'type' 
# with actual metadata column names that exist in your colData(dds)
# Assuming 'Cultivar' is used for color and 'Tissue' for shape
pcaData <- plotPCA(vst_data, intgroup=c("Cultivar", "tissue"), returnData=TRUE)

# Percent variance explained by the PCs
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Plot the PCA with tissue type included
pcaPlot <- ggplot(pcaData, aes(PC1, PC2, color=Cultivar, shape=tissue)) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(title="PCA of VST Data") +
  scale_shape_manual(values=1:length(unique(pcaData$tissue))) + # Optional: Customize shapes if needed
  theme_minimal()

# Print the PCA plot
print(pcaPlot)


