library(DESeq2)
library(tidyverse)


original_metadata_path <- "/Users/victorenglof/Downloads/deseq2/Metadata.csv"
count_files_path <- "/Users/victorenglof/Downloads/deseq2/Musan King"

#loads the metadata
metadata <- read.csv(original_metadata_path, header = TRUE, row.names = 1)



#count data
count_files <- list.files(
  path = count_files_path,
  pattern = "*.txt",
  full.names = TRUE
)

#check the count files so they are correct
if (length(count_files) == 0) {
  stop("No count files found in the specified directory.")
}
print("Loaded count files:")
print(count_files)

#reads the count data
count_data_list <- lapply(count_files, function(file) {

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
# Run DESeq analysis
dds <- DESeq(dds)

res <- results(dds)

write.csv(as.data.frame(res), file.path("01_DESeq2_results_mk.csv"))


res_df <- as.data.frame(res)  

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
  ggtitle('Tissue comparison in Musang King')


ggsave("aril_comparison_multiplecultivar_Volcano.png", plot = volcano_plot, width = 10, height = 8,path="/Users/victorenglof/Downloads/deseq2")
print(volcano_plot)


