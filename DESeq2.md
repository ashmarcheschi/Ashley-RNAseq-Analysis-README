# Running differential expression analysis with DESeq2

## Running DESeq2
This script utilizes LFC shrinkage, which adjusts extreme, noisy fold-change estimates with low read count/high variability - can be run without. The inputs are the counts matrix and the sample information table. It also creates files with ordered DEGs, and pulls those that are significant by your cutoff:
```
#RUNNING DESeq2 with LFC shrink
library(DESeq2)

# Load data
nutrientAllLeaves_count_data <- read.csv("L_combined_counts_matrix.tsv", sep = "\t", row.names = 1)
nutrientAllLeaves_coldata <- read.csv("L_combined_coldata.tsv", sep = "\t", row.names = 1)

# Format
nutrientAllLeaves_count_data <- round(nutrientAllLeaves_count_data)
nutrientAllLeaves_coldata$condition <- factor(nutrientAllLeaves_coldata$condition)

# DESeq2 object
nutrientAllLeaves_dds <- DESeqDataSetFromMatrix(
  countData = nutrientAllLeaves_count_data,
  colData = nutrientAllLeaves_coldata,
  design = ~ condition
)

# Run DESeq2
nutrientAllLeaves_dds <- DESeq(nutrientAllLeaves_dds)
conditions_leaves <- levels(nutrientAllLeaves_dds$condition)
control <- "Ctrl"
leaves_results_list <- list()
leaves_sig_degs_list <- list()
leaves_results_ordered_list <- list()

# Pairwise DE + export
for (cond in conditions_leaves) {
  if (cond != control) {
    res <- results(nutrientAllLeaves_dds, contrast = c("condition", cond, control))
    resLFC <- lfcShrink(nutrientAllLeaves_dds, contrast = c("condition", cond, control), res = res, type = "ashr")

    leaves_results_list[[cond]] <- resLFC
    write.csv(as.data.frame(resLFC), paste0("DESeq2_", cond, "_vs_Ctrl_leaves.csv"))
    
    # Annotate with transcript + gene info
    resLFC_df <- as.data.frame(resLFC) %>%
    rownames_to_column(var = "gene_id")
    leaf_degs_annotated <- resLFC_df %>%
    inner_join(mapping, by = "gene_id") %>%
    select(transcript_id, gene_id, gene_name, product, everything())
    write.csv(leaf_degs_annotated, paste0("DESeq2_", cond, "_vs_Ctrl_leaves_annotated.csv"), row.names = FALSE)

    res_sorted <- resLFC[order(resLFC$padj, na.last = TRUE), ]
    leaves_results_ordered_list[[cond]] <- res_sorted
    write.csv(res_sorted, paste0("DESeq2_", cond, "_vs_Ctrl_leaves_ordered.csv"))

    sig <- resLFC[which(resLFC$padj < 0.05 & abs(resLFC$log2FoldChange) > 1), ]
    leaves_sig_degs_list[[cond]] <- sig
    write.csv(as.data.frame(sig), paste0("sig_DEGs_", cond, "_vs_Ctrl_leaves.csv"))
  }
}
```

