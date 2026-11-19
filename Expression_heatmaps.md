# How to make expression heatmaps

## Significant DEG Heatmap
This plots genes that are significant DEGs for at least one treatment (leaves and roots together - produces very big heatmap!! beware)
```
##LARGE EXPRESSION HEATMAP FOR ALL LEAF AND ROOT SAMPLES
# COMBINED DEG HEATMAP (Leaves + Roots)
# Same file-reading logic as your barplot script, but instead of collapsing
# each comparison to an Up/Down count, this builds a genes x (nutrient+tissue)
# matrix of log2FC values for any gene significant in at least one comparison.

library(dplyr)
library(pheatmap)

padj_cutoff <- 0.05
lfc_cutoff  <- 1

# -----------------------------
# 1. Load leaves results into one long table (gene, comparison, log2FC, padj)
# -----------------------------
leaves_files <- list.files(pattern = "^DESeq2_.*_vs_Ctrl_leaves\\.csv$", ignore.case = TRUE)
leaves_files <- leaves_files[!grepl("ordered|sig_DEGs", leaves_files)]

leaves_long <- data.frame()
for (file in leaves_files) {
  res <- read.csv(file, row.names = 1)
  comparison <- gsub("DESeq2_|_vs_Ctrl_leaves\\.csv", "", file)

  res$gene <- rownames(res)
  res$Comparison <- paste0(comparison, "_Leaves")
  leaves_long <- rbind(leaves_long, res[, c("gene", "Comparison", "log2FoldChange", "padj")])
}

# -----------------------------
# 2. Load roots results the same way
# -----------------------------
roots_files <- list.files(pattern = "^DESeq2_.*_vs_Ctrl_roots\\.csv$", ignore.case = TRUE)
roots_files <- roots_files[!grepl("ordered|sig_DEGs", roots_files)]

roots_long <- data.frame()
for (file in roots_files) {
  res <- read.csv(file, row.names = 1)
  comparison <- gsub("DESeq2_|_vs_Ctrl_roots\\.csv", "", file)

  res$gene <- rownames(res)
  res$Comparison <- paste0(comparison, "_Roots")
  roots_long <- rbind(roots_long, res[, c("gene", "Comparison", "log2FoldChange", "padj")])
}

# -----------------------------
# 3. Combine leaves + roots into one long table
# -----------------------------
all_long <- rbind(leaves_long, roots_long)

# -----------------------------
# 4. Identify genes significant in AT LEAST ONE comparison
# -----------------------------
sig_genes <- all_long %>%
  filter(!is.na(padj) & padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff) %>%
  pull(gene) %>%
  unique()

cat("Number of genes significant in at least one comparison:", length(sig_genes), "\n")

# -----------------------------
# 5. Pivot to a wide gene x comparison matrix of log2FC values
# -----------------------------
# Genes not tested/present in a given comparison, or not significant there,
# will show their actual log2FC value (even if not significant) so you can see
# direction/magnitude everywhere the gene was measured. Adjust the filter below
# if you'd rather show NA for non-significant cells instead.

heatmap_long <- all_long %>% filter(gene %in% sig_genes)

heatmap_wide <- heatmap_long %>%
  select(gene, Comparison, log2FoldChange) %>%
  tidyr::pivot_wider(names_from = Comparison, values_from = log2FoldChange)

heatmap_matrix <- as.data.frame(heatmap_wide)
rownames(heatmap_matrix) <- heatmap_matrix$gene
heatmap_matrix$gene <- NULL
heatmap_matrix <- as.matrix(heatmap_matrix)

# Genes missing from a given comparison (e.g., filtered out entirely) will be
# NA - replace with 0 (no change) so pheatmap can cluster properly, OR leave as
# NA if you'd rather visually flag "not tested" cells (pheatmap shows NA as grey)
heatmap_matrix[is.na(heatmap_matrix)] <- 0

cat("Heatmap matrix dimensions:", nrow(heatmap_matrix), "genes x",
    ncol(heatmap_matrix), "comparisons\n")

# -----------------------------
# 6. Optional: if too many genes to display legibly, subset to top N by
# max absolute log2FC across all comparisons (uncomment to use)
# -----------------------------
# top_n <- 200
# gene_max_abs_fc <- apply(abs(heatmap_matrix), 1, max, na.rm = TRUE)
# top_genes <- names(sort(gene_max_abs_fc, decreasing = TRUE))[1:top_n]
# heatmap_matrix <- heatmap_matrix[top_genes, ]

# -----------------------------
# 7. Generate the heatmap
# -----------------------------
# Diverging palette centered at 0: blue = down, white = no change, red/orange = up
#
# IMPORTANT: cap the color scale at a reasonable log2FC range instead of using
# the raw max/min. A handful of extreme outlier genes (often low-count genes
# with inflated fold changes) will otherwise stretch the scale so far that the
# bulk of genes with typical, biologically meaningful fold changes all look
# washed out near white. Capping does NOT change your underlying data - it
# just means genes beyond the cap display as fully saturated color instead of
# an even more extreme shade.

fc_cap <- 5   # adjust based on your data - check summary(as.vector(heatmap_matrix))
              # to see where the bulk of your values actually sit before picking this

heatmap_matrix_capped <- heatmap_matrix
heatmap_matrix_capped[heatmap_matrix_capped > fc_cap] <- fc_cap
heatmap_matrix_capped[heatmap_matrix_capped < -fc_cap] <- -fc_cap

color_palette <- colorRampPalette(c("#4477AA", "white", "#EE6677"))(101)
breaks <- seq(-fc_cap, fc_cap, length.out = 101)

pdf("combined_leaves_roots_DEG_heatmap.pdf", width = 10, height = 14)
pheatmap(heatmap_matrix_capped,
         color = color_palette,
         breaks = breaks,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "average",
         show_rownames = FALSE,   # set TRUE if you want gene IDs labeled (only legible for small gene sets)
         fontsize_col = 10,
         border_color = NA)
dev.off()

# ---- Notes ----
# - show_rownames is set to FALSE by default since with hundreds/thousands of
#   genes, labels become unreadable. Set to TRUE only if you've subset to a
#   small curated gene list (e.g., your 64 overlapping DEGs, or the top_n
#   subset in step 6).
# - Genes are clustered (rows) so similar expression patterns across nutrients
#   group together - this is often the most informative part of this plot,
#   showing you clusters of genes that respond similarly across the dataset.
# - If the plot still feels too dense/tall with all significant genes included,
#   uncomment the top_n subsetting step above, or consider making this heatmap
#   specifically for your 64 overlapping DEGs instead of the full significant
#   gene set - that version would use overlap_gene_ids to filter heatmap_matrix
#   instead of using ALL significant genes across all comparisons.
```

This one plots genes that are significant DEGs for at least one treatment, leaves only (change to fit your needs)
```
# Combined correlation heatmap for all treatments (Leaves)
library(DESeq2)
library(pheatmap)

# Use the DESeq2 object you already created
dds_leaves <- nutrientAllLeaves_dds  # same object from your DESeq2 analysis

# Apply variance stabilizing transformation
vsd_leaves <- vst(dds_leaves, blind = TRUE)
vsd_leaves_mat <- assay(vsd_leaves)

# Compute Pearson correlation across all samples
cor_matrix <- cor(vsd_leaves_mat, method = "pearson")

# Define consistent color palette and breaks
breaks <- seq(0.90, 1.00, length.out = 101)  # fixed range for all plots

# Generate heatmap
png(paste0("new_Leaves_Combined_Correlation_Heatmap.png"), width = 3000, height = 2800, res = 300)
pheatmap(cor_matrix,
         main = "Sample Correlation - Leaves",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "average",
         breaks = breaks,
         border_color = NA)
dev.off()
```
