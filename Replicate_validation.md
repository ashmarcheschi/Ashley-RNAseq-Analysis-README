# How to validate replicates

## PCA Plots
Plots variance between samples, this script VST transforms them and plots them in two different PCA plots due to sample size
```
library(DESeq2)
library(ggplot2)
library(ggrepel)

# Create output directory
dir.create("PCA_Plots_Leaves", showWarnings = FALSE)

# VST transform
vsd_leaves <- vst(nutrientAllLeaves_dds, blind = FALSE)
pcaData <- plotPCA(vsd_leaves, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Get all conditions
all_conditions <- sort(unique(pcaData$condition))
ctrl_condition <- "Ctrl"
nutrient_conditions <- setdiff(all_conditions, ctrl_condition)

# Set factor levels so "Ctrl" is first
pcaData$condition <- factor(pcaData$condition, levels = c(ctrl_condition, nutrient_conditions))
# Add replicate info to pcaData
pcaData$replicate <- colData(nutrientAllLeaves_dds)$replicate

# Now split into two sets
set1_conditions <- nutrient_conditions[1:6]
set2_conditions <- nutrient_conditions[7:12]
plot1_conditions <- c(ctrl_condition, set1_conditions)
plot2_conditions <- c(ctrl_condition, set2_conditions)

# ✅ Use a truly colorblind-friendly palette (Okabe-Ito)
colorblind_palette_7 <- c(
  "black",        # Ctrl
  "#4477aa",      # orange
  "#EE6677",      # sky blue
  "#228833",      # bluish green
  "#CCBB44",      # yellow
  "#66CCEE",      # blue
  "#AA3377"       # vermillion (orange-red)
)

# Use 7 visually distinct shapes
shape_palette <- c(0, 1, 2, 3, 4, 5, 6)

# --- Plot 1 ---
pca_subset1 <- pcaData[pcaData$condition %in% plot1_conditions, , drop = FALSE]
pca_plot_file1 <- file.path("PCA_Plots_Leaves", "leaves_PCA_plot_color_shape_set1.pdf")

pdf(pca_plot_file1, width = 8, height = 6)
ggplot(pca_subset1, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(aes(label = replicate), size = 3, max.overlaps = 100, show.legend = FALSE) +
  scale_shape_manual(values = shape_palette[1:length(plot1_conditions)]) +
  scale_color_manual(values = colorblind_palette_7[1:length(plot1_conditions)]) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA Plot: Leaves (Set 1)",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    shape = "Condition",
    color = "Condition"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold")
  )

dev.off()

# --- Plot 2 ---
pca_subset2 <- pcaData[pcaData$condition %in% plot2_conditions, , drop = FALSE]
pca_plot_file2 <- file.path("PCA_Plots_Leaves", "leaves_PCA_plot_color_shape_set2.pdf")

pdf(pca_plot_file2, width = 8, height = 6)
ggplot(pca_subset2, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(aes(label = replicate), size = 3, max.overlaps = 100, show.legend = FALSE) +
  scale_shape_manual(values = shape_palette[1:length(plot2_conditions)]) +
  scale_color_manual(values = colorblind_palette_7[1:length(plot2_conditions)]) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA Plot: Leaves (Set 2)",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    shape = "Condition",
    color = "Condition"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold")
  )

dev.off()

```

## Sample Correlations
Here is code to make one big giant combined correlation heatmap:
```
#one big happy correlation heatmap with leaf and root samples combined 

# Combined correlation heatmap for all treatments (Leaves + Roots together)
library(DESeq2)
library(pheatmap)

# ---- 1. Combine raw counts from both DESeq2 objects ----
# This assumes nutrientAllLeaves_dds and nutrientAllRoots_dds were built from
# the same gene set (same rows/gene order). If they weren't, intersect the
# rownames first - see note at the bottom of the script.

heatmap_counts_leaves <- counts(nutrientAllLeaves_dds)
heatmap_counts_roots  <- counts(nutrientAllRoots_dds)

stopifnot(identical(rownames(heatmap_counts_leaves), rownames(heatmap_counts_roots)))

combined_counts <- cbind(heatmap_counts_leaves, heatmap_counts_roots)

# ---- 2. Build combined colData with a tissue label ----
coldata_leaves <- as.data.frame(colData(nutrientAllLeaves_dds))
coldata_roots  <- as.data.frame(colData(nutrientAllRoots_dds))

coldata_leaves$tissue <- "Leaf"
coldata_roots$tissue  <- "Root"

# keep only shared columns so rbind doesn't complain, then add tissue back
shared_cols <- intersect(colnames(coldata_leaves), colnames(coldata_roots))
combined_coldata <- rbind(coldata_leaves[, shared_cols], coldata_roots[, shared_cols])

# sanity check: sample order must match combined_counts columns
stopifnot(identical(rownames(combined_coldata), colnames(combined_counts)))

# ---- 3. Build a combined DESeqDataSet and run VST ----
dds_combined <- DESeqDataSetFromMatrix(
  countData = combined_counts,
  colData   = combined_coldata,
  design    = ~1   # design doesn't matter here since we're just using vst() for visualization
)

vsd_combined <- vst(dds_combined, blind = TRUE)
vsd_combined_mat <- assay(vsd_combined)

# ---- 4. Compute Pearson correlation across all samples (leaf + root) ----
cor_matrix <- cor(vsd_combined_mat, method = "pearson")

# ---- 5. Set up annotation bar so leaf/root samples are visible on the heatmap ----
annotation_col <- data.frame(Tissue = combined_coldata$tissue)
rownames(annotation_col) <- rownames(combined_coldata)

ann_colors <- list(Tissue = c(Leaf = "#4DAF4A", Root = "#A65628"))

# ---- 6. Define consistent color palette and breaks ----
breaks <- seq(0.60, 1.00, length.out = 101)  # fixed range, same as your original plots

# ---- 7. Generate heatmap ----
png("new_Combined_Leaf_Root_Correlation_Heatmap_unclustered.png", width = 4500, height = 4200, res = 300)
pheatmap(cor_matrix,
         main = "Sample Correlation - Leaves & Roots",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         breaks = breaks,
         border_color = NA,
         annotation_col = annotation_col,
         annotation_colors = ann_colors)
dev.off()

```

## Leave-one-out analysis 
This is what I did to do a leave-one-out (LOO) analysis to see if there were any replicates that were outliers in my data that could disproportionately influence the expression data. Wouldn't say this is necessary, I simply ran it to err on the side of caution!
```
BOTH tissues (Leaves + Roots), producing one master summary table.
#
# WARNING: this reruns DESeq2 multiple times per condition (1 full run + 1 per
# replicate dropped). With 12 nutrients x 2 tissues x ~5 fits each, this is
# computationally heavy - expect a long runtime. Consider running this on your
# supercomputer rather than locally, and/or submitting as a batch job.

library(DESeq2)
library(dplyr)

# ---- Config ----
condition_col <- "condition"
ctrl_label    <- "Ctrl"
padj_cutoff_for_tiers <- c(0.05, 0.01, 0.001)  # confidence tiers to test correlation within
lfc_drift_threshold   <- 0.5  # used in the "why did genes lose significance" diagnostic

# The two full dds objects (all conditions, all reps) for each tissue
tissue_dds_list <- list(
  Leaves = nutrientAllLeaves_dds,
  Roots  = nutrientAllRoots_dds
)

# ---- Function: run the full LOO log2FC correlation analysis for ONE condition + ONE tissue ----
# IMPORTANT FIX: rather than subsetting down to just Ctrl + target condition
# (which loses the pooled dispersion information from all other conditions and
# can badly underpower weaker-responding nutrients), this version drops
# individual replicates from the FULL multi-condition dds object and refits
# DESeq2 on the full design each time, then pulls the same Ctrl-vs-target
# contrast. This matches how your original pipeline (and its 46 S-leaf DEGs)
# was actually generated.

run_loo_for_condition <- function(dds_full, condition_col, ctrl_label, target_label, tissue_name) {
  # NOTE: dds_full here is expected to already be a FITTED DESeq2 object
  # (i.e., DESeq() already run on it) - fit once per tissue outside this
  # function and reuse across all nutrient contrasts, since refitting the
  # identical full model for every nutrient would be redundant.

  cat("\n=====================================================\n")
  cat("Processing:", tissue_name, "-", target_label, "vs", ctrl_label, "\n")
  cat("=====================================================\n")

  target_sample_names <- colnames(dds_full)[colData(dds_full)[[condition_col]] == target_label]
  cat("Replicates found for", target_label, ":", paste(target_sample_names, collapse = ", "), "\n")

  # ---- Pull the contrast from the already-fitted full model ----
  res_full <- results(dds_full, contrast = c(condition_col, target_label, ctrl_label))
  res_full <- lfcShrink(dds_full, contrast = c(condition_col, target_label, ctrl_label),
                         res = res_full, type = "ashr")

  full_df <- as.data.frame(res_full)
  full_df$gene <- rownames(full_df)
  full_df <- full_df %>% select(gene, log2FC_full = log2FoldChange, padj_full = padj)

  # ---- LOO runs: drop each replicate from the FULL multi-condition dataset,
  # refit DESeq2 fresh each time (this refit is necessary since the sample set
  # actually changes), then pull the same contrast ----
  loo_results <- list()

  for (dropped_sample in target_sample_names) {
    keep_this_round <- colnames(dds_full) != dropped_sample
    dds_loo <- dds_full[, keep_this_round]
    colData(dds_loo)[[condition_col]] <- droplevels(colData(dds_loo)[[condition_col]])

    dds_loo_fit <- DESeq(dds_loo)
    res_loo <- results(dds_loo_fit, contrast = c(condition_col, target_label, ctrl_label))
    res_loo <- lfcShrink(dds_loo_fit, contrast = c(condition_col, target_label, ctrl_label),
                          res = res_loo, type = "ashr")

    loo_results[[dropped_sample]] <- res_loo
  }

  # ---- Merge full + all LOO log2FC values into one table ----
  merged_lfc <- full_df
  for (dropped_sample in names(loo_results)) {
    loo_df <- as.data.frame(loo_results[[dropped_sample]])
    loo_df$gene <- rownames(loo_df)
    lfc_col  <- paste0("log2FC_without_", dropped_sample)
    padj_col <- paste0("padj_without_", dropped_sample)
    loo_df_sub <- loo_df %>% select(gene, !!lfc_col := log2FoldChange, !!padj_col := padj)
    merged_lfc <- merge(merged_lfc, loo_df_sub, by = "gene", all = TRUE)
  }

  lfc_loo_cols <- grep("^log2FC_without_", colnames(merged_lfc), value = TRUE)

  # ---- Correlations: all genes + each confidence tier ----
  condition_results <- data.frame()

  # All genes (no padj filter)
  for (col in lfc_loo_cols) {
    dropped_label <- gsub("log2FC_without_", "", col)
    valid <- !is.na(merged_lfc$log2FC_full) & !is.na(merged_lfc[[col]])
    if (sum(valid) >= 3) {
      r <- cor(merged_lfc$log2FC_full[valid], merged_lfc[[col]][valid], method = "pearson")
      condition_results <- rbind(condition_results, data.frame(
        tissue = tissue_name, condition = target_label, dropped_sample = dropped_label,
        tier = "all_genes", n_genes = sum(valid), pearson_r = r
      ))
    }
  }

  # Confidence tiers
  for (cutoff in padj_cutoff_for_tiers) {
    tier_genes <- merged_lfc %>% filter(!is.na(padj_full) & padj_full < cutoff)
    tier_name <- paste0("padj_lt_", cutoff)

    for (col in lfc_loo_cols) {
      dropped_label <- gsub("log2FC_without_", "", col)
      valid <- !is.na(tier_genes$log2FC_full) & !is.na(tier_genes[[col]])
      if (sum(valid) >= 3) {
        r <- cor(tier_genes$log2FC_full[valid], tier_genes[[col]][valid], method = "pearson")
        condition_results <- rbind(condition_results, data.frame(
          tissue = tissue_name, condition = target_label, dropped_sample = dropped_label,
          tier = tier_name, n_genes = sum(valid), pearson_r = r
        ))
      }
    }
  }

  return(condition_results)
}

# ---- Loop across all tissues and all non-control conditions ----
all_loo_summary <- data.frame()

for (tissue_name in names(tissue_dds_list)) {
  cat("\n\n### Fitting full multi-condition model for", tissue_name, "(once, reused for all nutrients) ###\n")
  dds_full_fitted <- DESeq(tissue_dds_list[[tissue_name]])

  all_conditions <- levels(colData(dds_full_fitted)[[condition_col]])
  nutrient_conditions <- setdiff(all_conditions, ctrl_label)

  for (nutrient in nutrient_conditions) {
    result <- tryCatch({
      run_loo_for_condition(dds_full_fitted, condition_col, ctrl_label, nutrient, tissue_name)
    }, error = function(e) {
      cat("ERROR processing", tissue_name, nutrient, ":", conditionMessage(e), "\n")
      return(data.frame())
    })
    all_loo_summary <- rbind(all_loo_summary, result)
  }
}

# ---- Save master summary table ----
write.csv(all_loo_summary, "All_Treatments_LOO_log2FC_correlation_summary.csv", row.names = FALSE)

# ---- Quick overview: average correlation per condition/tissue/tier ----
overview <- all_loo_summary %>%
  group_by(tissue, condition, tier) %>%
  summarise(mean_r = mean(pearson_r), min_r = min(pearson_r), max_r = max(pearson_r), .groups = "drop")

write.csv(overview, "All_Treatments_LOO_log2FC_correlation_overview.csv", row.names = FALSE)

cat("\n\n=== Overview: mean/min/max Pearson R per condition/tissue/tier ===\n")
print(overview, n = Inf)

```

I also have this script that makes a heatmap with the LOO data - maybe double check this to make sure it connects and everything is named correctly:
```
# Heatmap of LOO log2FC correlation overview across all treatments and tissues
# Reads the overview CSV produced by the prior script and visualizes mean
# Pearson R (Full vs LOO runs) as a heatmap: rows = condition+tissue, columns = tier

library(dplyr)
library(tidyr)
library(pheatmap)

# ---- 1. Load the overview table ----
overview <- read.csv("All_Treatments_LOO_log2FC_correlation_overview.csv")

# ---- 2. Create a combined row label (condition + tissue) ----
overview$condition_tissue <- paste0(overview$condition, "_", overview$tissue)

# ---- 3. Set a sensible tier order (rather than alphabetical) ----
tier_order <- c("all_genes", "padj_lt_0.05", "padj_lt_0.01", "padj_lt_0.001")
overview$tier <- factor(overview$tier, levels = tier_order)

# ---- 4. Pivot to wide: rows = condition_tissue, columns = tier, values = mean_r ----
heatmap_wide <- overview %>%
  select(condition_tissue, tier, mean_r) %>%
  pivot_wider(names_from = tier, values_from = mean_r)

heatmap_matrix <- as.data.frame(heatmap_wide)
rownames(heatmap_matrix) <- heatmap_matrix$condition_tissue
heatmap_matrix$condition_tissue <- NULL
heatmap_matrix <- heatmap_matrix[, tier_order[tier_order %in% colnames(heatmap_matrix)]]
heatmap_matrix <- as.matrix(heatmap_matrix)

# ---- 5. Build a tissue annotation bar (Leaf vs Root) ----
annotation_row <- overview %>%
  distinct(condition_tissue, tissue) %>%
  as.data.frame()
rownames(annotation_row) <- annotation_row$condition_tissue
annotation_row$condition_tissue <- NULL

ann_colors <- list(tissue = c(Leaves = "#4DAF4A", Roots = "#A65628"))

# ---- 6. Color scale: correlation is bounded 0-1, so fix the scale accordingly ----
# Adjust the lower bound if your correlations dip below 0.5 - check
# range(heatmap_matrix, na.rm = TRUE) first
color_palette <- colorRampPalette(c("#4477AA", "white", "#EE6677"))(101)
breaks <- seq(0.5, 1.0, length.out = 101)

# ---- 7. Generate the heatmap ----
png("All_Treatments_LOO_log2FC_correlation_heatmap.png", width = 3200, height = 3600, res = 300)
pheatmap(heatmap_matrix,
         color = color_palette,
         breaks = breaks,
         main = "Mean Log2FC Correlation (Full vs LOO) by Treatment, Tissue, and Confidence Tier",
         cluster_rows = FALSE,
         cluster_cols = FALSE,   # keep tier order fixed (all_genes -> most stringent)
         clustering_distance_rows = "euclidean",
         clustering_method = "average",
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         border_color = NA,
         fontsize_row = 9,
         fontsize_col = 10,
         display_numbers = TRUE,
         number_format = "%.2f")
dev.off()

```

