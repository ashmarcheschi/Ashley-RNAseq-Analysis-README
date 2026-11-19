# How to make various RNAseq analysis plots/figures

## DEG Barplot (from DESeq2 results)
This barplot is made from all combined data based on a DEG significance cutoff of p-adj < 0.05 and |log2FC| >= 1.
```
#DEG BARPLOT FOR MANUSCRIPT - BIGGER!!!
#COMBINED DEG BARPLOT (enlarged for visibility)
library(dplyr)
library(ggplot2)
# -----------------------------
# 1. Load and summarize leaves
# -----------------------------
leaves_files <- list.files(pattern = "^DESeq2_.*_vs_Ctrl_leaves\\.csv$", ignore.case = TRUE)
leaves_files <- leaves_files[!grepl("ordered|sig_DEGs", leaves_files)]
leaves_summary <- data.frame()
for (file in leaves_files) {
  res <- read.csv(file, row.names = 1)
  comparison <- gsub("DESeq2_|_vs_Ctrl_leaves\\.csv", "", file)
  res_filtered <- res[!is.na(res$padj), ]
  up <- sum(res_filtered$log2FoldChange > 1 & res_filtered$padj < 0.05)
  down <- sum(res_filtered$log2FoldChange < -1 & res_filtered$padj < 0.05)
  leaves_summary <- rbind(leaves_summary,
                          data.frame(Comparison = comparison, Regulation = "Up", Count = up, Tissue = "Leaves"),
                          data.frame(Comparison = comparison, Regulation = "Down", Count = down, Tissue = "Leaves"))
}
# -----------------------------
# 2. Load and summarize roots
# -----------------------------
roots_files <- list.files(pattern = "^DESeq2_.*_vs_Ctrl_roots\\.csv$", ignore.case = TRUE)
roots_files <- roots_files[!grepl("ordered|sig_DEGs", roots_files)]
roots_summary <- data.frame()
for (file in roots_files) {
  res <- read.csv(file, row.names = 1)
  comparison <- gsub("DESeq2_|_vs_Ctrl_roots\\.csv", "", file)
  res_filtered <- res[!is.na(res$padj), ]
  up <- sum(res_filtered$log2FoldChange > 1 & res_filtered$padj < 0.05)
  down <- sum(res_filtered$log2FoldChange < -1 & res_filtered$padj < 0.05)
  roots_summary <- rbind(roots_summary,
                         data.frame(Comparison = comparison, Regulation = "Up", Count = up, Tissue = "Roots"),
                         data.frame(Comparison = comparison, Regulation = "Down", Count = down, Tissue = "Roots"))
}
# -----------------------------
# 3. Combine leaves and roots
# -----------------------------
# Combine leaves and roots summaries as before
deg_summary <- rbind(leaves_summary, roots_summary)
# Create a new combined category for Tissue + Regulation
deg_summary$Category <- paste(deg_summary$Tissue, deg_summary$Regulation, sep = "_")
# Define colors for the four categories
category_colors <- c(
  "Leaves_Up" = "#AAAA00",
  "Leaves_Down" = "#EEDD88",
  "Roots_Up" = "#99DDFF",
  "Roots_Down" = "#77AADD"
)
# Plot: x-axis = Comparison (nutrient), bars = Category
p <- ggplot(deg_summary, aes(x = Comparison, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.95), width = 0.9) +
  scale_fill_manual(values = category_colors) +
  labs(
       x = "Nutrient Deficiency",
       y = "Number of DEGs",
       fill = "Category") +
  theme_minimal(base_size = 20) +   # increased base_size scales most text elements up
  theme(
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 20),
    axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5),   # larger nutrient labels
    axis.text.y = element_text(size = 16),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.2, "lines"),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  geom_text(aes(label = Count),
            position = position_dodge(width = 0.95),
            vjust = -0.5, size = 5) +   # larger DEG count labels above bars
  scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(
    breaks = seq(0, 3000, by = 250),
    expand = expansion(mult = c(0, 0.08)))  # slightly more headroom so labels aren't clipped

# Larger overall canvas so the bigger text/bars have room to breathe
ggsave("bigger_combined_leaves_roots_DEG_barplot.pdf", plot = p, width = 20, height = 10)
```
## Volcano Plots
Plots the log2FC against the significance
```
library(DESeq2)
library(ggplot2)

custom_palette_volcano <- c(
  "Not significant" = "grey70",
  "padj < 0.05 & log2FoldChange > 1" = "#AAAA00",
  "padj < 0.05 & log2FoldChange < -1" = "#77AADD"
)

# Compute global max Y axis across all results
global_max_y <- max(sapply(leaves_results_list, function(res) {
  max(-log10(res$padj), na.rm = TRUE)
}))
global_max_y <- ceiling(global_max_y) + 1

dir.create("Volcano_Plots_Leaves", showWarnings = FALSE)

for (cond in names(leaves_results_list)) {
  res <- as.data.frame(leaves_results_list[[cond]])
  res$gene <- rownames(res)
  res <- res[!is.na(res$padj), ]

  res$threshold <- "Not significant"
  res$threshold[res$padj < 0.05 & res$log2FoldChange > 1] <- "padj < 0.05 & log2FoldChange > 1"
  res$threshold[res$padj < 0.05 & res$log2FoldChange < -1] <- "padj < 0.05 & log2FoldChange < -1"
  res$threshold <- factor(res$threshold,
                          levels = c("padj < 0.05 & log2FoldChange < -1", 
                                     "Not significant",
                                     "padj < 0.05 & log2FoldChange > 1"))

  volcano_file <- file.path("Volcano_Plots_Leaves", paste0("volcano_", cond, "_vs_Ctrl_leaves.png"))
  png(volcano_file, width = 7, height = 6, units = "in", res = 300)

  print(
    ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = threshold), alpha = 0.8, size = 1.6) +
      scale_color_manual(values = custom_palette_volcano) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
      labs(
        title = paste(cond, "Leaves"),
        x = "log2 Fold Change",
        y = "-log10(p-adj)"
      ) +
      ylim(0, global_max_y) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "none"   # ⬅️ remove legend
      )
  )

  dev.off()
}
```
## MA Plots
Plots the log2FC against the mean expression
```
#MA Plots for each condition

library(ggplot2)

# Define a publication-friendly color palette
pal <- c(
  "Not significant" = "grey70",
  "padj < 0.05" = "#AA3377",         # blue, colorblind-friendly
  "|log2FoldChange| > 1 & padj < 0.05" = "#4477aa"  # green, colorblind-friendly
)

for (cond in names(leaves_results_ordered_list)) {
  df <- as.data.frame(leaves_results_ordered_list[[cond]])
  df$gene <- rownames(df)
  df <- df[!is.na(df$padj) & df$baseMean > 0, ]
  
  # Define significance thresholds
  df$threshold <- "Not significant"
  df$threshold[df$padj < 0.05] <- "padj < 0.05"
  df$threshold[df$padj < 0.05 & abs(df$log2FoldChange) > 1] <- "|log2FoldChange| > 1 & padj < 0.05"
  
  pdf(paste0("MA_plot_", cond, "_vs_Ctrl_leaves.pdf"), width = 7, height = 6)
  
  p <- ggplot(df, aes(x = log10(baseMean), y = log2FoldChange)) +
    geom_point(aes(color = threshold), size = 1.5, alpha = 0.6) +
    scale_color_manual(values = pal) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(size = 12),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste("MA Plot:", cond, "vs Ctrl - Leaves"),
      x = expression(log[10](baseMean)),
      y = expression(log[2](Fold~Change)),
      color = "Significance"
    )
  
  print(p)
  dev.off()
}
```
## UpSet Plot
Like a Venn diagram but for >5 categories
```
# UpSet plot
library(ggplot2)
library(dplyr)
library(tidyr)
library(UpSetR)

control <- "Ctrl"

# Get all CSV files
leaves_all_deg_files <- list.files(
  path = ".",
  pattern = "^sig_DEGs_.*_vs_Ctrl_leaves\\.csv$",
  full.names = TRUE
)

# Read DEG lists
sig_degs_leaves <- list()

for (file in leaves_all_deg_files) {
  cond <- sub("^sig_DEGs_(.*)_vs_Ctrl_leaves\\.csv$", "\\1", basename(file))
  df <- read.csv(file, row.names = 1)
  sig_degs_leaves[[cond]] <- rownames(df)
}

# Long format
deg_long <- tibble(
  gene = unlist(sig_degs_leaves),
  treatment = rep(names(sig_degs_leaves), lengths(sig_degs_leaves))
)

# Wide binary matrix
deg_wide <- deg_long %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = treatment,
    values_from = present,
    values_fill = list(present = 0),
    values_fn = max   # 🔑 THIS FIXES IT
  )

input_matrix <- as.data.frame(deg_wide)
rownames(input_matrix) <- input_matrix$gene
input_matrix$gene <- NULL

# Compute the intersection pattern for each row
row_pattern <- apply(input_matrix, 1, paste0, collapse = "")
intersection_counts <- table(row_pattern)

# Keep only rows that belong to intersections ≥10 
valid_patterns <- names(intersection_counts[intersection_counts >= 10])

# Keep patterns where only ONE treatment is present
single_patterns <- names(intersection_counts[
  rowSums(do.call(rbind, strsplit(names(intersection_counts), "")) == "1") == 1
])

valid_patterns <- union(valid_patterns, single_patterns)

input_matrix <- input_matrix[row_pattern %in% valid_patterns, ]

# UpSet plot
#pdf(UpSet_DEGs_All_Treatments_Leaves.pdf", width = 16, height = 12)
png("UpSet_DEGs_All_Treatments_Leaves.png", width = 3000, height = 2000, res = 300)

print(
  UpSetR::upset(
    input_matrix,
    nsets = ncol(input_matrix),
    nintersects = NA,
    mb.ratio = c(0.6, 0.4),
    order.by = "freq",    # sort by intersection size
    decreasing = TRUE      # largest overlaps first
  )
)
dev.off()
```
And how I was able to extract DEGs that intersected between specific sets in UpSet plot:
```
#to pull out all DEGs that intersect between N, K, Mg, and Zn leaf treatments from UpSet plot

required_sets <- c("N", "K", "Mg", "Zn")
other_sets <- setdiff(colnames(input_matrix), required_sets)

NKMgZn_treatments_sigDEGs_leaves <- rownames(input_matrix)[
  rowSums(input_matrix[, required_sets, drop = FALSE]) == length(required_sets) &
  rowSums(input_matrix[, other_sets, drop = FALSE]) == 0
]

write.table(
  data.frame(gene = NKMgZn_treatments_sigDEGs_leaves),
  file = "sigDEGs_NKMgZn_Treatments_Overlap_Leaves.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
```




