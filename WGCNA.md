# How to run WGCNA

## Running WGCNA
IMPORTANT: this script will produce a soft power threshold plot. It is important to pick the soft power value at which R^2 > 0.8. Then you need to assign that soft power value in the script for the WGCNA to run. Parameters here might need to be tweaked; I believe it is dependent on the experiment. Also important, this WGCNA is signed, you could also run unsigned.
```
#WGCNA

#==============================
# WGCNA for nutrient RNA-seq (Leaves)
# Using same sample names as PLS-DA
# Each deficiency as a separate trait
#==============================

library(WGCNA)
library(DESeq2)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()  # enable multi-threading

#-----------------------------
# 1. Load VST data from DESeq2
#-----------------------------
vsd_leaves <- vst(nutrientAllLeaves_dds, blind = FALSE)
expr_data <- assay(vsd_leaves)  # genes x samples

#-----------------------------
# 2. Filter low-variance genes (optional)
#-----------------------------
var_genes <- apply(expr_data, 1, var)
expr_data_filtered <- expr_data[var_genes > quantile(var_genes, 0.25), ]  # top 75% variable

#-----------------------------
# 3. Prepare trait matrix
#-----------------------------
sample_info <- colData(nutrientAllLeaves_dds)
conditions <- levels(sample_info$condition)
conditions <- conditions[conditions != "Ctrl"]  # exclude control

# Initialize binary trait matrix
trait_matrix <- matrix(0, nrow = ncol(expr_data_filtered), ncol = length(conditions))
rownames(trait_matrix) <- colnames(expr_data_filtered)
colnames(trait_matrix) <- conditions

for (cond in conditions) {
  trait_matrix[which(sample_info$condition == cond), cond] <- 1
}
trait_df <- as.data.frame(trait_matrix)

#-----------------------------
# 4. Transpose for WGCNA
#-----------------------------
datExpr <- t(expr_data_filtered)  # samples x genes

# Check for good samples and genes
gsg <- goodSamplesGenes(datExpr, verbose = 3)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

#-----------------------------
# 5. Pick soft-threshold
#-----------------------------
powers <- c(1:10, seq(12,20,2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot to choose power
plot(sft$fitIndices[,1], sft$fitIndices[,2], type="n",
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2")
text(sft$fitIndices[,1], sft$fitIndices[,2], labels=powers)

#-----------------------------
# 6. Build network and detect modules
#-----------------------------
cor <- WGCNA::cor
softPower <- 5  # choose based on previous plot
net <- blockwiseModules(datExpr,
                        power = softPower,
                        TOMType = "signed",
                        minModuleSize = 30,
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        verbose = 3)

# Module eigengenes
MEs <- net$MEs

#-----------------------------
# 7. Correlate modules with nutrient traits
#-----------------------------
# Calculate correlations and p-values
module_trait_cor <- cor(MEs, trait_df, use = "p")
module_trait_p <- corPvalueStudent(module_trait_cor, nrow(datExpr))

# Combine both into a single text matrix
text_matrix <- paste0(
  signif(module_trait_cor, 2), "\n(",
  signif(module_trait_p, 1), ")"
)

dim(text_matrix) <- dim(module_trait_cor)

# Save the heatmap as PNG
png("WGCNA_Module_Trait_Heatmap_Leaves_signed.png", width = 15, height = 15, units = "in", res = 300)


labeledHeatmap(
  Matrix = module_trait_cor,
  xLabels = colnames(trait_df),
  yLabels = names(MEs),
  colorLabels = TRUE,
  colors = blueWhiteRed(50),
  textMatrix = text_matrix,
  setStdMargins = FALSE,
  cex.text = 0.7,
  zlim = c(-1, 1),
  main = "Module–Trait Relationships (Leaves)"
)

dev.off()

```

Here is a bunch of code to extract various bits of information from the WGCNA data. I usually will just run the few lines I need at a time depending on what I want:
```
#to extract genes from WGCNA plot

#-----------------------------
# Extract genes from a specific module
#-----------------------------

# First, let's see what modules you have
#table(net$colors)  # shows how many genes are in each module

# The module colors are numeric (0, 1, 2, etc.)
# Module 0 is typically the "grey" module (unassigned genes)

# Example: Extract genes from module 1
module_of_interest <- 0

# Get gene names for this module
module_genes <- colnames(datExpr)[net$colors == module_of_interest]

# View the first few gene names
#head(module_genes)

# Get all gene names as a vector
#length(module_genes)  # how many genes in this module

# Save to a file
# write.table(module_genes, 
#             file = paste0("WGCNA_Module_", module_of_interest, "_NPositive_Leaves_genes.txt"),
#             row.names = FALSE, 
#             col.names = FALSE, 
#             quote = FALSE)

# Convert to dataframe for joining
module_genes_df <- data.frame(gene_id = module_genes)

# Join to mapping to get transcript IDs
module_transcripts <- module_genes_df %>%
  left_join(mapping %>% select(gene_id, transcript_id), by = "gene_id")

#save transcripot names to a file
write.table(
  module_transcripts$transcript_id,
  file = paste0("WGCNA_Module_", module_of_interest, "_NPositive_Leaves_transcripts.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

#-----------------------------
# Extract genes from ALL modules
#-----------------------------

# Create a dataframe with gene names and their module assignments
gene_module_df <- data.frame(
  Gene = colnames(datExpr),
  Module = net$colors
)

# View the first few rows
head(gene_module_df)

# Save complete list
write.csv(gene_module_df, 
          "WGCNA_All_Genes_Module_Assignments_Leaves.csv",
          row.names = FALSE)

#-----------------------------
# Extract genes from module with HIGHEST correlation to a specific trait
#-----------------------------

# Example: Find module most correlated with Nitrogen deficiency
trait_of_interest <- "N"  # change to your trait name

# Get correlation values for this trait
trait_correlations <- module_trait_cor[, trait_of_interest]

# Find module with highest absolute correlation
best_module <- which.max(abs(trait_correlations))
best_module_number <- as.numeric(gsub("ME", "", names(best_module)))

cat("Module", best_module_number, "has highest correlation with", trait_of_interest, 
    "\nCorrelation:", trait_correlations[best_module],
    "\nP-value:", module_trait_p[best_module, trait_of_interest], "\n")

# Extract genes from this module
best_module_genes <- colnames(datExpr)[net$colors == best_module_number]
length(best_module_genes)

#-----------------------------
# Extract genes from modules with SIGNIFICANT correlations (p < 0.05)
#-----------------------------

# For a specific trait
trait_of_interest <- "Zn"

# Find significant modules
sig_modules <- which(module_trait_p[, trait_of_interest] < 0.05)

# Extract module numbers (remove "ME" prefix)
sig_module_numbers <- as.numeric(gsub("ME", "", names(sig_modules)))

# Create a list to store genes from each significant module
sig_module_genes_list <- list()

for (mod_num in sig_module_numbers) {
  module_name <- paste0("Module_", mod_num)
  sig_module_genes_list[[module_name]] <- colnames(datExpr)[net$colors == mod_num]
  
  cat("Module", mod_num, ":", length(sig_module_genes_list[[module_name]]), "genes\n")
  cat("Correlation:", module_trait_cor[paste0("ME", mod_num), trait_of_interest], "\n")
  cat("P-value:", module_trait_p[paste0("ME", mod_num), trait_of_interest], "\n\n")
}

# Save each significant module's genes to separate files
for (mod_name in names(sig_module_genes_list)) {

  module_genes <- sig_module_genes_list[[mod_name]]

  annotated_module <- data.frame(gene_id = module_genes) %>%
    left_join(mapping, by = "gene_id")

  write.table(
    annotated_module,
    file = paste0("WGCNA_", trait_of_interest, "_", mod_name, "_genes_annotated_leaves.txt"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

#-----------------------------
# Get hub genes for a module
#-----------------------------

module_of_interest <- 24
ME_name <- paste0("ME", module_of_interest)
ME_of_interest <- MEs[, ME_name]

# Genes in this module
module_genes <- colnames(datExpr)[net$colors == module_of_interest]

# Calculate module membership (kME)
gene_module_membership <- cor(
  datExpr[, module_genes],
  ME_of_interest,
  use = "p"
)

# Build hub gene table
hub_gene_df <- data.frame(gene_id = module_genes, ModuleMembership = gene_module_membership )

# Annotate
hub_gene_df <- hub_gene_df %>%
  left_join(mapping, by = "gene_id")

# Sort by absolute kME (hub strength)
hub_gene_df <- hub_gene_df[order(-abs(hub_gene_df$ModuleMembership)), ]

# Optional: add module info
hub_gene_df$module <- module_of_interest
hub_gene_df$module_eigengene <- ME_name

# View top 10 hub genes
top_hub_genes <- head(hub_gene_df, 10)
print(top_hub_genes)

# Save full ranked hub list
write.csv(
  hub_gene_df,
  paste0("WGCNA_Module_", module_of_interest, "_HubGenes_Annotated_Ranked_leaves.csv"),
  row.names = FALSE
)

```

## Preparing data for Cytoscape (to make category netplots)
Also, I made category netplots with the WGCNA data using Cytoscape. This is the code I used to make the files necessary for this:
```
#to get hub gene lists for cytoCNA category netplots for WGCNA

module_of_interest <- 24
ME_name <- paste0("ME", module_of_interest)
ME_of_interest <- MEs[, ME_name]

module_genes <- colnames(datExpr)[net$colors == module_of_interest]

# kME (module membership)
kME <- cor(datExpr[, module_genes], ME_of_interest, use = "p")

hub_df <- data.frame(
  gene_id = module_genes,
  kME = kME
) %>%
  arrange(desc(abs(kME)))

top_hub_genes <- hub_df$gene_id[1:20]

datExpr_hubs <- datExpr[, top_hub_genes]

TOM_hubs <- TOMsimilarityFromExpr(
  datExpr_hubs,
  power = softPower
)

# CRITICAL: preserve gene names
colnames(TOM_hubs) <- top_hub_genes
rownames(TOM_hubs) <- top_hub_genes

library(dplyr)
library(reshape2)

TOM_df <- melt(TOM_hubs) %>%
  rename(
    source = Var1,
    target = Var2,
    weight = value
  ) %>%
  filter(source != target, weight > 0.1)

gene_name_lookup <- mapping %>%
  select(gene_id, gene_name) %>%
  distinct() %>%
  mutate(
    gene_name = ifelse(
      is.na(gene_name) | gene_name == "" | gene_name == "character(0)",
      gene_id,        # fallback to gene_id
      gene_name
    )
  )

TOM_named <- TOM_df %>%
  left_join(gene_name_lookup, by = c("source" = "gene_id")) %>%
  rename(source_name = gene_name) %>%
  left_join(gene_name_lookup, by = c("target" = "gene_id")) %>%
  rename(target_name = gene_name) %>%
  select(source = source_name, target = target_name, weight)

write.table(
  TOM_named,
  file = paste0(
    "Module_", module_of_interest,
    "_Top20HubGenes_Leaves_signed.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
```
I looked at core WGNCA transcripts by doing functional enrichment of GO terms. This process is outlined in the functional enrichment - TRAPID page!
