# How to run PLS-DA

## PLS-DA on all treatments combined against the control (stress vs. no stress)
This main script will provide you with plots of each condition replicates agains the control replicates, a table of 50 genes with the highest VIP scores, a variable plot of this top 50 highest scoring genes, and a plot showing the separation of all samples. Is also VST transformed, like the PCA plots
```
#==============================
# PLS-DA (Leaves)
#==============================

library(DESeq2)
library(mixOmics)
library(tidyverse)

# VST transform
vsd_leaves <- vst(nutrientAllLeaves_dds, blind = FALSE)

# Extract the matrix
expr_mat_leaves <- assay(vsd_leaves)

# Create stress group label
leaf_groups <- ifelse(colData(nutrientAllLeaves_dds)$condition == "Ctrl",
                      "Control", "Stress")
leaf_groups <- factor(leaf_groups)

# Run PLS-DA
plsda_res_leaves <- plsda(t(expr_mat), leaf_groups, ncomp = 2)

#---------------------------------
# Plot 1: Individuals (Stress vs Control)
#---------------------------------
# png("PLSDA_Leaves_Plots/PLSDA_Stress_vs_Control_allTreatments.png", width = 7, height = 6, units = "in", res = 300)
# plotIndiv(plsda_res_leaves,
#           group = leaf_groups,
#           ind.names = FALSE,
#           ellipse = TRUE,
#           legend = TRUE,
#           title = "PLS-DA: Leaves (Stress vs Control)")
# dev.off()

#---------------------------------
# Cross-validation
#---------------------------------
set.seed(123)
cv <- perf(plsda_res_leaves, validation = "Mfold", folds = 4, nrepeat = 50, progressBar = TRUE)
cv$error.rate
cv$Q2.total  # model predictive strength

#---------------------------------
# VIP Scores
#---------------------------------
vip_scores_leaves <- vip(plsda_res_leaves)
top_vip_leaves <- sort(vip_scores_leaves[,1], decreasing = TRUE)[1:50]
top_vip_leaves

write.csv(data.frame(Gene = names(top_vip_leaves),
                     VIP = top_vip_leaves),
          "PLSDA_Leaves_Plots/PLSDA_top50_genes_leaves.csv",
          row.names = FALSE)

#---------------------------------
# Plot 2: Variable Plot (Top 50 VIP genes only)
#---------------------------------
top_genes_leaves <- names(top_vip_leaves)
plsda_sub_leaves <- plsda(t(expr_mat_leaves[top_genes_leaves, ]), leaf_groups, ncomp = 2)

png("PLSDA_Leaves_Plots/PLSDA_Top50_Variables_leaves.png", width = 7, height = 6, units = "in", res = 300)
plotVar(plsda_sub_leaves, comp = c(1,2),
        title = "Variable Plot: Top 50 VIP Genes",
        var.names = FALSE)
dev.off()

#---------------------------------
# Plot 3: Individuals Colored by Nutrient Deficiency
#---------------------------------
nutrient_labels <- colData(nutrientAllLeaves_dds)$condition
png("PLSDA_Leaves_Plots/PLSDA_By_Nutrient_leaves.png", width = 7, height = 6, units = "in", res = 300)
plotIndiv(plsda_res_leaves, group = nutrient_labels, legend = TRUE,
          title = "PLS-DA: Leaves by Nutrient Deficiency")
dev.off()

```

I call this next script an "extension" of the PLS-DA, which takes all the information and annotates it based on the mapping element/GFF. It provides transcript IDs and gene names. It also incorporates one line that I added in that will extract the putative gene name from an annotation if it doesn't have a name already (lots of times I found that the name was stored in the description, so this helps more genes get labeled instead of having to use their gene ID).
```
library(tidyverse)
#extension of PLSDA

extract_putative_name <- function(product) {
  if (is.na(product)) return(NA_character_)

  if (grepl("hypothetical protein", product, ignore.case = TRUE)) {
    return(NA_character_)
  }

  tokens <- unlist(strsplit(product, "[ ,;()]+"))
  if (length(tokens) == 0) return(NA_character_)

  last_token <- tail(tokens, 1)

  bad_tokens <- c(
    "vacuolar", "mitochondrial", "chloroplastic", "chloroplast",
    "cytosolic", "peroxisomal", "nuclear", "extracellular",
    "kinase", "decarboxylase", "homolog", "fragment", "oxidase", "reductase", "protein", "Kirola"
  )

  if (tolower(last_token) %in% bad_tokens) {
    return(NA_character_)
  }

  if (grepl("^[A-Za-z][A-Za-z0-9._-]{1,15}$", last_token)) {
    return(last_token)
  }

  NA_character_
}

gene_name_lookup <- mapping %>%
  dplyr::select(gene_id, gene_name, product) %>%
  distinct() %>%
  rowwise() %>%
  mutate(
    putative_name = extract_putative_name(product),
    gene_name_final = case_when(
      !is.na(gene_name) & gene_name != "" & gene_name != "character(0)" ~ gene_name,
      !is.na(putative_name) ~ putative_name,
      TRUE ~ gene_id
    )
  ) %>%
  ungroup() %>% 
  dplyr::select(gene_id, gene_name = gene_name_final)



# Convert top_vip named vector to data frame
top_vip_df_leaves <- data.frame(
  gene_id = names(top_vip_leaves),
  VIP = as.numeric(top_vip_leaves)
)

# Join with annotation table
top_vip_annotated_leaves <- top_vip_df_leaves %>%
  left_join(gene_name_lookup, by = "gene_id") %>%
  arrange(desc(VIP)) %>%
  slice_head(n = 50)

# Save
write_csv(top_vip_annotated_leaves, "PLSDA_Leaves_Plots/PLSDA_top50_genes_allTreatments_annotated_leaves.csv")

#extract transcript IDs from top PLSDA VIPs for leaves and roots

library(tidyverse)
library(Biostrings)

# 2. Define PLSDA file paths
plsda_paths <- c(
  "/Users/ashmarch/Nutrient/PLSDA_Leaves_Plots/PLSDA_top50_genes_leaves_allTreatments_annotated.csv",
  "/Users/ashmarch/Nutrient/PLSDA_Roots_Plots/PLSDA_top50_genes_roots_all_Treatments_annotated.csv"
)

# 3. Loop over PLSDA files
for (file in plsda_paths) {
  # Read PLSDA top genes
  plsda_genes <- read_csv(file)
  
  # Assuming gene IDs are in the first column (adjust if needed)
  gene_ids <- plsda_genes[[1]]
  gene_names <- plsda_genes[[2]]
  
  # Map gene → transcript
  transcripts <- mapping %>%
    filter(gene_id %in% gene_ids) %>%
    pull(transcript_id)
  
  # Save transcript list
  out_file <- gsub("\\.csv$", "_transcriptIDs.txt", file)
  writeLines(transcripts, out_file)
  
  #Save gene name list
  out_file <- gsub("\\.csv$", "_geneNames.txt", file)
  writeLines(gene_names, out_file)
}

```
## PLS-DA on individual treatments vs. control 

Fundamentally the same as above, but this loops over all treatments and runs them against the control instead of lumping all stress treatments together
```
# Treatment specific PLSDA (Leaves)

all_conditions_leaves <- unique(colData(nutrientAllLeaves_dds)$condition)
treatments <- setdiff(all_conditions_leaves, "Ctrl")

library(DESeq2)
library(mixOmics)
library(tidyverse)

for (trt in treatments) {

  message("Running PLS-DA for treatment: ", trt)

  ## ----------------------------
  ## Subset DDS: Ctrl + treatment
  ## ----------------------------
  dds_sub_leaves <- nutrientAllLeaves_dds[
    ,
    colData(nutrientAllLeaves_dds)$condition %in% c("Ctrl", trt)
  ]

  dds_sub_leaves$condition <- droplevels(dds_sub_leaves$condition)
  message("Sanity check for ", trt)

  dds_sub_leaves$group <- factor(
    ifelse(
      colData(dds_sub_leaves)$condition == "Ctrl",
      "Control",
      trt
    )
  )

  groups <- dds_sub_leaves$group
  print(colnames(dds_sub_leaves))
  print(table(colData(dds_sub_leaves)$condition))

  ## ----------------------------
  ## VST transform
  ## ----------------------------
  vsd_sub_leaves <- vst(dds_sub_leaves, blind = FALSE)
  expr_mat_sub_leaves <- assay(vsd_sub_leaves)

  ## ----------------------------
  ## Run PLS-DA
  ## ----------------------------
  plsda_res_sub_leaves <- plsda(t(expr_mat_sub_leaves), groups, ncomp = 2)

  ## ----------------------------
  ## Plot 1: Individuals
  ## ----------------------------
  png(
    file.path("PLSDA_Leaves_Plots",
              paste0("PLSDA_", trt, "_vs_Control.png")),
    width = 7, height = 6, units = "in", res = 300
  )

  plotIndiv(
    plsda_res_sub_leaves,
    group = groups,
    ind.names = FALSE,
    ellipse = TRUE,
    legend = TRUE,
    title = "PLS-DA - Leaves"
  )

  dev.off()

  ## ----------------------------
  ## Cross-validation
  ## ----------------------------
  #set.seed(123)
  # cv <- perf(
  #   plsda_res_sub_leaves,
  #   validation = "Mfold",
  #   folds = 4,
  #   nrepeat = 50,
  #   progressBar = FALSE
  # )

  ## ----------------------------
  ## VIP Scores
  ## ----------------------------
  # vip_scores_sub_leaves <- vip(plsda_res_sub_leaves)
  # top_vip_sub_leaves <- sort(vip_scores_sub_leaves[, 1], decreasing = TRUE)[1:50]
  # 
  # write.csv(
  #   data.frame(
  #     gene_id = names(top_vip_sub_leaves),
  #     VIP = top_vip_sub_leaves
  #   ),
  #   file.path("PLSDA_Leaves_Plots",
  #             paste0("PLSDA_top50_VIP_", trt, "_leaves.csv")),
  #   row.names = FALSE
  # )

  # #to extract a txt file with transcript names of genes with VIP scores >= 1 for enrichment
  # vip_scores_sub_leaves <- vip(plsda_res_sub_leaves)
  # # Extract component 1 VIP scores
  # vip_df <- data.frame(
  #   gene_id = rownames(vip_scores_sub_leaves),
  #   VIP = vip_scores_sub_leaves[, 1]
  # )
  # 
  # # Keep genes with VIP >= 1.5
  # vip_genes <- vip_df %>%
  #   filter(VIP >= 1.5)
  # 
  # ## ----------------------------
  # ## Map genes -> transcripts
  # ## ----------------------------
  # vip_transcripts <- vip_genes %>%
  #   left_join(mapping, by = "gene_id") %>%
  #   dplyr::select(transcript_id) %>%
  #   distinct() %>%
  #   drop_na()
  # 
  # ## ----------------------------
  # ## Write transcript list
  # ## ----------------------------
  # writeLines(
  #   vip_transcripts$transcript_id,
  #   file.path(
  #     "PLSDA_Leaves_Plots",
  #     paste0("PLSDA_VIP_transcripts_", trt, "_leaves.txt")
  #   )
  # )
  
  # ## ----------------------------
  # ## Plot 2: Variable plot (Top 50)
  # ## ----------------------------
  # plsda_sub_leaves <- plsda(
  #   t(expr_mat_sub_leaves[names(top_vip_sub_leaves), ]),
  #   groups,
  #   ncomp = 2
  # )
  # 
  # png(
  #   file.path("PLSDA_Leaves_Plots",
  #             paste0("PLSDA_", trt, "_Top50_Variables.png")),
  #   width = 7, height = 6, units = "in", res = 300
  # )
  # 
  # plotVar(
  #   plsda_sub_leaves,
  #   comp = c(1, 2),
  #   title = paste("Top 50 VIP Genes:", trt),
  #   var.names = FALSE
  # )

  #dev.off()
}

```

And then to annotate the produced PLS-DA data:
```
# treatment specific PLSDA - annotating VIP lists (Leaves)

library(tidyverse)

extract_putative_name <- function(product) {
  if (is.na(product)) return(NA_character_)

  if (grepl("hypothetical protein", product, ignore.case = TRUE)) {
    return(NA_character_)
  }

  tokens <- unlist(strsplit(product, "[ ,;()]+"))
  if (length(tokens) == 0) return(NA_character_)

  last_token <- tail(tokens, 1)

  bad_tokens <- c(
    "vacuolar", "mitochondrial", "chloroplastic", "chloroplast",
    "cytosolic", "peroxisomal", "nuclear", "extracellular",
    "kinase", "decarboxylase", "homolog", "fragment",
    "oxidase", "reductase", "protein", "Kirola", "transposase"
  )

  if (tolower(last_token) %in% bad_tokens) {
    return(NA_character_)
  }

  if (grepl("^[A-Za-z][A-Za-z0-9._-]{1,15}$", last_token)) {
    return(last_token)
  }

  NA_character_
}

gene_name_lookup <- mapping %>%
  dplyr::select(transcript_id, gene_id, gene_name, product) %>%
  distinct() %>%
  rowwise() %>%
  mutate(
    putative_name = extract_putative_name(product),
    gene_name_final = case_when(
      !is.na(gene_name) & gene_name != "" & gene_name != "character(0)" ~ gene_name,
      !is.na(putative_name) ~ putative_name,
      TRUE ~ gene_id
    )
  ) %>%
  ungroup() %>%
  dplyr::select(gene_id, transcript_id, gene_name = gene_name_final)


vip_files_leaves <- list.files(
  "PLSDA_Leaves_Plots",
  pattern = "^PLSDA_top50_VIP_.*_leaves\\.csv$",
  full.names = TRUE
)

# loop
for (vip_file in vip_files_leaves) {

  message("Annotating: ", basename(vip_file))

  vip_df <- read_csv(vip_file, show_col_types = FALSE)

  vip_annotated <- vip_df %>%
    left_join(gene_name_lookup, by = "gene_id") %>%
    arrange(desc(VIP))

  out_file <- sub(
    "_leaves\\.csv$",
    "_annotated_leaves.csv",
    vip_file
  )

  write_csv(vip_annotated, out_file)
  
  vip_transcript_list <- vip_annotated %>%
  filter(!is.na(transcript_id)) %>%
  distinct(transcript_id)
  
  transcript_outfile <- sub(
  "_leaves\\.csv$",
  "_VIP_transcripts.txt",
  vip_file
)

write.table(
  vip_transcript_list,
  file = transcript_outfile,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
}

```
I looked at transcripts with the highest VIP scores by running functional enrichment. To do this, I extracted the transcript IDs and uploaded them to TRAPID. This process is outlined in the functional enrichment - TRAPID page!
