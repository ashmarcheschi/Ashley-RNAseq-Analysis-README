# How to work with GFF/mapping element in R

For a lot of this work, I had to use a "mapping" element, which was information from the annotation GFF. Mostly, it was transcript to gene mapping, as DESeq2 provides data on transcripts, not genes. I have various lines of code to make and use this mapping element across various analyses.
## Extracting transcript IDs from lists of significant DEGs
```
library(tidyverse)
library(Biostrings)

# 1. Read CDS fasta headers (longest isoforms only)
cds <- readDNAStringSet("quinoa_braker_functional_longest_isoform.cds")  # <-- update with your file
headers <- names(cds)

# Build mapping: transcript_id ↔ gene_id
mapping <- tibble(
  transcript_id = sub(" .*", "", headers),
  gene_id = sub(".*gene=([^ ]+).*", "\\1", headers)
)

# 2. List all DEG result files
deg_files <- list.files(
  path = "~/Nutrient/", 
  pattern = "^sig_DEGs_.*_vs_Ctrl_.*\\.csv$", 
  full.names = TRUE
)

# 3. Loop over DEG files
for (file in deg_files) {
  # Read DEGs
  degs <- read_csv(file)
  gene_ids <- degs[[1]]
  
  # Map gene → transcript (1:1 since longest isoform only)
  transcripts <- mapping %>%
    filter(gene_id %in% gene_ids) %>%
    pull(transcript_id)
  
  # Save plain transcript list
  out_file <- gsub("\\.csv$", "_transcriptsID.txt", file)
  writeLines(transcripts, out_file)
}
```
## Extracting gene IDs from lists of significant DEGs
```
library(tidyverse)

# --- 1️⃣ Define your input folder ---
deg_path <- "~/Nutrient/"   # <-- change if needed

# --- 2️⃣ Find all matching DEG result files ---
deg_files <- list.files(
  path = deg_path,
  pattern = "^sig_DEGs_.*_(up|down)regulated\\.tsv$",
  full.names = TRUE
)

# --- 3️⃣ Loop through each file ---
for (file in deg_files) {
  
  # Read TSV (tab-separated)
  degs <- read_tsv(file, show_col_types = FALSE)
  
  # Make sure 'gene_id' exists
  if (!"gene_id" %in% colnames(degs)) {
    warning("Skipping file (no gene_id column): ", basename(file))
    next
  }
  
  # Extract only the gene_id column, drop NAs, unique
  gene_ids <- degs %>%
    filter(!is.na(gene_id)) %>%
    distinct(gene_id) %>%
    pull(gene_id)
  
  # --- 4️⃣ Build output filename ---
  base_name <- gsub("\\.tsv$", "", basename(file))
  out_file <- file.path(deg_path, paste0(base_name, "_geneIDs.txt"))
  
  # --- 5️⃣ Write to file ---
  writeLines(gene_ids, out_file)
  
  message("✅ Extracted ", length(gene_ids), " gene IDs from: ", basename(file))
}

message("🎉 Done! Gene ID lists created for all DEG files.")
```
## Extracting Uniprot IDs from lists of significant DEGs
```
library(tidyverse)
library(rtracklayer)

# Load GFF annotation
gff <- import("/Users/ashmarch/Nutrient/qq74_quinoa_braker_functional_longest_isoform.gff3") %>% as.data.frame()

mapping <- gff %>%
  filter(type %in% c("mRNA","transcript")) %>%
  mutate(gene_id = as.character(Parent)) %>%
  select(gene_id, uniprot_id) %>%
  filter(!is.na(uniprot_id)) %>%
  mutate(gene_id = str_trim(gene_id))

# DEG files
deg_path <- "/Users/ashmarch/Nutrient/"
deg_files <- list.files(
  path = deg_path,
  pattern = "^sig_DEGs_.*_(up|down)regulated\\.tsv$",
  full.names = TRUE
)

for (file in deg_files) {
  
  degs <- read_tsv(file, show_col_types = FALSE) %>%
    mutate(gene_id = str_trim(as.character(gene_id)))
  
  uniprot_ids <- mapping %>%
    filter(gene_id %in% degs$gene_id) %>%
    pull(uniprot_id) %>%
    unique()
  
  base_name <- gsub("\\.tsv$", "", basename(file))
  out_file <- file.path(deg_path, paste0(base_name, "_uniprotIDs.txt"))
  writeLines(uniprot_ids, out_file)
  
  message("✅ Processed ", basename(file), " → ", length(uniprot_ids), " UniProt IDs")
}

```
