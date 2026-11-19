# How to do functional enrichment with GSEA

I used two different methods for my nutrient deficiency project, one of which was GSEA. GSEA uses a ranked list of unfiltered DESeq2 results to do random walk! I used clusterProfiler to do GSEA on GO terms and KEGG pathways.

## GSEA with clusterProfiler - GO terms
Before you run: this script is optimized to use quinoa transcript GO term information from TRAPID. When you upload transcript info (cds) to TRAPID, it does its own annotation with eggnog. These GO terms are needed to build the TERM2NAME and TERM2GENE elements. I believe I had to do this because there was not an existing database for quinoa that clusterProfiler could natively use, there could be another way around this!

This script is loaded with some nested loops. It outputs lists of enriched GO terms (also sorted by up/downregulated based on the NES score), and various plots: dotplots, ridgeplots, category netplots. It also prints out GO term coverage (which percentage of transcripts had GO terms annotated). I really only ended up using the dotplots. 
```
library(clusterProfiler)
library(tidyverse)
library(enrichplot)

## -----------------------------
## A. Load TRAPID transcript→GO
## -----------------------------


trapid_anno <- read.delim(
  "transcripts_go_exp7891.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  fill = TRUE
)


# Load GFF annotation
#gff <- import("/Users/ashmarch/Nutrient/qq74_quinoa_braker_functional_longest_isoform.gff3") %>% as.data.frame()

# Extract mapping: gene_id → Uniprot ID
mapping_uniprot <- gff %>%
  filter(type %in% c("mRNA", "transcript")) %>%
  mutate(gene_id = as.character(Parent)) %>%
  select(gene_id, uniprot_id) %>%
  filter(!is.na(uniprot_id)) %>%
  mutate(gene_id = str_trim(gene_id)) %>%
  distinct()

# Keep only informative terms (optional but recommended):
# TRAPID's is_hidden==1 tends to mark very broad GO like 'cell', 'organelle'.
trapid_anno <- trapid_anno %>%
  mutate(go = trimws(go),
         transcript_id = trimws(transcript_id)) %>%
  filter(!is.na(go), !is.na(transcript_id)) %>%
  filter(is_hidden == 0)         # <-- remove overly generic terms; set to TRUE to include all

## -----------------------------
## C. Join to get gene-level GO
## -----------------------------
trapid_gene_go <- trapid_anno %>%
  inner_join(mapping, by = "transcript_id") %>%  # only transcripts we can map to genes
  transmute(go_id = go, gene_id = gene_id, description = description) %>%
  distinct()

## -----------------------------
## D. Build TERM2GENE & TERM2NAME
## -----------------------------
cp_TERM2GENE <- trapid_gene_go %>%
  select(go_id, gene_id) %>%
  distinct()

cp_TERM2NAME <- trapid_gene_go %>%
  select(go_id, description) %>%
  mutate(description = if_else(is.na(description) | description == "",
                               go_id,  # fallback to GO ID text
                               description)) %>%
  distinct(go_id, .keep_all = TRUE)

# -----------------------------
# 2. Define directory and patterns
# -----------------------------
cp_data_dir <- "~/Nutrient/"
cp_treatments <- c("Ca", "Fe", "K", "Mg", "Mn", "N", "P", "S", "Zn", "B", "Cu", "Mo")
cp_tissues <- c("leaves", "roots")

# Create output folder for plots and results
cp_plot_dir <- "clusterProfiler"
if (!dir.exists(cp_plot_dir)) {
  dir.create(cp_plot_dir)
}

# -----------------------------
# 3. Loop through treatments and tissues
# -----------------------------
for (cp_treatment in cp_treatments) {
  for (cp_tissue in cp_tissues) {
    
    # Build file path
    cp_pattern <- paste0("DESeq2_", cp_treatment, "_vs_Ctrl_", cp_tissue, ".csv")
    cp_file_path <- file.path(cp_data_dir, cp_pattern)
    
    # Check if file exists
    if (!file.exists(cp_file_path)) {
      message("File not found: ", cp_file_path)
      next
    }
    
    # -----------------------------
    # 4. Read DESeq output
    # -----------------------------
    cp_degs <- read_csv(cp_file_path)
    
    # Expect columns: gene_id, log2FoldChange
    if (!all(c("gene_id", "log2FoldChange") %in% colnames(cp_degs))) {
      stop("File ", cp_file_path, " does not have expected columns.")
    }
    
    # Create ranked gene list (full list)
    cp_gene_list <- cp_degs %>%
      select(gene_id, log2FoldChange) %>%
      drop_na() %>%
      arrange(desc(log2FoldChange)) %>%
      deframe()  # named numeric vector
    
      message("Processing: ", cp_treatment, "-", cp_tissue)
      message("DESeq rows: ", nrow(cp_degs))
      message("GO coverage: ", mean(cp_degs$gene_id %in% cp_TERM2GENE$gene_id) * 100, "%")

    # -----------------------------
    # 5. Run GSEA
    # -----------------------------
    cp_gsea_res <- GSEA(
      geneList = cp_gene_list,
      TERM2GENE = cp_TERM2GENE,
      TERM2NAME = cp_TERM2NAME,
      minGSSize = 1,
      maxGSSize = 50000,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH"
    )
    
    cp_gsea_df <- as.data.frame(cp_gsea_res)
    
    # -----------------------------
    # 6. Save full results
    # -----------------------------
    cp_out_file <- file.path(cp_plot_dir, paste0("GSEA_", cp_treatment, "_", cp_tissue, "_all.tsv"))
    write.table(cp_gsea_df, file = cp_out_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    # -----------------------------
    # 7. NES filtering: up vs down
    # -----------------------------
    cp_up_terms <- cp_gsea_df %>% filter(NES > 0)
    cp_down_terms <- cp_gsea_df %>% filter(NES < 0)
    
    write.table(cp_up_terms, file = file.path(cp_plot_dir, paste0("GSEA_", cp_treatment, "_", cp_tissue, "_upregulated.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(cp_down_terms, file = file.path(cp_plot_dir, paste0("GSEA_", cp_treatment, "_", cp_tissue, "_downregulated.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # -----------------------------
    # 8. Save plots (only if results exist)
    # -----------------------------
    title_text <- paste("Enriched GO Terms —", cp_treatment, cp_tissue)
    if (nrow(cp_gsea_df) > 0) {
      png(
        file.path(
          cp_plot_dir,
          paste0("GSEA_", cp_treatment, "_", cp_tissue, "_dotplot.png")
        ),
        width = 1600,
        height = 1800,
        res = 150
      )
      print(
        dotplot(
          cp_gsea_res,
          showCategory = 10,
          title = title_text,
          split = ".sign"
        ) + facet_grid(. ~ .sign)
      )
      dev.off()
    } else {
      print(paste("No enriched terms for", cp_treatment, "-", cp_tissue))

    }
    # Enrichment map plots - with error handling
    if (nrow(cp_gsea_df) > 0) {
      # emapplot needs at least 2 terms to calculate similarity
      if (nrow(cp_gsea_df) >= 2) {
        tryCatch({
          # pairwise_termsim calculates GO term similarity first
          cp_gsea_res_sim <- pairwise_termsim(cp_gsea_res)
          
          png(
            file.path(
              cp_plot_dir,
              paste0(
                "GSEA_",
                cp_treatment,
                "_",
                cp_tissue,
                "_enrichmentMap.png"
              )
            ),
            width = 1600,
            height = 1800,
            res = 150
          )
          print(
            emapplot(cp_gsea_res_sim, showCategory = 25) +
              ggtitle(paste("GSEA Enrichment Map —", cp_treatment, cp_tissue))
            )
          dev.off()
        }, error = function(e) {
          print(paste("Could not create enrichment map for", cp_treatment, cp_tissue, ":", e$message))

        })
      } else {
        print(paste("Not enough terms (< 2) for enrichment map:", cp_treatment, cp_tissue))

      }
    }
    # 
    # # -----------------------------
    # # 9. Ridgeplot (all terms)
    # # -----------------------------
    # if (nrow(cp_gsea_df) >= 2) {
    #   png(
    #     file.path(
    #       cp_plot_dir,
    #       paste0("GSEA_", cp_treatment, "_", cp_tissue, "_ridgeplot.png")
    #     ),
    #     width = 1600,
    #     height = 1800,
    #     res = 150
    #   )
    #   
    #   print(
    #     ridgeplot(cp_gsea_res, showCategory=15) +
    #       labs(x = "Enrichment distribution") +
    #       ggtitle(paste("GSEA Ridgeplot —", cp_treatment, cp_tissue))
    #   )
    #   
    #   dev.off()
    # } else {
    #   message("Not enough terms for ridgeplot: ", cp_treatment, "-", cp_tissue)
    # }
    # 
    # # -----------------------------
    # # 10. Category network plot (cnetplot)
    # # -----------------------------
    # if (nrow(cp_gsea_df) >= 2) {
    #   png(
    #     file.path(
    #       cp_plot_dir,
    #       paste0("GSEA_", cp_treatment, "_", cp_tissue, "_cnetplot.png")
    #     ),
    #     width = 1800,
    #     height = 1800,
    #     res = 150
    #   )
    #   
    #   print(
    #     cnetplot(
    #       cp_gsea_res,
    #       foldChange = cp_gene_list,
    #       showCategory = 5
    #     ) +
    #       ggtitle(paste("GO Term–Gene Network —", cp_treatment, cp_tissue))
    #   )
    #   
    #   dev.off()
    # } else {
    #   message("Not enough terms for cnetplot: ", cp_treatment, "-", cp_tissue)
    # }
    
    print(paste("Completed GSEA for", cp_treatment, "-", cp_tissue))

  }
}

```

## Enriched GO term clustering with rrvgo
In order to better group/reduce GO terms, I used rrvgo to cluster the enriched terms following GSEA. This will produce scatterplots and treemaps - I only really used the treemaps. Make sure you change the ontology accordingly (MF, BP, CC). One thing to keep in mind is that you have to pick an organism from orgdb, I couldn't get the quinoa to work, so I had to use Arabidopsis, which limits the GO term assignment!
```
##rrvgo code to condense and cluster enriched GO terms - GSEA GO terms 

##GSEA GO ENRICHMENT
library(rrvgo)
#library(tools)
library(org.At.tair.db)


GSEA_dir <- path.expand("~/Nutrient/clusterProfiler")
GSEA_GO_files <- list.files(
  path = GSEA_dir,
  pattern = "GSEA_.*_all\\.tsv",
  full.names = TRUE,
)
GSEA_GO_files <- GSEA_GO_files[!grepl("_KEGG_|_KO_", GSEA_GO_files)]

for (f in GSEA_GO_files) {

  message("Processing: ", f)
  GSEA_GO <- read.delim(f, comment.char = "")
  file_base <- tools::file_path_sans_ext(basename(f))
  
    scores <- setNames(abs(GSEA_GO$NES), GSEA_GO$ID)

    simMatrix <- calculateSimMatrix(
  names(scores),
  orgdb = 'org.At.tair.db',  # no database
  ont = "MF",
  method = "Rel"
)

if (is.null(simMatrix) || !is.matrix(simMatrix) && !is.data.frame(simMatrix) || nrow(simMatrix) < 2) {
  message("  Skipping (simMatrix is NULL, not a matrix, or too small)")
  next
}


    # Ensure scores match simMatrix rownames
    scores <- scores[rownames(simMatrix)]

    reducedTerms <- reduceSimMatrix(
      simMatrix,
      scores,
      orgdb = 'org.At.tair.db',
      threshold = 0.7
    )

    ## ---- Scatter plot ----
    # if (nrow(simMatrix) < 3) {
    #   message("  Skipping scatter plot (only ", nrow(simMatrix), " GO terms)")
    #   next
    # }
    # 
    # png(
    #   filename = file.path(GSEA_dir, paste0(file_base, "_BP_scatter.png")),
    #   width = 1800, height = 1600, res = 300
    # )
    # print(scatterPlot(simMatrix, reducedTerms))
    # dev.off()

    # ---- Treemap plot ----
    png(
      filename = file.path(GSEA_dir, paste0(file_base, "_MF_treemap.png")),
      width = 2000, height = 1600, res = 300
    )
    treemapPlot(reducedTerms)
    dev.off()
}

```

## KEGG enrichment with clusterProfiler
This one is also a bit complex. KEGG has a quinoa database ('cqi'), however, it was build using the QQ74 v1 annotation, and I could not get it to work except for doing it this way. This is how I wrote about it in my thesis: "Quinoa protein sequences were aligned to the RefSeq quinoa genome annotation (GCF_001683475.1) using BLASTP, and the best hit per transcript was retained based on bitscore. RefSeq protein accessions were mapped to NCBI GeneIDs using GFF annotations and subsequently converted to quinoa KEGG gene identifiers (organism code: cqi) to enable species-specific pathway analyses. Custom gene IDs were joined to this mapping, retaining the highest-scoring hit per gene where multiple transcripts mapped to the same NCBI GeneID."

The output of this script is a tsv of enriched KEGG pathways, as well as dotplots. There is an addendum that does a Pathview rendering of each enriched KEGG pathway. Would recommend keeping/looking into, as you can see the log2FC of the genes involved in the enriched pathways! Does take some time, and the really huge pathways do not end up rendering.
```
library(tidyverse)
library(clusterProfiler)
library(KEGGREST)
library(pathview)

refseq_gff_path <- "GCF_001683475.1_ASM168347v1_genomic.gff"       # the KEGG/RefSeq quinoa GFF you showed
blast_tsv_path  <- "ourQQ74_v_KEGGcqi.blast"       # your BLASTP output (CquiM… → XP_…)


library(tidyverse)
library(stringr)

# Read the RefSeq/KEGG quinoa GFF (9 columns)
refseq_gff <- readr::read_tsv(refseq_gff_path, comment = "#",
                              col_names = FALSE, show_col_types = FALSE)
colnames(refseq_gff) <- c("seqid","source","type","start","end","score","strand","phase","attributes")

# mRNA rows: capture mRNA ID and NCBI GeneID (number)
mrna_df <- refseq_gff %>%
  filter(type == "mRNA") %>%
  mutate(
    mrna_id = str_extract(attributes, "(?<=ID=)[^;]+"),          # e.g., rna-XM_021904519.1
    gene_id_num = str_extract(attributes, "(?<=GeneID:)\\d+")    # e.g., 110725042
  ) %>%
  select(mrna_id, gene_id_num) %>%
  distinct()

# CDS rows: capture protein accession (XP_...) and parent mRNA ID
cds_df <- refseq_gff %>%
  filter(type == "CDS") %>%
  mutate(
    mrna_id = str_extract(attributes, "(?<=Parent=)[^;]+"),                                # rna-XM_...
    protein_id = coalesce(
      str_extract(attributes, "(?<=protein_id=)XP_\\d+\\.\\d+"),                           # protein_id=XP_...
      str_extract(attributes, "(?<=GenBank:)XP_\\d+\\.\\d+")                               # fallback if in Dbxref
    )
  ) %>%
  filter(!is.na(mrna_id), !is.na(protein_id)) %>%
  select(mrna_id, protein_id) %>%
  distinct()

# XP -> GeneID -> cqi
xp2cqi <- cds_df %>%
  inner_join(mrna_df, by = "mrna_id") %>%
  transmute(
    protein_id,                          # XP_021734789.1
    ncbi_geneid = gene_id_num,           # "110725042"
    cqi_id      = paste0("cqi:", gene_id_num)
  ) %>%
  distinct()

message("XP→cqi mappings: ", nrow(xp2cqi),
        " | unique XP: ", n_distinct(xp2cqi$protein_id),
        " | unique GeneIDs: ", n_distinct(xp2cqi$ncbi_geneid))


blast <- readr::read_tsv(blast_tsv_path, col_names = FALSE, show_col_types = FALSE)
colnames(blast) <- c("query","subject","pident","length","mismatch","gapopen",
                     "qstart","qend","sstart","send","evalue","bitscore")


# Best hit per transcript (by bitscore)
blast_best <- blast %>%
  group_by(query) %>%
  slice_max(order_by = bitscore, n = 1, with_ties = FALSE) %>%
  ungroup()

message("BLAST queries: ", n_distinct(blast$query),
        " | best-hits retained: ", nrow(blast_best))


# Join BLAST best hits to XP→cqi
blast_cqi <- blast_best %>%
  left_join(xp2cqi, by = c("subject" = "protein_id")) %>%
  filter(!is.na(cqi_id))     # keep only those with a cqi mapping

# Add your custom gene_id via `mapping` (transcript_id → gene_id)

blast_cqi_gene <- blast_cqi %>%
  # Make an explicit transcript_id column robustly
  mutate(
    transcript_id = stringr::str_trim(as.character(query)),
    subject       = as.character(subject),
    cqi_id        = as.character(cqi_id),
    ncbi_geneid   = as.character(ncbi_geneid)
  ) %>%
  # Clean mapping keys and join
  left_join(
    mapping %>%
      mutate(
        transcript_id = stringr::str_trim(as.character(transcript_id)),
        gene_id       = stringr::str_trim(as.character(gene_id))
      ) %>%
      select(transcript_id, gene_id) %>%
      distinct(),
    by = "transcript_id"
  ) %>%
  filter(!is.na(gene_id)) %>%
  select(gene_id, transcript_id, subject, pident, bitscore, evalue, cqi_id, ncbi_geneid) %>%
  distinct()


# Collapse to one cqi per custom gene (choose best hit across its transcripts)
gene2cqi <- blast_cqi_gene %>%
  group_by(gene_id) %>%
  slice_max(order_by = bitscore, n = 1, with_ties = FALSE) %>%   # or by lowest evalue
  ungroup() %>%
  select(gene_id, cqi_id, ncbi_geneid)

message("Custom genes with cqi mapping: ", nrow(gene2cqi))


mapping_with_cqi <- mapping %>%
  left_join(gene2cqi, by = "gene_id")

# Quick diagnostics
sum(is.na(mapping_with_cqi$cqi_id))  # how many mapping rows lack cqi
length(unique(na.omit(mapping_with_cqi$cqi_id)))  # how many distinct cqi IDs assigned

# Save for reuse
readr::write_tsv(mapping_with_cqi, "mapping_with_cqi.tsv")



# Paths
cp_data_dir <- "~/Nutrient/"
cp_plot_dir <- "clusterProfiler"
if (!dir.exists(cp_plot_dir)) dir.create(cp_plot_dir)

# Treatments/tissues
cp_treatments <- c("Ca", "Fe", "K", "Mg", "Mn", "N", "P", "S", "Zn", "B", "Cu", "Mo")
cp_tissues    <- c("leaves", "roots")

# Load gene→cqi (and ncbi_geneid) produced earlier
mapping_with_cqi <- readr::read_tsv("mapping_with_cqi.tsv", show_col_types = FALSE)

stopifnot(all(c("gene_id","cqi_id","ncbi_geneid") %in% names(mapping_with_cqi)))

# Keep a clean key; drop NAs
cqi_keys <- mapping_with_cqi %>%
  mutate(
    gene_id     = str_trim(as.character(gene_id)),
    cqi_id      = str_trim(as.character(cqi_id)),
    ncbi_geneid = str_trim(as.character(ncbi_geneid))
  ) %>%
  filter(!is.na(ncbi_geneid), ncbi_geneid != "") %>%
  distinct(gene_id, cqi_id, ncbi_geneid)


for (cp_treatment in cp_treatments) {
  for (cp_tissue in cp_tissues) {

    #safety checks 
      # Safety: close any dangling graphics devices from prior iterations
      while (!is.null(dev.list())) dev.off()

      # Log that we started this contrast
      message(sprintf("[START] %s - %s  @ %s", cp_treatment, cp_tissue, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      
      # Report sinks & suppression state
      message(sprintf("message sinks: %s", sink.number(type = "message")))

    
    cp_file_path <- file.path(cp_data_dir, paste0("DESeq2_", cp_treatment, "_vs_Ctrl_", cp_tissue, ".csv"))
    if (!file.exists(cp_file_path)) {
      print(paste0("File not found: ", cp_file_path))
      next
    }

    # Must contain gene_id and log2FoldChange
    cp_degs <- readr::read_csv(cp_file_path, show_col_types = FALSE)
    if (!all(c("gene_id", "log2FoldChange") %in% names(cp_degs))) {
      stop("File ", cp_file_path, " missing columns 'gene_id' and/or 'log2FoldChange'.")
    }

    # Rank by log2FC (your current choice); you could switch to stat/padj-weighted ranks later
    cp_gene_list <- cp_degs %>%
      select(gene_id, log2FoldChange) %>%
      drop_na() %>%
      arrange(desc(log2FoldChange)) %>%
      deframe()

    print(paste0("Processing: ", cp_treatment, "-", cp_tissue,
            " | DE rows: ", nrow(cp_degs),
            " | unique genes ranked: ", length(cp_gene_list)))

    # Map your custom gene_ids (CquiG...) → KEGG quinoa gene IDs
    degs_with_cqi <- tibble(
        gene_id = names(cp_gene_list),
        score   = as.numeric(cp_gene_list)
      ) %>%
      inner_join(cqi_keys, by = "gene_id") %>%
      # If a custom gene maps to multiple cqi IDs (shouldn't, given your dedup), keep max score
      group_by(ncbi_geneid) %>%
      summarize(score = max(score, na.rm = TRUE), .groups = "drop")

    covered <- nrow(degs_with_cqi)
    print(paste0("Genes with KEGG cqi mapping in this contrast: ", covered, " / ", length(cp_gene_list),
            " (", round(100 * covered / length(cp_gene_list), 1), "% coverage)"))

    if (covered < 10) {
      print(paste0("Coverage is low; results may be underpowered. Consider adding more mappings or checking namespace/whitespace."))
    }

    ## use numeric NCBI GeneID with keyType = 'ncbi-geneid' ---
    kegg_gene_list <- degs_with_cqi$score
    names(kegg_gene_list) <- degs_with_cqi$ncbi_geneid
    kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)

   
    marker <- file.path(cp_plot_dir, paste0("MARKER_ran_", cp_treatment, "_", cp_tissue, ".txt"))
    writeLines(sprintf("Started: %s %s at %s", cp_treatment, cp_tissue, Sys.time()), marker)

  
    
     # Run GSEA
    cqi_gsea <- tryCatch({
      gseKEGG(
        geneList      = kegg_gene_list,
        organism      = "cqi",         # quinoa’s KEGG code
        keyType       = "ncbi-geneid", # explicit numeric GeneIDs
        minGSSize     = 5,
        maxGSSize     = 50000,
        pvalueCutoff  = 0.05,
        pAdjustMethod = "BH",
        verbose       = FALSE,
        seed          = TRUE
      )
    
      }, error = function(e) {
        print(paste0("KEGG GSEA failed for ", cp_treatment, "-", cp_tissue, ": ", e$message))
        NULL
      })
      
      if (is.null(cqi_gsea)) {
        writeLines(sprintf("GSEA_NULL: %s %s at %s", cp_treatment, cp_tissue, Sys.time()), marker, sep="\n", useBytes=TRUE)
      } else {
        writeLines(sprintf("GSEA_ROWS=%d", nrow(as.data.frame(cqi_gsea))), marker, sep="\n", useBytes=TRUE)
      }


    # Save results if non-empty
    KEGG_title_text <- paste("Enriched KEGG Pathways —", cp_treatment, cp_tissue)   
    if (!is.null(cqi_gsea) && nrow(as.data.frame(cqi_gsea)) > 0) {
      cqi_df  <- as.data.frame(cqi_gsea)

      # Add pathway category info (metabolism vs genetic info, etc.)
      cqi_df_cat <- tryCatch({
        append_kegg_category(cqi_gsea) %>% as.data.frame()
      }, error = function(e) cqi_df)

      out_tsv <- file.path(cp_plot_dir, paste0("GSEA_KEGG_cqi_", cp_treatment, "_", cp_tissue, "_all.tsv"))
      out_png <- file.path(cp_plot_dir, paste0("GSEA_KEGG_cqi_", cp_treatment, "_", cp_tissue, "_dotplot.png"))

      write.table(cqi_df_cat, file = out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

      png(out_png, width = 1600, height = 1800, res = 150)
      on.exit(dev.off(), add = TRUE)
      print(dotplot(cqi_gsea, showCategory = 15, title = KEGG_title_text , split=".sign") +
              facet_grid(.~.sign))
      dev.off()

      print(paste0("Wrote GSEA results: ", basename(out_tsv), " & ", basename(out_png)))
      
      
 # ## --------------------- NEW: Pathview rendering --------------------- ##
 #      # Make a subfolder for this contrast's pathway images
 #      pv_dir <- file.path(cp_plot_dir, paste0("pathview_", cp_treatment, "_", cp_tissue))
 #      if (!dir.exists(pv_dir)) dir.create(pv_dir)
 # 
 #      # Filter enriched pathways to render (adjust thresholds as needed)
 #      pv_hits <- cqi_df %>%
 #        dplyr::filter(p.adjust < 0.05) %>%          # FDR threshold
 #        dplyr::arrange(p.adjust)
 # 
 #      if (nrow(pv_hits) == 0) {
 #        message("No pathways pass FDR < 0.05 for Pathview in ", cp_treatment, "-", cp_tissue)
 #      } else {
 #        # Extract numeric pathway IDs (e.g., "cqi04130" -> "04130")
 #        pv_ids <- gsub("^cqi", "", pv_hits$ID)
 # 
 #        message("Rendering Pathview PNGs for ", length(pv_ids), " pathways in ", cp_treatment, "-", cp_tissue)
 # 
 #        # Loop pathways and render
 #        for (pid in pv_ids) {
 #          # out.suffix to keep files distinguishable per contrast
 #          out_suffix <- paste0(cp_treatment, "_", cp_tissue)
 #          # Native KEGG PNG
 #          try({
 #            pv.out <- pathview(
 #              gene.data   = kegg_gene_list,  # named numeric vector: names = NCBI GeneIDs
 #              pathway.id  = pid,             # numeric ID like "04130"
 #              species     = "cqi",           # quinoa KEGG code
 #              gene.idtype = "entrez",        # Entrez/NCBI GeneID
 #              kegg.native = TRUE,            # PNG native KEGG view
 #              out.suffix  = out_suffix,
 #              kegg.dir    = pv_dir,
 #              # Optional visual tuning:
 #              low   = list(gene = "blue"),   # down-regulated
 #              mid   = list(gene = "gray"),
 #              high  = list(gene = "red"),    # up-regulated
 #              limit = list(gene = 2)         # color limit ~2 (log2FC-like); adjust if using other ranks
 #            )
 #          }, silent = TRUE)
 # 
 #        }
 #      }
 #      ## ------------------- END: Pathview rendering block ------------------ ##
    } else {
      print(paste0("No significant KEGG GSEA for ", cp_treatment, "-", cp_tissue, " (cqi)."))
    }
  }
}

```
