# How to do functional enrichment with TRAPID (GO terms and/or protein domains)

Aside from GSEA, I also did functional enrichment via overrepresentation analysis with TRAPID. This requires you to have a cds from your annotation, which must be processed and re-annotated on their website. It allows you to explore functional enrichment fairly easily, and I found that it produced a pretty high GO term annotation rate for quinoa. It also can find enriched protein domains which can assist in your interpretation. I used this for all functional enrichment outside of the core GSEA for each nutrient; unique DEGs for each deficiency, transcripts found from WGCNA and PLSDA, etc. I thought that its application was much more broad!

For the TRAPID annotation I used the PLAZA dicots database - you'll have to choose this very early on in the upload stage. To do functional enrichment, you'll have to upload lists of transcript IDs (only!!!) from your genes of interest (DEGs, etc) as subsets. Then on the sidebar for enrichment, you'll just pick a subset, then run for GO terms or protein domains. Then, I usually just exported the tsv files and used them in Rstudio! It was fun to play around with and look at overlaps between subsets too. 

Then, I usually ended up clustering the GO terms the same way as for GSEA, using rrvgo. The process is pretty much the same, but here is some general code I used: 
```
##rrvgo code to condense and cluster enriched GO terms - TRAPID overrepresented GO terms 

##TRAPID GO ENRICHMENT
library(rrvgo)
#library(tools)
library(org.At.tair.db)


trapid_dir <- path.expand("~/Nutrient/unique_transcripts")
trapid_unique_GO_files <- list.files(
  path = trapid_dir,
  pattern = "^go_enrichment_.*_0.05\\.tsv$",
  full.names = TRUE,
)

for (f in trapid_unique_GO_files) {

  message("Processing: ", f)
  trapid <- read.delim(f, comment.char = "")
  colnames(trapid)[1] <- "Aspect"
  file_base <- tools::file_path_sans_ext(basename(f))
  
  trapid <- subset(trapid, Aspect == "MF") #change if you run a different ontology
  # Remove NA / duplicate GO terms
  trapid <- trapid[!is.na(trapid$GO) & trapid$GO != "", ]
  trapid <- trapid[!duplicated(trapid$GO), ]
  
  # Skip if too few terms to cluster
  if (nrow(trapid) < 3) {
    message("  Skipping: fewer than 3 GO terms")
    next
  }
  
  scores <- setNames(-log10(trapid$q.value), trapid$GO)
  # Remove invalid scores
  scores <- scores[is.finite(scores)]
  
  if (length(scores) < 3) {
    message("  Skipping: insufficient valid scores")
    next
  }

  simMatrix <- calculateSimMatrix(
  names(scores),
  orgdb = 'org.At.tair.db',
  ont = "MF",
  method = "Rel"
)

if (is.null(simMatrix) || !is.matrix(simMatrix) && !is.data.frame(simMatrix) || nrow(simMatrix) < 2) {
  message("  Skipping ont (simMatrix is NULL, not a matrix, or too small)")
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
    # png(
    #   filename = file.path(trapid_dir, paste0(file_base, "_BP_scatter.png")),
    #   width = 1800, height = 1600, res = 300
    # )
    # print(scatterPlot(simMatrix, reducedTerms))
    # dev.off()

    # ## ---- Treemap plot ----
     png(
       filename = file.path(trapid_dir, paste0(file_base, "_MF_treemap.png")),
       width = 2000, height = 1600, res = 300
     )
     treemapPlot(reducedTerms)
     dev.off()
}

```
I used this process for PLS-DA, WGCNA, unique transcripts for each nutrient deficiency. Originally I also ran it on all DEG transcripts themselves, but ended up using GSEA instead for that. 
