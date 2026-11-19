# How to calculate TPM counts for eFP browser
To make the eFP browser for Nick Provart, I had to transform my counts into TPM counts. This is how I did it:
```
#getting TPM for gene level counts for efpBrowser- Leaves

library(GenomicRanges)

exons <- gff[gff$type == "exon", ]
exon_tx <- exons$Parent

exon_gene <- mapping$gene_id[
  match(exon_tx, mapping$transcript_id)
]
valid <- !is.na(exon_gene)
exons <- exons[valid, ]
exon_gene <- exon_gene[valid]

exon_gr <- GRanges(
  seqnames = exons$seqnames,
  ranges = IRanges(start = exons$start, end = exons$end),
  strand = exons$strand
)

exons_by_gene <- split(exon_gr, exon_gene)

exons_by_gene_reduced <- reduce(exons_by_gene)

gene_lengths <- sapply(exons_by_gene_reduced, function(x) sum(width(x)))

# gene_lengths is now a named numeric vector: names = gene IDs, values = length in bp

# Counts matrix (remove first column with gene IDs)
leaves_counts <- as.matrix(nutrientAllLeaves_count_data)

# Match lengths to counts
leaves_lengths_kb <- gene_lengths[rownames(leaves_counts)] / 1000

# Compute RPK
leaves_rpk <- leaves_counts / leaves_lengths_kb

# Scaling factor per sample
leaves_scaling_factors <- colSums(leaves_rpk)

# TPM
leaves_allReps_counts_tpm <- sweep(leaves_rpk, 2, leaves_scaling_factors, "/") * 1e6

```

And then to combine all of my data into one table, this is what I did:
```
#to combine leaf and root TPM tables together for eFP browser

# combine TPM matrices (columns)
combined_nutrient_tpm <- cbind(leaves_allReps_counts_tpm, roots_allReps_counts_tpm)

# export to CSV
write.csv(combined_nutrient_tpm, "tpm_leaves_and_root_quinoa_nutrient_eFP.csv", row.names = TRUE)

```
