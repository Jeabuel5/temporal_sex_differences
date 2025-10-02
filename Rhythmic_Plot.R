library(ggplot2)
library(dplyr)
library(biomaRt)

# -----------------------------
# 1. Prepare expression matrix
# -----------------------------
expr_mat <- as.data.frame(t(voom_log2_CPM))
expr_mat$SampleID <- rownames(expr_mat)

# -----------------------------
# 2. Map Ensembl IDs to gene symbols
# -----------------------------
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(voom_log2_CPM),
  mart = ensembl
)

# Clean Ensembl IDs in mapping table
gene_map$ensembl_gene_id <- sub("\\..*", "", gene_map$ensembl_gene_id)

# -----------------------------
# 3. Merge expression with metadata
# -----------------------------
metadata$SampleID <- rownames(metadata)
expr_merged <- merge(metadata, expr_mat, by = "SampleID")

# Clean column names
colnames(expr_merged) <- sub("\\..*", "", colnames(expr_merged))

# -----------------------------
# 4. Define your genes of interest
# -----------------------------
# 1. Clean rownames (remove version numbers)
expr_ids <- sub("\\..*", "", rownames(voom_log2_CPM))

# 2. Your genes of interest (Ensembl IDs from biomaRt)
genes_of_interest <- c("ENSMUSG00000027981", "ENSMUSG00000028248", "ENSMUSG00000035545",
                       "ENSMUSG00000071337", "ENSMUSG00000074024", "ENSMUSG00000074415",
                       "ENSMUSG00000080316", "ENSMUSG00000086370", "ENSMUSG00000090086",
                       "ENSMUSG00000091712", "ENSMUSG00000092274", "ENSMUSG00000092341",
                       "ENSMUSG00000096929", "ENSMUSG00000097156", "ENSMUSG00000097545",
                       "ENSMUSG00000097767", "ENSMUSG00000098202", "ENSMUSG00000102657",
                       "ENSMUSG00000102854", "ENSMUSG00000109321")

# 3. Check which genes are present
present_genes <- genes_of_interest[genes_of_interest %in% expr_ids]
missing_genes <- setdiff(genes_of_interest, expr_ids)

cat("Genes present in the expression matrix:\n")
print(present_genes)

cat("\nGenes NOT found in the expression matrix:\n")
print(missing_genes)






















library(ggplot2)
library(dplyr)
library(tidyr)
library(org.Mm.eg.db)
library(edgeR)

install.packages("circacompare")
library(circacompare)


# -----------------------------
# 1. Clean the expression matrix first
# -----------------------------
voom_log2_CPM_clean <- voom_log2_CPM
voom_log2_CPM_clean[is.na(voom_log2_CPM_clean)] <- 0

# Strip version numbers from rownames
rownames(voom_log2_CPM_clean) <- sub("\\..*", "", rownames(voom_log2_CPM_clean))

# -----------------------------
# 2. Map gene symbols to Ensembl IDs
# -----------------------------
genes_of_interest <- c("Meg3", "Rnpc3", "Pnisr", "Leng8", "Snhg11", "Tia1",
                       "Mir100hg", "4632427E13Rik", "Spaca6", "Ftx",
                       "Al480526", "Sec14l5", "Neat1", "Malat1", "A330023F24Rik",
                       "Gm3764", "Mir124a-1hg", "Miat", "B830012L14Rik", "Kcnq1ot1",
                       "Gm37899", "C130023A14Rik", "A330076H08Rik")

ens_ids <- mapIds(org.Mm.eg.db,
                  keys = genes_of_interest,
                  column = "ENSEMBL",
                  keytype = "SYMBOL",
                  multiVals = "first")

# Keep only IDs present in expression matrix
present_ensembl <- ens_ids[ens_ids %in% rownames(voom_log2_CPM_clean)]
present_symbols <- names(present_ensembl)

cat("Genes found in expression matrix:\n")
print(present_symbols)
cat("Genes missing:\n")
print(setdiff(genes_of_interest, present_symbols))

# -----------------------------
# 3. Prepare data for CircaCompare
# -----------------------------
# Ensure metadata matches expression samples
stopifnot(all(colnames(voom_log2_CPM_clean) %in% rownames(metadata)))

# Force metadata rows to use the same IDs as expression matrix
metadata$SampleID <- colnames(voom_log2_CPM_clean)
rownames(metadata) <- metadata$SampleID

# Assign sex from prefix
metadata$Sex <- ifelse(grepl("^A", metadata$SampleID), "Male", "Female")
metadata$Sex <- factor(metadata$Sex, levels = c("Male", "Female"))

# Convert Timepoint to numeric for CircaCompare
print("Original Timepoint values:")
print(unique(metadata$Timepoint))

# Method 1: Extract numbers from timepoint strings
metadata$Timepoint <- as.numeric(gsub("[^0-9]", "", as.character(metadata$Timepoint)))

# If Method 1 doesn't work, use manual mapping:
# timepoint_mapping <- c("your_timepoint1" = hours1, "your_timepoint2" = hours2, ...)
# metadata$Timepoint <- timepoint_mapping[as.character(metadata$Timepoint)]

print("Converted Timepoint values:")
print(unique(metadata$Timepoint))
print(paste("Timepoint class:", class(metadata$Timepoint)))

# Remove any samples with NA timepoints
valid_samples <- !is.na(metadata$Timepoint)
metadata <- metadata[valid_samples, ]
voom_log2_CPM_clean <- voom_log2_CPM_clean[, valid_samples]


# -----------------------------
# Simplified Rhythmicity analysis followed by CircaCompare
# -----------------------------
circacompare_results <- list()
rhythmicity_results <- list()

for(gene in present_symbols) {
  ens_id <- present_ensembl[gene]
  
  # Create data frame
  circ_data <- data.frame(
    time = metadata$Timepoint,
    measure = as.numeric(voom_log2_CPM_clean[ens_id, ]),
    group = metadata$Sex
  )
  
  circ_data <- circ_data[complete.cases(circ_data), ]
  
  # Try CircaCompare directly with error handling
  tryCatch({
    cc_result <- circacompare(x = circ_data, 
                              col_time = "time", 
                              col_group = "group", 
                              col_outcome = "measure",
                              alpha_threshold = 0.05)
    
    summary_df <- as.data.frame(cc_result$summary)
    summary_df$Gene <- gene
    summary_df$Parameter <- rownames(summary_df)
    
    circacompare_results[[gene]] <- summary_df
    cat("✓ CircaCompare completed for:", gene, "\n")
    
  }, error = function(e) {
    # Store the error information
    error_msg <- e$message
    cat("✗ CircaCompare failed for:", gene, "-", error_msg, "\n")
    
    # Classify the error type
    if(grepl("arrhythmic", error_msg, ignore.case = TRUE)) {
      rhythmicity_results[[gene]] <- data.frame(
        Gene = gene,
        Status = "Arrhythmic",
        Group = ifelse(grepl("Female", error_msg), "Female", 
                       ifelse(grepl("Male", error_msg), "Male", "Unknown")),
        Message = error_msg
      )
    } else {
      rhythmicity_results[[gene]] <- data.frame(
        Gene = gene,
        Status = "Error",
        Group = "Unknown", 
        Message = error_msg
      )
    }
  })
}

# Check what we got
if(length(circacompare_results) > 0) {
  successful_genes <- names(circacompare_results)
  cat("\nSuccessfully analyzed genes:", length(successful_genes), "\n")
  print(successful_genes)
  
  all_results <- do.call(rbind, circacompare_results)
  print(head(all_results))
}

if(length(rhythmicity_results) > 0) {
  error_summary <- do.call(rbind, rhythmicity_results)
  cat("\nGenes with issues:", nrow(error_summary), "\n")
  print(error_summary)
}




# 1. Extract and format the successful results
if(length(circacompare_results) > 0) {
  final_results <- do.call(rbind, circacompare_results)
  
  # Clean up the parameter names
  final_results$Parameter <- rownames(final_results)
  rownames(final_results) <- NULL
  
  # Extract key comparisons
  key_params <- c("Mesor difference estimate", "P-value for mesor difference",
                  "Amplitude difference estimate", "P-value for amplitude difference", 
                  "Phase difference estimate", "P-value for phase difference")
  
  key_results <- final_results[final_results$parameter %in% key_params, ]
  print("Key comparison results:")
  print(key_results)
}

# 2. Analyze the arrhythmic genes pattern
if(length(rhythmicity_results) > 0) {
  arrhythmic_genes <- do.call(rbind, rhythmicity_results)
  
  cat("\n=== ARRHYTHMIC GENE SUMMARY ===\n")
  cat("Total genes with arrhythmic female data:", nrow(arrhythmic_genes), "\n")
  cat("Percentage of tested genes:", round(nrow(arrhythmic_genes)/length(present_symbols)*100, 1), "%\n")
  
  # These genes represent sex-specific circadian regulation
  sex_specific_genes <- arrhythmic_genes$Gene
  cat("\nGenes with male-specific rhythmicity (arrhythmic in females):\n")
  print(sex_specific_genes)
}

# 3. Visualize some examples - FIXED VERSION
library(ggplot2)
library(gridExtra)  # Needed for arranging multiple ggplot plots

plot_gene_comparison <- function(gene1, gene2) {
  # Successful gene
  ens_id1 <- present_ensembl[gene1]
  data1 <- data.frame(
    time = metadata$Timepoint,
    expression = as.numeric(voom_log2_CPM_clean[ens_id1, ]),
    sex = metadata$Sex
  )
  
  p1 <- ggplot(data1, aes(x = time, y = expression, color = sex)) +
    geom_point(size = 2) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
    scale_color_manual(values = c("Male" = "blue", "Female" = "red")) +
    ggtitle(paste("Rhythmic in both sexes:", gene1)) +
    labs(x = "Time", y = "Expression (log2 CPM)") +
    theme_minimal() +
    theme(legend.position = "top")
  
  # Arrhythmic gene
  ens_id2 <- present_ensembl[gene2]
  data2 <- data.frame(
    time = metadata$Timepoint,
    expression = as.numeric(voom_log2_CPM_clean[ens_id2, ]),
    sex = metadata$Sex
  )
  
  p2 <- ggplot(data2, aes(x = time, y = expression, color = sex)) +
    geom_point(size = 2) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
    scale_color_manual(values = c("Male" = "blue", "Female" = "red")) +
    ggtitle(paste("Arrhythmic in females:", gene2)) +
    labs(x = "Time", y = "Expression (log2 CPM)") +
    theme_minimal() +
    theme(legend.position = "top")
  
  # Arrange plots side by side
  grid.arrange(p1, p2, ncol = 2)
}

# Compare one successful vs one arrhythmic gene
plot_gene_comparison("A330023F24Rik", "Neat1")






### More plots here!
# Compare multiple genes in one plot
compare_multiple_genes <- function(genes) {
  plot_data <- data.frame()
  
  for(gene in genes) {
    ens_id <- present_ensembl[gene]
    gene_data <- data.frame(
      time = metadata$Timepoint,
      expression = as.numeric(voom_log2_CPM_clean[ens_id, ]),
      sex = metadata$Sex,
      gene = gene
    )
    plot_data <- rbind(plot_data, gene_data)
  }
  
  ggplot(plot_data, aes(x = time, y = expression, color = sex)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
    scale_color_manual(values = c("Male" = "blue", "Female" = "red")) +
    facet_wrap(~ gene, scales = "free_y") +
    ggtitle("Comparison of Gene Expression Patterns") +
    labs(x = "Time", y = "Expression (log2 CPM)") +
    theme_minimal()
}

# Compare successful vs arrhythmic genes
successful_genes <- c("A330023F24Rik", "Miat")
arrhythmic_genes <- c("Neat1", "Malat1")
compare_multiple_genes(c(successful_genes, arrhythmic_genes))





# --------------------------------------------------------------------------------------------
### Create a result table identifying specific genes that are rhythmic / arrhythmic in male vs female or both.
summary_table <- data.frame(
  Gene = character(),
  Status = character(),
  Male_Rhythmic = logical(),
  Female_Rhythmic = logical(),
  Mesor_Diff_Pval = numeric(),
  stringsAsFactors = FALSE
)

for(gene in present_symbols) {
  status <- ifelse(gene %in% names(circacompare_results), "Both Rhythmic", "Female Arrhythmic")
  
  if(gene %in% names(circacompare_results)) {
    result <- circacompare_results[[gene]]
    mesor_pval <- result[result$parameter == "P-value for mesor difference", "value"]
  } else {
    mesor_pval <- NA
  }
  
  summary_table <- rbind(summary_table, data.frame(
    Gene = gene,
    Status = status,
    Male_Rhythmic = TRUE,  # Since CircaCompare requires male to be rhythmic
    Female_Rhythmic = gene %in% names(circacompare_results),
    Mesor_Diff_Pval = mesor_pval
  ))
}

print("Summary of all genes:")
print(summary_table)





# ------------------------------------------------------------------------------
## Functional Enrichment Analysis for Male-Specific Genes
# Define your gene groups
male_specific_genes <- c("Meg3", "Rnpc3", "Pnisr", "Leng8", "Snhg11", "Tia1",
                         "Mir100hg", "Spaca6", "Ftx", "Sec14l5", "Neat1", 
                         "Malat1", "Kcnq1ot1", "A330076H08Rik")

both_rhythmic_genes <- c("A330023F24Rik", "Gm3764", "Mir124a-1hg", "Miat")

# Install required packages if needed
if (!require("clusterProfiler")) BiocManager::install("clusterProfiler")
if (!require("org.Mm.eg.db")) BiocManager::install("org.Mm.eg.db")
if (!require("enrichplot")) BiocManager::install("enrichplot")
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

# Convert gene symbols to Entrez IDs for enrichment analysis
male_specific_entrez <- mapIds(org.Mm.eg.db,
                               keys = male_specific_genes,
                               column = "ENTREZID",
                               keytype = "SYMBOL",
                               multiVals = "first")

both_rhythmic_entrez <- mapIds(org.Mm.eg.db,
                               keys = both_rhythmic_genes,
                               column = "ENTREZID",
                               keytype = "SYMBOL",
                               multiVals = "first")

# Remove NA values
male_specific_entrez <- na.omit(male_specific_entrez)
both_rhythmic_entrez <- na.omit(both_rhythmic_entrez)

# Run GO enrichment analysis for male-specific genes
go_enrichment_male <- enrichGO(gene          = male_specific_entrez,
                               OrgDb         = org.Mm.eg.db,
                               keyType       = "ENTREZID",
                               ont           = "BP",  # Biological Process
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.10,
                               readable      = TRUE)

# Run KEGG pathway enrichment
kegg_enrichment_male <- enrichKEGG(gene         = male_specific_entrez,
                                   organism     = "mmu",
                                   pvalueCutoff = 0.05)

# Display results
cat("=== GO ENRICHMENT - MALE-SPECIFIC GENES ===\n")
if(!is.null(go_enrichment_male) && nrow(go_enrichment_male) > 0) {
  print(head(go_enrichment_male, 10))
  
  # Plot top enriched terms
  dotplot(go_enrichment_male, showCategory=15, title="GO Biological Process - Male-Specific Genes")
} else {
  cat("No significant GO enrichment found\n")
}

cat("\n=== KEGG PATHWAYS - MALE-SPECIFIC GENES ===\n")
if(!is.null(kegg_enrichment_male) && nrow(kegg_enrichment_male) > 0) {
  print(head(kegg_enrichment_male, 10))
  dotplot(kegg_enrichment_male, title="KEGG Pathways - Male-Specific Genes")
} else {
  cat("No significant KEGG enrichment found\n")
}

### Most of these genes are lncRNAs like Neat1, Malat1, Meg3, Kcnq1ot1 which do not code for proteins.
### Does not best represented in KEGG which mainly focuses on protein-coding pathways
### these genes function as regulators rather than pathway components and could be involved in non-cononical
### mechanisms such as chromatin and nuclear organization as well as post-transcriptional regulation.



# 1. GO Biological Process (might work better than KEGG)
go_bp_male <- enrichGO(gene = male_specific_entrez,
                       OrgDb = org.Mm.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05)

if(nrow(go_bp_male) > 0) {
  dotplot(go_bp_male, title="GO BP - Male-Specific Circadian Genes")
} else {
  cat("No GO BP enrichment either - suggests highly specialized functions\n")
}

# 2. Molecular Function and Cellular Component
go_mf_male <- enrichGO(gene = male_specific_entrez,
                       OrgDb = org.Mm.eg.db,
                       ont = "MF",
                       pvalueCutoff = 0.05)

go_cc_male <- enrichGO(gene = male_specific_entrez,
                       OrgDb = org.Mm.eg.db,
                       ont = "CC",
                       pvalueCutoff = 0.05)


# Add plotting for MF
if(nrow(go_mf_male) > 0) {
  dotplot(go_mf_male, title="GO MF - Male-Specific Circadian Genes")
} else {
  cat("No GO MF enrichment found\n")
}

# Add plotting for CC  
if(nrow(go_cc_male) > 0) {
  dotplot(go_cc_male, title="GO CC - Male-Specific Circadian Genes")
} else {
  cat("No GO CC enrichment found\n")
}



# Run comprehensive GO analysis
go_analysis <- enrichGO(gene = male_specific_entrez,
                        OrgDb = org.Mm.eg.db,
                        ont = "ALL",  # Try ALL ontologies
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.10,
                        readable = TRUE)

if(nrow(go_analysis) > 0) {
  cat("=== GO ENRICHMENT FOUND! ===\n")
  
  # Separate by ontology
  bp_terms <- go_analysis[go_analysis$ONTOLOGY == "BP", ]
  mf_terms <- go_analysis[go_analysis$ONTOLOGY == "MF", ]
  cc_terms <- go_analysis[go_analysis$ONTOLOGY == "CC", ]
  
  cat("Biological Process terms:", nrow(bp_terms), "\n")
  cat("Molecular Function terms:", nrow(mf_terms), "\n") 
  cat("Cellular Component terms:", nrow(cc_terms), "\n")
  
  # Show top terms
  if(nrow(bp_terms) > 0) {
    print(head(bp_terms[, c("Description", "pvalue", "Count")]))
  }
  
} else {
  cat("=== NO GO ENRICHMENT - HIGHLY SPECIALIZED FUNCTIONS ===\n")
  
  # Even if no enrichment, we can still analyze manually
  manual_analysis()
}



# Which genes are involved in genomic imprinting?
imprinting_genes <- c("Meg3", "Kcnq1ot1")  # Both are well-known imprinted lncRNAs
cat("Genomic imprinting genes:", paste(imprinting_genes, collapse=", "), "\n")

# Splicing-related genes in your set
splicing_genes <- c("Rnpc3", "Pnisr", "Tia1", "Malat1")
cat("RNA splicing genes:", paste(splicing_genes, collapse=", "), "\n")



# Get the complete results to see all enriched terms
cat("=== COMPLETE GO ANALYSIS ===\n")

# Biological Process terms
if(nrow(bp_terms) > 0) {
  cat("\n--- ALL BIOLOGICAL PROCESS TERMS ---\n")
  print(bp_terms[, c("Description", "pvalue", "Count", "geneID")])
}

# Molecular Function terms
if(nrow(mf_terms) > 0) {
  cat("\n--- MOLECULAR FUNCTION TERMS ---\n")
  print(mf_terms[, c("Description", "pvalue", "Count", "geneID")])
}

# Cellular Component terms  
if(nrow(cc_terms) > 0) {
  cat("\n--- CELLULAR COMPONENT TERMS ---\n")
  print(cc_terms[, c("Description", "pvalue", "Count", "geneID")])
}

# Create a manual functional categorization
functional_categories <- data.frame(
  Gene = male_specific_genes,
  Known_Function = c(
    "Meg3" = "Imprinted lncRNA, cell cycle regulator",
    "Rnpc3" = "RNA binding protein, splicing factor", 
    "Pnisr" = "Splicing factor",
    "Leng8" = "tRNA modification",
    "Snhg11" = "Small nucleolar RNA host gene",
    "Tia1" = "RNA binding protein, stress granules",
    "Mir100hg" = "MicroRNA host gene",
    "Spaca6" = "Sperm-associated antigen",
    "Ftx" = "X-inactivation regulator",
    "Sec14l5" = "Lipid binding protein",
    "Neat1" = "Architectural lncRNA, paraspeckles",
    "Malat1" = "Nuclear speckle lncRNA, splicing regulator",
    "Kcnq1ot1" = "Imprinted lncRNA, chromatin silencing",
    "A330076H08Rik" = "Uncharacterized"
  )
)

print(functional_categories)





### Rain Plot
library(ggplot2)
library(dplyr)
library(circular)  # For circular statistics

goi_data <- read.csv("GOI_rhythmicity_results_annotated.csv")

head(goi_data)
str(goi_data)

# Make sure your data is properly cleaned
goi_data_clean <- goi_data %>%
  filter(!is.na(phase_hours), !is.na(p_value), !is.na(Sex)) %>%
  mutate(
    phase_hours_adj = phase_hours %% 24,
    log_pval = -log10(p_value)
  )

rain_plot <- ggplot(goi_data_clean, 
                    aes(x = phase_hours_adj, 
                        y = log_pval,
                        color = Sex)) +
  geom_point(alpha = 0.6, size = 2) +  # THIS LINE ADDS THE DATA POINTS
  scale_x_continuous(limits = c(0, 24), 
                     breaks = c(0, 6, 12, 18, 24),
                     labels = c("ZT0", "ZT6", "ZT12", "ZT18", "ZT24")) +
  labs(x = "Zeitgeber Time", 
       y = "-log10(p-value)",  # CORRECTED LABEL
       color = "Sex") +
  theme_minimal() +
  scale_color_manual(values = c("Male" = "blue", "Female" = "red"))
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    legend.position = "top"
  )

print(rain_plot)

# Enhanced version with your actual ZT timepoints
enhanced_rain <- ggplot(goi_data, 
                        aes(x = phase_hours, 
                            y = -log10(p_value))) +
  geom_point(aes(color = Sex, size = amplitude), alpha = 0.6) +
  scale_x_continuous(limits = c(0, 24), 
                     breaks = c(0, 3, 6, 9, 12, 15, 18, 21, 24),
                     labels = c("ZT0", "ZT3", "ZT6", "ZT9", "ZT12", "ZT15", "ZT18", "ZT21", "ZT24")) +
  labs(x = "Zeitgeber Time", 
       y = "-log10(p-value)",
       title = "Circadian Phase Distribution",
       color = "Sex", size = "Amplitude") +
  theme_minimal() +
  facet_wrap(~ Sex) +
  scale_color_manual(values = c("Male" = "blue", "Female" = "red"))

print(enhanced_rain)

# Phase histogram with ZT labeling
phase_histogram <- ggplot(goi_data, aes(x = phase_hours, fill = Sex)) +
  geom_histogram(binwidth = 2, alpha = 0.7, position = "identity") +
  scale_x_continuous(limits = c(0, 24), 
                     breaks = c(0, 3, 6, 9, 12, 15, 18, 21, 24),
                     labels = c("ZT0", "ZT3", "ZT6", "ZT9", "ZT12", "ZT15", "ZT18", "ZT21", "ZT24")) +
  labs(x = "Zeitgeber Time", 
       y = "Number of Genes",
       title = "Phase Distribution Histogram") +
  theme_minimal() +
  scale_fill_manual(values = c("Male" = "blue", "Female" = "red")) +
  facet_wrap(~ Sex, ncol = 1)

print(phase_histogram)

# Statistical analysis (same as before)
goi_data$phase_circular <- goi_data$phase_hours %% 24

male_phases <- circular(goi_data$phase_circular[goi_data$Sex == "Male"] * (2*pi/24))
female_phases <- circular(goi_data$phase_circular[goi_data$Sex == "Female"] * (2*pi/24))

# Rayleigh test
rayleigh_male <- rayleigh.test(male_phases)
rayleigh_female <- rayleigh.test(female_phases)

cat("Rayleigh test for non-uniform distribution:\n")
cat("Male p-value:", rayleigh_male$p.value, "\n")
cat("Female p-value:", rayleigh_female$p.value, "\n")

# Watson-Wheeler test
if(length(male_phases) > 0 & length(female_phases) > 0) {
  watson_test <- watson.wheeler.test(list(male_phases, female_phases))
  cat("\nWatson-Wheeler test for phase distribution differences:\n")
  cat("p-value:", watson_test$p.value, "\n")
}

# Summary statistics
phase_summary <- goi_data %>%
  group_by(Sex) %>%
  summarise(
    n_genes = n(),
    mean_phase = mean(phase_circular),
    sd_phase = sd(phase_circular),
    mean_amplitude = mean(amplitude),
    mean_pvalue = mean(p_value),
    n_sig = sum(p_value < 0.05)
  )

print(phase_summary)
