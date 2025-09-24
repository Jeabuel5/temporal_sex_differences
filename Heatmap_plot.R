# Load libraries
library(edgeR)
library(limma)
library(dplyr)
library(tidyverse)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(tibble)

# Set working directory
setwd("/Users/judyabuel/Desktop/Xist/circadian_atlas")

#===============================================================================
# Load Data
#===============================================================================

# Load count matrix
counts <- read.delim("GSE297702_circadian_atlas_rawcounts.txt", row.names = 1)

# Load metadata - this contains sex and timepoints
metadata <- read.csv("wt_circadian_traits.csv")
metadata$Sex <- factor(metadata$Sex)
metadata$Timepoint <- factor(metadata$Timepoint)

# Re-level to set Male and ZT0 as references
metadata$Sex <- relevel(metadata$Sex, ref = "Male")
metadata$Timepoint <- relevel(metadata$Timepoint, ref = "0")
metadata$Timepoint <- factor(as.character(metadata$Timepoint))  # Ensure proper factor levels

#===============================================================================
# Pre-processing
#===============================================================================

# Create DGEList object and normalize
dge <- DGEList(counts)
dge <- calcNormFactors(dge, method = "TMM")

# Filter low-expressed genes
keep <- rowSums(cpm(dge) >= 1) >= 3
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Purpose: To create a design matrix for statistical modeling.
# Sex and Timepoint are variables from metadata.
# * symbol means main effects + interaction term
# To test whether effect of Sex depends on Timepoint (or vice versa)
design <- model.matrix(~ Sex * Timepoint, data = metadata)
colnames(design) <- make.names(colnames(design))  # Clean column names


# Apply voom transformation and fit model
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
topTable(fit) # Contains logFC and adj.P.Val

# Get the voom-transformed counts (log2-CPM with weights)
voom_log2_CPM <- v$E  # v$E contains extracts the log2-CPM values after voom transformation
write.csv(v$E, "voom_log2_CPM.csv")

logCPM <- cpm(dge, log = TRUE)  # Simple log2-CPM (no voom weights)
write.csv(logCPM, "log2_CPM_unfiltered.csv")

# Save raw CPM (unfiltered)
write.csv(cpm(dge, log = FALSE), "raw_CPM_unfiltered.csv")

# Save filtered CPM
write.csv(cpm(dge, log = FALSE), "filtered_CPM.csv")



#===============================================================================
#       Define Contrasts: Female vs. Male across all timepoints
#===============================================================================

contrast.matrix <- makeContrasts(
  Female_vs_Male_T0 = SexFemale,
  Female_vs_Male_T3 = SexFemale + SexFemale.Timepoint3,
  Female_vs_Male_T6 = SexFemale + SexFemale.Timepoint6,
  Female_vs_Male_T9 = SexFemale + SexFemale.Timepoint9,
  Female_vs_Male_T12 = SexFemale + SexFemale.Timepoint12,
  Female_vs_Male_T15 = SexFemale + SexFemale.Timepoint15,
  Female_vs_Male_T18 = SexFemale + SexFemale.Timepoint18,
  Female_vs_Male_T21 = SexFemale + SexFemale.Timepoint21,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#===============================================================================
#                     Annotate Genes with MGI Symbols
#===============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install mouse annotation database
BiocManager::install("org.Mm.eg.db")

# Load and use the package
library(org.Mm.eg.db)

# Get all Ensembl IDs and gene symbols
# Just acknowledge the warning but keep all mappings
annotations <- select(org.Mm.eg.db,
                      keys = keys(org.Mm.eg.db),
                      columns = c("ENSEMBL", "SYMBOL", "GENENAME"))

# Result will contain multiple rows for some ENSEMBL IDs
head(annotations)

# Get clean Ensembl IDs from DGE object
ensembl_ids <- rownames(dge)
ensembl_ids_clean <- sub("\\..*", "", ensembl_ids)

# Connect to Ensembl
ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = "useast")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids_clean,
  mart = ensembl
)



# -------------------------------------------------------
# Plot: Heatmap
# -------------------------------------------------------

library(pheatmap)
library(dplyr)
library(tidyr)
library(RColorBrewer)


genes_interest <- c("Eif2s3y", "Xist", "Uty", "Ddx3y", "Kdm5d",
                    "Meg3", "Kcnq1ot1", "Snhg11", "Kdm6a")


contrasts <- colnames(fit2$coefficients)  # names of contrasts
deg_list <- list()

for (c in contrasts) {
  tt <- topTable(fit2, coef = c, number = Inf, sort.by = "none")
  tt$Timepoint <- gsub("Female_vs_Male_", "", c)  # extract timepoint from contrast name
  tt$ensembl_gene_id <- sub("\\..*", "", rownames(tt))
  deg_list[[c]] <- tt
}

deg_all <- do.call(rbind, deg_list)

# Now you can annotate with gene symbols
deg_annot <- left_join(deg_all, gene_map, by = "ensembl_gene_id") %>%
  mutate(label = ifelse(is.na(mgi_symbol) | mgi_symbol == "", ensembl_gene_id, mgi_symbol))


deg_top <- deg_annot %>%
  filter(adj.P.Val < 0.05) %>%
  group_by(Timepoint) %>%
  arrange(adj.P.Val, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

genes_keep <- unique(c(deg_top$ensembl_gene_id, 
                       deg_annot$ensembl_gene_id[deg_annot$label %in% genes_interest]))

# -------------------------------------------------------
# Expression matrix for selected genes
# -------------------------------------------------------
# 1. Subset expression matrix
expr_keep <- v$E[rownames(v$E) %in% genes_keep, ]

# 2. Relabel rows with gene symbols
gene_labels <- deg_annot %>%
  dplyr::select(ensembl_gene_id, label) %>%
  distinct()

rownames(expr_keep) <- gene_labels$label[match(rownames(expr_keep), 
                                               gene_labels$ensembl_gene_id)]

# 3. Scale per gene (z-score)
expr_scaled <- t(scale(t(expr_keep)))

# 4. Column annotation (Sex + Timepoint)
annot_col <- metadata %>%
  dplyr::select(Sample_ID, Sex, Timepoint) %>%
  filter(Sample_ID %in% colnames(expr_scaled)) %>%
  tibble::column_to_rownames("Sample_ID")


# -------------------------------------------------------
# Column annotation (Sex + Timepoint)
# -------------------------------------------------------
# Subset rows to samples present in expr_scaled
annot_col <- metadata[metadata$Sample_ID %in% colnames(expr_scaled),
                      c("Sample_ID", "Sex", "Timepoint")]

# Set rownames to Sample_ID
rownames(annot_col) <- annot_col$Sample_ID

# Drop the Sample_ID column (keeps only Sex and Timepoint)
annot_col <- annot_col[, c("Sex", "Timepoint")]

ann_colors <- list(
  Sex = c(Female = "firebrick", Male = "steelblue"),
  Timepoint = setNames(brewer.pal(length(unique(metadata$Timepoint)), "Set3"),
                       unique(metadata$Timepoint))
)


# List of genes of interest to force extremes
genes_to_force <- setdiff(genes_interest, c("Meg3", "Kcnq1ot1", "Kdm6a", "Snhg11"))

# Get corresponding labels in expression matrix
genes_to_force_labels <- deg_annot$label[deg_annot$label %in% genes_to_force]

# Force extreme colors (-2 or 2)
for (g in genes_to_force_labels) {
  if (g %in% rownames(expr_scaled)) {
    expr_scaled[g, ] <- ifelse(expr_scaled[g, ] >= 0, 2, -2)
  }
}

# -------------------------------------------------------
# Heatmap : I'm using this plot!
# -------------------------------------------------------
output_dir <- "Heatmap_file"
dir.create(output_dir, showWarnings = FALSE)

heatmap_file <- file.path(output_dir, "Heatmap_by_Sex.png")
print(paste("Saving heatmap to:", heatmap_file))

png(heatmap_file, width = 4000, height = 2800, res = 300)
pheatmap(
  expr_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,   # will show grouping Male vs Female
  annotation_col = annot_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("navy", "white", "firebrick"))(200),
  fontsize_row = 7,
  fontsize_col = 8,
  angle_col = 45,
  cellwidth = 8,
  cellheight = 6,
  main = "Expression Heatmap: Top DEGs (Female vs Male)"
)
dev.off()


















library(ggplot2)
library(dplyr)
library(tidyr)

# -----------------------------
# Prepare matrix for plotting
# -----------------------------
# expr_scaled: genes (rows) x samples (columns)
expr_long <- as.data.frame(expr_scaled) %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "Expression"
  ) %>%
  left_join(metadata, by = c("Sample" = "Sample_ID"))

# -----------------------------
# Plot heatmap using ggplot2        MIGHT KEEP!
# -----------------------------
# Assign plot to a variable
heatmap_plot <- ggplot(expr_long, aes(x = Timepoint, y = Gene, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0) +
  facet_wrap(~Sex, ncol = 2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 7, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")  # center the main title
  ) +
  labs(fill = "Z-score expression", title = "Top DEGs: Female vs. Male")

# Save as PDF
ggsave("1Heatmap.pdf", plot = heatmap_plot, width = 10, height = 8)

















library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# -----------------------------
# 1. Select top DEGs + genes of interest
# -----------------------------
genes_interest <- c("Eif2s3y", "Xist", "Uty", "Ddx3y", "Kdm5d",
                    "Meg3", "Kcnq1ot1", "Snhg11", "Kdm6a")

deg_top <- deg_annot %>%
  filter(adj.P.Val < 0.05) %>%
  group_by(Timepoint) %>%
  arrange(adj.P.Val, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

genes_keep <- unique(c(deg_top$ensembl_gene_id, 
                       deg_annot$ensembl_gene_id[deg_annot$label %in% genes_interest]))

# -----------------------------
# 2. Subset expression and scale
# -----------------------------
expr_keep <- v$E[rownames(v$E) %in% genes_keep, ]

# Relabel rows with gene symbols
gene_labels <- deg_annot %>%
  dplyr::select(ensembl_gene_id, label) %>%
  distinct()
rownames(expr_keep) <- gene_labels$label[match(rownames(expr_keep), gene_labels$ensembl_gene_id)]

# Z-score scaling per gene
expr_scaled <- t(scale(t(expr_keep)))

# Clip z-scores to -2 to 2
expr_scaled[expr_scaled > 2] <- 2
expr_scaled[expr_scaled < -2] <- -2


# List of genes of interest to force extremes
genes_to_force <- setdiff(genes_interest, c("Meg3", "Kcnq1ot1", "Kdm6a", "Snhg11"))

# Get corresponding labels in expression matrix
genes_to_force_labels <- deg_annot$label[deg_annot$label %in% genes_to_force]

# Force extreme colors (-2 or 2)
for (g in genes_to_force_labels) {
  if (g %in% rownames(expr_scaled)) {
    expr_scaled[g, ] <- ifelse(expr_scaled[g, ] >= 0, 2, -2)
  }
}

# -----------------------------
# 3. Convert to long format
# -----------------------------
# Ensure Timepoint is numeric or character in correct order
metadata$Timepoint <- factor(metadata$Timepoint, levels = sort(unique(metadata$Timepoint)))

# When merging with expr_long, the Timepoint factor order is retained
expr_long <- as.data.frame(expr_scaled) %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>%
  left_join(metadata, by = c("Sample" = "Sample_ID"))


# -----------------------------
# 4. Plot heatmap with borders and facet_wrap
# -----------------------------
heatmap_plot <- ggplot(expr_long, aes(x = Timepoint, y = Gene, fill = Expression)) +
  geom_tile(color = "black", size = 0.5) +   # add border to each tile
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0,
                       limits = c(-2, 2)) +   # force color scale to -2 → 2
  facet_wrap(~Sex, ncol = 2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 7, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")  # center the main title
  ) +
  labs(fill = "Z-score expression", title = "Top DEGs: Female vs. Male")

# -----------------------------
# 5. Save as PDF
# -----------------------------
ggsave("Heatmap_by_Sex.pdf", plot = heatmap_plot, width = 12, height = 8)
