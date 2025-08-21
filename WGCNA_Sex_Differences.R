
#===============================================================================
# PURPOSE: DEG Analysis on SEX DIFFERENCES
# Circadian Data
#===============================================================================

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


# Define Contrasts: Female vs. Male across all timepoints
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

expr <- v$E

#===============================================================================
# Differential Co-Expression (Gene Network Topology Changes)
# Purpose: To analyze how co-expression network differ between sexes.
# Perform differential network analysis using Weighted Gene Co-expression
#                                             Network Analysis (WGCNA)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install missing WGCNA dependencies
BiocManager::install(c("impute", "preprocessCore"))

# Install WGCNA
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  install.packages("WGCNA")
}

# Load the package
library(WGCNA)


exprMale_t <- t(exprMale)
exprFemale_t <- t(exprFemale)

# Prepare male and female expression matrices
exprMale <- expr[, metadata$Sex == "Male"]
exprFemale <- expr[, metadata$Sex == "Female"]

# Subset metadata
male_samples <- metadata$Sex == "Male"
female_samples <- metadata$Sex == "Female"

# Subset expression matrix
exprMale <- expr[, male_samples]
exprFemale <- expr[, female_samples]

# Create networks in Males & plot the modules
netMale <- blockwiseModules(
  datExpr = exprMale_t,
  power = 6,
  TOMType = "signed", # directionally of correlation you want to see
  minModuleSize = 30, # minimum amount of genes 
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

plotDendroAndColors(
  netMale$dendrograms[[1]],           # the dendrogram for the first block
  netMale$colors[netMale$blockGenes[[1]]],  # module colors for that block
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)


# Create networks in Females & Plot module 
netFemale <- blockwiseModules(
  datExpr = exprFemale_t,
  power = 6,
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

plotDendroAndColors(
  netFemale$dendrograms[[1]],           # the dendrogram for the first block
  netFemale$colors[netFemale$blockGenes[[1]]],  # module colors for that block
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)





#===============================================================================
# Assess module preservation
multiExpr <- list(Male = list(data = exprMale_t),
                  Female = list(data = exprFemale_t))

# Create module color vector from Male network
multiColor <- list(Male = labels2colors(netMale$colors))

# Run module preservation
# This test whether those module in Male samples exist (are preserved) in Female samples.
preservation <- modulePreservation(
  multiExpr,
  multiColor,
  referenceNetworks = 1,  # Male is the reference
  nPermutations = 200,    # increase to 1000+ for publication
  verbose = 3
)

# View preservation stats
presZ <- preservation$preservation$Z$ref.Male
head(presZ)


# Visual Comparison of Dendrograms between sexes
moduleColorsMale <- labels2colors(netMale$colors)
moduleColorsFemale <- labels2colors(netFemale$colors)

# Plot side by side
# Save to a PDF file (adjust width and height as needed)
pdf("dendrogram_male_vs_female.pdf", width = 12, height = 6)

par(mfrow = c(1, 2))  # 1 row, 2 columns

plotDendroAndColors(
  netMale$dendrograms[[1]],
  moduleColorsMale[netMale$blockGenes[[1]]],
  groupLabels = "Male Modules",
  dendroLabels = FALSE,
  hang = 0.03,
  guideHang = 0.05
)

plotDendroAndColors(
  netFemale$dendrograms[[1]],
  moduleColorsFemale[netFemale$blockGenes[[1]]],
  groupLabels = "Female Modules",
  dendroLabels = FALSE,
  hang = 0.03,
  guideHang = 0.05
)

# Close the PDF device
dev.off()






#===============================================================================
# To analyze module-trait in Males

# Subset metadata to match Male samples
traitsMale <- metadata[metadata$Sex == "Male", ]


MEsMale <- moduleEigengenes(t(exprMale), colors = moduleColorsMale)$eigengenes

# Make a copy to avoid modifying original
traitsMale_numeric <- traitsMale
traitsMale_numeric <- traitsMale_numeric[match(colnames(exprMale), traitsMale$Sample_ID), ]

# Step 1: Remove non-numeric identifier columns
traitsMale_numeric <- traitsMale_numeric[, !(names(traitsMale_numeric) %in% c("Sample_ID", "Sample_Name"))]

# Step 2: Convert factor columns (like Timepoint and Sex) to numeric
# This will convert factor levels to 1, 2, 3, etc.
traitsMale_numeric$Timepoint <- as.numeric(as.character(traitsMale$Timepoint))  # or just as.numeric(traitsMale$Timepoint)
traitsMale_numeric$Sex <- as.numeric(traitsMale$Sex)  # e.g., Male = 1, Female = 2

# Step 3: Ensure the entire dataframe is numeric
traitsMale_numeric <- data.frame(lapply(traitsMale_numeric, as.numeric))

# Step 4: Optional - convert to matrix for heatmap
traitsMale_matrix <- as.matrix(traitsMale_numeric)

# Correlation and p-value matrices
moduleTraitCor_Male <- cor(MEsMale, traitsMale_numeric, use = "p") # p means use pairwise comparison

nSamples <- nrow(traitsMale_numeric)  # or nrow(exprMale_t)
moduleTraitPvalue_Male <- corPvalueStudent(moduleTraitCor_Male, nSamples)


traitsMale <- metadata[metadata$Sex == "Male", ]
traitsMale_numeric$Sex
traitsMale_numeric <- traitsMale_numeric[, sapply(traitsMale_numeric, function(x) length(unique(x)) > 1)]


# Text matrix: show correlations + p-values in heatmap
textMatrix_Male <- paste(signif(moduleTraitCor_Male, 2), "\n(",
                         signif(moduleTraitPvalue_Male, 1), ")", sep = "")

# Format correlation values to 2 digits
textMatrix_Male <- paste(signif(moduleTraitCor_Male, 2))


# Create a text matrix with stars only based on significance thresholds
textMatrix_Male <- matrix("", nrow = nrow(moduleTraitPvalue_Male), ncol = ncol(moduleTraitPvalue_Male))
textMatrix_Male[moduleTraitPvalue_Male < 0.001] <- "***"
textMatrix_Male[moduleTraitPvalue_Male < 0.01 & moduleTraitPvalue_Male >= 0.001] <- "**"
textMatrix_Male[moduleTraitPvalue_Male < 0.05 & moduleTraitPvalue_Male >= 0.01] <- "*"

# Create heatmap
pdf("module_trait_heatmap_male.pdf", width = 20, height = 18)

par(mar = c(10, 10, 4, 2))  # Adjust margins
par(family = "sans", font = 2)  # Use bold font

labeledHeatmap(
  Matrix = moduleTraitCor_Male,
  xLabels = colnames(moduleTraitCor_Male),
  yLabels = names(MEsMale),
  ySymbols = names(MEsMale),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_Male,  # only stars here
  setStdMargins = FALSE,
  cex.text = 1.0,
  zlim = c(-1, 1),
  main = "Module–Trait Relationships (Male)"
)

dev.off()






#===============================================================================
# To analyze module-trait in FEMALES

# Subset metadata to match Female samples
traitsFemale <- metadata[metadata$Sex == "Female", ]
traitsFemale <- traitsFemale[match(colnames(exprFemale), traitsFemale$Sample_ID), ]

# Calculate module eigengenes
MEsFemale <- moduleEigengenes(t(exprFemale), colors = moduleColorsFemale)$eigengenes

# Prepare numeric trait matrix
traitsFemale_numeric <- traitsFemale
traitsFemale_numeric <- traitsFemale_numeric[, !(names(traitsFemale_numeric) %in% c("Sample_ID", "Sample_Name"))]

traitsFemale_numeric$Timepoint <- as.numeric(as.character(traitsFemale$Timepoint))
traitsFemale_numeric$Sex <- as.numeric(traitsFemale$Sex)
traitsFemale_numeric <- data.frame(lapply(traitsFemale_numeric, as.numeric))

# Remove traits with no variation
traitsFemale_numeric <- traitsFemale_numeric[, sapply(traitsFemale_numeric, function(x) length(unique(x)) > 1)]

# Correlation and p-values
moduleTraitCor_Female <- cor(MEsFemale, traitsFemale_numeric, use = "p")
nSamples <- nrow(traitsFemale_numeric)  # or nrow(exprMale_t)
moduleTraitPvalue_Female <- corPvalueStudent(moduleTraitCor_Female, nSamples)


# Text matrix (stars only)
textMatrix_Female <- matrix("", nrow = nrow(moduleTraitPvalue_Female), ncol = ncol(moduleTraitPvalue_Female))
textMatrix_Female[moduleTraitPvalue_Female < 0.001] <- "***"
textMatrix_Female[moduleTraitPvalue_Female < 0.01 & moduleTraitPvalue_Female >= 0.001] <- "**"
textMatrix_Female[moduleTraitPvalue_Female < 0.05 & moduleTraitPvalue_Female >= 0.01] <- "*"

# Plot
pdf("module_trait_heatmap_female.pdf", width = 20, height = 18)

par(mar = c(10, 10, 4, 2))
par(family = "sans", font = 2)

labeledHeatmap(
  Matrix = moduleTraitCor_Female,
  xLabels = colnames(moduleTraitCor_Female),
  yLabels = names(MEsFemale),
  ySymbols = names(MEsFemale),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_Female,
  setStdMargins = FALSE,
  cex.text = 1.0,
  zlim = c(-1, 1),
  main = "Module–Trait Relationships (Female)"
)
dev.off()






#===============================================================================
# To analyze module-trait for Both sexes

# Prepare full metadata
traitsAll <- metadata[match(colnames(expr), metadata$Sample_ID), ]

# Create new module colors (optional — or reuse netMale if using male modules)
# moduleColorsAll <- labels2colors(netAll$colors)  # if available

# Eigengenes — using Male module colors (for consistency)
MEsAll <- moduleEigengenes(t(expr), colors = moduleColorsMale)$eigengenes

# Convert traits to numeric
traitsAll_numeric <- traitsAll[, !(names(traitsAll) %in% c("Sample_ID", "Sample_Name"))]
traitsAll_numeric$Timepoint <- as.numeric(as.character(traitsAll$Timepoint))
traitsAll_numeric$Sex <- as.numeric(traitsAll$Sex)
traitsAll_numeric <- data.frame(lapply(traitsAll_numeric, as.numeric))
traitsAll_numeric <- traitsAll_numeric[, sapply(traitsAll_numeric, function(x) length(unique(x)) > 1)]

# Correlation and p-values
moduleTraitCor_All <- cor(MEsAll, traitsAll_numeric, use = "p")
nSamples <- nrow(traitsAll_numeric)
moduleTraitPvalue_All <- corPvalueStudent(moduleTraitCor_All, nSamples)


# Text matrix (stars only)
textMatrix_All <- matrix("", nrow = nrow(moduleTraitPvalue_All), ncol = ncol(moduleTraitPvalue_All))
textMatrix_All[moduleTraitPvalue_All < 0.001] <- "***"
textMatrix_All[moduleTraitPvalue_All < 0.01 & moduleTraitPvalue_All >= 0.001] <- "**"
textMatrix_All[moduleTraitPvalue_All < 0.05 & moduleTraitPvalue_All >= 0.01] <- "*"

# Plot
pdf("module_trait_heatmap_all.pdf", width = 20, height = 18)

par(mar = c(10, 10, 4, 2))
par(family = "sans", font = 2)

labeledHeatmap(
  Matrix = moduleTraitCor_All,
  xLabels = colnames(moduleTraitCor_All),
  yLabels = names(MEsAll),
  ySymbols = names(MEsAll),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_All,
  setStdMargins = FALSE,
  cex.text = 1.0,
  zlim = c(-1, 1),
  main = "Module–Trait Relationships (All Samples)"
)
dev.off()





#===============================================================================
# To compare two correlations using Fisher's Z-Test
# Here we are comparing correlaton coefficients between two independent groups (Male vs. Female)

# z1 = atanh(r1)  # Fisher transform of correlation in males
# z2 = atanh(r2)  # Fisher transform of correlation in females
# SE = sqrt(1/(n1 - 3) + 1/(n2 - 3))  # Standard error
# Z = (z1 - z2)/SE                    # Z-score
# p = 2 * pnorm(-abs(Z))              # two-tailed p-value



traitsMale_numeric <- traitsMale_numeric[match(colnames(exprMale), traitsMale$Sample_ID), ]
traitsFemale_numeric <- traitsFemale_numeric[match(colnames(exprFemale), traitsFemale$Sample_ID), ]

traitsMale_numeric <- traitsMale_numeric[, sapply(traitsMale_numeric, function(x) {
  is.numeric(x) && sum(!is.na(x)) > 1 && length(unique(x[!is.na(x)])) > 1
})]

traitsFemale_numeric <- traitsFemale_numeric[, sapply(traitsFemale_numeric, function(x) {
  is.numeric(x) && sum(!is.na(x)) > 1 && length(unique(x[!is.na(x)])) > 1
})]

# Align trait columns between Male and Female correlation matrices
common_traits <- intersect(colnames(moduleTraitCor_Male), colnames(moduleTraitCor_Female))
moduleTraitCor_Male_aligned <- moduleTraitCor_Male[, common_traits, drop = FALSE]
moduleTraitCor_Female_aligned <- moduleTraitCor_Female[, common_traits, drop = FALSE]

# Ensure matrices are aligned
stopifnot(identical(dim(moduleTraitCor_Male_aligned), dim(moduleTraitCor_Female_aligned)))
stopifnot(identical(colnames(moduleTraitCor_Male_aligned), colnames(moduleTraitCor_Female_aligned)))
stopifnot(identical(rownames(moduleTraitCor_Male_aligned), rownames(moduleTraitCor_Female_aligned)))

# Get correct sample sizes (number of samples)
nMale <- sum(metadata$Sex == "Male")
nFemale <- sum(metadata$Sex == "Female")

# Initialize p-value matrix
diffPval <- matrix(NA, nrow = nrow(moduleTraitCor_Male_aligned), ncol = ncol(moduleTraitCor_Male_aligned))
rownames(diffPval) <- rownames(moduleTraitCor_Male_aligned)
colnames(diffPval) <- colnames(moduleTraitCor_Male_aligned)


diffPval_adj <- matrix(p.adjust(as.vector(diffPval), method = "fdr"), # this corrects for multiple testing
                       nrow = nrow(diffPval), 
                       dimnames = dimnames(diffPval))

# Fisher Z-test for each module–trait pair
for (i in 1:nrow(diffPval)) {
  for (j in 1:ncol(diffPval)) {
    r1 <- moduleTraitCor_Male_aligned[i, j]
    r2 <- moduleTraitCor_Female_aligned[i, j]
    
    if (!is.na(r1) && !is.na(r2) && abs(r1) < 1 && abs(r2) < 1) {
      z1 <- atanh(r1)
      z2 <- atanh(r2)
      SE <- sqrt(1 / (nMale - 3) + 1 / (nFemale - 3))
      z_diff <- (z1 - z2) / SE
      diffPval[i, j] <- 2 * pnorm(-abs(z_diff))  # two-tailed
    }
  }
}

# FDR correction: adjust all p-values in the matrix
diffPval_adj <- matrix(
  p.adjust(as.vector(diffPval), method = "fdr"),
  nrow = nrow(diffPval),
  dimnames = dimnames(diffPval)
)

# Transform to -log10(p) for heatmap
heatmapMatrix <- -log10(diffPval_adj) # transformation enhances visual interpretation
heatmapMatrix[!is.finite(heatmapMatrix)] <- 0
heatmapMatrix[heatmapMatrix > 10] <- 10  # cap extreme values for better color scale

# Create significance stars based on FDR-adjusted p-values
textMatrix <- matrix("", nrow = nrow(diffPval_adj), ncol = ncol(diffPval_adj))
for (i in 1:nrow(diffPval_adj)) {
  for (j in 1:ncol(diffPval_adj)) {
    p <- diffPval_adj[i, j]
    if (!is.na(p) && is.finite(p)) {
      textMatrix[i, j] <- ifelse(p < 0.001, "***",
                                 ifelse(p < 0.01, "**",
                                        ifelse(p < 0.05, "*", "")))
    }
  }
}

# Save heatmap to PDF
pdf("Difference_in_Module_Trait_M_vs_F.pdf", width = 20, height = 18)

par(mar = c(10, 10, 4, 2))
par(family = "sans", font = 2)  # bold font for axis labels

labeledHeatmap(
  Matrix = heatmapMatrix,
  xLabels = colnames(heatmapMatrix),
  yLabels = rownames(heatmapMatrix),
  ySymbols = rownames(heatmapMatrix),
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  colorLabels = FALSE,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(0, max(heatmapMatrix, na.rm = TRUE)),
  main = expression(paste("-log"[10], "(FDR-adjusted p) of Difference in Module–Trait Correlation (M vs F)"))
)

dev.off()



# This is a heatmap visualizing the significance of differences in module–trait correlations 
# between males and females (based on FDR-adjusted p-values from Fisher’s Z-test).

#===============================================================================

library(lattice) # for levelplot (optional)
library(ggplot2)  # for heatmap.2 (optional)

# Flatten the p-value matrix for sorting
pvals_flat <- as.data.frame(as.table(diffPval_adj))
colnames(pvals_flat) <- c("Module", "Trait", "FDR_pval")

# Remove NA and infinite values
pvals_flat <- pvals_flat[!is.na(pvals_flat$FDR_pval) & is.finite(pvals_flat$FDR_pval), ]

# Order by ascending FDR-adjusted p-value (most significant first)
pvals_flat <- pvals_flat[order(pvals_flat$FDR_pval), ]

# Select top N significant pairs
N <- 15
top_pairs <- head(pvals_flat, N)

# Extract row and column indices
top_modules <- unique(top_pairs$Module)
top_traits <- unique(top_pairs$Trait)

# Subset heatmapMatrix (-log10 FDR p-values) for these pairs
summaryHeatmapMatrix <- heatmapMatrix[top_modules, top_traits, drop = FALSE]

# Subset original adjusted p-values for stars
summaryPvals <- diffPval_adj[top_modules, top_traits, drop = FALSE]

# Create star matrix
summaryStars <- matrix("", nrow = length(top_modules), ncol = length(top_traits),
                       dimnames = list(top_modules, top_traits))
for (i in seq_along(top_modules)) {
  for (j in seq_along(top_traits)) {
    p <- summaryPvals[i, j]
    if (!is.na(p) && is.finite(p)) {
      summaryStars[i, j] <- ifelse(p < 0.001, "***",
                                   ifelse(p < 0.01, "**",
                                          ifelse(p < 0.05, "*", "")))
    }
  }
}

# Plot the heatmap
pdf("Top_Significant_Module_Trait_Differences_Heatmap.pdf", width = 12, height = 10)

par(mar = c(12, 12, 6, 3))  # bottom, left, top, right margins

labeledHeatmap(
  Matrix = summaryHeatmapMatrix,
  xLabels = colnames(summaryHeatmapMatrix),
  yLabels = rownames(summaryHeatmapMatrix),
  ySymbols = rownames(summaryHeatmapMatrix),
  colors = blueWhiteRed(50),
  textMatrix = summaryStars,
  colorLabels = FALSE,
  setStdMargins = FALSE,
  cex.text = 1.5,
  cex.lab = 1.3,
  zlim = c(0, max(summaryHeatmapMatrix, na.rm = TRUE)),
  main = paste0("Top ", " Significant Module-Trait Correlation Differences (M vs F)")
)

dev.off()



# Order modules by mean significance (across traits)
row_order <- rownames(summaryHeatmapMatrix)[
  order(rowMeans(summaryHeatmapMatrix, na.rm = TRUE), decreasing = TRUE)
]

# Order traits by mean significance (across modules)
col_order <- colnames(summaryHeatmapMatrix)[
  order(colMeans(summaryHeatmapMatrix, na.rm = TRUE), decreasing = TRUE)
]
# Reorder the matrix and star matrix
summaryHeatmapMatrix_ordered <- summaryHeatmapMatrix[row_order, col_order, drop = FALSE]
summaryStars_ordered <- summaryStars[row_order, col_order, drop = FALSE]

pheatmap(
  summaryHeatmapMatrix_ordered,
  display_numbers = summaryStars_ordered,
  number_color = "black",
  fontsize_number = 8,    # smaller font inside tiles
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = paste0("Top ", " Significant Module-Trait Correlation Differences (M vs F)"),
  cellwidth = 20,    # smaller cell width
  cellheight = 15,   # smaller cell height
  legend = TRUE,
  angle_col = 45,
  cluster_rows = FALSE,   # no clustering, respect order
  cluster_cols = FALSE
)




#===============================================================================
# To identify the number of genes in each module:
table(moduleColorsMale)
table(moduleColorsFemale)


names(moduleColorsMale) <- rownames(exprMale)  # Ensembl IDs as names
exprMEpink <- exprMale[genesInMEpink_Male, ]

geneIDs <- rownames(exprMale)
names(moduleColorsMale) <- geneIDs

exprMEpink_df <- data.frame(ensembl_gene_id = rownames(exprMEpink), exprMEpink)
exprMEpink_annotated <- merge(mapping, exprMEpink_df, by = "ensembl_gene_id", all.x = TRUE)

# For Male MEpink:
genesInMEpink_Male <- names(moduleColorsMale[moduleColorsMale == "pink"])

library(biomaRt)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mapping <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = genesInMEpink_Male,
  mart = mart
)

head(mapping)

write.csv(exprMEpink_annotated, file = "MEpink_genes_male.csv", row.names = FALSE)

#---------------------------------------------


# Extract MEmagenta genes from moduleColorsFemale
names(moduleColorsFemale) <- rownames(exprFemale)
genesInMEmagenta_Female <- names(moduleColorsFemale[moduleColorsFemale == "magenta"])

# Load biomaRt and query gene symbols
library(biomaRt)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

mapping_Female <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = genesInMEmagenta_Female,
  mart = mart
)

head(mapping_Female)

exprMEmagenta_Female <- exprFemale[genesInMEmagenta_Female, ]

exprMEmagenta_df_Female <- data.frame(
  ensembl_gene_id = rownames(exprMEmagenta_Female),
  exprMEmagenta_Female
)

exprMEmagenta_annotated_Female <- merge(
  mapping_Female,
  exprMEmagenta_df_Female,
  by = "ensembl_gene_id",
  all.x = TRUE
)

write.csv(exprMEmagenta_annotated_Female, file = "MEmagenta_genes_female.csv", row.names = FALSE)






