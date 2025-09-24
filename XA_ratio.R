

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
voom_log2_CPM <- v$E 


# Set Contrast for Female vs Male
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

# Extract results from Timepoint 6
tt_ZT6 <- topTable(fit2, coef = "Female_vs_Male_T6", number = Inf, sort.by = "P")

# Annotate genes
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install mouse annotation database
BiocManager::install("org.Mm.eg.db")

# Load and use the package
library(org.Mm.eg.db)

library(biomaRt)

# Connect to Ensembl (mouse)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Get chromosome info
gene_annot <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "external_gene_name"),
                    filters = "ensembl_gene_id",
                    values = rownames(voom_log2_CPM),
                    mart = mart)

# Clean: keep only standard chromosomes
gene_annot <- gene_annot[gene_annot$chromosome_name %in% c(1:19, "X", "Y"), ]


# Define Gene sets
auto_genes <- gene_annot$ensembl_gene_id[gene_annot$chromosome_name %in% as.character(1:19)]
x_genes    <- gene_annot$ensembl_gene_id[gene_annot$chromosome_name == "X"]
y_genes    <- gene_annot$ensembl_gene_id[gene_annot$chromosome_name == "Y"]

# Intersect with expression data
x_use    <- intersect(rownames(voom_log2_CPM), x_genes)
auto_use <- intersect(rownames(voom_log2_CPM), auto_genes)
y_use    <- intersect(rownames(voom_log2_CPM), y_genes)

xa_ratio <- colMeans(voom_log2_CPM[x_use, , drop = FALSE], na.rm = TRUE) /
  colMeans(voom_log2_CPM[auto_use, , drop = FALSE], na.rm = TRUE)

metadata$XA_ratio <- xa_ratio
summary(xa_ratio)  # check finite values

# Plot
boxplot(XA_ratio ~ Sex, data = metadata,
        main = "X:A dosage ratio by Sex")


ggplot(metadata, aes(x = Sex, y = XA_ratio, fill = Sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  scale_fill_manual(values = c("Male" = "#1F77B4", "Female" = "#E377C2")) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  theme_bw() +
  labs(title = "X:A dosage ratio by Sex", y = "X:A expression ratio")






ggplot(metadata, aes(x = Sex, y = XA_ratio, fill = Sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  scale_fill_manual(values = c("Male" = "#1F77B4", "Female" = "#E377C2")) +
  stat_compare_means(method = "t.test", label = "p.signif") +  # t-test on plot
  theme_bw() +
  labs(title = "X:A dosage ratio by Sex", y = "X:A expression ratio") +
  theme(
    panel.background = element_blank(),   # remove panel background
    plot.background = element_blank(),    # remove overall background
    panel.grid = element_blank(),         # remove grid lines
    axis.line = element_line(color = "black"), # add axis lines
    axis.ticks = element_line(color = "black")
  )





#-------------------KEEP THIS!--------------------------------------
library(ggplot2)
library(ggpubr)

ggplot(metadata, aes(x = Sex, y = XA_ratio, fill = Sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  scale_fill_manual(values = c("Male" = "#1F77B4", "Female" = "#E377C2")) +
  stat_compare_means(
    method = "t.test",
    label = "p.format",   # shows exact p-value instead of ns or stars
    label.y = max(metadata$XA_ratio, na.rm = TRUE) * 1.05  # position above boxes
  ) +
  labs(title = "X:A dosage ratio by Sex", y = "X:A expression ratio") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )




# Base R boxplot
boxplot(XA_ratio ~ Sex, data = metadata,
        main = "X:A dosage ratio by Sex",
        ylab = "X:A expression ratio",
        col = c("#1F77B4", "#E377C2"))

# Add t-test
t_test_res <- t.test(XA_ratio ~ Sex, data = metadata)
t_test_res

# Optional: add p-value to the plot
pval <- signif(t_test_res$p.value, 3)
mtext(paste("p =", pval), side = 3, line = 0.5, cex = 0.9)



# Quality Check
# The purpose of this code is to confirm that the samples declared as "Male" in the 
# metadata actually express genes found on the Y chromosome, while samples declared as 
# "Female" do not (or express them at very low levels). 
# This is a common and effective sanity check because the Y chromosome is male-specific.

y_interest <- c("Eif2s3y","Uty","Ddx3y","Kdm5d")
y_map <- gene_annot$ensembl_gene_id[gene_annot$external_gene_name %in% y_interest]

expr_y <- voom_log2_CPM[y_map, ]
boxplot(t(expr_y) ~ metadata$Sex, las=2, main="Y-linked gene expression")


# Boxplot grouped by Sex and Timepoint
# This will give you a side-by-side comparison at each timepoint, 
#     with males and females next to each other.
# Standard Error - error bars

boxplot(XA_ratio ~ interaction(Timepoint, Sex), 
        data = metadata, las = 2,
        col = c(Male = "lightblue", Female = "pink"),
        main = "X:A ratio across Timepoints and Sex",
        ylab = "X:A expression ratio")



#-----------------------------KEEP THIS!----------------------------------------
# Define all desired timepoint
all_timepoints <- c(0, 3, 6, 9, 12, 15, 18, 21)

# Convert Timepoint to factor with all levels
metadata$Timepoint <- factor(metadata$Timepoint, levels = all_timepoints)

# Updated plot
ggplot(metadata, aes(x = Timepoint, y = XA_ratio, fill = Sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Sex), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
              alpha = 0.5) +
  stat_compare_means(
    method = "t.test",
    label = "p.format",   # shows exact p-value instead of ns or stars
    label.y = max(metadata$XA_ratio, na.rm = TRUE) * 1.05) +
  labs(title = "X:A ratio by Sex and Timepoint",
       x = "Timepoint (ZT)",
       y = "X:A expression ratio") +
  scale_fill_manual(values = c("Male" = "#1F77B4", "Female" = "#E377C2")) +
  scale_color_manual(values = c("Male" = "#1F77B4", "Female" = "#E377C2")) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )





install.packages("ggpubr")   # only once
library(ggpubr)


ggplot(metadata, aes(x = factor(Timepoint, levels = c(0,3,6,9,12,15,18,21)), 
                     y = XA_ratio, fill = Sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Sex), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
              alpha = 0.6, size = 2) +
  labs(title = "X:A ratio by Sex and Timepoint",
       x = "Timepoint (ZT)",
       y = "X:A expression ratio") +
  scale_fill_manual(values = c("Male" = "#1F77B4", "Female" = "#E377C2")) +
  scale_color_manual(values = c("Male" = "#1F77B4", "Female" = "#E377C2")) +
  theme_bw() +
  stat_compare_means(
    aes(group = Sex),           # compare Male vs Female
    method = "t.test",          # t-test
    label = "p.signif"          # show significance stars
  )


  
  
# Fit linear model for XA_ratio
xa_lm <- lm(XA_ratio ~ Sex * Timepoint, data = metadata)

# ANOVA table
anova(xa_lm)

# Check estimated effects
summary(xa_lm)







