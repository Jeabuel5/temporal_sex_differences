
#===============================================================================
#         Rhythmic Analysis: Core Circadian Detection
#===============================================================================

# Purpose: 
# Identify genes that are rhythmic in both sexes, only in male or only in females.
# Compare amplitude, phase, period, and baseline expression of rhythmic genes across sexes.


# Install necessary packages if not already installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "cosinor", "cosinor2", "circacompare", "tidyverse", "ggplot2", "ggpubr"))

# Load libraries
library(limma)
library(edgeR)
library(cosinor2)
library(circacompare)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(dplyr)

# Set working directory
setwd("/Users/judyabuel/Desktop/Xist/circadian_atlas")

# Load count matrix & metadata
count_matrix <- read.delim("GSE297702_circadian_atlas_rawcounts.txt", row.names = 1)
metadata <- read.csv("wt_circadian_traits.csv")


#-------------------------------------------------------------------------------
###                 Organize your data
#-------------------------------------------------------------------------------

# Create metadata from your sample names
sample_names <- colnames(count_matrix)

# Align both data to match information
# First, ensure metadata uses consistent naming
metadata <- metadata %>% 
  mutate(Sample = Sample_ID)  # Create a clean 'Sample' column matching count matrix

# Find ACTUAL overlapping samples
common_samples <- intersect(colnames(count_matrix), metadata$Sample)

# Filter BOTH datasets to only these samples
count_matrix <- count_matrix[, common_samples, drop=FALSE]
metadata <- metadata %>% filter(Sample %in% common_samples)

# Verify
print(data.frame(
  CountMatrix = colnames(count_matrix),
  Metadata = metadata$Sample[match(colnames(count_matrix), metadata$Sample)]
))

# Check timepoints and sexes
table(metadata$Timepoint)
table(metadata$Sex)

# Ensure metadata rows match count matrix columns exactly
metadata <- metadata[match(colnames(count_matrix), metadata$Sample), ]

###-----------------------------------------------------------------------------
# Load required packages
library(edgeR)
library(rain)
library(VennDiagram)

### 1. Data Preparation -------------------------------------------------
# Read count matrices (genes as rows, samples as columns)
count_matrix <- read.delim("GSE297702_circadian_atlas_rawcounts.txt", row.names = 1)
metadata <- read.csv("wt_circadian_traits.csv")

### 2. Data Cleaning ---------------------------------------------------
# Remove samples with zero counts
clean_data <- function(count_matrix, metadata) {
  valid_samples <- colSums(count_matrix) > 0
  list(
    counts = count_matrix[, valid_samples],
    timepoints = metadata$Timepoint[match(colnames(count_matrix), metadata$Sample_ID)][valid_samples]
  )
}

male_clean <- clean_data(male_data, metadata)
female_clean <- clean_data(female_data, metadata)

### 3. Normalization ---------------------------------------------------
normalize_data <- function(count_matrix) {
  dge <- DGEList(counts = count_matrix)
  dge <- calcNormFactors(dge)
  cpm(dge, log=TRUE)
}

male_norm <- normalize_data(male_clean$counts)
female_norm <- normalize_data(female_clean$counts)

### 4. Quality Control - PCA --------------------------------------------
plot_pca <- function(norm_data, timepoints, title) {
  time_colors <- rainbow(length(unique(timepoints)))[as.factor(timepoints)]
  pca <- prcomp(t(norm_data), scale.=TRUE)
  
  plot(pca$x[,1], pca$x[,2],
       col=time_colors, pch=16, cex=1.5,
       main=title,
       xlab=paste0("PC1 (", round(100*pca$sdev[1]^2/sum(pca$sdev^2),1), "%"),
       ylab=paste0("PC2 (", round(100*pca$sdev[2]^2/sum(pca$sdev^2),1), "%"))
  
  legend("topright", legend=unique(timepoints),
         col=rainbow(length(unique(timepoints))), pch=16,
         title="Timepoints")
}

par(mfrow=c(1,2))
plot_pca(male_norm, male_clean$timepoints, "Male Samples")
plot_pca(female_norm, female_clean$timepoints, "Female Samples")
par(mfrow=c(1,1))


### 5. Visualization: Heatmap Plot --------------------------------------------
library(cosinor)

results <- list()


for(sex in c("Male", "Female")) {
  sex_samples <- metadata %>% filter(Sex == sex)
  sex_expr <- logCPM[, sex_samples$Sample, drop = FALSE]  # Ensure matrix structure
  
  results[[sex]] <- apply(sex_expr, 1, function(y) {
    model_data <- data.frame(
      y = y,
      time = sex_samples$Timepoint
    )
    
    tryCatch({
      fit <- cosinor.lm(y ~ time(time), data = model_data, period = 24)
      data.frame(
        p_value = cosinor.PR(fit)$p,
        amplitude = fit$coefficients[2],
        phase = fit$coefficients[3],
        mesor = fit$coefficients[1]
      )
    }, error = function(e) {
      data.frame(
        p_value = NA,
        amplitude = NA,
        phase = NA,
        mesor = NA
      )
    })
  }) %>% bind_rows(.id = "gene")
}

results$Male   # Rhythm parameters for male-specific genes
results$Female # Rhythm parameters for female-specific genes


### Phase comparison:
# Convert phase to hours
results$Male$phase_hours <- results$Male$phase * 24/(2*pi)

### Amplitude Filtering
sig_genes <- results$Female %>% filter(p_value < 0.05 & amplitude > 1)



# 1. Prepare annotation data
time_annotation <- data.frame(
  Timepoint = factor(
    metadata$Timepoint, 
    levels = sort(unique(metadata$Timepoint)),
    labels = paste0("ZT", sort(unique(metadata$Timepoint)))
  ),
  row.names = metadata$Sample
)

# 2. Create matching color scheme
time_colors <- list(
  Timepoint = setNames(
    viridis::viridis(length(unique(metadata$Timepoint))),
    paste0("ZT", sort(unique(metadata$Timepoint)))
  )
)

library(ComplexHeatmap)
library(circlize)

# Create heatmap
hm <- Heatmap(
  matrix = logCPM[sig_genes$gene, ],
  name = "Expression",
  col = colorRamp2(c(-3, 0, 3), c("navy", "white", "firebrick3")),
  column_title = "Circadian Expression Patterns\nby Sex and Timepoint",
  show_row_names = ifelse(nrow(sig_genes) <= 50, TRUE, FALSE),
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8)
)

# Save as PDF
pdf("acircadian_heatmap.pdf", width = 12, height = 9)
draw(hm, newpage = FALSE)  # Crucial: newpage = FALSE
dev.off()












# 3. Identify rhythmic genes (FDR < 0.05)
rhythmic_genes <- list(
  Male = results$Male %>% 
    mutate(p_adj = p.adjust(p_value, "fdr")) %>% 
    filter(p_adj < 0.05),
  Female = results$Female %>% 
    mutate(p_adj = p.adjust(p_value, "fdr")) %>% 
    filter(p_adj < 0.05)
)

# 4. Classify genes by sex-specificity
gene_categories <- data.frame(
  gene = unique(c(rhythmic_genes$Male$gene, rhythmic_genes$Female$gene)),
  Male = FALSE,
  Female = FALSE
)

gene_categories$Male[gene_categories$gene %in% rhythmic_genes$Male$gene] <- TRUE
gene_categories$Female[gene_categories$gene %in% rhythmic_genes$Female$gene] <- TRUE

gene_categories <- gene_categories %>%
  mutate(category = case_when(
    Male & Female ~ "Both",
    Male ~ "Male_only",
    Female ~ "Female_only"
  ))

gene_categories <- gene_categories %>% 
  left_join(results$Male, by = "gene") %>% 
  mutate(
    phase_hours = (phase / (2 * pi)) * 21
  ) %>% 
  filter(
    !is.na(phase_hours),
    phase_hours >= 0 & phase_hours <= 21
  )



#-------------------------------------------------------------------------------
# Statistical Test: Watson-Wheeler Test (Non-parametric, for comparing 2 groups)
#-------------------------------------------------------------------------------

# Install and load libraries
install.packages("circular")  # For circular statistics
install.packages("CircStats") # Additional circular tools

library(circular)
library(CircStats)

gene_categories <- gene_categories %>%
  mutate(
    phase_radians = (phase_hours / 21) * 2 * pi  # Convert hours to radians (0–21h → 0–2π)
  )

# Split data by category
both_phases <- gene_categories %>% filter(category == "Both") %>% pull(phase_radians)
female_phases <- gene_categories %>% filter(category == "Female_only") %>% pull(phase_radians)
male_phases <- gene_categories %>% filter(category == "Male_only") %>% pull(phase_radians)

# Run Statistical Test (Compare All Groups)
watson.wheeler.test(list(
  "Both" = circular(both_phases),
  "Female" = circular(female_phases),
  "Male" = circular(male_phases)
))

# Extract p-value
watson_p_value <- ww_test$p.value


#-------------------------------------------------------------------------------
###                           PLOTTING DATA
#-------------------------------------------------------------------------------

# This plot cluster genes by phase shift between sexes across timepoint
# 1. First, create the base plot without adding it to itself
phase_plot <- ggplot(gene_categories, 
                     aes(x = phase_hours, fill = category)) + 
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c(
    "Both" = "purple", 
    "Female_only" = "red", 
    "Male_only" = "blue"
  )) +
  scale_x_continuous(
    limits = c(0, 21),
    breaks = seq(0, 21, by = 3),
    labels = seq(0, 21, by = 3)
  ) +
  labs(
    title = "Circadian Phase Distribution by Gene Category (0-21h)",
    x = "Circadian Phase (hours)",
    y = "Density",
    fill = "Category"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.position = "top"
  )

# Format p-value to avoid "0"
p_label <- ifelse(
  watson_p_value < 1e-16, 
  "Watson-Wheeler p < 2e-16", 
  paste("Watson-Wheeler p =", signif(watson_p_value, 3))
)

# 2. Add the Watson-Wheeler p-value as a subtitle
phase_plot <- phase_plot + 
  labs(subtitle = paste("Watson-Wheeler p =", format.pval(watson_p_value, eps = 1e-16))) +
  theme(plot.subtitle = element_text(size = 10, hjust = 0.5))

# 3. Save the plot
ggsave(
  filename = "Phase_Distribution_by_Category.pdf", 
  plot = phase_plot,
  device = "pdf",
  width = 8,
  height = 5,
  bg = "white"
)



#### Rose plot for circular data------------------------------------------------


ggplot(gene_categories, aes(x = phase_hours, fill = category)) +
  geom_histogram(alpha = 0.5, binwidth = 1) +
  coord_polar() +  
  scale_x_continuous(limits = c(0, 21), breaks = seq(0, 21, by = 3))


gene_categories %>%
  group_by(category) %>%
  summarise(mean_phase_h = mean(phase_hours))


set.seed(123)
gene_categories %>%
  filter(category == "Both") %>%
  sample_n(5000)  # Downsample for balanced comparisons


#-------------------------------------------------------------------------------
# To analyze phase shifts between male and female


# Example data (replace with your actual data)
phase_shift <- gene_categories(
  Time = rep(c(0, 3, 6, 9, 12, 15, 18, 21), each = 4),
  Expression = c(rnorm(32, mean = -1.0, sd = 0.2)), 
  Sex = rep(c("Female", "Male"), each = 16)
)


# Fit a model comparing sexes
result <- circacompare(
  x = phase_shift, 
  col_time = "Time", 
  col_group = "Sex", 
  col_outcome = "Expression",
  period = 24,  # Circadian period (hours)
  alpha_threshold = 0.05  # Significance cutoff
)

# View phase difference (in hours)
print(result$summary)





### Create and save the category counts plot: BOX PLOT--------------------------


category_plot <- ggplot(gene_categories, aes(x=category, fill=category)) +
  geom_bar() +
  scale_fill_manual(values=c("Both"="purple", "Male_only"="blue", "Female_only"="red")) +
  labs(title="Rhythmic Genes by Category", x="", y="Count")

ggsave("Rhythmic_Gene_Categories.pdf", 
       plot = category_plot,
       width = 8, 
       height = 6,
       device = "pdf")


### Create a Rain Plot----------------------------------------------------------

# Install required packages
install.packages("ggdist")  # For raincloud plots
library(ggdist)


# Create the plot
rain_plot_landscape <- ggplot(gene_categories, aes(x = phase_hours, y = category, fill = category)) +
  # Half-density (mirrored)
  stat_halfeye(
    adjust = 0.5,
    width = 0.7,
    .width = 0,
    justification = -0.1,
    orientation = "y"  # Critical for landscape
  ) +
  # Boxplot (horizontal)
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.5,
    orientation = "y"  # Critical for landscape
  ) +
  # Raw data (optional for large datasets)
  geom_point(
    size = 1.3,
    alpha = 0.05,  # Reduced opacity for visibility
    color = "black",
    position = position_jitter(height = 0.1, width = 0)
  ) +
  # Custom colors
  scale_fill_manual(values = c("Both" = "purple", "Female_only" = "red", "Male_only" = "blue")) +
  # Axes/labels
  scale_x_continuous(
    limits = c(0, 21),
    breaks = seq(0, 21, by = 3)
  ) +
  labs(
    title = "Circadian Phase Distribution by Gene Category (0-21h)",
    x = "Circadian Phase (hours)",
    y = "Category",
    fill = "Category"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5), # bold axis lines
    panel.grid = element_blank(), # Remove all grid lines
    axis.line = element_line(colour = "black", linewidth = 1),
    axis.text = element_text(face = "bold"), # Bold axis labels
    axis.title = element_text(face = "bold", size = 10),
    legend.position = "none",
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )

# Save as landscape PDF
ggsave(
  filename = "Phase_Distribution_RainPlot_Landscape.pdf",
  plot = rain_plot_landscape,
  device = "pdf",
  width = 10,  # Wider for landscape
  height = 6,
  bg = "white"
)





# Load required packages
library(clusterProfiler)
library(org.Mm.eg.db)
library(biomaRt)
library(ggplot2)
library(AnnotationDbi)
library(GenomicRanges)
library(multiMiR)

# 1. Enhanced gene ID conversion function
convert_to_entrez <- function(ensembl_ids) {
  # Input validation
  if (length(ensembl_ids) == 0 || all(is.na(ensembl_ids))) {
    warning("Empty or invalid input gene list")
    return(character(0))
  }
  
  # Clean IDs
  clean_ids <- unique(sub("\\..*", "", ensembl_ids))
  clean_ids <- clean_ids[!is.na(clean_ids) & nchar(clean_ids) > 0]
  if (length(clean_ids) == 0) {
    warning("No valid Ensembl IDs after cleaning")
    return(character(0))
  }
  
  # Try biomaRt first
  mart <- tryCatch({
    useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  }, error = function(e) NULL)
  
  if (!is.null(mart)) {
    id_map <- tryCatch({
      getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
            filters = "ensembl_gene_id",
            values = clean_ids,
            mart = mart)
    }, error = function(e) data.frame())
    
    if (nrow(id_map) > 0 && "entrezgene_id" %in% colnames(id_map)) {
      return(na.omit(unique(id_map$entrezgene_id)))
    }
  }
  
  # Fallback to AnnotationDbi
  suppressWarnings({
    mapped <- tryCatch({
      mapIds(org.Mm.eg.db,
             keys = clean_ids,
             column = "ENTREZID",
             keytype = "ENSEMBL",
             multiVals = "first")
    }, error = function(e) character(0))
  })
  
  return(na.omit(unique(as.character(mapped))))
}

# 2. Extract and validate male active-phase genes
get_phase_genes <- function(data, sex, phase_start, phase_end) {
  genes <- subset(data, 
                  category == paste0(sex, "_only") &
                    phase_hours >= phase_start & 
                    phase_hours < phase_end)$gene
  
  if (length(genes) == 0) {
    genes <- subset(data, category == paste0(sex, "_only"))$gene
    message("No genes in ZT", phase_start, "-", phase_end, 
            " window. Using all ", sex, "-specific genes: ", length(genes))
  }
  
  # Clean IDs
  genes <- sub("\\..*", "", genes)
  genes <- genes[grepl("^ENSMUSG", genes)]
  
  if (length(genes) == 0) {
    stop("No valid ", sex, " genes found after filtering.")
  }
  return(genes)
}

# 3. Get gene biotypes
get_gene_biotypes <- function(ensembl_ids) {
  mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  getBM(attributes = c("ensembl_gene_id", "gene_biotype"),
        filters = "ensembl_gene_id",
        values = ensembl_ids,
        mart = mart)
}

# 4. Main analysis pipeline
analyze_sex_phase_genes <- function(data, sex, phase_start, phase_end) {
  # Get genes
  genes <- get_phase_genes(data, sex, phase_start, phase_end)
  
  # Get biotypes
  biotypes <- get_gene_biotypes(genes)
  
  # Separate coding and non-coding
  coding_genes <- biotypes[biotypes$gene_biotype == "protein_coding", "ensembl_gene_id"]
  noncoding_genes <- setdiff(genes, coding_genes)
  
  # Convert coding genes to Entrez
  entrez_ids <- convert_to_entrez(coding_genes)
  
  # Functional enrichment
  if (length(entrez_ids) >= 5) {
    ego <- enrichGO(entrez_ids, OrgDb = org.Mm.eg.db, ont = "BP",
                    pAdjustMethod = "BH", pvalueCutoff = 0.05)
    kk <- enrichKEGG(entrez_ids, organism = "mmu", pvalueCutoff = 0.05)
  } else {
    warning("Insufficient ", sex, " genes for enrichment (", length(entrez_ids), ")")
    ego <- NULL
    kk <- NULL
  }
  
  # Non-coding analysis
  if (length(noncoding_genes) > 0) {
    nc_biotypes <- table(biotypes$gene_biotype[biotypes$ensembl_gene_id %in% noncoding_genes])
    message("\nNon-coding gene biotypes (", sex, "):")
    print(nc_biotypes)
  }
  
  return(list(
    genes = genes,
    coding_genes = coding_genes,
    noncoding_genes = noncoding_genes,
    go_results = ego,
    kegg_results = kk
  ))
}

# 5. Run analyses
female_results <- analyze_sex_phase_genes(gene_categories, "Female", 0, 6)
male_results <- analyze_sex_phase_genes(gene_categories, "Male", 12, 18)

# 6. Visualization
plot_enrichment <- function(enrich_result, title) {
  if (!is.null(enrich_result)) {
    dotplot(enrich_result, showCategory = 15, font.size = 10) +
      ggtitle(title) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold", size = 12))
  }
}

female_go_plot <- plot_enrichment(female_results$go_results, "Female Dawn Genes (ZT0-6)")
male_go_plot <- plot_enrichment(male_results$go_results, "Male Active-Phase Genes (ZT12-18)")

# 7. Save results
save(female_results, male_results, file = "circadian_enrichment_results.RData")

if (!is.null(female_go_plot)) {
  ggsave("female_go_enrichment.pdf", female_go_plot, width = 10, height = 6)
}

if (!is.null(male_go_plot)) {
  ggsave("male_go_enrichment.pdf", male_go_plot, width = 10, height = 6)
}

# 8. Pseudogene miRNA interaction analysis
if (length(male_results$noncoding_genes) > 0) {
  # Get pseudogene symbols
  mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  pseudogenes <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = c("ensembl_gene_id", "biotype"),
    values = list(male_results$noncoding_genes, 
                  c("processed_pseudogene", "unprocessed_pseudogene")),
    mart = mart
  )
  
  if (nrow(pseudogenes) > 0) {
    # Run multiMiR analysis
    results <- get_multimir(
      org = "mmu",
      target = pseudogenes$external_gene_name,
      table = "all"
    )
    
    # Process results
    mirna_results <- as.data.frame(results@data)
    write.csv(mirna_results, "pseudogene_miRNA_interactions.csv", row.names = FALSE)
    
    # Filter significant results
    if ("support_type" %in% colnames(mirna_results)) {
      sig_results <- subset(mirna_results, 
                            support_type == "validated" |
                              (support_type == "predicted" & 
                                 as.numeric(score) >= 50))
    } else {
      sig_results <- head(mirna_results, 100)  # Fallback
    }
    
    write.csv(sig_results, "significant_miRNA_interactions.csv", row.names = FALSE)
  }
}








### Enhanced Specific Gene plot-------------------------------------------------
plot_gene <- function(gene_id) {
  plot_data <- data.frame(
    Expression = logCPM[gene_id, ],
    Time = metadata$Timepoint,
    Sex = metadata$Sex
  )
  
  ggplot(plot_data, aes(x=Time, y=Expression, color=Sex, group=Sex)) +
    geom_point(
      position = position_jitter(width = 0.3),
      alpha = 0.5,
      size = 3
    ) +
    geom_smooth(
      method = "loess", 
      se = FALSE, 
      color = "gray50",      # Medium gray line
      linewidth = 0.8,       # Thinner line
      linetype = "solid"     # Clean line style
    ) +
    scale_x_continuous(
      breaks = seq(0, 21, by = 3),
      labels = paste0("ZT", seq(0, 21, by = 3))
    ) +
    scale_color_manual(
      values = c("Male" = "#1f77b4", "Female" = "#ff7f0e")  # Maintain color scheme
    ) +
    labs(
      title = "Xist Rhythmic Gene Analysis in both Sexes",
      y = "Expression (logCPM)",
      x = "Zeitgeber Time"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),         # Remove all grid lines
      panel.background = element_blank(),   # Clean white background
      plot.background = element_blank(),    # No plot border
      legend.position = "top",              # Legend on top
      plot.title = element_text(hjust = 0.5, face = "bold"),  # Centered bold title
      axis.line = element_line(color = "black")  # Add axis lines for orientation
    )
}

# SAVE the plot
enhanced_plot <- plot_gene("ENSMUSG00000086503") 
ggsave("Xist_Rhythmic_Analysis.pdf", 
       plot = enhanced_plot,
       width = 8, 
       height = 6,
       bg = "white")  # White background for PDF




### Plotting Core Circadian Genes (Expected in both Sexes)----------------------
plot_gene <- function(gene_id, custom_title = NULL) {
  plot_data <- data.frame(
    Expression = logCPM[gene_id, ],
    Time = metadata$Timepoint,
    Sex = metadata$Sex
  )
  
  # Use custom title if provided, otherwise default
  plot_title <- ifelse(is.null(custom_title), 
                       paste("Rhythmic Analysis:", gene_id),
                       paste("Rhythmic Analysis:", custom_title))
  
  ggplot(plot_data, aes(x=Time, y=Expression, color=Sex, group=Sex)) +
    geom_point(
      position = position_jitter(width = 0.3),
      alpha = 0.5,
      size = 3
    ) +
    geom_smooth(
      method = "loess", 
      se = FALSE, 
      color = "gray50",
      linewidth = 0.8,
      linetype = "solid"
    ) +
    scale_x_continuous(
      breaks = seq(0, 21, by = 3),
      labels = paste0("ZT", seq(0, 21, by = 3))
    ) +
    scale_color_manual(
      values = c("Male" = "#1f77b4", "Female" = "#ff7f0e")
    ) +
    labs(
      title = plot_title,
      y = "Expression (logCPM)",
      x = "Zeitgeber Time"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.line = element_line(color = "black")
    )
}

# Gene list with custom titles
test_genes <- c(
  "ENSMUSG00000064842" = "Sry (Y-linked)",            # Male-biased
  "ENSMUSG00000086503" = "Xist (X-inactivation)",     # Female-biased
  "ENSMUSG00000031138" = "Androgen Receptor",         # Male-biased
  "ENSMUSG00000037432" = "Cyp17a1 (Steroidogenesis)", # Male-biased
  "ENSMUSG00000024406" = "Aromatase (Cyp19a1)"        # Female-biased
)

# Generate plots
for(i in seq_along(test_genes)) {
  gene_id <- names(test_genes)[i]
  if(gene_id %in% rownames(logCPM)) {
    plot_gene(gene_id, test_genes[i])  # Now passing both arguments correctly
    ggsave(paste0(gsub("/", "_", test_genes[i]), "_expression.pdf"), 
           width = 6, height = 4, bg = "white")
  } else {
    message(gene_id, " (", test_genes[i], ") not found")
  }
}


### PLOT lncRNAs----------------------------------------------------------------
# Define gene list with Ensembl IDs and symbols
lncRNAs <- c(
  "ENSMUSG00000040583" = "Meg3",
  "ENSMUSG00000037974" = "Snhg11",  
  "ENSMUSG00000092341" = "Malat1",
  "ENSMUSG00000055526" = "Kcnq1ot1",   # Not available
  "ENSMUSG00000089157" = "Ftx"  # Mouse FTX (X-linked)
)

# Verify which genes exist in your data
available_genes <- lncRNAs[names(lncRNAs) %in% rownames(logCPM)]
print(paste("Found", length(available_genes), "of", length(lncRNAs), "genes:"))
print(available_genes)


plot_lncRNA <- function(gene_id, gene_name) {
  plot_data <- data.frame(
    Expression = logCPM[gene_id, ],
    Time = metadata$Timepoint,
    Sex = metadata$Sex
  )
  
  ggplot(plot_data, aes(x=Time, y=Expression, color=Sex, group=Sex)) +
    geom_point(
      position = position_jitter(width = 0.3),
      alpha = 0.5,
      size = 3
    ) +
    geom_smooth(
      method = "loess", 
      se = FALSE, 
      color = "gray40",
      linewidth = 0.8
    ) +
    scale_x_continuous(
      breaks = seq(0, 21, by = 3),
      labels = paste0("ZT", seq(0, 21, by = 3))
    ) +
    scale_color_manual(
      values = c("Male" = "#3b7b9e", "Female" = "#e67c63")
    ) +
    labs(
      title = paste("Circadian Pattern of", gene_name, "\n", gene_id),
      y = "Expression (logCPM)",
      x = "Zeitgeber Time"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.line = element_line(color = "black")
    )
}

pdf("lncRNA_plots.pdf", width=7, height=5)  # Single PDF with multiple pages

for(i in seq_along(available_genes)) {
  gene_id <- names(available_genes)[i]
  gene_name <- available_genes[i]
  print(plot_lncRNA(gene_id, gene_name))
}

dev.off()  # Important - closes the PDF device



