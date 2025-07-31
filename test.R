#===============================================================================
# PURPOSE: DEG Analysis on Sex Differences
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
}

#===============================================================================
#               Generate Volcano Plots for Each Timepoint
#===============================================================================

# Create output directories with consistent naming
dir.create("volcano_plots(1)", showWarnings = FALSE)
dir.create("deg_results(1)", showWarnings = FALSE)

# Initialize list to store all DEG results
all_deg_results <- list()

# Define color scheme
deg_colors <- c("Upregulated" = "#E41A1C",  # Red
                "Downregulated" = "#377EB8", # Blue
                "Not Significant" = "grey80") # Light grey

# Process all contrasts
for (contrast_name in colnames(contrast.matrix)) {
  tryCatch({
    message("\nProcessing ", contrast_name, "...")
    
    # Get DEG results with robust settings
    top_table <- topTable(fit2, coef = contrast_name, number = Inf, 
                          sort.by = "p", adjust.method = "BH")
    
    # Skip if no results
    if (nrow(top_table) == 0) {
      warning("No DEGs found for contrast: ", contrast_name)
      next
    }
    
    # Extract and annotate genes
    top_table$ensembl_gene_id <- sub("\\..*", "", rownames(top_table))
    top_table_annotated <- merge(top_table, gene_map, 
                                 by = "ensembl_gene_id", 
                                 all.x = TRUE) %>%
      tibble::as_tibble()
    
    # Create labels - prioritize MGI symbol
    top_table_annotated <- top_table_annotated %>%
      dplyr::mutate(
        label = dplyr::coalesce(mgi_symbol, ensembl_gene_id),
        Timepoint = stringr::str_extract(contrast_name, "T\\d+$") %>% 
          stringr::str_replace("T", "ZT"),
        Expression = dplyr::case_when(
          adj.P.Val < 0.05 & logFC > 1  ~ "Upregulated",
          adj.P.Val < 0.05 & logFC < -1 ~ "Downregulated",
          TRUE ~ "Not Significant"
        ),
        Significance = ifelse(adj.P.Val < 0.05, "FDR < 0.05", "Not Sig")
      )
    
    # Store results
    all_deg_results[[contrast_name]] <- top_table_annotated
    
    # Get significant genes for labeling
    sig_genes <- top_table_annotated %>%
      dplyr::filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>%
      dplyr::arrange(adj.P.Val) %>%
      dplyr::distinct(label, .keep_all = TRUE) %>%
      dplyr::group_by(sign(logFC)) %>%  # Balance up/down regulated labels
      dplyr::slice_head(n = 15)         # Top 15 per direction
    
    # Create volcano plot
    volcano_plot <- ggplot(top_table_annotated, 
                           aes(x = logFC, 
                               y = -log10(adj.P.Val),
                               color = Expression,
                               alpha = Significance)) +
      geom_point(size = 1.8) +
      scale_color_manual(values = deg_colors) +
      scale_alpha_manual(values = c("FDR < 0.05" = 0.8, "Not Sig" = 0.4)) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
      labs(
        title = paste("Female vs Male at", unique(top_table_annotated$Timepoint)),
        subtitle = "Differentially Expressed Genes (logFC > 1, FDR < 0.05)",
        x = expression(Log[2]~"Fold Change (Female/Male)"),
        y = expression(-Log[10]~"Adjusted P-value"),
        color = "Expression"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
      ) +
      geom_text_repel(
        data = sig_genes,
        aes(label = label),
        size = 2.8,
        max.overlaps = 30,
        min.segment.length = 0.2,
        box.padding = 0.4,
        force = 1,
        show.legend = FALSE
      ) +
      coord_cartesian(xlim = c(-max(abs(top_table_annotated$logFC))-0.5, 
                               max(abs(top_table_annotated$logFC))+0.5))
    
    # Save plot in multiple formats
    plot_name <- paste0("volcano_plots(1)/", make.names(contrast_name))
    ggsave(paste0(plot_name, ".pdf"), volcano_plot, width = 9, height = 7)
    ggsave(paste0(plot_name, ".png"), volcano_plot, width = 9, height = 7, dpi = 300)
    
    message("Saved plots for ", contrast_name)
    
  }, error = function(e) {
    message("Error processing ", contrast_name, ": ", e$message)
  })
}

# Combine and save all results with better formatting
combined_deg <- dplyr::bind_rows(all_deg_results) %>%
  dplyr::select(Contrast = Timepoint,
                Ensembl_ID = ensembl_gene_id,
                Gene_Symbol = mgi_symbol,
                LogFC = logFC,
                PValue = P.Value,
                FDR = adj.P.Val,
                Expression = Expression,
                AveExpr,
                dplyr::everything()) %>%
  dplyr::arrange(Contrast, FDR)

readr::write_csv(combined_deg, "deg_results(1)/combined_DEG_results.csv")





#===============================================================================
#                     Generate Heatmap across all Timepoint
#===============================================================================

library(pheatmap)
library(dplyr)
library(RColorBrewer)

install.packages("RColorBrewer")

# Load DEG CSV and gene_map
deg_all <- read.csv("DEG_Female_vs_Male_by_Timepoint.csv")
deg_all <- deg_all %>%
  mutate(ensembl_gene_id = sub(".*\\.", "", X))  # Extract Ensembl ID

# Join with gene_map to add gene symbols
deg_annot <- left_join(deg_all, gene_map, by = "ensembl_gene_id")

# Use gene symbol for labels; fall back to Ensembl if missing
deg_annot <- deg_annot %>%
  mutate(label = ifelse(is.na(mgi_symbol) | mgi_symbol == "", ensembl_gene_id, mgi_symbol))

# Get top 15 DEGs per timepoint (adj.P.Val < 0.05)
deg_top <- deg_annot %>%
  filter(adj.P.Val < 0.05) %>%
  group_by(Timepoint) %>%
  arrange(adj.P.Val, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

# Get unique top genes
top_genes <- unique(deg_top$ensembl_gene_id)

# Subset voom-normalized expression matrix
expr_top <- v$E[rownames(v$E) %in% top_genes, ]

# Match rownames to gene symbols
gene_labels <- deg_top %>%
  dplyr::select(ensembl_gene_id, label) %>%
  distinct()


rownames(expr_top) <- gene_labels$label[match(rownames(expr_top), gene_labels$ensembl_gene_id)]

# Optional: scale per gene (z-score)
expr_scaled <- t(scale(t(expr_top)))

# Annotation for columns (samples)
ann_colors <- list(Sex = c(Female = "firebrick", Male = "steelblue"))
annot_col <- metadata %>%
  dplyr::select(Sample_ID, Sex, Timepoint) %>%
  filter(Sample_ID %in% colnames(expr_scaled)) %>%
  tibble::column_to_rownames("Sample_ID")


# Save heatmap
output_dir <- "heatmaps"
dir.create(output_dir, showWarnings = FALSE)

heatmap_file <- file.path(output_dir, "Heatmap_Top10_DEGs_AllTimepoints.png")
print(paste("Saving heatmap to:", heatmap_file))

# Plot heatmap
png("heatmaps/Heatmap_DEGs_AllTimepoints.png", width = 2200, height = 1600, res = 300)
pheatmap(
  expr_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annot_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("navy", "white", "firebrick"))(100),
  fontsize_row = 6,
  fontsize_col = 6,
  angle_col = 45,
  main = "Top DEGs at Each Timepoint: Female vs Male",
  fontsize = 10
)

dev.off()


#===============================================================================
#              Line graph of the Top 10 genes across timepoints
#===============================================================================

# Load required packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggforce)  # For pagination

# Read data and select genes
data <- read.csv("Heatmap_Top_DEGs_logFC_matrix.csv", row.names = 1)
specific_genes <- c("Uty", "Ddx3y", "Kdm5d", "Eif2s3y", "Xist", "Eif2s3x", "Kdm6a", "Glp2r")  # Your genes

# Reshape data
data_long <- data %>%
  filter(rownames(.) %in% specific_genes) %>%
  mutate(Gene = factor(rownames(.))) %>%  # Ensure proper ordering
  pivot_longer(cols = -Gene, names_to = "Timepoint", values_to = "logFC")

# Order timepoints
data_long$Timepoint <- factor(data_long$Timepoint, levels = paste0("ZT", sprintf("%02d", seq(0, 21, 3))))

# Create a function to plot in batches
plot_genes_paginated <- function(data, ncol = 3, nrow = 2, per_page = ncol * nrow) {
  n_pages <- ceiling(length(unique(data$Gene)) / per_page)
  
  for (i in 1:n_pages) {
    p <- ggplot(data, aes(x = Timepoint, y = logFC, group = Gene, color = Gene)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      ggforce::facet_wrap_paginate(~ Gene, ncol = ncol, nrow = nrow, page = i) +
      scale_y_continuous(limits = c(-10, 15), breaks = seq(-10, 15, by = 5)) +
      labs(title = paste("Changes in Gene Expression per Timepoint"), x = "Timepoint", y = "log2 Fold Change") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            strip.background = element_rect(fill = "lightgrey"))
    
    ggsave(
      filename = paste0("Selected_Genes_Page_", i, ".pdf"),
      plot = p,
      device = "pdf",
      width = 10,
      height = 7,
      dpi = 300
    )
  }
}

# Run the function (adjust nrow/ncol as needed)
plot_genes_paginated(data_long, ncol = 3, nrow = 2)  # 6 genes per page (3 columns x 2 rows)



#===============================================================================
#                         Generate UpSet Plot
#===============================================================================

# Load libraries
install.packages("ComplexUpset")
install.packages("ggplot2")  # needed for plotting

library(ComplexUpset)
library(ggplot2)


# Define timepoints
all_tp <- c("ZT00", "ZT03", "ZT06", "ZT09", "ZT12", "ZT15", "ZT18", "ZT21")

# Basic UpSet plot
upset_plot <- upset(
  deg_df,
  intersect = all_tp,
  name = "DEG Timepoints",
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(size = 3)
    )
  ),
  set_sizes = upset_set_size(),
  sort_intersections_by = "cardinality",
  n_intersections = 15,
  min_size = 1,
  width_ratio = 0.3,
  themes = upset_modify_themes(list(
    'intersect_size' = theme(text = element_text(size = 10)),
    'overall_sizes' = theme(axis.text.y = element_text(size = 10)),
    'intersections_matrix' = theme(axis.text.x = element_text(size = 8, angle = 90))
  ))
) +
  labs(
    title = "Top DEG Intersections Across Circadian Timepoints",
    subtitle = NULL
  )

# Save to PNG
ggsave("UpSet_DEGs_Top15_Basic.png", plot = upset_plot, width = 10, height = 6, dpi = 300)








#===============================================================================
#                 Generate MDS plot Sex:Timepoint
#===============================================================================

mds <- plotMDS(d, 
               col = ifelse(metadata$Sex == "Male", "blue", "pink"),  # Colors by sex
               pch = ifelse(metadata$Timepoint == 0, 16,               # Shapes by timepoint
                            ifelse(metadata$Timepoint == 3, 17,
                                   ifelse(metadata$Timepoint == 6, 15, 18))),
               cex = 1.5,  # Point size
               main = "MDS Plot: Sex and Timepoint",
               xlab = "Dimension 1",
               ylab = "Dimension 2")

# Add legend for Sex
legend("topright", 
       legend = c("Male", "Female"), 
       col = c("blue", "pink"), 
       pch = 16, 
       title = "Sex")

# Add legend for Timepoint
legend("bottomright", 
       legend = metadata$Timepoint, 
       #pch = c(16, 17, 15, 18, 8, 9, 10, 12), 
       title = "Timepoint")


#-------MDS PLOT between Sexes-------------
# Convert MDS coordinates to data frame
mds_data <- data.frame(
  Dim1 = mds$x,
  Dim2 = mds$y,
  Sex = metadata$Sex
)

# Plot
ggplot(mds_data, aes(Dim1, Dim2, color = Sex)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Male" = "blue", "Female" = "pink")) +
  labs(x = "Dimension 1", y = "Dimension 2", title = "MDS Plot: Male vs. Female") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove gray background
    axis.line = element_line(color = "black")  # Keep axis lines
  )

