
#===============================================================================
#         Volcano Plot: Female vs. Male across all timepoints
#===============================================================================

# Load packages
library(limma)
library(edgeR)

# Load data
counts <- read.delim("GSE297702_circadian_atlas_rawcounts.txt", row.names = 1)
metadata <- read.csv("wt_circadian_traits.csv")

# Clean factors
metadata$Sex <- factor(metadata$Sex, levels = c("Male", "Female"))
metadata$Timepoint <- factor(metadata$Timepoint, 
                             levels = c("0", "3", "6", "9", "12", "15", "18", "21"))

# Create DGEList and filter
dge <- DGEList(counts)
keep <- rowSums(cpm(dge) >= 1) >= 3
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

# Create design matrix
design <- model.matrix(~ Sex * Timepoint, data = metadata)
colnames(design) <- make.names(colnames(design))

# Limma/Voom transformation
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)

# Define contrasts
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

# Fit contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Test output
topTable(fit2, coef = "Female_vs_Male_T0", number = 5)


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
#                   Volcano Plot: Separate Male and Female
#===============================================================================

# Load your data
deg_data <- read.csv("/Users/judyabuel/Desktop/Xist/circadian_atlas/volcano_plots(1)/deg_results(1)/combined_DEG_results.csv")

# Define DEG thresholds (adjust as needed)
logFC_threshold <- 1    # Absolute log2 fold change
pval_threshold <- 0.05  # Adjusted p-value (FDR)

# Female DEGs (upregulated in Female = positive LogFC)
female_degs <- subset(deg_data, 
                      LogFC > logFC_threshold & FDR < pval_threshold)

# Male DEGs (upregulated in Male = negative LogFC)
male_degs <- subset(deg_data, 
                    LogFC < -logFC_threshold & FDR < pval_threshold)

# Non-DEGs (for background in volcano plot)
non_degs <- subset(deg_data, 
                   abs(LogFC) <= logFC_threshold | FDR >= pval_threshold)


write.csv(female_degs, "/Users/judyabuel/Desktop/Xist/circadian_atlas/Female_DEGs.csv", row.names = FALSE)
write.csv(male_degs, "/Users/judyabuel/Desktop/Xist/circadian_atlas/Male_DEGs.csv", row.names = FALSE)




#===============================================================================
#                       Female DEGs Volcano Plot
#===============================================================================


library(ggplot2)
library(dplyr)

# Create output directory
dir.create("female_vplots(1)", showWarnings = FALSE)

# Define gene lists
female_specific_genes <- c("Xist", "Eif2s3x", "Kdm6a", "Kdm5c")
female_associated <- c("Xist", "Malat1", "Meg3", "Neat1", "Kdm6a", "Eif2s3x", "Ftx", "Kcnq1ot1", "Snhg11")  # Updated list
male_specific_genes <- c("Uty", "Ddx3y", "Kdm5d", "Eif2s3y")  # Genes to remove

# Remove male-specific genes and classify remaining genes
deg_data <- deg_data %>%
  filter(!Gene_Symbol %in% male_specific_genes) %>%
  mutate(
    # First define gene_type
    gene_type = case_when(
      Gene_Symbol %in% female_specific_genes ~ "Female-specific",
      Gene_Symbol %in% female_associated ~ "Female-associated",
      FDR < 0.05 & abs(LogFC) > 1 ~ "Significant (both sexes)",
      TRUE ~ "Non-significant"
    ),
    # Then create to_label USING gene_type
    to_label = ifelse(
      Gene_Symbol == "Xist" |  # Force label Xist if present
        gene_type %in% c("Female-specific", "Female-associated", "Significant (both sexes)"),
      Gene_Symbol, 
      NA
    )
  )

# Volcano plot function
plot_volcano <- function(subset_data, timepoint) {
  ggplot(subset_data, aes(x = LogFC, y = -log10(PValue))) +
    # Non-significant points (no labels)
    geom_point(
      data = subset(subset_data, gene_type == "Non-significant"),
      color = "gray70", alpha = 0.5, size = 2
    ) +
    # Significant points
    geom_point(
      data = subset(subset_data, gene_type != "Non-significant"),
      aes(color = gene_type), alpha = 0.8, size = 2.5
    ) +
    scale_color_manual(
      values = c(
        "Female-specific" = "red",
        "Female-associated" = "magenta",  # Magenta for emphasized genes
        "Significant (both sexes)" = "green",
        "Non-significant" = "gray70"
      ),
      name = "Gene Label"
    ) +
    # Label all significant and female-associated genes
    geom_text_repel(
      aes(label = to_label),
      color = "black",
      size = 3,
      segment.color = NA,
      max.overlaps = Inf,
      na.rm = TRUE
    ) +
    theme_minimal() +
    labs(
      title = paste("Volcano Plot (", timepoint, ")", sep = ""),
      x = "Log2 Fold Change",
      y = "-log10(p-value)"
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.3) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
}

# Generate and save plots for each timepoint
for (timepoint in unique(deg_data$Contrast)) {
  subset_data <- deg_data %>% filter(Contrast == timepoint)
  p <- plot_volcano(subset_data, timepoint)
  
  ggsave(
    filename = paste0("female_vplots(1)/Volcano_", timepoint, ".pdf"),
    plot = p,
    width = 10,
    height = 7
  )
}











#===============================================================================
#                       Male DEGs Volcano Plot
#===============================================================================

library(ggplot2)
library(dplyr)
library(ggrepel)


# Define gene lists
male_specific_genes <- c("Uty", "Ddx3y", "Kdm5d", "Eif2s3y")
male_associated <- c("Ddx3y", "Uty", "Eif2s3y", "Kdm5d", "Ftx", "Kcnq1ot1", "Snhg11")
female_specific_genes <- c("Xist", "Eif2s3x", "Kdm6a", "Kdm5c")

# Process data
deg_data <- deg_data %>%
  filter(!Gene_Symbol %in% female_specific_genes) %>%
  mutate(
    gene_type = case_when(
      Gene_Symbol %in% male_specific_genes ~ "Male-specific",
      Gene_Symbol %in% male_associated ~ "Male-associated",
      FDR < 0.05 & LogFC < -1 ~ "Downregulated (male)",
      FDR < 0.05 & LogFC > 1 ~ "Upregulated (female)",
      TRUE ~ "Non-significant"
    ),
    to_label = ifelse(
      gene_type %in% c("Male-specific", "Male-associated", "Downregulated (male)"),
      Gene_Symbol, 
      NA
    )
  )

# Create output directory
dir.create("male_vplots", showWarnings = FALSE)

# Corrected volcano plot function
plot_volcano <- function(subset_data, timepoint) {
  ggplot(subset_data, aes(x = LogFC, y = -log10(FDR))) +  # Removed color from here
    geom_point(
      data = subset(subset_data, gene_type == "Non-significant"),
      color = "gray80", alpha = 0.6, size = 2.5
    ) +
    geom_point(
      data = subset(subset_data, gene_type != "Non-significant"),
      aes(color = gene_type), alpha = 0.8, size = 2.5
    ) +
    scale_color_manual(
      values = c(
        "Male-specific" = "blue",
        "Male-associated" = "darkorange",
        "Downregulated (male)" = "green",
        "Upregulated (female)" = "gray60"
      ),
      name = "Gene Label"
    ) +
    geom_text_repel(
      aes(label = to_label),
      color = "black",
      size = 3,
      segment.color = NA,
      max.overlaps = Inf,
      na.rm = TRUE
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.4) +
    labs(
      title = paste("Volcano Plot (", timepoint, ")"),
      x = "log2 Fold Change",
      y = "-log10(FDR)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "bottom"
    )
}

# Generate and save plots (now properly outside the function)
for (timepoint in unique(deg_data$Contrast)) {
  subset_data <- deg_data %>% filter(Contrast == timepoint)
  p <- plot_volcano(subset_data, timepoint)
  ggsave(
    filename = paste0("male_vplots/Volcano_", timepoint, ".pdf"),
    plot = p,
    width = 10,
    height = 8
  )
  print(paste("Saved plot for", timepoint))  # Verification message
}









