
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


#===============================================================================
#               Generate Volcano Plots for Each Timepoint
#===============================================================================

# Create output directories with consistent naming
dir.create("volcano_plots", showWarnings = FALSE)
dir.create("deg_results", showWarnings = FALSE)

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
    
    # Define genes to highlight: Brown4_modules gene network interacting with imprinted & xist
    genes_of_interest <- c("kdm6a", "Eif2s3y", "Xist", "Uty", "Ddx3y", "Kdm5d", "Meg3", "Kcnq1ot1", "Snhg11")
    
    # Inside your volcano plotting loop (replace the sig_genes block with this):
    highlight_genes <- top_table_annotated %>%
      dplyr::filter(label %in% genes_of_interest)
    
    # Volcano plot
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
        #subtitle = "Candidate Genes Highlighted",
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
        data = highlight_genes,
        aes(label = label),
        color ="black",
        size = 3.0,
        max.overlaps = 30,
        min.segment.length = 0.2,
        box.padding = 0.4,
        force = 1,
        show.legend = FALSE
      ) +
      coord_cartesian(xlim = c(-max(abs(top_table_annotated$logFC))-0.5, 
                               max(abs(top_table_annotated$logFC))+0.5))
    
    
    # Save plot in multiple formats
    plot_name <- paste0("volcano_plots/", make.names(contrast_name))
    ggsave(paste0(plot_name, ".pdf"), volcano_plot, width = 9, height = 7)
    #ggsave(paste0(plot_name, ".png"), volcano_plot, width = 9, height = 7, dpi = 300)
    
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

readr::write_csv(combined_deg, "deg_results/combined_DEG_results.csv")



# Count DEGs per timepoint
deg_summary <- combined_deg %>%
  dplyr::filter(FDR < 0.05 & abs(LogFC) > 1) %>%   # same cutoff as volcano
  dplyr::group_by(Contrast, Expression) %>%
  dplyr::summarize(DEG_count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Expression, values_from = DEG_count, values_fill = 0) %>%
  dplyr::mutate(Total_DEGs = Upregulated + Downregulated)

deg_summary

dir.create("deg_results", showWarnings = FALSE)
readr::write_csv(deg_summary, "deg_results/DEG_counts_per_timepoint.csv")



### What are the specific DEGs in Males and Females or Both.
### Our designn is Female vs Male (logFC = Female / Male)
### logFC > 0 --> higher in females
### logFC < 0 --> higher in Males

# Let's filter the Significant DEGs
sig_deg <- combined_deg %>%
  dplyr::filter(FDR < 0.05 & abs(LogFC) > 1) %>%
  dplyr::mutate(Regulation = ifelse(LogFC > 0, "Female_up", "Male_up"))


# Check unique vs shared per timepoint
library(purrr)

unique_deg <- sig_deg %>%
  split(.$Contrast) %>%
  map(~ {
    female_genes <- .x %>% filter(Regulation == "Female_up") %>% pull(Gene_Symbol)
    male_genes   <- .x %>% filter(Regulation == "Male_up") %>% pull(Gene_Symbol)
    
    tibble(
      Timepoint   = unique(.x$Contrast),
      Female_only = list(setdiff(female_genes, male_genes)),
      Male_only   = list(setdiff(male_genes, female_genes)),
      Shared      = list(intersect(female_genes, male_genes))
    )
  }) %>%
  bind_rows()

unique_deg$Female_only[[1]]   # all female-only genes at ZT0

# Summarize counts
deg_counts <- sig_deg %>%
  group_by(Contrast, Regulation) %>%
  summarise(DEG_count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Regulation, values_from = DEG_count, values_fill = 0)

deg_counts <- deg_counts %>%
  mutate(Total = Female_up + Male_up)

# Write unique/shared gene lists
saveRDS(unique_deg, "deg_results/unique_DEGs_by_timepoint.rds")

# Save summary counts
readr::write_csv(deg_counts, "deg_results/DEG_counts_unique.csv")






#--------------------------------------------------------------------------------
### Let's make a Venn Diagram per timepoint comparing Female_up and Male_up DEGs
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}

library(VennDiagram)

# Create output folder
dir.create("venn_plots", showWarnings = FALSE)

# Loop over each timepoint
for(tp in unique(sig_deg$Contrast)) {
  
  sub <- sig_deg %>% filter(Contrast == tp)
  
  female_genes <- sub %>% filter(Regulation == "Female_up") %>% pull(Gene_Symbol)
  male_genes   <- sub %>% filter(Regulation == "Male_up") %>% pull(Gene_Symbol)
  
  # Save as image
  # png(filename = paste0("venn_plots/Venn_", tp, ".png"), width = 1200, height = 1000, res = 150)
  pdf("Venn_ZT0.pdf", width = 6, height = 6
      
  draw.pairwise.venn(
    area1 = length(female_genes),
    area2 = length(male_genes),
    cross.area = length(intersect(female_genes, male_genes)),
    category = c("Female_DEGs", "Male_DEGs"),
    fill = c("#E41A1C", "#377EB8"),
    alpha = 0.5,
    cex = 1.5,
    cat.cex = 1.5,
    cat.pos = c(-20, 20),
    cat.dist = 0.05,
    main = paste("ZT", tp, "- Female vs Male Upregulated DEGs")
  )
  
  dev.off()
}





# Combine multiple VENN DIAGRAMS in one figure =================================
if (!requireNamespace("ggVennDiagram", quietly = TRUE)) install.packages("ggVennDiagram")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")

library(ggVennDiagram)
library(dplyr)
library(purrr)
library(patchwork)

# Set the order of timepoints
timepoint_order <- c("ZT0","ZT3","ZT6","ZT9","ZT12","ZT15","ZT18","ZT21")

# Split DEGs per timepoint into Female_up vs Male_up
venn_list <- sig_deg %>%
  split(.$Contrast) %>%
  map(~ list(
    Female = .x %>% filter(Regulation == "Female_up") %>% pull(Gene_Symbol),
    Male   = .x %>% filter(Regulation == "Male_up") %>% pull(Gene_Symbol)
  ))

venn_plots <- lapply(timepoint_order, function(tp) {
  if(!tp %in% names(venn_list)) return(NULL)
  
  ggVennDiagram(
    venn_list[[tp]],
    label_alpha = 0,                 # hide numbers inside circles
    edge_size = 0.5,                 # keep subtle edges
    set_size = 0,                    # show set labels with smaller size
    category.names = c("", "")
  ) +
    ggtitle(tp) +
    scale_fill_gradient(low = "white", high = "lightblue", guide = "none") +
    #scale_color_manual(values = c("#FFB6C1", "#ADD8E6")) +  # Circle border colors
    theme(
      plot.title = element_text(hjust=0.5, size=14, face="bold"),
      legend.position = "none"
    )
})

# Remove NULLs
venn_plots <- venn_plots[!sapply(venn_plots, is.null)]

# Combine plots
combined_venn <- wrap_plots(venn_plots) +
  plot_annotation(
    title = "Venn Diagrams of Female vs Male Upregulated DEGs per Timepoint",
    theme = theme(plot.title = element_text(size=16, face="bold", hjust=0.5))
  )

ggsave("1venn_plots_all_timepoints.pdf", combined_venn, width=12, height=8)






# working on this! -------------------------------------------------------------
venn_plots <- lapply(timepoint_order, function(tp) {
  if(!tp %in% names(venn_list)) return(NULL)
  
  ggVennDiagram(
    venn_list[[tp]],
    label_alpha = 0,                 # hide numbers inside circles
    edge_size = 0.5,                 # keep subtle edges
    set_size = 0,                    # show set labels with smaller size
    category.names = c("", ""),
    set_color = c("#87CEEB", "#FFB6C1")  # Top (Male)=Blue, Bottom (Female)=Pink
  ) +
    ggtitle(tp) +
    scale_fill_gradient(low = "white", high = "lightblue", guide = "none") +
    theme(
      plot.title = element_text(hjust=0.5, size=14, face="bold"),
      legend.position = "none"
    )
})

# Remove NULLs
venn_plots <- venn_plots[!sapply(venn_plots, is.null)]

# Combine plots
combined_venn <- wrap_plots(venn_plots) +
  plot_annotation(
    title = "Venn Diagrams of Female vs Male Upregulated DEGs per Timepoint",
    theme = theme(plot.title = element_text(size=16, face="bold", hjust=0.5))
  )

ggsave("1venn_plots_all_timepoints.pdf", combined_venn, width=12, height=8)















#-------------------------------------------------------------------------------
# PLOT individual DEGs per group (Male_up, Female_up, Shared) at each timepoint
# Filter only significant DEGs and keep needed columns
# This is to create Box plots

library(ggplot2)
library(dplyr)
library(patchwork)

# Calculate counts for each group at each timepoint
deg_counts <- sig_deg %>%
  filter(Regulation %in% c("Female_up", "Male_up")) %>%
  mutate(Group = case_when(
    Regulation == "Female_up" ~ "Female",
    Regulation == "Male_up" ~ "Male"
  )) %>%
  group_by(Contrast, Group) %>%
  summarise(Count = n(), .groups = "drop")

# Calculate shared DEGs (intersection of Female and Male upregulated)
shared_counts <- map_df(timepoint_order, function(tp) {
  female_genes <- sig_deg %>% 
    filter(Contrast == tp, Regulation == "Female_up") %>% 
    pull(Gene_Symbol)
  
  male_genes <- sig_deg %>% 
    filter(Contrast == tp, Regulation == "Male_up") %>% 
    pull(Gene_Symbol)
  
  shared_genes <- intersect(female_genes, male_genes)
  
  data.frame(
    Contrast = tp,
    Group = "Shared",
    Count = length(shared_genes)
  )
})

# Combine both datasets
all_counts <- bind_rows(deg_counts, shared_counts) %>%
  mutate(Contrast = factor(Contrast, levels = timepoint_order),
         Group = factor(Group, levels = c("Female", "Male", "Shared")))

# Create the bar plot
ggplot(all_counts, aes(x = Group, y = Count, fill = Group)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = Count), vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Female" = "#FFB6C1", "Male" = "#ADD8E6", "Shared" = "#D8BFD8")) +
  facet_wrap(~ Contrast, ncol = 4) +
  labs(title = "DEG Counts by Group at Each Timepoint", 
       x = NULL, y = "Number of DEGs") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    #panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Save the plot
ggsave("1deg_counts_by_group_timepoints.pdf", width = 12, height = 8)

