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
library(org.Mm.eg.db)

# Set working directory
setwd("/Users/judyabuel/Desktop/Xist/circadian_atlas")

# Load Data ====================================================================
# Load count matrix
counts <- read.delim("GSE297702_circadian_atlas_rawcounts.txt", row.names = 1)

# Load metadata - this contains sex and timepoints
metadata <- read.csv("wt_circadian_traits.csv")
metadata$Sex <- factor(metadata$Sex)
metadata$Timepoint <- factor(metadata$Timepoint)

# Re-level to set Male and ZT0 as references
metadata$Sex <- relevel(metadata$Sex, ref = "Male")
metadata$Timepoint <- relevel(metadata$Timepoint, ref = "0")
metadata$Timepoint <- factor(as.character(metadata$Timepoint))

# Pre-processing ===============================================================
# Create DGEList object and normalize
dge <- DGEList(counts)
dge <- calcNormFactors(dge, method = "TMM")

# Filter low-expressed genes
keep <- rowSums(cpm(dge) >= 1) >= 3
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Create design matrix
design <- model.matrix(~ Sex * Timepoint, data = metadata)
colnames(design) <- make.names(colnames(design))  # Clean column names

# Apply voom transformation and fit model
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
topTable(fit)

# Get the voom-transformed counts (log2-CPM with weights)
voom_log2_CPM <- v$E
logCPM <- cpm(dge, log = TRUE)

# Define Contrasts: Female vs. Male across timepoints ==========================
contrast.matrix <- makeContrasts(
  Female_vs_Male_T0 = SexFemale,
  Female_vs_Male_T6 = SexFemale + SexFemale.Timepoint6,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Annotate Genes with MGI Symbols ==============================================
annotations <- select(org.Mm.eg.db,
                      keys = keys(org.Mm.eg.db),
                      columns = c("ENSEMBL", "SYMBOL", "GENENAME"))

# Get clean Ensembl IDs from DGE object
ensembl_ids <- rownames(dge)
ensembl_ids_clean <- sub("\\..*", "", ensembl_ids)

# Connect to Ensembl (safer call without mirror)
ensembl <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids_clean,
  mart = ensembl
)

# Volcano Plot =================================================================
all_deg_results <- list()

# Define color scheme
deg_colors <- c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")

# Loop through contrasts and generate volcano plots
for (contrast_name in colnames(contrast.matrix)) {
  tryCatch({
    
    message("Processing ", contrast_name, "...")
    
    # Fit model for this contrast
    fit2 <- contrasts.fit(fit, contrast.matrix[, contrast_name])
    fit2 <- eBayes(fit2)
    top_table <- topTable(fit2, number = Inf, sort.by = "P")
    
    # Annotate with gene symbols
    top_table_annotated <- top_table %>%
      rownames_to_column("ensembl_gene_id") %>%
      left_join(gene_map, by = "ensembl_gene_id")
    
    
    # Add DE labels
    top_table_annotated <- top_table_annotated %>%
      mutate(
        Expression = case_when(
          logFC > 1 & adj.P.Val < 0.05 ~ "Up",
          logFC < -1 & adj.P.Val < 0.05 ~ "Down",
          TRUE ~ "Not Sig"
        ),
        Significance = ifelse(adj.P.Val < 0.05, "FDR < 0.05", "Not Sig")
      )
    
    # Pick genes of interest present in results
    highlight_genes <- top_table_annotated %>%
      filter(mgi_symbol %in% genes_of_interest)
    
    # Volcano plot
    volcano_plot <- ggplot(top_table_annotated,
                           aes(x = logFC, y = -log10(adj.P.Val),
                               color = Expression, alpha = Significance)) +
      geom_point(size = 1.8) +
      scale_color_manual(values = deg_colors) +
      scale_alpha_manual(values = c("FDR < 0.05" = 0.8, "Not Sig" = 0.4)) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
      labs(
        title = paste("vplot of brown4 module at ZT6", unique(top_table_annotated$Timepoint)),
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
        aes(label = mgi_symbol),
        color = "black",
        size = 3.0,
        max.overlaps = Inf,
        box.padding = 0.4,
        force = 1,
        show.legend = FALSE
      ) +
      coord_cartesian(
        xlim = c(-max(abs(top_table_annotated$logFC)) - 0.5,
                 max(abs(top_table_annotated$logFC)) + 0.5)
      )
    
    # Save plot directly in current working directory
    plot_name <- paste0(make.names(contrast_name), ".pdf")
    ggsave(plot_name, volcano_plot, width = 9, height = 7)
    

# ------------------------------------------------------------------------------  
    # Check your gene_map columns
    colnames(gene_map)
    # Usually it has: "ensembl_gene_id", "mgi_symbol"
    
    # Your genes of interest
    #genes_of_interest <- c("Meg3", "Rnpc3", "Pnisr", "Leng8", "Snhg11", "Tia1",
                           #"Mir100hg", "4632427E13Rik", "Spaca6", "Ftx",
                           #"AI480526", "Sec14l5", "Neat1", "Malat1",
                           #"A330023F24Rik", "Gm3764", "Mir124a-1hg", "Miat",
                           #"B830012L14Rik", "Kcnq1ot1", "Gm37899",
                           #"C130023A14Rik", "A330076H08Rik")

  gene_map_of_interest <- c("ENSMUSG00000021268", "ENSMUSG00000027981", "ENSMUSG00000028248",
                            "ENSMUSG00000035545", "ENSMUSG00000044349", "ENSMUSG00000071337",
                            "ENSMUSG00000074024", "ENSMUSG00000074415", "ENSMUSG00000080316",
                            "ENSMUSG00000086370", "ENSMUSG00000090086", "ENSMUSG00000091712",
                            "ENSMUSG00000092274", "ENSMUSG00000092341", "ENSMUSG00000096929",
                            "ENSMUSG00000097156", "ENSMUSG00000097545", "ENSMUSG00000097767",
                            "ENSMUSG00000098202", "ENSMUSG00000101609", "ENSMUSG00000102657",
                            "ENSMUSG00000102854", "ENSMUSG00000109321")  
    
    
        
    # Filter gene_map for these genes
    gene_map_of_interest <- gene_map %>%
      filter(mgi_symbol %in% genes_of_interest)
    
    gene_map_of_interest
    

    
    
    
    
    
    
    
    
    
#-------------------------------------------------------------------------------
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(org.Mm.eg.db)
    
    # 1. Map gene symbols to Ensembl IDs (already done)
    # gene_map_present contains only the 23 genes found
    
    # 2. Subset expression matrix for present genes
    expr_ids_clean <- sub("\\..*", "", rownames(voom_log2_CPM))
    present_ensembl <- gene_map_present$ensembl_gene_id
    expr_subset <- voom_log2_CPM[present_ensembl, , drop = FALSE]
    rownames(expr_subset) <- gene_map_present$mgi_symbol  # rename to symbols
    
    # 3. Merge with metadata
    expr_df <- as.data.frame(t(expr_subset))
    expr_df$SampleID <- rownames(expr_df)
    metadata$SampleID <- rownames(metadata)
    expr_merged <- merge(metadata, expr_df, by = "SampleID")
    
    # 4. Pivot longer for ggplot
    expr_long <- expr_merged %>%
      pivot_longer(
        cols = gene_map_present$mgi_symbol,   # only present genes
        names_to = "Gene",
        values_to = "Expression"
      )
    
    nrow(expr_long)  # should now be >0
    
    
    # 5. Plot rhythmic expression
    ggplot(expr_long, aes(x = as.numeric(as.character(Timepoint)), y = Expression, color = Sex)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_line(aes(group = interaction(Sex, Gene)), alpha = 0.7) +
      facet_wrap(~ Gene, scales = "free_y") +
      scale_x_continuous(breaks = seq(0, 21, by = 3), labels = paste0("ZT", seq(0, 21, by = 3))) +
      labs(
        x = "Zeitgeber Time (ZT)",
        y = "Log2-CPM Expression",
        color = "Sex",
        title = "Rhythmic expression of genes of interest"
      ) +
      theme_minimal(base_size = 12)

    
    
    
    
    
#-------------------------------------------------------------------------------

    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(org.Mm.eg.db)
    
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
    # 3. Create the plot data
    # -----------------------------
    # Ensure metadata matches expression samples
    stopifnot(all(colnames(voom_log2_CPM_clean) %in% rownames(metadata)))
    
    # Force metadata rows to use the same IDs as expression matrix
    metadata$SampleID <- colnames(voom_log2_CPM_clean)
    
    # If metadata was originally in the same order as expression columns:
    rownames(metadata) <- metadata$SampleID
    
    # Assign sex from prefix
    metadata$Sex <- ifelse(grepl("^A", metadata$SampleID), "Male", "Female")
    metadata$Sex <- factor(metadata$Sex, levels = c("Male", "Female"))
    
    
    plot_data <- lapply(present_symbols, function(gene) {
      ens_id <- present_ensembl[gene]
      expr_values <- voom_log2_CPM_clean[ens_id, ]
      data.frame(
        SampleID = names(expr_values),
        Gene = gene,
        Expression = as.numeric(expr_values)
      )
    }) %>% bind_rows()
    
    plot_data <- merge(plot_data, metadata, by = "SampleID")
    
    # -----------------------------
    # 4. Plot with individual points + mean lines
    # -----------------------------
    mean_data <- plot_data %>%
      group_by(Gene, Timepoint, Sex) %>%
      summarize(Mean_Expression = mean(Expression, na.rm = TRUE), .groups = "drop")
    
    ggplot(plot_data, aes(x = Timepoint, y = Expression, color = Sex)) +
      geom_point(alpha = 0.3, size = 2) +
      geom_line(data = mean_data, 
                aes(x = Timepoint, y = Mean_Expression, group = Sex),
                size = 1, alpha = 0.8) +
      facet_wrap(~ Gene, scales = "free_y") +
      scale_x_continuous(breaks = seq(0, 21, by = 3), 
                         labels = paste0("ZT", seq(0, 21, by = 3))) +
      scale_color_manual(values = c("Male" = "blue", "Female" = "red")) +
      labs(
        x = "Zeitgeber Time (ZT)",
        y = "Log2-CPM Expression",
        color = "Sex",
        title = "Rhythmic expression of genes of interest"
      ) +
      theme_minimal(base_size = 12)
    
    
    
    

 # Plot eaach genes individually------------------------------------------------
    
    for (g in unique(plot_data$Gene)) {
      gene_data <- plot_data %>% filter(Gene == g)
      
      p <- ggplot(gene_data, aes(x = Timepoint, y = Expression, color = Sex)) +
        geom_point(alpha = 0.3, size = 2) +
        geom_smooth(aes(group = Sex), method = "loess", se = TRUE, size = 1, alpha = 0.2) +
        scale_color_manual(values = c("Male" = "blue", "Female" = "red")) +
        scale_fill_manual(values = c("Male" = "blue", "Female" = "red")) +
        scale_x_continuous(breaks = seq(0, 21, by = 3), 
                           labels = paste0("ZT", seq(0, 21, by = 3))) +
        labs(
          x = "Zeitgeber Time (ZT)",
          y = "Log2-CPM Expression",
          color = "Sex",
          title = paste("Rhythmic expression of", g)
        ) +
        theme_minimal(base_size = 12) +
        theme(panel.grid = element_blank(),
              axis.line = element_line(color = "black"))
      
      ggsave(filename = paste0("Gene_", g, ".pdf"), plot = p, width = 6, height = 4)
    }

    
    
    
    

    
###=============================================================================    
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    BiocManager::install("rain")
    
    library(rain)
    library(ggplot2)
    library(dplyr)
    library(biomaRt)
    library(tidyr)
    

    expr_mat <- as.data.frame(t(voom_log2_CPM))
    expr_mat$Sample_ID <- rownames(expr_mat)
    

    ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    gene_map <- getBM(
      attributes = c("ensembl_gene_id", "mgi_symbol"),
      filters = "ensembl_gene_id",
      values = sub("\\..*", "", rownames(voom_log2_CPM)),  # strip version if present
      mart = ensembl
    )
    
    colnames(gene_map) <- c("Ensembl_ID", "Gene")
    

    metadata$Sample_ID <- rownames(metadata)
    
    expr_merged <- expr_mat %>%
      left_join(metadata, by = "Sample_ID")
    

    genes_of_interest <- c("ENSMUSG00000027981", "ENSMUSG00000028248", "ENSMUSG00000035545",
                           "ENSMUSG00000071337", "ENSMUSG00000074024", "ENSMUSG00000074415",
                           "ENSMUSG00000080316", "ENSMUSG00000086370", "ENSMUSG00000090086",
                           "ENSMUSG00000091712", "ENSMUSG00000092274", "ENSMUSG00000092341",
                           "ENSMUSG00000096929", "ENSMUSG00000097156", "ENSMUSG00000097545",
                           "ENSMUSG00000097767", "ENSMUSG00000098202", "ENSMUSG00000102657",
                           "ENSMUSG00000102854", "ENSMUSG00000109321")
    

    plot_data <- voom_log2_CPM %>%
      as.data.frame() %>%
      rownames_to_column("Ensembl_ID") %>%
      filter(sub("\\..*", "", Ensembl_ID) %in% genes_of_interest) %>%
      pivot_longer(-Ensembl_ID, names_to = "Sample_ID", values_to = "Expression") %>%
      left_join(metadata, by = "Sample_ID") %>%
      left_join(gene_map, by = c("Ensembl_ID"))  # now you have both Ensembl + symbol
    
    head(plot_data)
    
    
     # Assuming plot_data has columns: Gene, Timepoint, Expression, Sex
    
    results_list <- list()
    
    for (g in unique(plot_data$Gene)) {
      gene_data <- plot_data %>% filter(Gene == g)
      
      # Loop through Sex
      for (sex in unique(gene_data$Sex)) {
        sub_data <- gene_data %>% filter(Sex == sex) %>%
          arrange(Timepoint)
        
        # RAIN expects a matrix with timepoints as columns, replicates as rows
        # So we pivot wider to have timepoints as columns
        expr_mat <- sub_data %>%
          select(Timepoint, Expression) %>%
          group_by(Timepoint) %>%
          summarize(Expression = mean(Expression), .groups = "drop") %>%
          pivot_wider(names_from = Timepoint, values_from = Expression) %>%
          as.matrix()
        
        # Run RAIN
        rain_res <- rain(
          t(expr_mat),        # transpose so rows = timepoints
          deltat = diff(sort(unique(sub_data$Timepoint)))[1], # spacing between timepoints
          period = 24,
          method = "independent"
        )
        
        results_list[[paste(g, sex, sep = "_")]] <- data.frame(
          Gene = g,
          Sex = sex,
          pVal = rain_res$pVal,
          period = rain_res$period,
          peak = rain_res$peak
        )
      }
    }
    
    rain_results <- do.call(rbind, results_list)
    
    # Adjust p-values for multiple testing
    rain_results$adjP <- p.adjust(rain_results$pVal, method = "BH")
    
    rain_results %>% arrange(adjP) %>% head()
    
