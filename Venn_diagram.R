
#===============================================================================
# Set working directory and load libraries
setwd("/Users/judyabuel/Desktop/Xist/circadian_atlas")

library(edgeR)
library(limma)
library(ggVennDiagram)

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

# Create DGEList object and normalize
dge <- DGEList(counts)
dge <- calcNormFactors(dge, method = "TMM")

# Filter low-expressed genes
keep <- rowSums(cpm(dge) >= 1) >= 3
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Design matrix with main effects + interaction
design <- model.matrix(~ Sex * Timepoint, data = metadata)
colnames(design) <- make.names(colnames(design))  # Clean column names

# Voom transformation and fit linear model
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Voom log2-CPM matrix
voom_log2_CPM <- v$E 

#===============================================================================
# Set Contrast for Female vs Male at each timepoint
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
# Extract DEGs for each timepoint (adj.P.Val < 0.05)
deg_list <- lapply(colnames(contrast.matrix), function(contrast_name) {
  tt <- topTable(fit2, coef = contrast_name, number = Inf, adjust.method = "BH")
  degs <- rownames(tt)[tt$adj.P.Val < 0.05]
  return(degs)
})
names(deg_list) <- colnames(contrast.matrix)

#===============================================================================
# Combine all DEGs by sex
# Since contrasts are Female_vs_Male, these DEGs are "Female DEGs" per timepoint
deg_female_all <- unique(unlist(deg_list))

# Optional: Get Male DEGs if you have contrasts Male_vs_Female
# Or define Male DEGs as genes with logFC < 0
deg_male_all <- unique(unlist(lapply(deg_list, function(gene_set) {
  tt <- topTable(fit2, coef = 1, number = Inf)  # same tt
  rownames(tt)[tt$adj.P.Val < 0.05 & tt$logFC < 0]  # down in female = up in male
})))

#===============================================================================
library(edgeR)
library(limma)
library(ggVennDiagram)
library(gridExtra)  # for arranging multiple plots

#---------------------------
# List of timepoint contrasts
timepoints <- colnames(contrast.matrix)

# Create empty list to store plots
venn_plots <- list()

for(tp in timepoints) {
  
  # Extract top table for this contrast
  tt <- topTable(fit2, coef = tp, number = Inf, adjust.method = "BH")
  
  # Female DEGs (up in females)
  deg_female <- rownames(tt)[tt$adj.P.Val < 0.05 & tt$logFC > 0]
  
  # Male DEGs (up in males = down in females)
  deg_male <- rownames(tt)[tt$adj.P.Val < 0.05 & tt$logFC < 0]
  
  # Prepare Venn list
  venn_list <- list(
    Female = deg_female,
    Male   = deg_male
  )
  
  # Generate plot and store in list
  venn_plots[[tp]] <- ggVennDiagram(venn_list,
                                    category.names = c("Female", "Male"),
                                    edge_size = 1,
                                    set_size = 5,
                                    label_alpha = 0,
                                    fill_color = c("pink", "blue")) +
    ggtitle(paste("DEG Overlap at", tp)) +
    theme_void()
}


# Display all Venn diagrams in a grid
do.call(grid.arrange, c(venn_plots, ncol = 3))  # adjust ncol as desired










#-------------------------------------------------------------------------------
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


