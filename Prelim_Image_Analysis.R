
setwd("/Users/judyabuel/Desktop")

install.packages("tidyr")

# Load your data (modify the path as necessary)
df <- read.csv("prelim_fpws.csv")

# Load necessary libraries
library(dplyr)
library(readr)
library(tidyr)

  
# B dataset (B01 - B03 labeled as '12:12')
b_data <- df %>%
  filter(Well.Name %in% c("B01", "B02", "B03")) %>%
  mutate(Label = '12:12') %>%
  drop_na()
  
# C dataset (C01 - C03 labeled as '11:11')
c_data <- df %>%
  filter(Well.Name %in% c("C01", "C02", "C03")) %>%
  mutate(Label = '11:11') %>%
  drop_na()

# Save each dataset to CSV

write.csv(b_data, 'Deletion_data.csv', row.names = FALSE)
write.csv(c_data, 'Compensation_data.csv', row.names = FALSE)

# cleaned data after removing all NAs
a_data <- read.csv("WT_data.csv")
b_data <- read.csv("Deletion_data.csv")
c_data <- read.csv("Compensation_data.csv")

combined_data <- bind_rows(a_data, b_data, c_data)

write.csv(combined_data, 'Prelim_Analysis.csv', row.names = FALSE)


#------------------------------------------------------------------------------
# Sum of all the AREA per channel in all slides a = WT, b = 12:12, c = 11:11 
# in all 3 regions

wt_dapi_area <- sum(a_data$Dapi_Area)
wt_Meg3_area <- sum(a_data$Meg3_Area)
wt_Snhg14_area <- sum(a_data$Snhg14_Area)
wt_Xist_area <- sum(a_data$Xist_Area)

# PLOT the WT area
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)

a_area <- a_data %>%
  dplyr::select(Well.Name, Dapi_Area, Meg3_Area, Snhg14_Area, Xist_Area) %>%
  pivot_longer(cols = c(Dapi_Area, Meg3_Area, Snhg14_Area, Xist_Area),
               names_to = "Channel",
               values_to = "Area")

a_area$Channel <- gsub("\\_$", "", a_area$Channel)  # remove trailing dots
a_area$Channel <- gsub("_Area", "", a_area$Channel)

comparisons_to_dapi <- list(
  c("Dapi", "Meg3"),
  c("Dapi", "Snhg14"),
  c("Dapi", "Xist"),
)

ggplot(a_area, aes(x = Channel, y = Area, fill = Channel)) +
  geom_jitter(aes(color = Channel),
              width = 0.2, alpha = 0.4, size = 1.5) +  # dots first (background)
  geom_boxplot(alpha = 0.7, outlier.shape = NA, color = "black") +  # box in front
  stat_compare_means(
    comparisons = comparisons_to_dapi,
    method = "t.test",
    label = "p.signif",
    label.y = c(160, 170, 180)
  ) +
  stat_compare_means(label.y = 240) +
  scale_fill_manual(values = c("lightblue", "lightgreen", "pink", "cyan")) +
  scale_color_manual(values = c("lightblue", "lightgreen", "pink", "cyan")) +
  labs(
    title = " lncRNA Cloud Size WT",
    x = "lncRNA Channel",
    y = "Area (micron)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3)
  )





# --------------------------------------
# lncRNA Clous Size Deletion F 12:12

del_dapi_area <- sum(b_data$Dapi_Area)
del_Meg3_area <- sum(b_data$Meg3_Area)
del_Snhg14_area <- sum(b_data$Snhg14_Area)
del_Xist_area <- sum(b_data$Xist_Area)

b_area <- b_data %>%
  select(Well.Name, Dapi_Area, Meg3_Area, Snhg14_Area, Xist_Area) %>%
  pivot_longer(cols = c(Dapi_Area, Meg3_Area, Snhg14_Area, Xist_Area),
               names_to = "Channel",
               values_to = "Area")

b_area$Channel <- gsub("\\_$", "", b_area$Channel)  # remove trailing dots
b_area$Channel <- gsub("_Area", "", b_area$Channel)

comparisons_to_dapi <- list(
  c("Dapi", "Meg3"),
  c("Dapi", "Snhg14"),
  c("Dapi", "Xist")
)

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)


ggplot(b_area, aes(x = Channel, y = Area, fill = Channel)) +
  geom_jitter(aes(color = Channel),
              width = 0.2, alpha = 0.4, size = 1.5) +  # dots first (background)
  geom_boxplot(alpha = 0.7, outlier.shape = NA, color = "black") +  # box in front
  stat_compare_means(
    comparisons = comparisons_to_dapi,
    method = "t.test",
    label = "p.signif",
    label.y = c(160, 170, 180)
  ) +
  stat_compare_means(label.y = 240) +
  scale_fill_manual(values = c("lightblue", "lightgreen", "pink", "cyan")) +
  scale_color_manual(values = c("lightblue", "lightgreen", "pink", "cyan")) +
  labs(
    title = " lncRNA Cloud Size Deletion",
    x = "lncRNA Channel",
    y = "Area (micron)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3)
  )







# --------------------------------------
# lncRNA Cloud Size COMPENSATION F 11:11

comp_dapi_area <- sum(c_data$Dapi_Area)
comp_Meg3_area <- sum(c_data$Meg3_Area)
comp_Snhg14_area <- sum(c_data$Snhg14_Area)
comp_Xist_area <- sum(c_data$Xist_Area)

c_area <- c_data %>%
  select(Well.Name, Dapi_Area, Meg3_Area, Snhg14_Area, Xist_Area) %>%
  pivot_longer(cols = c(Dapi_Area, Meg3_Area, Snhg14_Area, Xist_Area),
               names_to = "Channel",
               values_to = "Area")

c_area$Channel <- gsub("\\_$", "", c_area$Channel)  # remove trailing dots
c_area$Channel <- gsub("_Area", "", c_area$Channel)

comparisons_to_dapi <- list(
  c("Dapi", "Meg3"),
  c("Dapi", "Snhg14"),
  c("Dapi", "Xist")
)

ggplot(c_area, aes(x = Channel, y = Area, fill = Channel)) +
  geom_jitter(aes(color = Channel),
              width = 0.2, alpha = 0.4, size = 1.5) +  # dots first (background)
  geom_boxplot(alpha = 0.7, outlier.shape = NA, color = "black") +  # box in front
  stat_compare_means(
    comparisons = comparisons_to_dapi,
    method = "t.test",
    label = "p.signif",
    label.y = c(160, 170, 180)
  ) +
  stat_compare_means(label.y = 240) +
  scale_fill_manual(values = c("lightblue", "lightgreen", "pink", "cyan")) +
  scale_color_manual(values = c("lightblue", "lightgreen", "pink", "cyan")) +
  labs(
    title = " lncRNA Cloud Size Compensation",
    x = "lncRNA Channel",
    y = "Area (micron)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3)
  )


# SANITY CHECK!
# Load necessary libraries
library(dplyr)
library(knitr)

# Summary statistics for a_data
summary_a <- a_data %>%
  summarize(
    wt_dapi_area = mean(Dapi_Area, na.rm = TRUE),
    wt_Meg3_area = mean(Meg3_Area, na.rm = TRUE),
    wt_Snhg14_area = mean(Snhg14_Area, na.rm = TRUE),
    wt_Xist_area = mean(Xist_Area, na.rm = TRUE)
  )

# Summary statistics for b_data
summary_b <- b_data %>%
  summarize(
    del_dapi_area = mean(Dapi_Area, na.rm = TRUE),
    del_Meg3_area = mean(Meg3_Area, na.rm = TRUE),
    del_Snhg14_area = mean(Snhg14_Area, na.rm = TRUE),
    del_Xist_area = mean(Xist_Area, na.rm = TRUE)
  )

# Summary statistics for c_data
summary_c <- c_data %>%
  summarize(
    comp_dapi_area = mean(Dapi_Area, na.rm = TRUE),
    comp_Meg3_area = mean(Meg3_Area, na.rm = TRUE),
    comp_Snhg14_area = mean(Snhg14_Area, na.rm = TRUE),
    comp_Xist_area = mean(Xist_Area, na.rm = TRUE)
  )


# Combine results into a single data frame for a clean table output
summary_stats <- data.frame(
  Channel = c("Dapi_Area", "Meg3_Area", "Snhg14_Area", "Xist_Area"),
  WT = unlist(summary_a),
  Twelve = unlist(summary_b),
  Eleven = unlist(summary_c)
)

# Print the table
kable(summary_stats, caption = "Summary Statistics of Areas for Different Datasets")


#-------------------------------------------------------------------
# Load ggplot2 if not already loaded
library(ggplot2)

# Combine datasets for plotting distributions
combined_data <- bind_rows(
  a_data %>% mutate(Type = "WT"),
  b_data %>% mutate(Type = "Deletion"),
  c_data %>% mutate(Type = "Compensation")
)

# Plot distributions for each area
ggplot(combined_data) +
  geom_density(aes(x = Dapi_Area, fill = Type), alpha = 0.5) +
  labs(title = "Distribution of Dapi Area", x = "Area", y = "Density") +
  theme_minimal() +  # Use a minimal theme
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_text(face = "bold"), # Bold title
    axis.line = element_line(color = "black", linewidth = 0.5) # add axis line
  )

  
ggplot(combined_data) +
  geom_density(aes(x = Meg3_Area, fill = Type), alpha = 0.5) +
  labs(title = "Distribution of Meg3 Area", x = "Area", y = "Density") +
  theme_minimal() +  # Use a minimal theme
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      plot.title = element_text(face = "bold"),  # Bold title
      axis.line = element_line(color = "black", linewidth = 0.5)
    )
  
  
ggplot(combined_data) +
  geom_density(aes(x = Snhg14_Area, fill = Type), alpha = 0.5) +
  labs(title = "Distribution of Snhg14 Area", x = "Area", y = "Density") +
  theme_minimal() +  # Use a minimal theme
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_text(face = "bold"),  # Bold title
    axis.line = element_line(color = "black", linewidth = 0.5)
  )

ggplot(combined_data) +
  geom_density(aes(x = Xist_Area, fill = Type), alpha = 0.5) +
  labs(title = "Distribution of Xist Area", x = "Area", y = "Density") +
  theme_minimal() +  # Use a minimal theme
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_text(face = "bold"),  # Bold title
    axis.line = element_line(color = "black", linewidth = 0.5)
  )



# ---------------------------------------------------

# Sum of all the fluorescence intensity per channel in all regions
a_Dapi_fluo <- sum(a_data$Dapi_.Integrated_Intensity)
a_Meg3_fluo <- sum(a_data$Meg3_Integrated_Intensity)
a_Snhg14_fluo <- sum(a_data$Snhg14_Integrated_Intensity)
a_Xist_fluo <- sum(a_data$Xist_Integrated_Intensity)

# Twelve : Twelve
b_Dapi_fluo <- sum(b_data$Dapi_Integrated_Intensity)
b_Meg3_fluo <- sum(b_data$Meg3_Integrated_Intensity)
b_Snhg14_fluo <- sum(b_data$Snhg14_Integrated_Intensity)
b_Xist_fluo <- sum(b_data$Xist_Integrated_Intensity)

# Eleven : Eleven
c_Dapi_fluo <- sum(c_data$Dapi_Integrated_Intensity)
c_Meg3_fluo <- sum(c_data$Meg3_Integrated_Intensity)
c_Snhg14_fluo <- sum(c_data$Snhg14_Integrated_Intensity)
c_Xist_fluo <- sum(c_data$Xist_Integrated_Intensity)

#----------------------------

# Co-localization Analysis: b_data = 12:12 Female DELETION sample
# Pearson Correlation between Snhg14 and Xist Co-localization
b_cleaned_data <- b_data %>%
  filter(!is.na(Snhg14_Integrated_Intensity) & !is.na(Xist_Integrated_Intensity))

# Pearson Correlation between Snhg14 and Meg3 Co-localization
b_cor_result <- cor.test(b_cleaned_data$Snhg14_Integrated_Intensity,
                         b_cleaned_data$Xist_Integrated_Intensity,
                         method = "pearson")

# Extract the correlation and R^2 values
r_value <- b_cor_result$estimate
R2_value <- r_value^2

# Create the plot
ggplot(b_cleaned_data, aes(x = Snhg14_Integrated_Intensity, y = Xist_Integrated_Intensity)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "red") +
  labs(
    title = "Pearson Correlation of Snhg14 vs Xist in F 12:12",
    subtitle = paste0(
      " R² = ", round(R2_value, 3),
      " | p < ", format(cor.test(b_cleaned_data$Snhg14_Integrated_Intensity, b_cleaned_data$Xist_Integrated_Intensity)$p.value, scientific = TRUE, digits = 2)
    ),
    x = "Snhg14 Integrated Intensity",
    y = "Xist Integrated Intensity"
  ) +
  theme(
    text = element_text(size = 10),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_text(face = "bold"),  # Bold title
    axis.line = element_line(color = "black", linewidth = 0.5)
  )

#-----------------------------
# Pearson Correlation between Snhg14 and Meg3 Co-localization
bb_cleaned_data <- b_data %>%
  filter(!is.na(Snhg14_Integrated_Intensity) & !is.na(Meg3_Integrated_Intensity))

# Pearson Correlation between Snhg14 and Meg3 Co-localization
bb_cor_result <- cor.test(bb_cleaned_data$Snhg14_Integrated_Intensity,
                         bb_cleaned_data$Meg3_Integrated_Intensity,
                         method = "pearson")

# Extract the correlation and R^2 values
r_value <- bb_cor_result$estimate
R2_value <- r_value^2

# Create the plot
ggplot(bb_cleaned_data, aes(x = Snhg14_Integrated_Intensity, y = Meg3_Integrated_Intensity)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "red") +
  labs(
    title = "Pearson Correlation of Snhg14 vs Meg3 in F 12:12",
    subtitle = paste0(
      "r = ", round(r_value, 3),
      " | R² = ", round(R2_value, 3),
      " | p < ", format(cor.test(b_cleaned_data$Snhg14_Integrated_Intensity, b_cleaned_data$Meg3_Integrated_Intensity)$p.value, scientific = TRUE, digits = 2)
    ),
    x = "Snhg14 Integrated Intensity",
    y = "Meg3 Integrated Intensity"
  ) +
  theme(
    text = element_text(size = 10),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_text(face = "bold"),  # Bold title
    axis.line = element_line(color = "black", linewidth = 0.5)
  )

#-----------------------------
# Co-localization Analysis: c_data = 11:11 Female COMPENSATION sample
# Pearson Correlation between Snhg14 and Xist Co-localization
c_cleaned_data <- c_data %>%
  filter(!is.na(Snhg14_Integrated_Intensity) & !is.na(Xist_Integrated_Intensity))

# Pearson Correlation between Snhg14 and Meg3 Co-localization
c_cor_result <- cor.test(c_cleaned_data$Snhg14_Integrated_Intensity,
                         c_cleaned_data$Xist_Integrated_Intensity,
                         method = "pearson")

# Extract the correlation and R^2 values
r_value <- c_cor_result$estimate
R2_value <- r_value^2

# Create the plot
ggplot(c_cleaned_data, aes(x = Snhg14_Integrated_Intensity, y = Xist_Integrated_Intensity)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "red") +
  labs(
    title = "Pearson Correlation of Snhg14 vs Xist in F 11:11",
    subtitle = paste0(
      "r = ", round(r_value, 3),
      " | R² = ", round(R2_value, 3),
      " | p < ", format(cor.test(c_cleaned_data$Snhg14_Integrated_Intensity, c_cleaned_data$Xist_Integrated_Intensity)$p.value, scientific = TRUE, digits = 2)
    ),
    x = "Snhg14 Integrated Intensity",
    y = "Xist Integrated Intensity"
  ) +
  theme(
    text = element_text(size = 10),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_text(face = "bold"),  # Bold title
    axis.line = element_line(color = "black", linewidth = 0.5)
  )


#--
# Pearson Correlation between Snhg14 and Meg3 Co-localization
cc_cleaned_data <- c_data %>%
  filter(!is.na(Snhg14_Integrated_Intensity) & !is.na(Meg3_Integrated_Intensity))

# Pearson Correlation between Snhg14 and Meg3 Co-localization
cc_cor_result <- cor.test(cc_cleaned_data$Snhg14_Integrated_Intensity,
                          cc_cleaned_data$Meg3_Integrated_Intensity,
                          method = "pearson")

# Extract the correlation and R^2 values
r_value <- cc_cor_result$estimate
R2_value <- r_value^2

# Create the plot
ggplot(cc_cleaned_data, aes(x = Snhg14_Integrated_Intensity, y = Meg3_Integrated_Intensity)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "red") +
  labs(
    title = "Pearson Correlation of Snhg14 vs Meg3 in F 11:11",
    subtitle = paste0(
      "r = ", round(r_value, 3),
      " | R² = ", round(R2_value, 3),
      " | p < ", format(cor.test(cc_cleaned_data$Snhg14_Integrated_Intensity, cc_cleaned_data$Meg3_Integrated_Intensity)$p.value, scientific = TRUE, digits = 2)
    ),
    x = "Snhg14 Integrated Intensity",
    y = "Meg3 Integrated Intensity"
  ) +
  theme(
    text = element_text(size = 10),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_text(face = "bold"),  # Bold title
    axis.line = element_line(color = "black", linewidth = 0.5)
  )


#-----------------------------
# Co-localization Analysis: a_data = WT
# Pearson Correlation between Snhg14 and Xist Co-localization
a_cleaned_data <- a_data %>%
  filter(!is.na(Snhg14_Integrated_Intensity) & !is.na(Xist_Integrated_Intensity))

# Pearson Correlation between Snhg14 and Meg3 Co-localization
a_cor_result <- cor.test(a_cleaned_data$Snhg14_Integrated_Intensity,
                         a_cleaned_data$Xist_Integrated_Intensity,
                         method = "pearson")

# Extract the correlation and R^2 values
r_value <- a_cor_result$estimate
R2_value <- r_value^2

# Create the plot
# Osman's version

ggplot(a_cleaned_data, aes(x = log2(Snhg14_Integrated_Intensity), y = log2(Xist_Integrated_Intensity))) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", color = "red") +
  labs(
    title = "Pearson Correlation of Snhg14 vs Xist in WT",
    subtitle = paste0(
      " R² = ", round(R2_value, 3),
      " | p < ", format(cor.test(a_cleaned_data$Snhg14_Integrated_Intensity, a_cleaned_data$Xist_Integrated_Intensity)$p.value, scientific = TRUE, digits = 2)
    ),
    x = "log2(Snhg14 Integrated Intensity)",
    y = "log2(Xist Integrated Intensity)"
  ) +
  theme(
    text = element_text(size = 10),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_text(face = "bold"),  # Bold title
    axis.line = element_line(color = "black", linewidth = 0.5)
  )
###


# Pearson Correlation between Snhg14 and Meg3 Co-localization
# Wild Type
aa_cleaned_data <- a_data %>%
  filter(!is.na(Snhg14_Integrated_Intensity) & !is.na(Meg3_Integrated_Intensity))

# Pearson Correlation between Snhg14 and Meg3 Co-localization
aa_cor_result <- cor.test(aa_cleaned_data$Snhg14_Integrated_Intensity,
                         aa_cleaned_data$Meg3_Integrated_Intensity,
                         method = "pearson")

# Extract the correlation and R^2 values
r_value <- aa_cor_result$estimate
R2_value <- r_value^2

# Create the plot
ggplot(aa_cleaned_data, aes(x = log2(Snhg14_Integrated_Intensity), y = log2(Meg3_Integrated_Intensity))) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", color = "red") +
  labs(
    title = "Pearson Correlation of Snhg14 vs Meg3 in WT",
    subtitle = paste0(
      " R² = ", round(R2_value, 3),
      " | p < ", format(cor.test(aa_cleaned_data$Snhg14_Integrated_Intensity, aa_cleaned_data$Meg3_Integrated_Intensity)$p.value, scientific = TRUE, digits = 2)
    ),
    x = "log2(Snhg14 Integrated Intensity)",
    y = "log2(Meg3 Integrated Intensity)"
  ) +
  theme(
    text = element_text(size = 10),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_text(face = "bold"),  # Bold title
    axis.line = element_line(color = "black", linewidth = 0.5)
  )








# ------------------------------------------------------------------------------
# MANDERS Correlation Analysis
# The Manders’ overlap coefficient (often M1 and M2) is a measure of co-localization 
# between two fluorescent signals (say, proteins or RNAs) in microscopy images.

# Define A and B for the first pair "Snhg14 - Xist"
A1 <- b_cleaned_data$Snhg14_Integrated_Intensity
B1 <- b_cleaned_data$Xist_Integrated_Intensity

# Ensure no NAs
A1[is.na(A1)] <- 0
B1[is.na(B1)] <- 0

# Compute Manders coefficients for "Snhg14 - Xist"
M1 <- sum(A1[B1 > 0], na.rm = TRUE) / sum(A1, na.rm = TRUE)
M2 <- sum(B1[A1 > 0], na.rm = TRUE) / sum(B1, na.rm = TRUE)

# Now define A and B for the second pair "Snhg14 - Meg3"
A2 <- b_cleaned_data$Snhg14_Integrated_Intensity
B2 <- b_cleaned_data$Meg3_Integrated_Intensity

# Ensure no NAs
A2[is.na(A2)] <- 0
B2[is.na(B2)] <- 0

# Compute Manders coefficients for "Snhg14 - Meg3"
M1_Snhg14_Meg3 <- sum(A2[B2 > 0], na.rm = TRUE) / sum(A2, na.rm = TRUE)
M2_Snhg14_Meg3 <- sum(B2[A2 > 0], na.rm = TRUE) / sum(B2, na.rm = TRUE)

# Now define A and B for the third pair "Xist - Meg3"
A3 <- b_cleaned_data$Xist_Integrated_Intensity
B3 <- b_cleaned_data$Meg3_Integrated_Intensity

# Ensure no NAs
A3[is.na(A3)] <- 0
B3[is.na(B3)] <- 0

# Compute Manders coefficients for "Xist - Meg3"
M1_Xist_Meg3 <- sum(A3[B3 > 0], na.rm = TRUE) / sum(A3, na.rm = TRUE)
M2_Xist_Meg3 <- sum(B3[A3 > 0], na.rm = TRUE) / sum(B3, na.rm = TRUE)

# Create data frame for Manders data
manders_df_b <- data.frame(
  Pair = c("Snhg14 - Xist", "Snhg14 - Meg3", "Xist - Meg3"),
  M1 = c(M1, M1_Snhg14_Meg3, M1_Xist_Meg3),
  M2 = c(M2, M2_Snhg14_Meg3, M2_Xist_Meg3)
)

# Check the contents of manders_df_b
print(manders_df_b)

# Convert to long format for plotting
manders_long_b <- manders_df_b %>%
  pivot_longer(cols = c(M1, M2), names_to = "Manders", values_to = "Coefficient") %>%
  mutate(SourceChannel = case_when(
    Manders == "M1" & Pair == "Snhg14 - Xist" ~ "Snhg14",
    Manders == "M2" & Pair == "Snhg14 - Xist" ~ "Xist",
    Manders == "M1" & Pair == "Snhg14 - Meg3" ~ "Snhg14",
    Manders == "M2" & Pair == "Snhg14 - Meg3" ~ "Meg3",
    Manders == "M1" & Pair == "Xist - Meg3" ~ "Xist",
    Manders == "M2" & Pair == "Xist - Meg3" ~ "Meg3"
  ))

# Define colors for channels
channel_colors <- c(
  "Xist" = "lightblue",
  "Meg3" = "lightgreen",
  "Snhg14" = "pink"
)

# Plotting the Manders Coefficients
ggplot(manders_long, aes(x = Pair, y = Coefficient, fill = SourceChannel)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
  geom_text(aes(label = round(Coefficient, 2)),
            position = position_dodge(width = 0.8),
            vjust = -0.6, size = 4) +
  #geom_text(aes(label = Label)) +
            #position = position_dodge(width = 0.8),
            #vjust = 2, size = 3, angle = 90, color = "black") +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 1)) +
  labs(
    title = "Manders Coefficients per lncRNA Pair in WT",
    x = "lncRNA Pair",
    y = "Manders Coefficient (%)",
    fill = "Channel"
  ) +
  scale_fill_manual(values = channel_colors) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3)
  )


#------------

# Calculate Manders coefficients for "Snhg14 - Meg3"
A1 <- b_cleaned_data$Snhg14_Integrated_Intensity
B1 <- b_cleaned_data$Meg3_Integrated_Intensity

# Ensure no NAs
A1[is.na(A1)] <- 0
B1[is.na(B1)] <- 0

# Compute Manders coefficients for "Snhg14 - Meg3"
M1_Snhg14_Meg3 <- sum(A1[B1 > 0], na.rm = TRUE) / sum(A1, na.rm = TRUE)
M2_Snhg14_Meg3 <- sum(B1[A1 > 0], na.rm = TRUE) / sum(B1, na.rm = TRUE)

# Calculate Manders coefficients for "Xist - Meg3"
A2 <- b_cleaned_data$Xist_Integrated_Intensity
B2 <- b_cleaned_data$Meg3_Integrated_Intensity

# Ensure no NAs
A2[is.na(A2)] <- 0
B2[is.na(B2)] <- 0

# Compute Manders coefficients for "Xist - Meg3"
M1_Xist_Meg3 <- sum(A2[B2 > 0], na.rm = TRUE) / sum(A2, na.rm = TRUE)
M2_Xist_Meg3 <- sum(B2[A2 > 0], na.rm = TRUE) / sum(B2, na.rm = TRUE)

# Create data frame for Manders data
manders_df_b <- data.frame(
  Pair = c("Snhg14 - Xist", "Snhg14 - Meg3", "Xist - Meg3"),
  M1 = c(M1, M1_Snhg14_Meg3, M1_Xist_Meg3),
  M2 = c(M2, M2_Snhg14_Meg3, M2_Xist_Meg3)
)





###-----------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyr)

# Reshape combined_data to long format for individual points
long_combined_data <- combined_data %>%
  pivot_longer(cols = c(Dapi_Area, Meg3_Area, Snhg14_Area, Xist_Area), 
               names_to = "Channel", 
               values_to = "Area")

# Summarizing the data to calculate mean and SE for each area type by genotype
summary_with_errors <- long_combined_data %>%
  group_by(Label) %>%
  summarise(
    Dapi_Mean = mean(Area[Channel == "Dapi_Area"], na.rm = TRUE),
    Dapi_SE = sd(Area[Channel == "Dapi_Area"], na.rm = TRUE) / sqrt(n()),
    Meg3_Mean = mean(Area[Channel == "Meg3_Area"], na.rm = TRUE),
    Meg3_SE = sd(Area[Channel == "Meg3_Area"], na.rm = TRUE) / sqrt(n()),
    Snhg14_Mean = mean(Area[Channel == "Snhg14_Area"], na.rm = TRUE),
    Snhg14_SE = sd(Area[Channel == "Snhg14_Area"], na.rm = TRUE) / sqrt(n()),
    Xist_Mean = mean(Area[Channel == "Xist_Area"], na.rm = TRUE),
    Xist_SE = sd(Area[Channel == "Xist_Area"], na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# Reshape summary for visualization
library(ggplot2)

ggplot(summary_long, aes(x = Channel, y = Mean, fill = Label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.6) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), 
                width = 0.2, position = position_dodge(width = 0.9), size = 1.2) +
  labs(title = "Mean Area by Genotype and Channel",
       x = "Channel", 
       y = "Mean Area") +
  geom_point(data = long_combined_data,
             aes(x = Channel, y = Area, color = Label),
             position = position_dodge(width = 0.9),
             shape = 21, size = 2.5, alpha = 0.6) +
  coord_cartesian(ylim = c(0, 200)) +  # Use coord_cartesian for limits
  theme_minimal() +
  scale_fill_brewer(palette = "Pastel1")



library(dplyr)
library(ggplot2)
library(tidyr)
install.packages("ggsignif")
library(ggsignif)

# Reshape combined_data to long format for individual points
long_combined_data <- combined_data %>%
  pivot_longer(cols = c(Dapi_Area, Meg3_Area, Snhg14_Area, Xist_Area), 
               names_to = "Channel", 
               values_to = "Area")

# Summarizing the data to calculate mean and SE for each area type by genotype
summary_with_errors <- long_combined_data %>%
  group_by(Channel, Label) %>%
  summarise(
    Mean = mean(Area, na.rm = TRUE),
    SE = sd(Area, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Plotting
ggplot(summary_with_errors, aes(x = Channel, y = Mean, fill = Label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.6) +
  geom_errorbar(aes(ymin = Mean - (SE * 2), ymax = Mean + (SE * 2)),  # Increased error bar range
                width = 0.2, position = position_dodge(width = 1), size = 1.5) +
  labs(title = "Cloud Size of lncRNAs per Genotype",
       x = "lncRNAs", 
       y = "Mean Area (micron)") +
  geom_jitter(data = long_combined_data,  # Use geom_jitter for individual cell data
              aes(x = Channel, y = Area, color = Label),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
              fill = "grey", shape = 21, size = 1, alpha = 0.2) +  # Adjusted jitter points
  coord_cartesian(ylim = c(0, 200)) +  # Use coord_cartesian for limits
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 1, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)) +
  scale_fill_brewer(palette = "Pastel1")







library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsignif)
library(broom)  # For tidy results
library(purrr)

# Reshape combined_data to long format for individual points
long_combined_data <- combined_data %>%
  pivot_longer(cols = c(Dapi_Area, Meg3_Area, Snhg14_Area, Xist_Area), 
               names_to = "Channel", 
               values_to = "Area")

# Summarizing the data to calculate mean and SE for each area type by genotype
summary_with_errors <- long_combined_data %>%
  group_by(Channel, Label) %>%
  summarise(
    Mean = mean(Area, na.rm = TRUE),
    SE = sd(Area, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


# Perform ANOVA and Tukey's test for each Channel
anova_results <- long_combined_data %>%
  group_by(Channel) %>%
  summarise(
    aov_result = list(aov(Area ~ Label, data = cur_data())),
    .groups = 'drop'
  )

# Calculate Tukey HSD for each ANOVA result
turkey_results <- anova_results %>%
  mutate(tukey = map(aov_result, ~TukeyHSD(.x))) %>%
  mutate(tukey_tidy = map(tukey, tidy)) %>%
  unnest(tukey_tidy)

# Extract significant comparisons for plotting
pairwise_comps <- turkey_results %>%
  filter(adj.p.value < 0.05) %>%
  mutate(comparison = contrast) %>%
  select(comparison)

# Prepare the comparisons list for geom_signif
comparison_list <- pairwise_comps$comparison %>% 
  strsplit("-") %>% 
  lapply(function(x) setNames(x, c("group1", "group2")))



# Plotting
p <- ggplot(summary_with_errors, aes(x = Channel, y = Mean, fill = Label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 2.5), alpha = 0.8) +
  geom_jitter(data = long_combined_data,  
              aes(x = Channel, y = Area),
              position = position_jitterdodge(jitter.width = .5, dodge.width = 2.5),
              color = "black", shape = 21, size = 0.2, alpha = 0.2) + 
  geom_errorbar(aes(ymin = Mean - (SE * 2), ymax = Mean + (SE * 2)),  
                width = 0.5, position = position_dodge(width = 2.5), linewidth = 1) +
  labs(title = "Cloud Size of lncRNAs per Genotype",
       x = "Genotype", 
       y = "Mean Area (micron)") +
  coord_cartesian(ylim = c(0, 200)) +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 1, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)) +
  scale_fill_brewer(palette = "Pastel1") +
  scale_x_discrete(labels = function(x) gsub("_Area", " ", x))



# Clean up the comparison_list if duplicates exist
comparison_list_cleaned <- lapply(comparison_list, function(x) {
  x[!duplicated(x)]  # Remove duplicates if any
})

# Plotting the graph with modified geom_signif
p + 
  geom_signif(comparisons = comparison_list,
              map_signif_level = TRUE,  # Automatically map significance levels
              position = position_dodge(width = 0.9),  # Ensure it's aligned with bars
              vjust = 0.5,
              size = 4) + 
  facet_wrap(~ Channel, scales = "free_x")  # Create separate panels for each channel








#--------------THIS WORKS, but it doesn't have statistical analysis---------------------
# Example of your raw data
data <- tibble(
  Channel = rep(c("Dapi_Area", "Meg3_Area", "Snhg14_Area", "Xist_Area"), each = 30),
  Label = rep(rep(c("WT", "11:11", "12:12"), each = 10), 4),
  Value = rnorm(120, mean = rep(c(80, 85, 99, 58, 63, 70, 50, 7, 7, 36, 14, 14), each = 10), sd = 5)
)

library(dplyr)
install.packages("rstatix")
library(rstatix)

# Perform one-way ANOVA or pairwise comparisons for each Channel
stats_results <- data %>%
  group_by(Channel) %>%
  anova_test(Value ~ Label) %>%
  get_anova_table()

# Or pairwise t-tests with correction
pairwise_results <- data %>%
  group_by(Channel) %>%
  pairwise_t_test(Value ~ Label, p.adjust.method = "bonferroni")


summary_with_p <- summary_with_errors %>%
  left_join(pairwise_results, by = "Channel")


summary_with_p %>%
  select(Channel, Label, Mean, SE, group1, group2, p, p.adj, p.signif)

write.csv(summary_with_p, "stat_cloudsize_genotypes.csv")

library(ggplot2)
library(ggpubr)

ggplot(data, aes(x = Label, y = Value, fill = Label)) +
  geom_boxplot() +
  facet_wrap(~ Channel, scales = "free_y") +
  stat_pvalue_manual(pairwise_results, label = "p.signif", y.position = 1.1 * max(data$Value)) +
  labs(title = "Cloud Size of lncRNAs per Genotype",
       x = "Genotype", 
       y = "Mean Area (µm²)") +
  coord_cartesian(ylim = c(0, 200)) +  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  facet_wrap(~ Channel, scales = "free_x") +
  scale_x_discrete(labels = function(x) gsub("_Area", " ", x))

#------------------------------------------------------

#This looks a better version of all the graphs I created. 
                   #This is the last version I created

library(ggplot2)
library(ggsignif)

library(dplyr)
install.packages("rstatix")
library(rstatix)

# Perform one-way ANOVA or pairwise comparisons for each Channel
stats_results <- data %>%
  group_by(Channel) %>%
  anova_test(Value ~ Label) %>%
  get_anova_table()

# Or pairwise t-tests with correction
pairwise_results <- data %>%
  group_by(Channel) %>%
  pairwise_t_test(Value ~ Label, p.adjust.method = "bonferroni")


summary_with_p <- summary_with_errors %>%
  left_join(pairwise_results, by = "Channel")


summary_with_p %>%
  select(Channel, Label, Mean, SE, group1, group2, p, p.adj, p.signif)
# Remove duplicates in your comparison list
comparison_list_cleaned <- lapply(comparison_list, unique)

summary_with_errors$Label <- recode(summary_with_errors$Label,
                                    "WT" = "Wild Type",
                                    "11:11" = "HET/TG",
                                    "12:12" = "HET/WT"
)

long_combined_data$Label <- recode(long_combined_data$Label,
                                   "WT" = "Wild Type",
                                   "11:11" = "HET/TG",
                                   "12:12" = "HET/WT"
)


p <- ggplot(summary_with_errors, aes(x = Label, y = Mean, fill = Label)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.9), 
           alpha = 0.8) +
  geom_jitter(data = long_combined_data,  
              aes(x = Label, y = Area),
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9),
              color = "black", shape = 21, size = 0.8, alpha = 0.3) + 
  geom_errorbar(aes(ymin = Mean - (SE * 2), ymax = Mean + (SE * 2)),  
                width = 0.3, 
                position = position_dodge(width = 0.9), 
                linewidth = 0.8) +
  labs(title = "Cloud Size of lncRNAs per Genotype",
       x = "Genotype", 
       y = "Mean Area (µm²)") +
  coord_cartesian(ylim = c(0, 100)) +  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  facet_wrap(~ Channel, scales = "free_x") +
  scale_x_discrete(labels = function(x) gsub("_Area", " ", x)) +
  
  # 🔹 Add significance
  geom_signif(
    comparisons = comparison_list_cleaned,
    map_signif_level = TRUE,
    test = "t.test",             # use t-test per pair
    y_position = c(150, 160, 170),  # adjust for your data range
    step_increase = 0.1,
    textsize = 3.5,
    tip_length = 0.01
  )

#----------------------------------------------------------------------End here




