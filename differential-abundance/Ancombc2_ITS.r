# -------------------------------
# ANCOM-BC2 Differential Abundance Script
# Author: ChatGPT
# Purpose: Identify differentially abundant species over time and temperature
# -------------------------------

# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ANCOMBC", quietly = TRUE)) BiocManager::install("ANCOMBC")
if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("FreddieWitherden/ancombc", subdir = "ANCOMBC2")
remove.packages("ANCOMBC")
devtools::install_github("FrederickHuangLin/ANCOMBC", subdir = "ANCOMBC2")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC")

library(ANCOMBC)

library(phyloseq)
library(ANCOMBC)
library(ggplot2)
#ls("package:ANCOMBC")
?ancombc2
# -------------------------------
# Load data
# -------------------------------

# Set file paths (edit these to match your setup)
otu_file <- "/Users/meyeanni/Desktop/git_sourdough/SourdoughFlow/LP4/shipping_analysis/shipping_data/ITS/OTUs/alpha_rarefaction/core-metrics-results-1000/20250511_relative_features_OTU_species_names_combined.csv"        # samples x species (relative or counts)
metadata_file <- "/Users/meyeanni/Desktop/git_sourdough/SourdoughFlow/LP4/shipping_analysis/shipping_data/20250513_ITS_shipping_general_metadata_all_samples_with_alpha_diversity_16S_ITS.csv"    # samples x metadata (sample_id as row names)

# Load and preprocess OTU table
otu_df <- read.csv(otu_file, row.names = 1, check.names = FALSE)
otu_mat <- t(as.matrix(otu_df))  # transpose to taxa x samples

# If your data is relative abundance, rescale to pseudo-counts (optional)
otu_mat <- round(otu_mat * 1000)

# Load metadata
meta_df <- read.csv(metadata_file, row.names = 1, check.names = FALSE)

# Create phyloseq object
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
META <- sample_data(meta_df)
ps <- phyloseq(OTU, META)
# Drop temperature 20 samples
ps_filtered <- subset_samples(ps, temperature != 20)

# (Optional) Drop taxa that are now all zero
ps_filtered <- prune_taxa(taxa_sums(ps_filtered) > 0, ps_filtered)

sample_data(ps_filtered)$temperature <- factor(
  sample_data(ps_filtered)$temperature, levels = c("4", "17", "30")
)

sample_data(ps_filtered)$day <- factor(
  sample_data(ps_filtered)$day, levels = c("1", "2", "3", "4", "5", "6", "7", "14", "21", "28")
)
table(sample_data(ps_filtered)$day)
levels(sample_data(ps_filtered)$day)


# -------------------------------
# Run ANCOM-BC2
# -------------------------------
# Construct the contrast matrix for increasing trend across 10 days
# Construct 10 × 10 square contrast matrix for increasing trend
# Build 9x9 matrix: day2 - day1 ≥ 0, day3 - day2 ≥ 0, ..., day10 - day9 ≥ 0
# Construct 9x10 contrast matrix for monotonic increase: day2 ≥ day1, ..., day10 ≥ day9
# Number of levels - 1 = 9
# 10 levels of day → 9 dummy variables
C <- matrix(0, nrow = 9, ncol = 9)
for (i in 1:8) {
  C[i, i] <- -1
  C[i, i + 1] <- 1
}
C[9, 9] <- 1  # to make it square





trend_contrast <- list(
  contrast = list(C),
  node = list(1)
)
res <- ancombc2(
  data = ps_filtered,
  assay_name = "counts",             # Use "counts" even for pseudo-counted rel. abundance
  fix_formula = "day + temperature", # Main + interaction effects
  p_adj_method = "fdr",
  lib_cut = 1000,
  group = "temperature", 
  pairwise = TRUE,# Optional grouping variable (e.g., for pairwise)
  global = TRUE,
  trend = FALSE,
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  #trend_control = list(
   # contrast = trend_contrast$contrast,
    #node = trend_contrast$node,
    #solver = "ECOS",
    #B = 100
  #)
)
#for trends over time: trend = TRUE and goup ='day'

# -------------------------------
# Save and inspect results
# -------------------------------
#res_df <- res$res
#write.csv(res_df, "/Users/meyeanni/Desktop/git_sourdough/SourdoughFlow/LP4/shipping_analysis/shipping_data/ancombc2_results.csv")
write.csv(res$res, file = "/Users/meyeanni/Desktop/git_sourdough/SourdoughFlow/LP4/shipping_analysis/shipping_data/Stats/ancombc2_main_results_temperature_group_ITS.csv")
write.csv(res$res_pair, file = "/Users/meyeanni/Desktop/git_sourdough/SourdoughFlow/LP4/shipping_analysis/shipping_data/Stats/ancombc2_pairwise_results_temp_group_ITS.csv")
write.csv(res$res_global, file = "/Users/meyeanni/Desktop/git_sourdough/SourdoughFlow/LP4/shipping_analysis/shipping_data/Stats/ancombc2_global_results_temp_group_ITS.csv")
write.csv(res$samp_frac, file = "/Users/meyeanni/Desktop/git_sourdough/SourdoughFlow/LP4/shipping_analysis/shipping_data/Stats/ancombc2_sample_fractions_temp_group_ITS.csv")


# View top differentially abundant taxa for a selected term
selected_term <- "temperature30:day14"
if (selected_term %in% colnames(res_df$q_val)) {
  sig_taxa <- res_df[res_df$q_val[[selected_term]] < 0.05, ]
  print(paste("Significant taxa for", selected_term))
  print(rownames(sig_taxa))
}

# -------------------------------
# Optional: basic volcano plot
# -------------------------------
if (selected_term %in% colnames(res_df$q_val)) {
  volcano_df <- data.frame(
    Taxon = rownames(res_df),
    logFC = res_df$lfc[[selected_term]],
    q = res_df$q_val[[selected_term]]
  )

  ggplot(volcano_df, aes(x = logFC, y = -log10(q))) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = paste("Volcano plot:", selected_term),
         x = "Log Fold Change",
         y = "-log10(q-value)")
  ggsave("volcano_plot_Temp30_Day14.pdf", width = 8, height = 6)
}


#get KEGG info:
install.packages("massdatabase")
library(massdatabase)

# Install devtools if you haven't already
install.packages("devtools")
library(devtools)
install.packages(c("cli", "glue", "tibble", "httr", "xml2", "RCurl", "rvest", "stringr", "dplyr", "purrr", "remotes"))

# Use devtools to install massdatabase from GitHub
devtools::install_github("tidymass/massdatabase")

library(massdatabase)

