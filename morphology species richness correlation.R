## Correlation Analysis of Phylogenetic Signal and Species Richness (Spearman and Pearson)

### Load required packages
library(dplyr)
library(openxlsx)
library(stats)

## Step 1: Load phylogenetic signal data from Excel
phylo_signal_file <- "D:/Akash Ajay/insect bite force project/Morphology/species richness/phylogenetic signal.xlsx"
phylo_signal_data <- read.xlsx(phylo_signal_file, sheet = 1)

## Step 2: Rename λ column to Lambda for consistency and coerce to numeric
colnames(phylo_signal_data)[colnames(phylo_signal_data) == "λ"] <- "Lambda"
phylo_signal_data$Lambda <- as.numeric(as.character(phylo_signal_data$Lambda))
phylo_signal_data$K <- as.numeric(as.character(phylo_signal_data$K))

## Step 3: Load species richness data from Excel
species_richness_file <- "D:/Akash Ajay/insect bite force project/Morphology/species richness/richness.xlsx"
species_richness_data <- read.xlsx(species_richness_file, sheet = 1)

## Step 4: Rename species richness column for consistency
colnames(species_richness_data)[colnames(species_richness_data) == "No.of.species"] <- "No_of_species"

## Step 5: Data quality checks before merging
cat("=== Data Quality Checks ===\n")
cat("Phylogenetic signal data - Unique orders:", length(unique(phylo_signal_data$Order)), "\n")
cat("Species richness data - Unique orders:", length(unique(species_richness_data$Order)), "\n")

# Check for duplicate Order+Trait combinations in phylogenetic data
duplicates <- phylo_signal_data %>%
  group_by(Order, Trait) %>%
  filter(n() > 1) %>%
  ungroup()

if(nrow(duplicates) > 0) {
  cat("Warning: Found", nrow(duplicates), "duplicate Order+Trait combinations:\n")
  print(duplicates %>% select(Order, Trait) %>% distinct())
  cat("These will be treated as separate data points in correlation analysis.\n")
}

## Step 6: Merge phylogenetic signal and species richness data
merged_data <- merge(phylo_signal_data, species_richness_data, by = "Order", all.x = FALSE)

cat("After merging - Total observations:", nrow(merged_data), "\n")
cat("Orders in merged data:", length(unique(merged_data$Order)), "\n")

## Step 7: Perform Spearman and Pearson correlations for each trait
traits <- unique(merged_data$Trait)
correlation_results <- data.frame(
  Trait = character(),
  N_Orders = numeric(),
  Spearman_Lambda_Cor = numeric(),
  Spearman_Lambda_P = numeric(),
  Spearman_K_Cor = numeric(),
  Spearman_K_P = numeric(),
  Pearson_Lambda_Cor = numeric(),
  Pearson_Lambda_P = numeric(),
  Pearson_K_Cor = numeric(),
  Pearson_K_P = numeric(),
  stringsAsFactors = FALSE
)

for (trait in traits) {
  cat("\n=== Correlation Analysis for Trait:", trait, "===\n")
  
  # Subset data for the current trait
  trait_data <- merged_data %>% filter(Trait == trait)
  
  # Remove rows with NA in Lambda, K, or No_of_species
  trait_data <- trait_data %>% filter(!is.na(Lambda) & !is.na(K) & !is.na(No_of_species))
  
  # Check if enough data points for correlation
  if (nrow(trait_data) < 3) {
    cat("Insufficient data for trait", trait, " (", nrow(trait_data), " orders). Skipping.\n")
    correlation_results <- rbind(correlation_results, data.frame(
      Trait = trait,
      N_Orders = nrow(trait_data),
      Spearman_Lambda_Cor = NA,
      Spearman_Lambda_P = NA,
      Spearman_K_Cor = NA,
      Spearman_K_P = NA,
      Pearson_Lambda_Cor = NA,
      Pearson_Lambda_P = NA,
      Pearson_K_Cor = NA,
      Pearson_K_P = NA
    ))
    next
  }
  
  cat("Analyzing", nrow(trait_data), "data points from", length(unique(trait_data$Order)), "orders\n")
  
  # Spearman correlation for Lambda
  spearman_lambda_cor <- tryCatch({
    cor.test(trait_data$Lambda, trait_data$No_of_species, method = "spearman", exact = FALSE)
  }, warning = function(w) {
    cat("Warning in Spearman Lambda correlation for", trait, ":", w$message, "\n")
    list(estimate = NA, p.value = NA)
  }, error = function(e) {
    cat("Error in Spearman Lambda correlation for", trait, ":", e$message, "\n")
    list(estimate = NA, p.value = NA)
  })
  
  # Spearman correlation for K
  spearman_k_cor <- tryCatch({
    cor.test(trait_data$K, trait_data$No_of_species, method = "spearman", exact = FALSE)
  }, warning = function(w) {
    cat("Warning in Spearman K correlation for", trait, ":", w$message, "\n")
    list(estimate = NA, p.value = NA)
  }, error = function(e) {
    cat("Error in Spearman K correlation for", trait, ":", e$message, "\n")
    list(estimate = NA, p.value = NA)
  })
  
  # Pearson correlation for Lambda
  pearson_lambda_cor <- tryCatch({
    cor.test(trait_data$Lambda, trait_data$No_of_species, method = "pearson")
  }, warning = function(w) {
    cat("Warning in Pearson Lambda correlation for", trait, ":", w$message, "\n")
    list(estimate = NA, p.value = NA)
  }, error = function(e) {
    cat("Error in Pearson Lambda correlation for", trait, ":", e$message, "\n")
    list(estimate = NA, p.value = NA)
  })
  
  # Pearson correlation for K
  pearson_k_cor <- tryCatch({
    cor.test(trait_data$K, trait_data$No_of_species, method = "pearson")
  }, warning = function(w) {
    cat("Warning in Pearson K correlation for", trait, ":", w$message, "\n")
    list(estimate = NA, p.value = NA)
  }, error = function(e) {
    cat("Error in Pearson K correlation for", trait, ":", e$message, "\n")
    list(estimate = NA, p.value = NA)
  })
  
  # Store results
  correlation_results <- rbind(correlation_results, data.frame(
    Trait = trait,
    N_Orders = nrow(trait_data),
    Spearman_Lambda_Cor = if (!is.na(spearman_lambda_cor$estimate)) spearman_lambda_cor$estimate else NA,
    Spearman_Lambda_P = if (!is.na(spearman_lambda_cor$p.value)) spearman_lambda_cor$p.value else NA,
    Spearman_K_Cor = if (!is.na(spearman_k_cor$estimate)) spearman_k_cor$estimate else NA,
    Spearman_K_P = if (!is.na(spearman_k_cor$p.value)) spearman_k_cor$p.value else NA,
    Pearson_Lambda_Cor = if (!is.na(pearson_lambda_cor$estimate)) pearson_lambda_cor$estimate else NA,
    Pearson_Lambda_P = if (!is.na(pearson_lambda_cor$p.value)) pearson_lambda_cor$p.value else NA,
    Pearson_K_Cor = if (!is.na(pearson_k_cor$estimate)) pearson_k_cor$estimate else NA,
    Pearson_K_P = if (!is.na(pearson_k_cor$p.value)) pearson_k_cor$p.value else NA
  ))
  
  # Print results
  cat("Spearman Lambda Correlation:", if (!is.na(spearman_lambda_cor$estimate)) round(spearman_lambda_cor$estimate, 4), 
      "p-value:", if (!is.na(spearman_lambda_cor$p.value)) round(spearman_lambda_cor$p.value, 4), "\n")
  cat("Spearman K Correlation:", if (!is.na(spearman_k_cor$estimate)) round(spearman_k_cor$estimate, 4), 
      "p-value:", if (!is.na(spearman_k_cor$p.value)) round(spearman_k_cor$p.value, 4), "\n")
  cat("Pearson Lambda Correlation:", if (!is.na(pearson_lambda_cor$estimate)) round(pearson_lambda_cor$estimate, 4), 
      "p-value:", if (!is.na(pearson_lambda_cor$p.value)) round(pearson_lambda_cor$p.value, 4), "\n")
  cat("Pearson K Correlation:", if (!is.na(pearson_k_cor$estimate)) round(pearson_k_cor$estimate, 4), 
      "p-value:", if (!is.na(pearson_k_cor$p.value)) round(pearson_k_cor$p.value, 4), "\n")
}

## Step 8: Write results to Excel
output_file <- "D:/Akash Ajay/insect bite force project/Morphology/species richness/order_corr.xlsx"
write.xlsx(correlation_results, output_file, sheetName = "Correlation_Results", rowNames = FALSE)

## Step 9: Print final summary
cat("\n===== Correlation Analysis Summary =====\n")
cat("Analyzed", length(traits), "traits across", length(unique(merged_data$Order)), "orders\n")
cat("Total observations in merged data:", nrow(merged_data), "\n")
cat("Output written to:", output_file, "\n")
print(correlation_results)