# =============================================================================
# TRAIT-SPECIFIC ANALYSIS: Evolutionary Models vs Species Richness
# Accounting for sample size effects and trait-specific patterns
# =============================================================================

# Load libraries
library(readxl)
library(dplyr)
library(openxlsx)

# 1. READ DATA
# =============================================================================
data_path <- "D:/Akash Ajay/insect bite force project/genome features/model_fits_output_with_EB_WN/all_aic_values_with_EB_WN.xlsx"
genomic_data <- read_excel(data_path)

# Remove skipped models
clean_data <- genomic_data[!genomic_data$Best_Model %in% c("Skipped", "NA", NA), ]

# Define species richness
species_rich <- c("Coleoptera", "Diptera", "Hymenoptera", "Lepidoptera", "Hemiptera")
species_poor <- c("Blattodea", "Odonata", "Orthoptera", "Phasmatodea", "Psocodea", "Trichoptera")
clean_data$Richness_Category <- ifelse(clean_data$Order %in% species_rich, "Rich", "Poor")

# 2. SAMPLE SIZE ANALYSIS BY TRAIT
# =============================================================================
cat("=== SAMPLE SIZE ANALYSIS BY TRAIT ===\n")

sample_size_by_trait <- clean_data %>%
  group_by(Trait, Richness_Category) %>%
  summarise(
    n_orders = n_distinct(Order),
    n_observations = n(),
    .groups = 'drop'
  )

print(sample_size_by_trait)

# Identify traits with adequate sample size (≥3 orders per category)
adequate_traits <- sample_size_by_trait %>%
  group_by(Trait) %>%
  filter(all(n_orders >= 3)) %>%
  pull(Trait) %>%
  unique()

cat("\nTraits with adequate sample size (≥3 orders per category):", paste(adequate_traits, collapse = ", "), "\n\n")

# 3. TRAIT-SPECIFIC ANALYSIS
# =============================================================================
cat("=== TRAIT-SPECIFIC MODEL DISTRIBUTIONS ===\n")

trait_results <- list()

for(trait in adequate_traits) {
  trait_data <- clean_data[clean_data$Trait == trait, ]
  
  if(nrow(trait_data) >= 6) { # Minimum threshold for meaningful analysis
    model_table <- table(trait_data$Richness_Category, trait_data$Best_Model)
    
    # Only run Fisher test if table has sufficient data
    if(all(dim(model_table) >= 2) && min(rowSums(model_table)) >= 3) {
      fisher_test <- tryCatch(
        fisher.test(model_table),
        error = function(e) list(p.value = NA)
      )
      
      trait_results[[trait]] <- list(
        table = model_table,
        p_value = fisher_test$p.value,
        rich_total = sum(model_table["Rich", ]),
        poor_total = sum(model_table["Poor", ])
      )
      
      cat("\n---", trait, "---\n")
      print(model_table)
      cat("Fisher p-value:", round(fisher_test$p.value, 5), "\n")
    }
  }
}

# 4. FOCUS ON HIGH-QUALITY TRAITS (excluding potential White Noise artifacts)
# =============================================================================
cat("\n=== CONSERVATIVE ANALYSIS: EXCLUDING WHITE NOISE ===\n")

# Remove WN models to reduce sample size artifacts
conservative_data <- clean_data[clean_data$Best_Model != "WN", ]

conservative_results <- list()
for(trait in adequate_traits) {
  trait_data <- conservative_data[conservative_data$Trait == trait, ]
  
  if(nrow(trait_data) >= 6) {
    model_table <- table(trait_data$Richness_Category, trait_data$Best_Model)
    
    if(all(dim(model_table) >= 2) && min(rowSums(model_table)) >= 3) {
      fisher_test <- tryCatch(
        fisher.test(model_table),
        error = function(e) list(p.value = NA)
      )
      
      conservative_results[[trait]] <- list(
        table = model_table,
        p_value = fisher_test$p.value
      )
      
      cat("\n---", trait, "(WN excluded) ---\n")
      print(model_table)
      cat("Fisher p-value:", round(fisher_test$p.value, 5), "\n")
    }
  }
}

# 5. KAPPA PATTERN BY TRAIT
# =============================================================================
cat("\n=== KAPPA MODEL PREVALENCE BY TRAIT ===\n")

kappa_analysis <- clean_data %>%
  group_by(Trait, Richness_Category) %>%
  summarise(
    total_models = n(),
    kappa_count = sum(Best_Model == "kappa"),
    kappa_proportion = kappa_count / total_models,
    .groups = 'drop'
  )

print(kappa_analysis)

# 6. SAVE DETAILED RESULTS
# =============================================================================
results_wb <- createWorkbook()

# Sample size summary
addWorksheet(results_wb, "Sample_Size")
writeData(results_wb, "Sample_Size", sample_size_by_trait)

# Trait-specific results
addWorksheet(results_wb, "Trait_Analysis")
trait_summary <- data.frame()
for(trait in names(trait_results)) {
  trait_summary <- rbind(trait_summary, data.frame(
    Trait = trait,
    P_Value = trait_results[[trait]]$p_value,
    Rich_Observations = trait_results[[trait]]$rich_total,
    Poor_Observations = trait_results[[trait]]$poor_total,
    Significance = ifelse(trait_results[[trait]]$p_value < 0.05, "SIGNIFICANT", "Non-significant")
  ))
}
writeData(results_wb, "Trait_Analysis", trait_summary)

# Kappa analysis
addWorksheet(results_wb, "Kappa_by_Trait")
writeData(results_wb, "Kappa_by_Trait", kappa_analysis)

# Conservative analysis (WN excluded)
addWorksheet(results_wb, "Conservative_Analysis")
conservative_summary <- data.frame()
for(trait in names(conservative_results)) {
  conservative_summary <- rbind(conservative_summary, data.frame(
    Trait = trait,
    P_Value = conservative_results[[trait]]$p_value,
    Significance = ifelse(conservative_results[[trait]]$p_value < 0.05, "SIGNIFICANT", "Non-significant")
  ))
}
writeData(results_wb, "Conservative_Analysis", conservative_summary)

saveWorkbook(results_wb, "Trait_Specific_Evolutionary_Analysis.xlsx", overwrite = TRUE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Detailed results saved to: Trait_Specific_Evolutionary_Analysis.xlsx\n")
cat("Focus on traits with adequate sample size for robust conclusions.\n")