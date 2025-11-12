## Family-Level Phylogenetic Signal Analysis with Excel Output - DEBUGGED VERSION

### Load required packages
library(ape)
library(geiger)
library(phytools)
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)

### Define paths and traits
data_file_path <- "D:/Akash Ajay/insect bite force project/Morphology Family/family_level_data_morphology.xlsx"
tree_base_path <- "D:/Akash Ajay/insect bite force project/Morphology Family"
output_file_path <- "D:/Akash Ajay/insect bite force project/Morphology Family/phylogenetic_signal_results1.xlsx"
selected_traits <- c("iBite", "head.w", "head.h", "head.l", "th.w", "wing.l", "body.l")
families <- c( "Mantidae", "Carabidae")

### Function to check and fix zero-length branches
check_and_fix_tree <- function(tree) {
  # Check for zero-length branches
  zero_branches <- which(tree$edge.length == 0)
  if (length(zero_branches) > 0) {
    cat("  Found", length(zero_branches), "zero-length branches. Adding small value (1e-8)...\n")
    # Add a very small value to zero-length branches
    tree$edge.length[tree$edge.length == 0] <- 1e-8
  }
  
  # Check for NA branch lengths
  na_branches <- which(is.na(tree$edge.length))
  if (length(na_branches) > 0) {
    cat("  Found", length(na_branches), "NA branch lengths. Replacing with mean branch length...\n")
    mean_branch <- mean(tree$edge.length, na.rm = TRUE)
    tree$edge.length[is.na(tree$edge.length)] <- mean_branch
  }
  
  return(tree)
}

### IMPROVED Function to calculate Pagel's lambda with proper error handling
calculate_pagel_lambda <- function(tree, trait_data) {
  result <- list(lambda = NA, pvalue = NA, converged = FALSE, error = NA)
  
  tryCatch({
    # First try phytools::phylosig
    lambda_test <- phytools::phylosig(tree, trait_data, method = "lambda", test = TRUE)
    
    result$lambda <- lambda_test$lambda
    result$pvalue <- lambda_test$P
    result$converged <- TRUE
    
  }, error = function(e) {
    result$error <- paste("phytools failed:", e$message)
    
    # Alternative method using geiger with bounds
    tryCatch({
      fit <- geiger::fitContinuous(tree, trait_data, model = "lambda", 
                                   bounds = list(lambda = c(0.0001, 1)),
                                   control = list(niter = 100))
      result$lambda <- fit$opt$lambda
      result$converged <- fit$opt$convergence == 0
      if (!is.null(fit$opt$pval)) {
        result$pvalue <- fit$opt$pval
      }
    }, error = function(e2) {
      result$error <- paste(result$error, "| geiger failed:", e2$message)
    })
  })
  
  return(result)
}

### IMPROVED Function to calculate Blomberg's K with proper error handling
calculate_blomberg_k <- function(tree, trait_data) {
  result <- list(K = NA, pvalue = NA, error = NA)
  
  # Check for constant traits (zero variance)
  if (var(trait_data, na.rm = TRUE) == 0) {
    result$error <- "Zero variance in trait data"
    return(result)
  }
  
  tryCatch({
    # Use fewer simulations for stability
    k_test <- phytools::phylosig(tree, trait_data, method = "K", test = TRUE, nsim = 500)
    result$K <- k_test$K
    result$pvalue <- k_test$P
  }, error = function(e) {
    result$error <- e$message
    
    # Try alternative calculation
    tryCatch({
      # Manual calculation of Blomberg's K as backup
      C <- vcv(tree)
      C <- C[names(trait_data), names(trait_data)]
      n <- length(trait_data)
      one <- matrix(1, n, 1)
      invC <- solve(C)
      
      # Calculate MSE for phylogenetic and non-phylogenetic cases
      MSE0 <- as.numeric((t(trait_data) %*% (diag(n) - one %*% solve(t(one) %*% one) %*% t(one)) %*% trait_data) / (n - 1))
      MSE <- as.numeric((t(trait_data) %*% (invC - invC %*% one %*% solve(t(one) %*% invC %*% one) %*% t(one) %*% invC) %*% trait_data) / (n - 1))
      
      result$K <- MSE0 / MSE
      result$pvalue <- NA  # Can't calculate p-value manually easily
    }, error = function(e2) {
      result$error <- paste(result$error, "| Manual calc failed:", e2$message)
    })
  })
  
  return(result)
}

### Step 1: Initialize results data frame
results <- data.frame(
  Family = character(),
  Trait = character(),
  Lambda = numeric(),
  Lambda_P = numeric(),
  K = numeric(),
  K_P = numeric(),
  Sample_Size = numeric(),
  Notes = character(),
  stringsAsFactors = FALSE
)

### Step 2: Check which families have tree files first
cat("Checking available tree files...\n")
available_families <- c()
for (family in families) {
  tree_file_path <- file.path(tree_base_path, family, paste0(family, ".nwk"))
  if (file.exists(tree_file_path)) {
    available_families <- c(available_families, family)
    cat("  Found tree for:", family, "\n")
  } else {
    cat("  MISSING tree for:", family, "-", tree_file_path, "\n")
  }
}

cat("\nFamilies with available trees:", paste(available_families, collapse = ", "), "\n\n")

### Step 3: Loop through each available family
for (family in available_families) {
  cat("\n", rep("=", 50), "\n", sep = "")
  cat("Starting analysis for", family, "\n")
  cat(rep("=", 50), "\n", sep = "")
  
  # Read morphology data from family-specific sheet
  family_data <- tryCatch({
    data <- read_excel(data_file_path, sheet = family)
    # Check for required columns
    if (!"ID" %in% colnames(data)) {
      cat("ERROR: ID column not found in", family, "data\n")
      next
    }
    data
  }, error = function(e) {
    cat("ERROR: Failed to read data from sheet", family, ":", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(family_data)) next
  
  n_species <- nrow(family_data)
  cat("Number of species in data:", n_species, "\n")
  
  # Read and prepare tree
  tree_file_path <- file.path(tree_base_path, family, paste0(family, ".nwk"))
  tree <- tryCatch({
    t <- read.tree(tree_file_path)
    # Check and fix tree structure
    t <- check_and_fix_tree(t)
    t
  }, error = function(e) {
    cat("ERROR: Failed to read tree file:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(tree)) next
  
  # Match species between tree and data
  species_in_tree <- tree$tip.label
  species_in_data <- family_data$ID
  common_species <- intersect(species_in_tree, species_in_data)
  
  if (length(common_species) < 5) {
    cat("WARNING: Only", length(common_species), "common species found (minimum 5 required). Skipping family.\n")
    next
  }
  
  # Prune tree and data to common species
  tree_pruned <- drop.tip(tree, setdiff(species_in_tree, common_species))
  family_data_pruned <- family_data %>% filter(ID %in% common_species)
  
  cat("Number of common species after pruning:", length(common_species), "\n")
  cat("Tree tips:", length(tree_pruned$tip.label), "\n")
  
  # Loop through each trait
  for (trait in selected_traits) {
    cat("\n--- Testing trait:", trait, "---\n")
    
    # Check if trait exists and has data
    if (!(trait %in% colnames(family_data_pruned))) {
      cat("  Trait not found in data\n")
      results <- rbind(results, data.frame(
        Family = family,
        Trait = trait,
        Lambda = NA,
        Lambda_P = NA,
        K = NA,
        K_P = NA,
        Sample_Size = 0,
        Notes = "Trait not found",
        stringsAsFactors = FALSE
      ))
      next
    }
    
    # Extract and validate trait data
    trait_data <- family_data_pruned[[trait]]
    names(trait_data) <- family_data_pruned$ID
    
    # Remove NA values and check for sufficient data
    valid_data <- !is.na(trait_data) & is.finite(trait_data)
    valid_species <- names(trait_data)[valid_data]
    trait_data_clean <- trait_data[valid_data]
    
    if (length(valid_species) < 5) {
      cat("  Insufficient data: only", length(valid_species), "valid observations (minimum 5 required)\n")
      results <- rbind(results, data.frame(
        Family = family,
        Trait = trait,
        Lambda = NA,
        Lambda_P = NA,
        K = NA,
        K_P = NA,
        Sample_Size = length(valid_species),
        Notes = paste("Insufficient data (n =", length(valid_species), ")"),
        stringsAsFactors = FALSE
      ))
      next
    }
    
    cat("  Valid data points:", length(valid_species), "\n")
    
    # Further prune tree for this specific trait
    tree_trait <- drop.tip(tree_pruned, setdiff(tree_pruned$tip.label, valid_species))
    
    # Ensure tree and data are in the same order
    trait_data_final <- trait_data_clean[tree_trait$tip.label]
    
    # Debug info
    cat("  Tree has", length(tree_trait$tip.label), "tips\n")
    cat("  Trait variance:", var(trait_data_final), "\n")
    
    # Calculate phylogenetic signals
    lambda_result <- calculate_pagel_lambda(tree_trait, trait_data_final)
    k_result <- calculate_blomberg_k(tree_trait, trait_data_final)
    
    # Prepare notes
    notes <- character()
    if (!is.na(lambda_result$error)) {
      notes <- c(notes, paste("Lambda:", lambda_result$error))
    }
    if (!is.na(k_result$error)) {
      notes <- c(notes, paste("K:", k_result$error))
    }
    if (length(valid_species) < 10) {
      notes <- c(notes, "Small sample size")
    }
    if (var(trait_data_final) == 0) {
      notes <- c(notes, "Zero variance in trait")
    }
    
    notes_text <- if (length(notes) > 0) paste(notes, collapse = "; ") else ""
    
    # SAFELY Store results with explicit values
    new_row <- data.frame(
      Family = family,
      Trait = trait,
      Lambda = ifelse(!is.na(lambda_result$lambda), as.numeric(lambda_result$lambda), NA),
      Lambda_P = ifelse(!is.na(lambda_result$pvalue), as.numeric(lambda_result$pvalue), NA),
      K = ifelse(!is.na(k_result$K), as.numeric(k_result$K), NA),
      K_P = ifelse(!is.na(k_result$pvalue), as.numeric(k_result$pvalue), NA),
      Sample_Size = as.integer(length(valid_species)),
      Notes = notes_text,
      stringsAsFactors = FALSE
    )
    
    # Ensure no empty rows
    if (nrow(new_row) == 1) {
      results <- rbind(results, new_row)
    } else {
      cat("  WARNING: Could not create result row for", trait, "\n")
    }
    
    # Print results
    cat("  Pagel's Lambda:", 
        ifelse(!is.na(lambda_result$lambda), round(lambda_result$lambda, 4), "NA"),
        "p-value:", 
        ifelse(!is.na(lambda_result$pvalue), round(lambda_result$pvalue, 4), "NA"), "\n")
    cat("  Blomberg's K:", 
        ifelse(!is.na(k_result$K), round(k_result$K, 4), "NA"),
        "p-value:", 
        ifelse(!is.na(k_result$pvalue), round(k_result$pvalue, 4), "NA"), "\n")
    if (notes_text != "") {
      cat("  Notes:", notes_text, "\n")
    }
  }
  
  # Save family-specific results
  family_results_file <- file.path(tree_base_path, family, paste0(family, "_phylogenetic_signal_results.rds"))
  family_results <- results %>% filter(Family == family)
  if (nrow(family_results) > 0) {
    saveRDS(family_results, family_results_file)
    cat("\nFamily results saved:", family_results_file, "\n")
  }
}

### Step 4: Write combined results to Excel
if (nrow(results) > 0) {
  # Create a workbook with multiple sheets
  wb <- createWorkbook()
  
  # Main results sheet
  addWorksheet(wb, "Phylogenetic_Signal")
  writeData(wb, "Phylogenetic_Signal", results)
  
  # Summary sheet
  summary_data <- results %>%
    group_by(Family) %>%
    summarise(
      n_traits = n(),
      n_significant_K = sum(K_P < 0.05, na.rm = TRUE),
      n_significant_lambda = sum(Lambda_P < 0.05, na.rm = TRUE),
      mean_sample_size = round(mean(Sample_Size, na.rm = TRUE), 1)
    )
  
  addWorksheet(wb, "Summary")
  writeData(wb, "Summary", summary_data)
  
  # Save workbook
  saveWorkbook(wb, output_file_path, overwrite = TRUE)
  
  cat("\n", rep("=", 50), "\n", sep = "")
  cat("ANALYSIS COMPLETED SUCCESSFULLY!\n")
  cat("Output written to:", output_file_path, "\n")
  cat("Families analyzed:", length(unique(results$Family)), "\n")
  cat("Total trait analyses:", nrow(results), "\n")
  cat(rep("=", 50), "\n", sep = "")
  
  # Print final summary
  print(summary_data)
  
} else {
  cat("\nERROR: No results were generated. Check your data and tree files.\n")
}

### Show any warnings
if (length(warnings()) > 0) {
  cat("\nWarnings encountered:\n")
  print(warnings()[1:min(10, length(warnings()))])  # Show first 10 warnings
}