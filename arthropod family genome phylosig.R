
## Family-Level Phylogenetic Signal Analysis for Genomic Traits

### Load required packages
library(ape)
library(phytools)
library(readxl)
library(openxlsx)
library(dplyr)

### Define file paths and traits
data_file_path <- "D:/Akash Ajay/insect bite force project/genome features/family/output_by_family.xlsx"
tree_base_path <- "D:/Akash Ajay/insect bite force project/genome features/family"
output_file_path <- "D:/Akash Ajay/insect bite force project/genome features/family/phylogenetic_signal_summary.xlsx"
selected_traits <- c("Size_Mb", "Size_bp", "GC_Content", "Chromosome_Number", "Coding_genes", "Non_coding_genes")

### Step 1: Get family names from Excel sheets
sheet_names <- excel_sheets(data_file_path)
cat("Families found in data file:", paste(sheet_names, collapse = ", "), "\n")

### Step 2: Initialize results data frame
results_df <- data.frame(
  Family = character(),
  Trait = character(),
  Valid_Species = integer(),
  Lambda = numeric(),
  P_lambda = numeric(),
  K = numeric(),
  P_K = numeric(),
  stringsAsFactors = FALSE
)

### Step 3: Function to process each family
process_family <- function(family_name, family_data, tree) {
  cat("\n\n===== Processing", family_name, "=====\n")
  
  ## Step 1: Prepare data
  family_data <- as.data.frame(family_data)
  
  # Check if Species column exists
  if (!"Species" %in% colnames(family_data)) {
    cat("Error: 'Species' column missing in", family_name, "data\n")
    cat("Available columns:", colnames(family_data), "\n")
    return(NULL)
  }
  
  # Check for NA values in Species
  if (any(is.na(family_data$Species))) {
    cat("Warning: NA values found in 'Species' for", family_name, "\n")
    family_data <- family_data[!is.na(family_data$Species), , drop = FALSE]
    if (nrow(family_data) == 0) {
      cat("Error: No valid rows remaining after removing NAs in", family_name, "\n")
      return(NULL)
    }
  }
  
  # Create proper row names by replacing spaces with underscores in Species
  tryCatch({
    family_data$Species <- trimws(as.character(family_data$Species))
    rownames(family_data) <- gsub(" ", "_", family_data$Species)
  }, error = function(e) {
    cat("Error setting row names for", family_name, ":", e$message, "\n")
    return(NULL)
  })
  
  # Check for selected traits
  missing_traits <- setdiff(selected_traits, colnames(family_data))
  if (length(missing_traits) > 0) {
    cat("Warning: Missing traits in", family_name, ":", missing_traits, "\n")
  }
  
  ## Step 2: Prepare tree
  tree <- tryCatch({
    multi2di(tree)  # Resolve polytomies
  }, error = function(e) {
    cat("Error resolving polytomies for", family_name, ":", e$message, "\n")
    return(NULL)
  })
  if (is.null(tree)) return(NULL)
  
  tree$node.label <- NULL
  min_branch <- 1e-6 * max(node.depth.edgelength(tree))
  tree$edge.length[tree$edge.length < min_branch] <- min_branch
  
  ## Step 3: Match data and tree
  common_species <- intersect(rownames(family_data), tree$tip.label)
  cat("Number of species in analysis:", length(common_species), "\n")
  
  if (length(common_species) < 10) {
    cat("Error: Fewer than 10 common species between tree and data for", family_name, ". Skipping analysis.\n")
    for (trait in selected_traits) {
      results_df <<- rbind(results_df, data.frame(
        Family = family_name,
        Trait = trait,
        Valid_Species = NA_integer_,
        Lambda = NA_real_,
        P_lambda = NA_real_,
        K = NA_real_,
        P_K = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    return(NULL)
  }
  
  # Subset data and tree to common species
  trait_data <- family_data[common_species, , drop = FALSE]
  tree <- keep.tip(tree, common_species)
  
  ## Step 4: Phylogenetic signal tests for each trait
  results <- list()
  analyzed_traits <- character()
  
  for (trait in selected_traits) {
    cat("\n=== Testing trait:", trait, "===\n")
    
    # Check if trait exists in data
    if (!trait %in% colnames(trait_data)) {
      cat("Trait", trait, "not found in data for", family_name, "\n")
      results_df <<- rbind(results_df, data.frame(
        Family = family_name,
        Trait = trait,
        Valid_Species = 0L,
        Lambda = NA_real_,
        P_lambda = NA_real_,
        K = NA_real_,
        P_K = NA_real_,
        stringsAsFactors = FALSE
      ))
      next
    }
    
    # Filter out NA values for this trait
    valid_data_trait <- trait_data[!is.na(trait_data[[trait]]) & is.numeric(trait_data[[trait]]), , drop = FALSE]
    valid_species_trait <- rownames(valid_data_trait)
    cat("Number of valid species for", trait, ":", length(valid_species_trait), "\n")
    
    # Add row with valid species count
    row_idx <- nrow(results_df) + 1
    results_df <<- rbind(results_df, data.frame(
      Family = family_name,
      Trait = trait,
      Valid_Species = length(valid_species_trait),
      Lambda = NA_real_,
      P_lambda = NA_real_,
      K = NA_real_,
      P_K = NA_real_,
      stringsAsFactors = FALSE
    ))
    
    # Check minimum sample size
    if (length(valid_species_trait) < 10) {
      cat("Skipping", trait, "due to insufficient valid species (< 10)\n")
      next
    }
    
    # Subset tree to valid species for this trait
    valid_tree <- keep.tip(tree, valid_species_trait)
    current_trait <- setNames(valid_data_trait[[trait]], valid_species_trait)
    
    # Remove Inf/-Inf values
    finite_trait <- current_trait[is.finite(current_trait)]
    if (length(finite_trait) < 10) {
      cat("Skipping", trait, "due to insufficient finite values (< 10)\n")
      results_df$Valid_Species[row_idx] <<- length(finite_trait)
      next
    }
    
    # Pagel's lambda
    lambda <- tryCatch({
      phylosig(valid_tree, finite_trait, method = "lambda", test = TRUE, nsim = 999)
    }, error = function(e) {
      cat("Error computing Pagel's lambda for", trait, ":", e$message, "\n")
      return(NULL)
    })
    if (!is.null(lambda)) {
      cat("Lambda:", round(lambda$lambda, 6), "p-value:", lambda$P, "\n")
      results_df$Lambda[row_idx] <<- lambda$lambda
      results_df$P_lambda[row_idx] <<- lambda$P
    }
    
    # Blomberg's K
    K <- tryCatch({
      phylosig(valid_tree, finite_trait, method = "K", test = TRUE, nsim = 999)
    }, error = function(e) {
      cat("Error computing Blomberg's K for", trait, ":", e$message, "\n")
      return(NULL)
    })
    if (!is.null(K)) {
      cat("K:", round(K$K, 6), "p-value:", K$P, "\n")
      results_df$K[row_idx] <<- K$K
      results_df$P_K[row_idx] <<- K$P
    }
    
    # Store results
    results[[trait]] <- list(lambda = lambda, K = K)
    analyzed_traits <- c(analyzed_traits, trait)
    
    # Plot trait evolution
    ylim_range <- range(finite_trait, na.rm = TRUE)
    ylim <- ylim_range + c(-0.05, 0.05) * diff(ylim_range)
    
    tryCatch({
      png(file.path(tree_base_path, family_name, paste0(family_name, "_", gsub("[^A-Za-z0-9]", "_", trait), "_trait_evolution.png")), 
          width = 10, height = 6, units = "in", res = 300)
      phenogram(valid_tree, finite_trait, main = paste(trait, "trait evolution in", family_name), ylim = ylim)
      dev.off()
    }, error = function(e) {
      cat("Error plotting phenogram for", trait, "in", family_name, ":", e$message, "\n")
      if (dev.cur() > 1) dev.off()
    })
  }
  
  # Save family-specific results
  saveRDS(results, file.path(tree_base_path, family_name, paste0(family_name, "_phylogenetic_signal_results.rds")))
  
  # Print summary
  cat("\n===== ", family_name, " Analysis Summary =====\n")
  cat("Analyzed", length(analyzed_traits), "traits across", length(common_species), "species\n")
  cat("Output files created:\n")
  cat("- ", family_name, "_phylogenetic_signal_results.rds\n")
  cat("- Trait-specific evolution plots (PNG)\n")
}

### Step 4: Loop through families
for (family in sheet_names) {
  # Read family data
  family_data <- tryCatch({
    read_excel(data_file_path, sheet = family)
  }, error = function(e) {
    cat("ERROR: Failed to read data from sheet", family, ":", e$message, "\n")
    next
  })
  
  # Clean column names to match selected_traits
  colnames(family_data) <- gsub(" ", "_", colnames(family_data))
  
  # Read family tree
  tree_file_path <- file.path(tree_base_path, family, paste0(family, "_species.nwk"))
  tree <- tryCatch({
    read.tree(tree_file_path)
  }, error = function(e) {
    cat("ERROR: Failed to read tree file for", family, ":", e$message, "\n")
    return(NULL)
  })
  if (is.null(tree)) next
  
  # Process the family
  process_family(family, family_data, tree)
}

### Step 5: Write results to Excel
write.xlsx(results_df, output_file_path, sheetName = "Phylogenetic_Signal", rowNames = FALSE)

### Step 6: Print final summary
cat("\n===== Phylogenetic Signal Analysis Summary =====\n")
cat("Analyzed", length(selected_traits), "traits across", length(sheet_names), "families\n")
cat("Results summary saved to:", output_file_path, "\n")
print(results_df)