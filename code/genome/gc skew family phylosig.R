## Family-Level Phylogenetic Signal Analysis for GC Skew with Duplicate Species Handling

### Load required packages
library(ape)
library(phytools)
library(readxl)
library(openxlsx)
library(dplyr)

### Define file paths and traits
original_data_path <- "D:/Akash Ajay/insect bite force project/genome features/family/output_by_family.xlsx"
gc_skew_data_path <- "D:/Akash Ajay/insect bite force project/sravya work/GC skew/gc_skew_family_wise.xlsx"
tree_base_path <- "D:/Akash Ajay/insect bite force project/genome features/family"
output_file_path <- "D:/Akash Ajay/insect bite force project/sravya work/GC skew/phylogenetic_signal_gc_skew_family.xlsx"
selected_traits <- c("GC_skew")

### Step 1: Get family names from GC skew Excel sheets
sheet_names <- excel_sheets(gc_skew_data_path)
cat("Families found in GC skew data file:", paste(sheet_names, collapse = ", "), "\n")

### Step 2: Initialize results data frames
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

unmatched_df <- data.frame(
  Family = character(),
  Accession_Number = character(),
  Source = character(),
  stringsAsFactors = FALSE
)

### Step 3: Function to process each family
process_family <- function(family_name, original_data, gc_skew_data, tree) {
  cat("\n\n===== Processing", family_name, "=====\n")
  
  ## Step 1: Prepare original data
  original_data <- as.data.frame(original_data)
  
  # Check if Species and Accession_Number columns exist
  if (!"Species" %in% colnames(original_data)) {
    cat("Error: 'Species' column missing in", family_name, "original data\n")
    cat("Available columns:", colnames(original_data), "\n")
    return(NULL)
  }
  if (!"Accession_Number" %in% colnames(original_data)) {
    cat("Error: 'Accession_Number' column missing in", family_name, "original data\n")
    cat("Available columns:", colnames(original_data), "\n")
    return(NULL)
  }
  
  # Check for NA values in Species
  if (any(is.na(original_data$Species))) {
    cat("Warning: NA values found in 'Species' for", family_name, "original data\n")
    original_data <- original_data[!is.na(original_data$Species), , drop = FALSE]
    if (nrow(original_data) == 0) {
      cat("Error: No valid rows remaining after removing NAs in", family_name, "\n")
      return(NULL)
    }
  }
  
  ## Step 2: Prepare GC_skew data and aggregate duplicates
  gc_skew_data <- as.data.frame(gc_skew_data)
  if ("Accession Number" %in% colnames(gc_skew_data) && "GC_skew" %in% colnames(gc_skew_data)) {
    gc_skew_data$GC_skew <- as.numeric(as.character(gc_skew_data$GC_skew))
    if (any(is.na(gc_skew_data$GC_skew) & !is.na(gc_skew_data$`GC_skew`))) {
      cat("Warning: NAs introduced by coercion for GC_skew in", family_name, "\n")
    }
    colnames(gc_skew_data) <- gsub(" ", "_", colnames(gc_skew_data))
    
    # Merge with original data
    merged_data <- merge(original_data[, c("Accession_Number", "Species")], 
                         gc_skew_data[, c("Accession_Number", "GC_skew")], 
                         by = "Accession_Number", all.x = FALSE)
    
    # Identify unmatched Accession Numbers
    unmatched_data1 <- original_data$Accession_Number[!original_data$Accession_Number %in% gc_skew_data$Accession_Number]
    unmatched_data2 <- gc_skew_data$Accession_Number[!gc_skew_data$Accession_Number %in% original_data$Accession_Number]
    if (length(unmatched_data1) > 0) {
      unmatched_df <<- rbind(unmatched_df, data.frame(
        Family = family_name,
        Accession_Number = unmatched_data1,
        Source = paste0(family_name, "_original_data"),
        stringsAsFactors = FALSE
      ))
    }
    if (length(unmatched_data2) > 0) {
      unmatched_df <<- rbind(unmatched_df, data.frame(
        Family = family_name,
        Accession_Number = unmatched_data2,
        Source = "gc_skew_data",
        stringsAsFactors = FALSE
      ))
    }
    
    if (nrow(merged_data) == 0) {
      cat("Warning: No matching accessions for GC_skew in", family_name, "\n")
      results_df <<- rbind(results_df, data.frame(
        Family = family_name,
        Trait = "GC_skew",
        Valid_Species = 0L,
        Lambda = NA_real_,
        P_lambda = NA_real_,
        K = NA_real_,
        P_K = NA_real_,
        stringsAsFactors = FALSE
      ))
      return(NULL)
    }
    
    # Aggregate duplicate species by averaging GC_skew
    merged_data <- merged_data %>%
      group_by(Species) %>%
      summarise(GC_skew = mean(GC_skew, na.rm = TRUE), .groups = "drop") %>%
      as.data.frame()
    
    # Set row names
    tryCatch({
      merged_data$Species <- trimws(as.character(merged_data$Species))
      rownames(merged_data) <- gsub(" ", "_", merged_data$Species)
    }, error = function(e) {
      cat("Error setting row names for merged data in", family_name, ":", e$message, "\n")
      return(NULL)
    })
    
    trait_data <- merged_data[, c("GC_skew"), drop = FALSE]
  } else {
    cat("Error: Required columns ('Accession Number', 'GC_skew') missing in GC_skew data for", family_name, "\n")
    results_df <<- rbind(results_df, data.frame(
      Family = family_name,
      Trait = "GC_skew",
      Valid_Species = 0L,
      Lambda = NA_real_,
      P_lambda = NA_real_,
      K = NA_real_,
      P_K = NA_real_,
      stringsAsFactors = FALSE
    ))
    return(NULL)
  }
  
  ## Step 3: Prepare tree
  tree <- tryCatch({
    multi2di(tree)
  }, error = function(e) {
    cat("Error resolving polytomies for", family_name, ":", e$message, "\n")
    return(NULL)
  })
  if (is.null(tree)) return(NULL)
  
  tree$node.label <- NULL
  min_branch <- 1e-6 * max(node.depth.edgelength(tree))
  tree$edge.length[tree$edge.length < min_branch] <- min_branch
  
  ## Step 4: Match data and tree
  common_species <- intersect(rownames(trait_data), tree$tip.label)
  cat("Number of species in analysis:", length(common_species), "\n")
  
  if (length(common_species) < 10) {
    cat("Error: Fewer than 10 common species between tree and data for", family_name, ". Skipping analysis.\n")
    results_df <<- rbind(results_df, data.frame(
      Family = family_name,
      Trait = "GC_skew",
      Valid_Species = NA_integer_,
      Lambda = NA_real_,
      P_lambda = NA_real_,
      K = NA_real_,
      P_K = NA_real_,
      stringsAsFactors = FALSE
    ))
    return(NULL)
  }
  
  # Subset data and tree to common species
  trait_data <- trait_data[common_species, , drop = FALSE]
  tree <- keep.tip(tree, common_species)
  
  ## Step 5: Phylogenetic signal test for GC_skew
  results <- list()
  trait <- "GC_skew"
  cat("\n=== Testing trait:", trait, "===\n")
  
  # Filter out NA values
  valid_data_trait <- trait_data[!is.na(trait_data[[trait]]) & is.numeric(trait_data[[trait]]), , drop = FALSE]
  valid_species_trait <- rownames(valid_data_trait)
  cat("Number of valid species for", trait, ":", length(valid_species_trait), "\n")
  
  # Add row to results
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
  
  if (length(valid_species_trait) < 10) {
    cat("Skipping", trait, "due to insufficient valid species (< 10)\n")
    return(NULL)
  }
  
  # Subset tree to valid species
  valid_tree <- keep.tip(tree, valid_species_trait)
  current_trait <- setNames(valid_data_trait[[trait]], valid_species_trait)
  
  # Remove Inf/-Inf values
  finite_trait <- current_trait[is.finite(current_trait)]
  if (length(finite_trait) < 10) {
    cat("Skipping", trait, "due to insufficient finite values (< 10)\n")
    results_df$Valid_Species[row_idx] <<- length(finite_trait)
    return(NULL)
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
  
  # Save family-specific results
  saveRDS(results, file.path(tree_base_path, family_name, paste0(family_name, "_phylogenetic_signal_gc_skew.rds")))
  
  # Print summary
  cat("\n===== ", family_name, " Analysis Summary =====\n")
  cat("Analyzed GC_skew across", length(common_species), "species\n")
  cat("Output files created:\n")
  cat("- ", family_name, "_phylogenetic_signal_gc_skew.rds\n")
  cat("- ", family_name, "_GC_skew_trait_evolution.png\n")
}

### Step 4: Loop through families with available trees
for (family in sheet_names) {
  # Check if tree file exists
  tree_file_path <- file.path(tree_base_path, family, paste0(family, "_species.nwk"))
  if (!file.exists(tree_file_path)) {
    cat("Tree file not found for", family, ":", tree_file_path, ". Skipping.\n")
    results_df <<- rbind(results_df, data.frame(
      Family = family,
      Trait = "GC_skew",
      Valid_Species = NA_integer_,
      Lambda = NA_real_,
      P_lambda = NA_real_,
      K = NA_real_,
      P_K = NA_real_,
      stringsAsFactors = FALSE
    ))
    next
  }
  
  # Read original family data
  original_data <- tryCatch({
    read_excel(original_data_path, sheet = family)
  }, error = function(e) {
    cat("ERROR: Failed to read original data from sheet", family, ":", e$message, "\n")
    next
  })
  colnames(original_data) <- gsub(" ", "_", colnames(original_data))
  
  # Read GC skew data
  gc_skew_data <- tryCatch({
    read_excel(gc_skew_data_path, sheet = family)
  }, error = function(e) {
    cat("ERROR: Failed to read GC skew data from sheet", family, ":", e$message, "\n")
    next
  })
  
  # Read family tree
  tree <- tryCatch({
    read.tree(tree_file_path)
  }, error = function(e) {
    cat("ERROR: Failed to read tree file for", family, ":", e$message, "\n")
    return(NULL)
  })
  if (is.null(tree)) next
  
  # Process the family
  process_family(family, original_data, gc_skew_data, tree)
}

### Step 5: Write results to Excel
write.xlsx(list(Results = results_df, Unmatched_Accessions = unmatched_df), 
           output_file_path, 
           sheetName = c("Results", "Unmatched_Accessions"), 
           rowNames = FALSE)

### Step 6: Print final summary
cat("\n===== Phylogenetic Signal Analysis Summary for GC Skew =====\n")
cat("Analyzed", length(selected_traits), "trait across", length(sheet_names), "families\n")
cat("Results summary saved to:", output_file_path, "\n")
cat("Unmatched Accession Numbers saved in the 'Unmatched_Accessions' sheet of", output_file_path, "\n")
print(results_df)

### Step 7: Create non-NA results data frame
results1_df <- results_df[complete.cases(results_df[, c("Valid_Species", "Lambda", "P_lambda", "K", "P_K")]), ]

### Step 8: Write results to Excel
write.xlsx(list(Results = results_df, 
                Unmatched_Accessions = unmatched_df, 
                Non_NA_Results = results1_df), 
           output_file_path, 
           sheetName = c("Results", "Unmatched_Accessions", "Non_NA_Results"), 
           rowNames = FALSE)

### Step 9: Print final summary
cat("\n===== Phylogenetic Signal Analysis Summary for GC Skew =====\n")
cat("Analyzed", length(selected_traits), "trait across", length(sheet_names), "families\n")
cat("Total results:", nrow(results_df), "rows\n")
cat("Non-NA results:", nrow(results1_df), "rows\n")
cat("Results summary saved to:", output_file_path, "\n")
cat("Unmatched Accession Numbers saved in the 'Unmatched_Accessions' sheet of", output_file_path, "\n")
cat("Non-NA results saved in the 'Non_NA_Results' sheet of", output_file_path, "\n")
print(results_df)
cat("\nNon-NA Results:\n")
print(results1_df)