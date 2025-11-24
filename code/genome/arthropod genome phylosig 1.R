library(ape)
library(caper)
library(phytools)
library(openxlsx)
library(readxl)

### Read Newick tree files
Blattodea_tree <- read.tree("D:/Akash Ajay/insect bite force project/genome features/Blattodea/Blattodea.NWK")
Coleoptera_tree <- read.tree("D:/Akash Ajay/insect bite force project/genome features/Coleoptera/Coleoptera.NWK")
Diptera_tree <- read.tree("D:/Akash Ajay/insect bite force project/genome features/Diptera/Diptera.NWK")
Hemiptera_tree <- read.tree("D:/Akash Ajay/insect bite force project/genome features/Hemiptera/Hemiptera.NWK")
Hymenoptera_tree <- read.tree("D:/Akash Ajay/insect bite force project/genome features/Hymenoptera/Hymenoptera.NWK")
Lepidoptera_tree <- read.tree("D:/Akash Ajay/insect bite force project/genome features/Lepidoptera/Lepidoptera.NWK")
Odonata_tree <- read.tree("D:/Akash Ajay/insect bite force project/genome features/Odonata/Odonata.NWK")
Orthoptera_tree <- read.tree("D:/Akash Ajay/insect bite force project/genome features/Orthoptera/Orthoptera.NWK")
Phasmatodea_tree <- read.tree("D:/Akash Ajay/insect bite force project/genome features/Phasmatodea/Phasmatodea.NWK")
Psocodea_tree <- read.tree("D:/Akash Ajay/insect bite force project/genome features/Psocodea/Psocodea.NWK")
Trichoptera_tree <- read.tree("D:/Akash Ajay/insect bite force project/genome features/Trichoptera/Trichoptera.NWK")

# Define the file path
file_path <- "D:/Akash Ajay/insect bite force project/genome features/output_by_order.xlsx"

# Get the names of all sheets in the Excel file
sheet_names <- excel_sheets(file_path)

# Read all sheets into a list of data frames
list_of_dfs <- lapply(sheet_names, function(sheet) {
  read_excel(file_path, sheet = sheet)
})

# Assign each sheet to a separate data frame in the global environment
names(list_of_dfs) <- sheet_names
list2env(list_of_dfs, envir = .GlobalEnv)

### Load required packages for analysis
library(picante)
library(adephylo)
library(ade4)
library(phylobase)
library(geiger)

# Initialize results dataframe
results_df <- data.frame(
  Order = character(),
  Trait = character(),
  Valid_Species = integer(),
  Lambda = numeric(),
  P_lambda = numeric(),
  K = numeric(),
  P_K = numeric(),
  stringsAsFactors = FALSE
)

# Function to process each order
process_order <- function(order_name, tree, data) {
  cat("\n\n===== Processing", order_name, "=====\n")
  
  ## Step 1: Prepare data
  data <- as.data.frame(data)
  
  # Check if Species column exists
  if (!"Species" %in% colnames(data)) {
    cat("Error: 'Species' column missing in", order_name, "data\n")
    cat("Available columns:", colnames(data), "\n")
    return(NULL)
  }
  
  # Check for NA values in Species
  if (any(is.na(data$Species))) {
    cat("Warning: NA values found in 'Species' for", order_name, "\n")
    # Remove rows with NA in Species
    data <- data[!is.na(data$Species), , drop = FALSE]
    if (nrow(data) == 0) {
      cat("Error: No valid rows remaining after removing NAs in", order_name, "\n")
      return(NULL)
    }
  }
  
  # Create proper row names by replacing spaces with underscores in Species
  tryCatch({
    data$Species <- trimws(as.character(data$Species))
    rownames(data) <- gsub(" ", "_", data$Species)
  }, error = function(e) {
    cat("Error setting row names for", order_name, ":", e$message, "\n")
    return(NULL)
  })
  
  # Select only the specified numeric traits (updated with exact column names including spaces and Size_bp)
  selected_traits <- c("Size_Mb", "Size_bp", "GC_Content", "Chromosome Number", "Coding genes", "Non - coding genes")
  missing_traits <- setdiff(selected_traits, colnames(data))
  if (length(missing_traits) > 0) {
    cat("Warning: Missing traits in", order_name, ":", missing_traits, "\n")
    selected_traits <- intersect(selected_traits, colnames(data))
    if (length(selected_traits) == 0) {
      cat("Error: No valid traits available for", order_name, "\n")
      return(NULL)
    }
  }
  
  # Coerce selected traits to numeric and log problematic values
  # Keep original data structure but create a clean trait subset
  trait_data <- data.frame(row.names = rownames(data))
  
  for (trait in selected_traits) {
    if (trait %in% colnames(data)) {
      original <- data[[trait]]
      numeric_trait <- as.numeric(as.character(data[[trait]]))
      
      if (any(is.na(numeric_trait) & !is.na(original))) {
        cat("Warning: NAs introduced by coercion for trait", trait, "in", order_name, "\n")
        cat("Problematic values:", original[is.na(numeric_trait) & !is.na(original)], "\n")
      }
      if (all(is.na(numeric_trait))) {
        cat("Warning: Trait", trait, "in", order_name, "contains only NA values after coercion\n")
      }
      # Add to trait_data, keeping species with NAs for this trait (will be handled per trait analysis)
      trait_data[[trait]] <- numeric_trait
    }
  }
  
  ## Step 2: Prepare tree
  tree <- multi2di(tree)
  tree$node.label <- NULL
  
  # Set minimum branch length to avoid singularity
  min_branch <- 1e-6 * max(node.depth.edgelength(tree))
  tree$edge.length[tree$edge.length < min_branch] <- min_branch
  
  ## Step 3: Match data and tree
  common_species <- intersect(rownames(trait_data), tree$tip.label)
  cat("Number of species in analysis:", length(common_species), "\n")
  
  if (length(common_species) < 10) {
    cat("Error: Fewer than 10 common species between tree and data for", order_name, ". Skipping analysis.\n")
    # Add placeholder rows for missing traits
    for (trait in selected_traits) {
      results_df <<- rbind(results_df, data.frame(
        Order = order_name,
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
  trait_data <- trait_data[common_species, , drop = FALSE]
  tree <- keep.tip(tree, common_species)
  
  ## Step 4: Filter species with complete data for Moran's I (across all available traits)
  available_traits <- colnames(trait_data)[colnames(trait_data) %in% selected_traits]
  if (length(available_traits) > 0) {
    valid_rows_moran <- complete.cases(trait_data[, available_traits, drop = FALSE])
    if (sum(valid_rows_moran) < 10) {
      cat("Warning: Fewer than 10 species with complete data for Moran's I in", order_name, ". Skipping Moran's I test.\n")
    } else {
      valid_species_moran <- rownames(trait_data)[valid_rows_moran]
      valid_data_moran <- trait_data[valid_species_moran, available_traits, drop = FALSE]
      valid_tree_moran <- keep.tip(tree, valid_species_moran)
      
      ## Step 5: Create phylo4d object and run Moran's I test
      phylotraits <- phylo4d(valid_tree_moran, valid_data_moran)
      moran <- tryCatch({
        abouheif.moran(phylotraits, method = "Abouheif")
      }, error = function(e) {
        cat("Error in Moran's I test for", order_name, ":", e$message, "\n")
        return(NULL)
      })
      
      if (!is.null(moran)) {
        # Plot function
        plot_moran <- function() {
          plot(moran, main = paste(order_name, "Moran's I Test"))
        }
        
        # Save plots
        png(paste0(order_name, "_moran_test.png"), width = 8, height = 6, units = "in", res = 300)
        plot_moran()
        dev.off()
        
        pdf(paste0(order_name, "_moran_test.pdf"), width = 8, height = 6)
        plot_moran()
        dev.off()
      }
    }
  }
  
  ## Step 6: Phylogenetic signal tests for each trait individually
  results <- list()
  analyzed_traits <- character()
  
  for (trait in selected_traits) {
    cat("\n=== Testing trait:", trait, "===\n")
    
    # Check if trait exists in data
    if (!trait %in% colnames(trait_data)) {
      cat("Trait", trait, "not found in data for", order_name, "\n")
      results_df <<- rbind(results_df, data.frame(
        Order = order_name,
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
    
    # Filter out NA values for THIS SPECIFIC TRAIT ONLY
    valid_data_trait <- trait_data[!is.na(trait_data[[trait]]), , drop = FALSE]
    valid_species_trait <- rownames(valid_data_trait)
    cat("Number of valid species for", trait, ":", length(valid_species_trait), "\n")
    
    # Add row with valid species count
    row_idx <- nrow(results_df) + 1
    results_df <<- rbind(results_df, data.frame(
      Order = order_name,
      Trait = trait,
      Valid_Species = length(valid_species_trait),
      Lambda = NA_real_,
      P_lambda = NA_real_,
      K = NA_real_,
      P_K = NA_real_,
      stringsAsFactors = FALSE
    ))
    
    # Check minimum sample size for this trait
    if (length(valid_species_trait) < 10) {
      cat("Skipping", trait, "due to insufficient valid species (< 10)\n")
      next
    }
    
    # Subset tree to valid species FOR THIS TRAIT
    valid_tree <- keep.tip(tree, valid_species_trait)
    current_trait <- setNames(valid_data_trait[[trait]], valid_species_trait)
    
    # Remove Inf/-Inf values for this trait
    finite_trait <- current_trait[is.finite(current_trait)]
    if (length(finite_trait) < 10) {
      cat("Skipping", trait, "due to insufficient finite values (< 10)\n")
      results_df$Valid_Species[row_idx] <<- length(finite_trait)
      next
    }
    
    # Update the row with actual values
    results_df$Lambda[row_idx] <<- NA_real_  # Will be set below
    results_df$P_lambda[row_idx] <<- NA_real_
    results_df$K[row_idx] <<- NA_real_
    results_df$P_K[row_idx] <<- NA_real_
    
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
    
    # Plot trait evolution with computed ylim
    ylim_range <- range(finite_trait, na.rm = TRUE)
    # Add small buffer to ylim to avoid clipping
    ylim <- ylim_range + c(-0.05, 0.05) * diff(ylim_range)
    
    tryCatch({
      png(paste0(order_name, "_", gsub("[^A-Za-z0-9]", "_", trait), "_trait_evolution.png"), 
          width = 10, height = 6, units = "in", res = 300)
      phenogram(valid_tree, finite_trait, main = paste(trait, "trait evolution in", order_name), ylim = ylim)
      dev.off()
    }, error = function(e) {
      cat("Error plotting phenogram for", trait, "in", order_name, ":", e$message, "\n")
      if (dev.cur() > 1) dev.off()
    })
  }
  
  # Save all results
  saveRDS(results, paste0(order_name, "_phylogenetic_signal_results.rds"))
  
  # Print final summary
  cat("\n\n===== ", order_name, " Analysis Summary =====\n")
  cat("Analyzed", length(analyzed_traits), "traits across", length(common_species), "species\n")
  cat("Output files created:\n")
  cat("- ", order_name, "_moran_test.png/pdf (if Moran's I was computed)\n")
  cat("- ", order_name, "_phylogenetic_signal_results.rds\n")
  cat("- Trait-specific evolution plots (PNG)\n")
}

# List of orders and their corresponding trees and data
orders <- list(
  list(name = "Blattodea", tree = Blattodea_tree, data = Blattodea),
  list(name = "Coleoptera", tree = Coleoptera_tree, data = Coleoptera),
  list(name = "Diptera", tree = Diptera_tree, data = Diptera),
  list(name = "Hemiptera", tree = Hemiptera_tree, data = Hemiptera),
  list(name = "Hymenoptera", tree = Hymenoptera_tree, data = Hymenoptera),
  list(name = "Lepidoptera", tree = Lepidoptera_tree, data = Lepidoptera),
  list(name = "Odonata", tree = Odonata_tree, data = Odonata),
  list(name = "Orthoptera", tree = Orthoptera_tree, data = Orthoptera),
  list(name = "Phasmatodea", tree = Phasmatodea_tree, data = Phasmatodea),
  list(name = "Psocodea", tree = Psocodea_tree, data = Psocodea),
  list(name = "Trichoptera", tree = Trichoptera_tree, data = Trichoptera)
)

# Process each order
for (order in orders) {
  process_order(order$name, order$tree, order$data)
}

# Write the results dataframe to Excel
write.xlsx(results_df, "D:/Akash Ajay/insect bite force project/genome features/phylogenetic_signal_summary.xlsx", rowNames = FALSE)

cat("\nAll analyses completed successfully!\n")
cat("Results summary saved to: phylogenetic_signal_summary.xlsx\n")