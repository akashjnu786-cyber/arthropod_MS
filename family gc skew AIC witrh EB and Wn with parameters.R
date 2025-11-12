### Family Trait Evolution models with AIC comparison for GC Skew
library(ape)
library(geiger)
library(readxl)
library(openxlsx)
library(dplyr)

# Define paths for GC SKEW analysis
file_path <- "D:/Akash Ajay/insect bite force project/sravya work/GC skew/gc_skew_family_wise.xlsx"
tree_base_path <- "D:/Akash Ajay/insect bite force project/genome features/family"

# Get family names from Excel sheets
families <- excel_sheets(file_path)
cat("Families found:", paste(families, collapse = ", "), "\n")

# Build tree paths for families
tree_paths <- list()
for (family in families) {
  # Construct the path: base_path/family/family_species.nwk
  tree_path <- file.path(tree_base_path, family, paste0(family, "_species.NWK"))
  
  # Check if file exists with case variations
  if (file.exists(tree_path)) {
    tree_paths[[family]] <- tree_path
  } else {
    # Try alternative extensions and naming
    possible_paths <- c(
      tree_path,
      file.path(tree_base_path, family, paste0(family, "_species.nwk")),
      file.path(tree_base_path, family, paste0(family, ".NWK")),
      file.path(tree_base_path, family, paste0(family, ".nwk")),
      file.path(tree_base_path, family, paste0(family, "_tree.NWK")),
      file.path(tree_base_path, family, paste0(family, "_tree.nwk"))
    )
    
    found <- FALSE
    for (path in possible_paths) {
      if (file.exists(path)) {
        tree_paths[[family]] <- path
        cat("Found tree for", family, "at:", path, "\n")
        found <- TRUE
        break
      }
    }
    
    if (!found) {
      cat("Warning: Tree file not found for", family, "in folder:", file.path(tree_base_path, family), "\n")
      # List files in the family directory to help debug
      family_dir <- file.path(tree_base_path, family)
      if (dir.exists(family_dir)) {
        cat("Files in", family, "directory:", list.files(family_dir), "\n")
      }
    }
  }
}

# GC Skew traits to analyze
selected_traits <- c("GC_skew")

# Create output directories
output_dir <- "D:/Akash Ajay/insect bite force project/sravya work/GC skew/model_fits_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Initialize data frames for results
all_aic <- data.frame(
  Family = character(),
  Trait = character(),
  BM_AIC = numeric(),
  OU_AIC = numeric(),
  Delta_AIC = numeric(),
  Kappa_AIC = numeric(),
  EB_AIC = numeric(),
  WN_AIC = numeric(),
  Best_Model = character(),
  Best_AIC = numeric(),
  Sample_Size = numeric(),
  stringsAsFactors = FALSE
)

# Initialize data frame for model parameters
all_parameters <- data.frame(
  Family = character(),
  Trait = character(),
  Model = character(),
  sigsq = numeric(),
  z0 = numeric(),
  alpha = numeric(),
  a = numeric(),
  kappa = numeric(),
  delta = numeric(),
  AIC = numeric(),
  Convergence = logical(),
  Sample_Size = numeric(),
  stringsAsFactors = FALSE
)

# Function to extract parameters from model fit
extract_parameters <- function(fit, model_name) {
  if (is.null(fit) || is.null(fit$opt)) {
    return(data.frame(
      sigsq = NA, z0 = NA, alpha = NA, a = NA, kappa = NA, delta = NA,
      AIC = NA, Convergence = FALSE
    ))
  }
  
  params <- data.frame(
    sigsq = ifelse(!is.null(fit$opt$sigsq), fit$opt$sigsq, NA),
    z0 = ifelse(!is.null(fit$opt$z0), fit$opt$z0, NA),
    alpha = ifelse(!is.null(fit$opt$alpha), fit$opt$alpha, NA),
    a = ifelse(!is.null(fit$opt$a), fit$opt$a, NA),
    kappa = ifelse(!is.null(fit$opt$kappa), fit$opt$kappa, NA),
    delta = ifelse(!is.null(fit$opt$delta), fit$opt$delta, NA),
    AIC = ifelse(!is.null(AIC(fit)), AIC(fit), NA),
    Convergence = ifelse(!is.null(fit$opt$convergence), fit$opt$convergence == 0, TRUE)
  )
  return(params)
}

# Process each family
for (family in families) {
  cat("PROCESSING FAMILY:", family, "\n")
  
  # Create family-specific output file
  output_file <- file.path(output_dir, paste0(family, "_gc_skew_model_fits.txt"))
  family_con <- file(output_file, "w")
  
  tryCatch({
    # Load data for this family
    data <- read_excel(file_path, sheet = family)
    data <- as.data.frame(data)
    
    # Check if Species column exists
    if (!"Species" %in% colnames(data)) {
      cat("ERROR: 'Species' column missing in", family, "\n", file = family_con)
      close(family_con)
      next
    }
    
    # Remove rows with NA in Species
    if (any(is.na(data$Species))) {
      cat("Warning: NA values found in 'Species' for", family, ". Removing affected rows.\n", file = family_con)
      data <- data[!is.na(data$Species), , drop = FALSE]
    }
    
    # Process Species names to match tree tip labels
    data$Species <- trimws(as.character(data$Species))
    data$Species <- gsub(" ", "_", data$Species)
    rownames(data) <- make.unique(data$Species)
    
    cat("Sample Species:", paste(head(data$Species, 3), collapse = ", "), "\n", file = family_con)
    cat("Total species in data:", length(data$Species), "\n", file = family_con)
    
    # Select available GC skew traits and convert to numeric
    available_traits <- intersect(selected_traits, colnames(data))
    cat("Available GC skew traits:", paste(available_traits, collapse = ", "), "\n", file = family_con)
    
    if (length(available_traits) == 0) {
      cat("ERROR: No GC skew traits available for", family, "\n", file = family_con)
      close(family_con)
      next
    }
    
    # Convert traits to numeric
    data_traits <- data[, available_traits, drop = FALSE]
    for (trait in available_traits) {
      original <- data_traits[[trait]]
      data_traits[[trait]] <- as.numeric(as.character(original))
      if (any(is.na(data_traits[[trait]]) & !is.na(original))) {
        cat("Warning: NAs introduced by coercion for GC skew trait", trait, "in", family, "\n", file = family_con)
      }
    }
    
    # Load tree for this family
    tree_path <- tree_paths[[family]]
    if (is.null(tree_path)) {
      cat("ERROR: No tree path found for", family, "\n", file = family_con)
      close(family_con)
      next
    }
    
    cat("Loading tree from:", tree_path, "\n", file = family_con)
    tree <- read.tree(tree_path)
    cat("Tree loaded:", length(tree$tip.label), "tips\n", file = family_con)
    cat("Tree tip labels sample:", paste(head(tree$tip.label, 3), collapse = ", "), "\n", file = family_con)
    
    # Basic tree processing
    tree <- multi2di(tree)
    tree$node.label <- NULL
    
    # Set minimum branch length
    if (length(tree$edge.length) > 0) {
      min_branch <- 1e-6 * max(node.depth.edgelength(tree))
      tree$edge.length[tree$edge.length < min_branch] <- min_branch
    }
    
    # Process each GC skew trait
    for (trait in available_traits) {
      cat("\n--- GC SKEW TRAIT:", trait, "---\n", file = family_con)
      
      # Select trait data
      trait_data <- data_traits[[trait]]
      names(trait_data) <- data$Species
      
      # Remove NAs and infinite values
      trait_data <- trait_data[!is.na(trait_data) & is.finite(trait_data)]
      
      # Match with tree
      common_species <- intersect(names(trait_data), tree$tip.label)
      trait_data <- trait_data[common_species]
      
      cat("Species in both data and tree:", length(common_species), "\n", file = family_con)
      cat("Sample common species:", paste(head(common_species, 3), collapse = ", "), "\n", file = family_con)
      
      if (length(trait_data) >= 4) {
        trait_tree <- keep.tip(tree, common_species)
        
        cat("Final data length:", length(trait_data), "\n", file = family_con)
        cat("Final tree tips:", length(trait_tree$tip.label), "\n", file = family_con)
        
        # Initialize AIC values
        aic_vals <- list(BM = NA, OU = NA, delta = NA, kappa = NA, EB = NA, WN = NA)
        
        # Fit models with parameter extraction
        models_to_fit <- c("BM", "OU", "delta", "kappa", "EB", "white")
        
        for (model in models_to_fit) {
          cat("Fitting", model, "model...\n", file = family_con)
          tryCatch({
            fit <- fitContinuous(trait_tree, trait_data, model = model)
            aic_vals[[ifelse(model == "white", "WN", model)]] <- AIC(fit)
            cat(model, "AIC:", round(AIC(fit), 3), "\n", file = family_con)
            
            # Extract parameters
            params <- extract_parameters(fit, model)
            if (!is.na(params$sigsq)) cat("  sigsq:", params$sigsq, "\n", file = family_con)
            if (!is.na(params$alpha)) cat("  alpha:", params$alpha, "\n", file = family_con)
            if (!is.na(params$kappa)) cat("  kappa:", params$kappa, "\n", file = family_con)
            if (!is.na(params$delta)) cat("  delta:", params$delta, "\n", file = family_con)
            if (!is.na(params$a)) cat("  a:", params$a, "\n", file = family_con)
            
            # Store parameters
            all_parameters <- rbind(all_parameters, data.frame(
              Family = family,
              Trait = trait,
              Model = model,
              sigsq = params$sigsq,
              z0 = params$z0,
              alpha = params$alpha,
              a = params$a,
              kappa = params$kappa,
              delta = params$delta,
              AIC = params$AIC,
              Convergence = params$Convergence,
              Sample_Size = length(trait_data),
              stringsAsFactors = FALSE
            ))
            
          }, error = function(e) {
            cat("❌", model, "failed:", e$message, "\n", file = family_con)
          })
        }
        
        # Determine best model
        valid_aics <- unlist(aic_vals)
        valid_aics <- valid_aics[!is.na(valid_aics)]
        
        if (length(valid_aics) > 0) {
          best_model <- names(which.min(valid_aics))
          best_aic <- min(valid_aics)
          cat("✅ BEST MODEL for GC Skew", trait, ":", best_model, "(AIC =", round(best_aic, 3), ")\n", file = family_con)
        } else {
          best_model <- "None"
          best_aic <- NA
          cat("❌ No models converged for GC Skew", trait, "\n", file = family_con)
        }
        
        # Add to AIC data frame
        all_aic <- rbind(all_aic, data.frame(
          Family = family,
          Trait = trait,
          BM_AIC = aic_vals$BM,
          OU_AIC = aic_vals$OU,
          Delta_AIC = aic_vals$delta,
          Kappa_AIC = aic_vals$kappa,
          EB_AIC = aic_vals$EB,
          WN_AIC = aic_vals$WN,
          Best_Model = best_model,
          Best_AIC = best_aic,
          Sample_Size = length(trait_data),
          stringsAsFactors = FALSE
        ))
        
      } else {
        cat("❌ SKIPPING - insufficient GC skew data (", length(trait_data), " values)\n", file = family_con)
        
        all_aic <- rbind(all_aic, data.frame(
          Family = family,
          Trait = trait,
          BM_AIC = NA,
          OU_AIC = NA,
          Delta_AIC = NA,
          Kappa_AIC = NA,
          EB_AIC = NA,
          WN_AIC = NA,
          Best_Model = "Skipped",
          Best_AIC = NA,
          Sample_Size = length(trait_data),
          stringsAsFactors = FALSE
        ))
      }
    }
    
  }, error = function(e) {
    cat("❌ FATAL ERROR processing", family, ":", e$message, "\n", file = family_con)
  })
  
  close(family_con)
  cat("GC Skew results saved to:", output_file, "\n")
}

# Save all GC skew results to Excel files
cat("Saving GC skew results to Excel files...\n")

excel_path_aic <- file.path(output_dir, "gc_skew_aic_results.xlsx")
excel_path_params <- file.path(output_dir, "gc_skew_model_parameters.xlsx")
excel_path_combined <- file.path(output_dir, "gc_skew_comprehensive_results.xlsx")

tryCatch({
  write.xlsx(all_aic, excel_path_aic, rowNames = FALSE)
  cat("GC Skew AIC results saved to:", excel_path_aic, "\n")
}, error = function(e) {
  cat("Error saving GC Skew AIC results:", e$message, "\n")
})

tryCatch({
  write.xlsx(all_parameters, excel_path_params, rowNames = FALSE)
  cat("GC Skew model parameters saved to:", excel_path_params, "\n")
}, error = function(e) {
  cat("Error saving GC Skew model parameters:", e$message, "\n")
})

# Create combined results for GC skew
tryCatch({
  combined_results <- list(
    GC_Skew_AIC_Results = all_aic,
    GC_Skew_Model_Parameters = all_parameters
  )
  write.xlsx(combined_results, excel_path_combined)
  cat("GC Skew combined results saved to:", excel_path_combined, "\n")
}, error = function(e) {
  cat("Error saving GC Skew combined results:", e$message, "\n")
})

# Final summary
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("GC SKEW FAMILY-LEVEL ANALYSIS COMPLETED!\n")
cat("GC Skew AIC results saved to:", excel_path_aic, "\n")
cat("GC Skew model parameters saved to:", excel_path_params, "\n")
cat("GC Skew combined results saved to:", excel_path_combined, "\n")
cat("Models fitted: BM, OU, Delta, Kappa, EB, White Noise\n")
cat(paste(rep("=", 60), collapse=""), "\n")