library(ape)
library(geiger)
library(readxl)

# Define paths
file_path <- "D:/Akash Ajay/insect bite force project/genome features/output_by_order.xlsx"
tree_paths <- list(
  Blattodea = "D:/Akash Ajay/insect bite force project/genome features/Blattodea/Blattodea.NWK",
  Coleoptera = "D:/Akash Ajay/insect bite force project/genome features/Coleoptera/Coleoptera.NWK",
  Diptera = "D:/Akash Ajay/insect bite force project/genome features/Diptera/Diptera.NWK",
  Hemiptera = "D:/Akash Ajay/insect bite force project/genome features/Hemiptera/Hemiptera.NWK",
  Hymenoptera = "D:/Akash Ajay/insect bite force project/genome features/Hymenoptera/Hymenoptera.NWK",
  Lepidoptera = "D:/Akash Ajay/insect bite force project/genome features/Lepidoptera/Lepidoptera.NWK",
  Odonata = "D:/Akash Ajay/insect bite force project/genome features/Odonata/Odonata.NWK",
  Orthoptera = "D:/Akash Ajay/insect bite force project/genome features/Orthoptera/Orthoptera.NWK",
  Phasmatodea = "D:/Akash Ajay/insect bite force project/genome features/Phasmatodea/Phasmatodea.NWK",
  Psocodea = "D:/Akash Ajay/insect bite force project/genome features/Psocodea/Psocodea.NWK",
  Trichoptera = "D:/Akash Ajay/insect bite force project/genome features/Trichoptera/Trichoptera.NWK"
)

# Traits to analyze
selected_traits <- c("Size_Mb", "GC_Content", "Chromosome Number", "Coding genes", "Non - coding genes")

# Models to fit
models_to_fit <- c("BM", "OU", "delta", "kappa")

# Orders to process
orders <- c("Blattodea", "Coleoptera", "Diptera", "Hemiptera", "Hymenoptera", 
            "Lepidoptera", "Odonata", "Orthoptera", "Phasmatodea", "Psocodea", "Trichoptera")

# Create output directory
output_dir <- "D:/Akash Ajay/insect bite force project/genome features/genomic_model_fits_output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Process each order
for (order in orders) {
  cat(paste0("\n", paste(rep("=", 60), collapse=""), "\n"))
  cat(paste0("PROCESSING ORDER: ", order, "\n"))
  cat(paste0(paste(rep("=", 60), collapse=""), "\n\n"))
  
  # Open file for this order
  output_file <- file.path(output_dir, paste0(order, "_genomic_model_fits.txt"))
  sink(output_file, split = TRUE)  # split=TRUE sends to both file and console
  
  cat(paste0("\n", paste(rep("=", 60), collapse=""), "\n"))
  cat(paste0("PROCESSING GENOMIC ORDER: ", order, "\n"))
  cat(paste0("Date: ", Sys.Date(), "\n"))
  cat(paste0(paste(rep("=", 60), collapse=""), "\n\n"))
  
  # Load data for this order
  data <- read_excel(file_path, sheet = order)
  
  # Check if Species column exists
  if (!"Species" %in% colnames(data)) {
    cat(paste0("ERROR: 'Species' column missing in ", order, "\n"))
    sink()
    next
  }
  
  # Remove rows with NA in Species
  if (any(is.na(data$Species))) {
    cat(paste0("Warning: NA values found in 'Species' for ", order, ". Removing affected rows.\n"))
    data <- data[!is.na(data$Species), , drop = FALSE]
  }
  
  # REPLACE WHITESPACE WITH UNDERSCORE IN Species COLUMN AND USE AS ROWNAMES
  data$Species <- trimws(as.character(data$Species))
  data$Species <- gsub(" ", "_", data$Species)  # Replace spaces with underscores in Species column
  rownames(data) <- make.unique(data$Species)   # Handle any duplicates
  
  cat("Sample Species column (modified):", paste(head(data$Species, 3), collapse = ", "), "\n")
  cat("Sample rownames created:", paste(head(rownames(data), 3), collapse = ", "), "\n")
  
  # Select available traits and convert to numeric
  available_traits <- intersect(selected_traits, colnames(data))
  cat("Available traits:", paste(available_traits, collapse = ", "), "\n")
  
  if (length(available_traits) == 0) {
    cat(paste0("ERROR: No valid traits available for ", order, "\n"))
    sink()
    next
  }
  
  # Convert data to regular data.frame (not tibble) to avoid rowname issues
  data <- as.data.frame(data)
  
  # Select only trait columns for analysis
  data_traits <- data[, available_traits, drop = FALSE]
  
  # Convert traits to numeric
  for (trait in available_traits) {
    original <- data_traits[[trait]]
    data_traits[[trait]] <- as.numeric(as.character(original))
    if (any(is.na(data_traits[[trait]]) & !is.na(original))) {
      cat(paste0("Warning: NAs introduced by coercion for trait ", trait, " in ", order, "\n"))
      cat("Problematic values:", paste(head(original[is.na(data_traits[[trait]]) & !is.na(original)], 3), collapse = ", "), "\n")
    }
    if (all(is.na(data_traits[[trait]]))) {
      cat(paste0("Warning: Trait ", trait, " in ", order, " contains only NA values after coercion\n"))
      available_traits <- setdiff(available_traits, trait)
    }
  }
  
  # Remove traits that are all NA
  data_traits <- data_traits[, available_traits, drop = FALSE]
  
  # Load tree for this order
  tree_path <- tree_paths[[order]]
  tree <- read.tree(tree_path)
  cat("Tree loaded:", length(tree$tip.label), "tips\n")
  cat("Sample tree tips:", paste(head(tree$tip.label, 3), collapse = ", "), "\n")
  
  tree <- multi2di(tree)
  tree$node.label <- NULL
  
  # Set minimum branch length to avoid singularity
  if (length(tree$edge.length) > 0) {
    min_branch <- 1e-6 * max(node.depth.edgelength(tree))
    tree$edge.length[tree$edge.length < min_branch] <- min_branch
  }
  
  cat("Data dimensions after cleaning:", nrow(data_traits), "rows x", ncol(data_traits), "traits\n")
  
  # Process each trait one by one
  for (trait in available_traits) {
    cat(paste0("\n--- TRAIT: ", trait, "---\n"))
    
    # Select trait data using modified Species names
    trait_data <- data_traits[[trait]]
    names(trait_data) <- data$Species  # Use the modified Species column with underscores
    
    # Remove NAs and infinite values
    trait_data <- trait_data[!is.na(trait_data) & is.finite(trait_data)]
    
    # Match with tree
    common_species <- intersect(names(trait_data), tree$tip.label)
    trait_data <- trait_data[common_species]
    trait_tree <- keep.tip(tree, common_species)
    
    cat("Data length:", length(trait_data), "\n")
    cat("Tree tips:", length(trait_tree$tip.label), "\n")
    cat("Names match:", identical(names(trait_data), trait_tree$tip.label), "\n")
    if (length(common_species) > 0) {
      cat("Sample matched names:", paste(head(common_species, 3), collapse = ", "), "\n")
    }
    
    # Only fit models if we have enough data
    if (length(trait_data) >= 4) {
      
      cat("\n=== FITTING BM MODEL ===\n")
      bm_fit <- fitContinuous(trait_tree, trait_data, model = "BM")
      cat("BM Summary:\n")
      print(bm_fit)
      cat("BM Parameters:", bm_fit$opt$par, "\n")
      cat("BM AIC:", AIC(bm_fit), "\n")
      
      cat("\n=== FITTING OU MODEL ===\n")
      ou_fit <- fitContinuous(trait_tree, trait_data, model = "OU")
      cat("OU Summary:\n")
      print(ou_fit)
      cat("OU Parameters:", ou_fit$opt$par, "\n")
      cat("OU AIC:", AIC(ou_fit), "\n")
      
      cat("\n=== FITTING DELTA MODEL ===\n")
      delta_fit <- fitContinuous(trait_tree, trait_data, model = "delta")
      cat("Delta Summary:\n")
      print(delta_fit)
      cat("Delta Parameters:", delta_fit$opt$par, "\n")
      cat("Delta AIC:", AIC(delta_fit), "\n")
      
      cat("\n=== FITTING KAPPA MODEL ===\n")
      kappa_fit <- fitContinuous(trait_tree, trait_data, model = "kappa")
      cat("Kappa Summary:\n")
      print(kappa_fit)
      cat("Kappa Parameters:", kappa_fit$opt$par, "\n")
      cat("Kappa AIC:", AIC(kappa_fit), "\n")
      
      # Show best model by AIC
      aic_values <- data.frame(
        Model = c("BM", "OU", "delta", "kappa"),
        AIC = c(AIC(bm_fit), AIC(ou_fit), AIC(delta_fit), AIC(kappa_fit))
      )
      best_model <- aic_values[which.min(aic_values$AIC), ]
      cat(paste0("\nBEST MODEL for ", trait, ": ", best_model$Model, " (AIC = ", round(best_model$AIC, 3), ")\n"))
      
    } else {
      cat(paste0("SKIPPING - insufficient data (", length(trait_data), " values)\n"))
    }
    
    cat(paste0("\n", paste(rep("-", 40), collapse=""), "\n"))
  }
  
  cat(paste0("\n", paste(rep("=", 60), collapse=""), "\n"))
  cat(paste0("COMPLETED GENOMIC ORDER: ", order, "\n"))
  cat(paste0("Analysis date: ", Sys.Date(), "\n"))
  cat(paste0(paste(rep("=", 60), collapse=""), "\n\n"))
  
  # Close file for this order
  sink()
  cat(paste0("Results saved to: ", output_file, "\n"))
}

cat(paste0("\n", paste(rep("=", 60), collapse=""), "\n"))
cat("ALL GENOMIC ORDERS COMPLETED!\n")
cat("Individual output files created in: ", output_dir, "\n")
cat(paste0("Files created for: ", paste(orders, collapse = ", "), "\n"))
cat(paste(rep("=", 60), collapse="") + "\n")