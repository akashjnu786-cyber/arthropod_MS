library(ape)
library(geiger)
library(readxl)

# Define paths
file_path <- "D:/Akash Ajay/insect bite force project/Morphology/Insect_Orders_Averaged.xlsx"
tree_paths <- list(
  Blattodea = "D:/Akash Ajay/insect bite force project/Morphology/Blattodea/Blattodea.NWK",
  Coleoptera = "D:/Akash Ajay/insect bite force project/Morphology/Coleoptera/Coleoptera.NWK",
  Hymenoptera = "D:/Akash Ajay/insect bite force project/Morphology/Hymenoptera/Hymenoptera.NWK",
  Mantodea = "D:/Akash Ajay/insect bite force project/Morphology/Mantodea/Mantodea.NWK",
  Odonata = "D:/Akash Ajay/insect bite force project/Morphology/Odonata/Odonata.NWK",
  Orthoptera = "D:/Akash Ajay/insect bite force project/Morphology/Orthoptera/Orthoptera.NWK"
)

# Traits to analyze
selected_traits <- c("iBite", "head.w", "head.h", "head.l", "th.w", 
                     "wing.l", "body.l")

# Orders to process
orders <- c("Blattodea", "Coleoptera", "Hymenoptera", "Mantodea", "Odonata", "Orthoptera")

# Create output directory
output_dir <- "D:/Akash Ajay/insect bite force project/Morphology/model_fits_output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Process each order
for (order in orders) {
  cat(paste0("\n", paste(rep("=", 60), collapse=""), "\n"))
  cat(paste0("PROCESSING ORDER: ", order, "\n"))
  cat(paste0(paste(rep("=", 60), collapse=""), "\n\n"))
  
  # Open file for this order
  output_file <- file.path(output_dir, paste0(order, "_model_fits.txt"))
  sink(output_file, split = TRUE)  # split=TRUE sends to both file and console
  
  cat(paste0("\n", paste(rep("=", 60), collapse=""), "\n"))
  cat(paste0("PROCESSING ORDER: ", order, "\n"))
  cat(paste0(paste(rep("=", 60), collapse=""), "\n\n"))
  
  # Load data for this order
  data <- read_excel(file_path, sheet = order)
  rownames(data) <- paste(data$genus, data$species, sep = "_")
  
  # Load tree for this order
  tree_path <- tree_paths[[order]]
  tree <- read.tree(tree_path)
  
  cat("Data loaded:", nrow(data), "rows\n")
  cat("Tree loaded:", length(tree$tip.label), "tips\n")
  
  # Process each trait one by one
  for (trait in selected_traits) {
    cat(paste0("\n--- TRAIT: ", trait, "---\n"))
    
    # Select trait data
    trait_data <- data[[trait]]
    names(trait_data) <- rownames(data)
    
    # Remove NAs
    trait_data <- trait_data[!is.na(trait_data)]
    
    # Match with tree
    common_species <- intersect(names(trait_data), tree$tip.label)
    trait_data <- trait_data[common_species]
    trait_tree <- keep.tip(tree, common_species)
    
    cat("Data length:", length(trait_data), "\n")
    cat("Tree tips:", length(trait_tree$tip.label), "\n")
    cat("Names match:", identical(names(trait_data), trait_tree$tip.label), "\n")
    
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
  cat(paste0("COMPLETED ORDER: ", order, "\n"))
  cat(paste0(paste(rep("=", 60), collapse=""), "\n\n"))
  
  # Close file for this order
  sink()
  cat(paste0("Results saved to: ", output_file, "\n"))
}

cat("\n" + paste(rep("=", 60), collapse="") + "\n")
cat("ALL ORDERS COMPLETED!\n")
cat("Individual output files created in: ", output_dir, "\n")
cat(paste(rep("=", 60), collapse="") + "\n")