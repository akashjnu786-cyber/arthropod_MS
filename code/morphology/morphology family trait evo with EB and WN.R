## Family-Level Model Fitting Analysis with Excel Output

### Load required packages
library(ape)
library(geiger)
library(readxl)
library(openxlsx)

# Define paths
file_path <- "D:/Akash Ajay/insect bite force project/Morphology Family/family_level_data_morphology.xlsx"
tree_paths <- list(
  Acrididae = "D:/Akash Ajay/insect bite force project/Morphology Family/Acrididae/Acrididae.NWK",
  Formicidae = "D:/Akash Ajay/insect bite force project/Morphology Family/Formicidae/Formicidae.NWK",
  Mantidae = "D:/Akash Ajay/insect bite force project/Morphology Family/Mantidae/Mantidae.NWK",
  Carabidae = "D:/Akash Ajay/insect bite force project/Morphology Family/Carabidae/Carabidae_species.nwk"
)

# Traits to analyze
selected_traits <- c("iBite", "head.w", "head.h", "head.l", "th.w", 
                     "wing.l", "body.l")

# Families to process
families <- c("Acrididae", "Formicidae", "Mantidae", "Carabidae")

# Create output directory
output_dir <- "D:/Akash Ajay/insect bite force project/Morphology Family/model_fits_output_with_EB_WN"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Initialize data frame for all AIC values
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
  stringsAsFactors = FALSE
)

# Process each family
for (family in families) {
  cat(paste0("\n", paste(rep("=", 60), collapse=""), "\n"))
  cat(paste0("PROCESSING FAMILY: ", family, "\n"))
  cat(paste0(paste(rep("=", 60), collapse=""), "\n\n"))
  
  # Open file for this family
  output_file <- file.path(output_dir, paste0(family, "_model_fits.txt"))
  sink(output_file, split = TRUE)  # split=TRUE sends to both file and console
  
  cat(paste0("\n", paste(rep("=", 60), collapse=""), "\n"))
  cat(paste0("PROCESSING FAMILY: ", family, "\n"))
  cat(paste0(paste(rep("=", 60), collapse=""), "\n\n"))
  
  # Load data for this family
  data <- read_excel(file_path, sheet = family)
  rownames(data) <- paste(data$genus, data$species, sep = "_")
  
  # Load tree for this family
  tree_path <- tree_paths[[family]]
  tree <- read.tree(tree_path)
  
  cat("Data loaded:", nrow(data), "rows\n")
  cat("Tree loaded:", length(tree$tip.label), "tips\n")
  
  # Process each trait one by one
  for (trait in selected_traits) {
    cat(paste0("\n--- TRAIT: ", trait, " ---\n"))
    
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
      
      cat("\n=== FITTING EB MODEL ===\n")
      eb_fit <- fitContinuous(trait_tree, trait_data, model = "EB")
      cat("EB Summary:\n")
      print(eb_fit)
      cat("EB Parameters:", eb_fit$opt$par, "\n")
      cat("EB AIC:", AIC(eb_fit), "\n")
      
      cat("\n=== FITTING WHITE NOISE MODEL ===\n")
      wn_fit <- fitContinuous(trait_tree, trait_data, model = "white")
      cat("White Noise Summary:\n")
      print(wn_fit)
      cat("White Noise Parameters (sigsq):", wn_fit$opt$sigsq, "\n")
      cat("White Noise AIC:", AIC(wn_fit), "\n")
      
      # Show best model by AIC
      aic_values <- data.frame(
        Model = c("BM", "OU", "delta", "kappa", "EB", "WN"),
        AIC = c(AIC(bm_fit), AIC(ou_fit), AIC(delta_fit), AIC(kappa_fit), 
                AIC(eb_fit), AIC(wn_fit))
      )
      best_model <- aic_values[which.min(aic_values$AIC), ]
      cat(paste0("\nBEST MODEL for ", trait, ": ", best_model$Model, " (AIC = ", round(best_model$AIC, 3), ")\n"))
      
      # Add to all_aic data frame
      all_aic <- rbind(all_aic, data.frame(
        Family = family,
        Trait = trait,
        BM_AIC = aic_values$AIC[aic_values$Model == "BM"],
        OU_AIC = aic_values$AIC[aic_values$Model == "OU"],
        Delta_AIC = aic_values$AIC[aic_values$Model == "delta"],
        Kappa_AIC = aic_values$AIC[aic_values$Model == "kappa"],
        EB_AIC = aic_values$AIC[aic_values$Model == "EB"],
        WN_AIC = aic_values$AIC[aic_values$Model == "WN"],
        Best_Model = best_model$Model,
        Best_AIC = best_model$AIC
      ))
      
    } else {
      cat(paste0("SKIPPING - insufficient data (", length(trait_data), " values)\n"))
      
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
        Best_AIC = NA
      ))
    }
    
    cat(paste0("\n", paste(rep("-", 40), collapse=""), "\n"))
  }
  
  cat(paste0("\n", paste(rep("=", 60), collapse=""), "\n"))
  cat(paste0("COMPLETED FAMILY: ", family, "\n"))
  cat(paste0(paste(rep("=", 60), collapse=""), "\n\n"))
  
  # Close file for this family
  sink()
  cat(paste0("Results saved to: ", output_file, "\n"))
}

# Write all AIC values to Excel
excel_path <- file.path(output_dir, "all_aic_values_with_EB_WN.xlsx")
write.xlsx(all_aic, excel_path, rowNames = FALSE)

cat("\n" + paste(rep("=", 60), collapse="") + "\n")
cat("ALL FAMILIES COMPLETED!\n")
cat("Individual output files created in: ", output_dir, "\n")
cat("AIC summary Excel saved to: ", excel_path, "\n")
cat(paste(rep("=", 60), collapse="") + "\n")