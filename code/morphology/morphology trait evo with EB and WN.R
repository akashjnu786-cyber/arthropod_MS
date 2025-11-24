library(ape)
library(geiger)
library(readxl)
library(openxlsx)

# Define paths
file_path <- "D:/Akash Ajay/insect bite force project/Morphology/Insect_Orders_Averaged.xlsx"
tree_paths <- list(
  Blattodea = "D:/Akash Ajay/insect bite force project/Morphology/Blattodea/Blattodea.NWK",
  Coleoptera = "D:/Akash Ajay/insect bite force project/Morphology/Coleoptera/Coleoptera.NWK",
  Hymenoptera = "D:/Akash Ajay/insect bite force project/Morphology/Hymenoptera/Hymenoptera.NWK",
  Mantodea = "D:/Akash Ajay/insect bite force project/Morphology/Mantodea/Mantodea.NWK",
  Odonata = "D:/Akash Ajay/insect bite force project/Morphology/Odonata/Odonata.NWK",
  Orthoptera = "D:/Akash Ajay/insect bite force project/Morphology/Orthoptera/Orthoptera.NWK",
  Phasmatodea = "D:/Akash Ajay/insect bite force project/Morphology/Phasmatodea/Phasmatodea.NWK"
)

# Traits to analyze
selected_traits <- c("iBite", "head.w", "head.h", "head.l", "th.w", 
                     "wing.l", "body.l")

# Orders to process
orders <- c("Blattodea", "Coleoptera", "Hymenoptera", "Mantodea", "Odonata", "Orthoptera","Phasmatodea")

# Create output directory (separate for this run with EB and WN)
output_dir <- "D:/Akash Ajay/insect bite force project/Morphology/model_fits_output_with_EB_WN"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Initialize data frame for all AIC values
all_aic <- data.frame(
  Order = character(),
  Trait = character(),
  BM_AIC = numeric(),
  OU_AIC = numeric(),
  Delta_AIC = numeric(),
  Kappa_AIC = numeric(),
  EB_AIC = numeric(),
  WN_AIC = numeric(),  # Added White Noise model
  Best_Model = character(),
  Best_AIC = numeric(),
  stringsAsFactors = FALSE
)

# Process each order
for (order in orders) {
  cat(paste0("\n", paste(rep("=", 60), collapse=""), "\n"))
  cat(paste0("PROCESSING ORDER: ", order, "\n"))
  cat(paste0(paste(rep("=", 60), collapse=""), "\n\n"))
  
  # Open file for this order
  output_file <- file.path(output_dir, paste0(order, "_model_fits.txt"))
  sink(output_file, split = TRUE)
  
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
    
    # Only fit models if we have enough data
    if (length(trait_data) >= 4) {
      
      cat("\n=== FITTING BM MODEL ===\n")
      bm_fit <- fitContinuous(trait_tree, trait_data, model = "BM")
      cat("BM AIC:", AIC(bm_fit), "\n")
      
      cat("\n=== FITTING OU MODEL ===\n")
      ou_fit <- fitContinuous(trait_tree, trait_data, model = "OU")
      cat("OU AIC:", AIC(ou_fit), "\n")
      
      cat("\n=== FITTING DELTA MODEL ===\n")
      delta_fit <- fitContinuous(trait_tree, trait_data, model = "delta")
      cat("Delta AIC:", AIC(delta_fit), "\n")
      
      cat("\n=== FITTING KAPPA MODEL ===\n")
      kappa_fit <- fitContinuous(trait_tree, trait_data, model = "kappa")
      cat("Kappa AIC:", AIC(kappa_fit), "\n")
      
      cat("\n=== FITTING EB MODEL ===\n")
      eb_fit <- fitContinuous(trait_tree, trait_data, model = "EB")
      cat("EB AIC:", AIC(eb_fit), "\n")
      
      # NEW: Fit White Noise model
      cat("\n=== FITTING WHITE NOISE MODEL ===\n")
      wn_fit <- fitContinuous(trait_tree, trait_data, model = "white")
      cat("White Noise AIC:", AIC(wn_fit), "\n")
      cat("White Noise Parameters (sigma2):", wn_fit$opt$sigsq, "\n")
      
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
        Order = order,
        Trait = trait,
        BM_AIC = aic_values$AIC[aic_values$Model == "BM"],
        OU_AIC = aic_values$AIC[aic_values$Model == "OU"],
        Delta_AIC = aic_values$AIC[aic_values$Model == "delta"],
        Kappa_AIC = aic_values$AIC[aic_values$Model == "kappa"],
        EB_AIC = aic_values$AIC[aic_values$Model == "EB"],
        WN_AIC = aic_values$AIC[aic_values$Model == "WN"],  # White Noise AIC
        Best_Model = best_model$Model,
        Best_AIC = best_model$AIC
      ))
      
    } else {
      cat(paste0("SKIPPING - insufficient data (", length(trait_data), " values)\n"))
      
      # Add row with NAs for skipped traits
      all_aic <- rbind(all_aic, data.frame(
        Order = order,
        Trait = trait,
        BM_AIC = NA,
        OU_AIC = NA,
        Delta_AIC = NA,
        Kappa_AIC = NA,
        EB_AIC = NA,
        WN_AIC = NA,  # White Noise AIC
        Best_Model = "Skipped",
        Best_AIC = NA
      ))
    }
    
    cat(paste0("\n", paste(rep("-", 40), collapse=""), "\n"))
  }
  
  cat(paste0("\n", paste(rep("=", 60), collapse=""), "\n"))
  cat(paste0("COMPLETED ORDER: ", order, "\n"))
  
  # Close file for this order
  sink()
  cat(paste0("Results saved to: ", output_file, "\n"))
}

# Write all AIC values to Excel
excel_path <- file.path(output_dir, "all_aic_values_with_WN.xlsx")
write.xlsx(all_aic, excel_path, rowNames = FALSE)

cat("\n" + paste(rep("=", 60), collapse="") + "\n")
cat("ALL ORDERS COMPLETED!\n")
cat("Individual output files created in: ", output_dir, "\n")
cat("AIC summary Excel saved to: ", excel_path, "\n")
cat("Models fitted: BM, OU, Delta, Kappa, EB, and White Noise\n")
cat(paste(rep("=", 60), collapse="") + "\n")