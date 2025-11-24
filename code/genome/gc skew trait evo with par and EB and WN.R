# Define paths
gc_skew_path <- "D:/Akash Ajay/insect bite force project/sravya work/GC skew/gc_skew_order_wise.xlsx"
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

# Load required libraries
library(ape)
library(geiger)
library(readxl)
library(openxlsx)

# Trait to analyze
selected_traits <- c("GC_skew")

# Orders to process
orders <- c("Blattodea", "Coleoptera", "Diptera", "Hemiptera", "Hymenoptera", 
            "Lepidoptera", "Odonata", "Orthoptera", "Phasmatodea", "Psocodea", "Trichoptera")

# Create output directory
output_dir <- "D:/Akash Ajay/insect bite force project/genome features/model_fits_output_GC_Skew2"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Verify file existence and list sheets for debugging
cat("Verifying Excel file:", gc_skew_path, "\n")
if (!file.exists(gc_skew_path)) {
  stop("ERROR: Excel file not found at ", gc_skew_path)
}
available_sheets <- excel_sheets(gc_skew_path)
cat("Available sheets in Excel file:", paste(available_sheets, collapse = ", "), "\n")

# Initialize data frames
all_aic <- data.frame(
  Order = character(),
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

# NEW: Data frame for best model parameters
best_model_params <- data.frame(
  Order = character(),
  Trait = character(),
  Best_Model = character(),
  Best_AIC = numeric(),
  Parameter = character(),
  Parameter_Value = numeric(),
  sigsq = numeric(),
  stringsAsFactors = FALSE
)

# Process each order
for (order in orders) {
  cat(paste0("\n", paste(rep("=", 60), collapse=""), "\n"))
  cat(paste0("PROCESSING ORDER: ", order, "\n"))
  cat(paste0(paste(rep("=", 60), collapse=""), "\n\n"))
  
  # Open file for this order
  output_file <- file.path(output_dir, paste0(order, "_GC_Skew_model_fits.txt"))
  sink(output_file, split = TRUE)
  
  cat(paste0("\n", paste(rep("=", 60), collapse=""), "\n"))
  cat(paste0("PROCESSING ORDER: ", order, "\n"))
  cat(paste0(paste(rep("=", 60), collapse=""), "\n\n"))
  
  # Check if sheet exists
  if (!order %in% available_sheets) {
    cat(paste0("ERROR: Sheet '", order, "' not found in ", gc_skew_path, "\n"))
    sink()
    next
  }
  
  # Load GC skew data for this order
  tryCatch({
    data <- read_excel(gc_skew_path, sheet = order)
    data <- as.data.frame(data)
  }, error = function(e) {
    cat(paste0("ERROR: Failed to read sheet '", order, "' from ", gc_skew_path, ": ", e$message, "\n"))
    sink()
    next
  })
  
  # Log column names for debugging
  cat("Column names in sheet:", paste(colnames(data), collapse = ", "), "\n")
  
  # Check if Species and GC_skew columns exist
  if (!"Species" %in% colnames(data)) {
    cat(paste0("ERROR: 'Species' column missing in ", order, " data\n"))
    sink()
    next
  }
  if (!"GC_skew" %in% colnames(data)) {
    cat(paste0("ERROR: 'GC_skew' column missing in ", order, " data\n"))
    sink()
    next
  }
  
  # Clean Species column
  data$Species <- trimws(as.character(data$Species))
  data$Species <- gsub(" ", "_", data$Species)
  
  # Remove rows with NA in Species
  if (any(is.na(data$Species))) {
    cat(paste0("Warning: NA values found in 'Species' for ", order, ". Removing affected rows.\n"))
    data <- data[!is.na(data$Species), , drop = FALSE]
  }
  
  # Use Species as rownames
  rownames(data) <- make.unique(data$Species)
  
  cat("Sample Species column (modified):", paste(head(data$Species, 3), collapse = ", "), "\n")
  cat("Sample rownames created:", paste(head(rownames(data), 3), collapse = ", "), "\n")
  
  # Select GC_skew trait and convert to numeric
  trait <- "GC_skew"
  data_traits <- data[, trait, drop = FALSE]
  original <- data_traits[[trait]]
  data_traits[[trait]] <- suppressWarnings(as.numeric(as.character(original)))
  if (any(is.na(data_traits[[trait]]) & !is.na(original))) {
    cat(paste0("Warning: NAs introduced by coercion for trait ", trait, " in ", order, "\n"))
    cat("Problematic values:", paste(head(original[is.na(data_traits[[trait]]) & !is.na(original)], 3), collapse = ", "), "\n")
  }
  if (all(is.na(data_traits[[trait]]))) {
    cat(paste0("ERROR: Trait ", trait, " in ", order, " contains only NA values after coercion\n"))
    sink()
    next
  }
  
  # Load tree for this order
  tryCatch({
    tree_path <- tree_paths[[order]]
    if (!file.exists(tree_path)) {
      stop(paste0("Tree file not found: ", tree_path))
    }
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
  }, error = function(e) {
    cat(paste0("ERROR: Failed to load tree for ", order, ": ", e$message, "\n"))
    sink()
    next
  })
  
  cat("Data dimensions after cleaning:", nrow(data_traits), "rows x", ncol(data_traits), "traits\n")
  
  # Process GC_skew
  cat(paste0("\n--- TRAIT: ", trait, " ---\n"))
  
  # Select trait data
  trait_data <- data_traits[[trait]]
  names(trait_data) <- data$Species
  
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
  
  # Only fit models if enough data
  if (length(trait_data) >= 4) {
    
    # Store all model fits to extract parameters later
    model_fits <- list()
    
    cat("\n=== FITTING BM MODEL ===\n")
    bm_fit <- tryCatch({
      fit <- fitContinuous(trait_tree, trait_data, model = "BM")
      model_fits[["BM"]] <- fit
      fit
    }, error = function(e) {
      cat(paste0("ERROR: BM model fitting failed for ", trait, " in ", order, ": ", e$message, "\n"))
      NULL
    })
    if (!is.null(bm_fit)) {
      cat("BM Summary:\n")
      print(bm_fit)
      cat("BM AIC:", AIC(bm_fit), "\n")
    }
    
    cat("\n=== FITTING OU MODEL ===\n")
    ou_fit <- tryCatch({
      fit <- fitContinuous(trait_tree, trait_data, model = "OU")
      model_fits[["OU"]] <- fit
      fit
    }, error = function(e) {
      cat(paste0("ERROR: OU model fitting failed for ", trait, " in ", order, ": ", e$message, "\n"))
      NULL
    })
    if (!is.null(ou_fit)) {
      cat("OU Summary:\n")
      print(ou_fit)
      cat("OU AIC:", AIC(ou_fit), "\n")
    }
    
    cat("\n=== FITTING DELTA MODEL ===\n")
    delta_fit <- tryCatch({
      fit <- fitContinuous(trait_tree, trait_data, model = "delta")
      model_fits[["delta"]] <- fit
      fit
    }, error = function(e) {
      cat(paste0("ERROR: Delta model fitting failed for ", trait, " in ", order, ": ", e$message, "\n"))
      NULL
    })
    if (!is.null(delta_fit)) {
      cat("Delta Summary:\n")
      print(delta_fit)
      cat("Delta AIC:", AIC(delta_fit), "\n")
    }
    
    cat("\n=== FITTING KAPPA MODEL ===\n")
    kappa_fit <- tryCatch({
      fit <- fitContinuous(trait_tree, trait_data, model = "kappa")
      model_fits[["kappa"]] <- fit
      fit
    }, error = function(e) {
      cat(paste0("ERROR: Kappa model fitting failed for ", trait, " in ", order, ": ", e$message, "\n"))
      NULL
    })
    if (!is.null(kappa_fit)) {
      cat("Kappa Summary:\n")
      print(kappa_fit)
      cat("Kappa AIC:", AIC(kappa_fit), "\n")
    }
    
    cat("\n=== FITTING EB MODEL ===\n")
    eb_fit <- tryCatch({
      fit <- fitContinuous(trait_tree, trait_data, model = "EB")
      model_fits[["EB"]] <- fit
      fit
    }, error = function(e) {
      cat(paste0("ERROR: EB model fitting failed for ", trait, " in ", order, ": ", e$message, "\n"))
      NULL
    })
    if (!is.null(eb_fit)) {
      cat("EB Summary:\n")
      print(eb_fit)
      cat("EB AIC:", AIC(eb_fit), "\n")
    }
    
    cat("\n=== FITTING WHITE NOISE MODEL ===\n")
    wn_fit <- tryCatch({
      fit <- fitContinuous(trait_tree, trait_data, model = "white")
      model_fits[["WN"]] <- fit
      fit
    }, error = function(e) {
      cat(paste0("ERROR: White Noise model fitting failed for ", trait, " in ", order, ": ", e$message, "\n"))
      NULL
    })
    if (!is.null(wn_fit)) {
      cat("White Noise Summary:\n")
      print(wn_fit)
      cat("White Noise AIC:", AIC(wn_fit), "\n")
    }
    
    # Collect AIC values for models that succeeded
    aic_values <- data.frame(
      Model = c("BM", "OU", "delta", "kappa", "EB", "WN"),
      AIC = c(
        if (!is.null(bm_fit)) AIC(bm_fit) else NA,
        if (!is.null(ou_fit)) AIC(ou_fit) else NA,
        if (!is.null(delta_fit)) AIC(delta_fit) else NA,
        if (!is.null(kappa_fit)) AIC(kappa_fit) else NA,
        if (!is.null(eb_fit)) AIC(eb_fit) else NA,
        if (!is.null(wn_fit)) AIC(wn_fit) else NA
      )
    )
    valid_aic <- aic_values[!is.na(aic_values$AIC), ]
    best_model <- if (nrow(valid_aic) > 0) valid_aic[which.min(valid_aic$AIC), ] else data.frame(Model = "None", AIC = NA)
    
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
      WN_AIC = aic_values$AIC[aic_values$Model == "WN"],
      Best_Model = best_model$Model,
      Best_AIC = best_model$AIC
    ))
    
    # NEW: Extract parameters for the best model
    if (best_model$Model != "None" && best_model$Model %in% names(model_fits)) {
      best_fit <- model_fits[[best_model$Model]]
      
      # Extract parameters based on model type
      if (best_model$Model == "BM") {
        # BM model: only sigsq parameter
        sigsq_value <- best_fit$opt$sigsq
        best_model_params <- rbind(best_model_params, data.frame(
          Order = order,
          Trait = trait,
          Best_Model = best_model$Model,
          Best_AIC = best_model$AIC,
          Parameter = "NA",
          Parameter_Value = NA,
          sigsq = sigsq_value
        ))
      }
      else if (best_model$Model == "OU") {
        # OU model: alpha and sigsq parameters
        alpha_value <- best_fit$opt$alpha
        sigsq_value <- best_fit$opt$sigsq
        best_model_params <- rbind(best_model_params, data.frame(
          Order = order,
          Trait = trait,
          Best_Model = best_model$Model,
          Best_AIC = best_model$AIC,
          Parameter = "alpha",
          Parameter_Value = alpha_value,
          sigsq = sigsq_value
        ))
      }
      else if (best_model$Model == "delta") {
        # Delta model: delta and sigsq parameters
        delta_value <- best_fit$opt$delta
        sigsq_value <- best_fit$opt$sigsq
        best_model_params <- rbind(best_model_params, data.frame(
          Order = order,
          Trait = trait,
          Best_Model = best_model$Model,
          Best_AIC = best_model$AIC,
          Parameter = "delta",
          Parameter_Value = delta_value,
          sigsq = sigsq_value
        ))
      }
      else if (best_model$Model == "kappa") {
        # Kappa model: kappa and sigsq parameters
        kappa_value <- best_fit$opt$kappa
        sigsq_value <- best_fit$opt$sigsq
        best_model_params <- rbind(best_model_params, data.frame(
          Order = order,
          Trait = trait,
          Best_Model = best_model$Model,
          Best_AIC = best_model$AIC,
          Parameter = "kappa",
          Parameter_Value = kappa_value,
          sigsq = sigsq_value
        ))
      }
      else if (best_model$Model == "EB") {
        # EB model: a (decay parameter) and sigsq parameters
        a_value <- best_fit$opt$a
        sigsq_value <- best_fit$opt$sigsq
        best_model_params <- rbind(best_model_params, data.frame(
          Order = order,
          Trait = trait,
          Best_Model = best_model$Model,
          Best_AIC = best_model$AIC,
          Parameter = "a",
          Parameter_Value = a_value,
          sigsq = sigsq_value
        ))
      }
      else if (best_model$Model == "WN") {
        # White Noise model: only sigsq parameter
        sigsq_value <- best_fit$opt$sigsq
        best_model_params <- rbind(best_model_params, data.frame(
          Order = order,
          Trait = trait,
          Best_Model = best_model$Model,
          Best_AIC = best_model$AIC,
          Parameter = "NA",
          Parameter_Value = NA,
          sigsq = sigsq_value
        ))
      }
      
      cat("Best model parameters extracted:\n")
      cat("  Parameter:", best_model_params[nrow(best_model_params), "Parameter"], "\n")
      cat("  Parameter Value:", best_model_params[nrow(best_model_params), "Parameter_Value"], "\n")
      cat("  sigsq:", best_model_params[nrow(best_model_params), "sigsq"], "\n")
    }
    
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
      WN_AIC = NA,
      Best_Model = "Skipped",
      Best_AIC = NA
    ))
    
    # Also add to best_model_params for skipped orders
    best_model_params <- rbind(best_model_params, data.frame(
      Order = order,
      Trait = trait,
      Best_Model = "Skipped",
      Best_AIC = NA,
      Parameter = "NA",
      Parameter_Value = NA,
      sigsq = NA
    ))
  }
  
  cat(paste0("\n", paste(rep("-", 40), collapse=""), "\n"))
  
  cat(paste0("\n", paste(rep("=", 60), collapse=""), "\n"))
  cat(paste0("COMPLETED ORDER: ", order, "\n"))
  
  # Close file for this order
  sink()
  cat(paste0("Results saved to: ", output_file, "\n\n"))
}

# Write all AIC values to Excel with multiple sheets
excel_path <- file.path(output_dir, "all_aic_values_GC_Skew.xlsx")

# Create a workbook
wb <- createWorkbook()

# Add all AIC values sheet
addWorksheet(wb, "All_AIC_Values")
writeData(wb, "All_AIC_Values", all_aic, rowNames = FALSE)

# NEW: Add best model parameters sheet
addWorksheet(wb, "Best_Model_Parameters")
writeData(wb, "Best_Model_Parameters", best_model_params, rowNames = FALSE)

# Save the workbook
saveWorkbook(wb, excel_path, overwrite = TRUE)

cat(paste0("\n", paste(rep("=", 60), collapse=""), "\n"))
cat("ALL ORDERS COMPLETED!\n")
cat(paste0("Individual output files created in: ", output_dir, "\n"))
cat(paste0("AIC summary Excel saved to: ", excel_path, "\n"))
cat("  - Sheet 'All_AIC_Values': All AIC values for all models\n")
cat("  - Sheet 'Best_Model_Parameters': Parameters for best models\n")
cat("Models fitted: BM, OU, Delta, Kappa, EB, and White Noise\n")
cat("Trait analyzed: GC_skew\n")
cat(paste0(paste(rep("=", 60), collapse=""), "\n"))

# Print summary of best models
cat("\nSUMMARY OF BEST MODELS:\n")
cat("=======================\n")
print(best_model_params[, c("Order", "Best_Model", "Best_AIC", "Parameter", "Parameter_Value", "sigsq")])