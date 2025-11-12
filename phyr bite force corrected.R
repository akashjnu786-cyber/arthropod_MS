library(ape)
library(readxl)
library(dplyr)
library(phyr)  # Install if needed: install.packages("phyr")
library(ggplot2)

# Define paths and parameters
file_path <- "D:/Akash Ajay/insect bite force project/Morphology/Insect_Orders_Averaged_Renamed.xlsx"
orders <- c("Blattodea", "Coleoptera", "Hymenoptera", "Mantodea", "Odonata", "Orthoptera","Phasmatodea")
tree_paths <- list(
  Blattodea = "D:/Akash Ajay/insect bite force project/Morphology/Blattodea/Blattodea.NWK",
  Coleoptera = "D:/Akash Ajay/insect bite force project/Morphology/Coleoptera/Coleoptera.NWK",
  Hymenoptera = "D:/Akash Ajay/insect bite force project/Morphology/Hymenoptera/Hymenoptera.NWK",
  Mantodea = "D:/Akash Ajay/insect bite force project/Morphology/Mantodea/Mantodea.NWK",
  Odonata = "D:/Akash Ajay/insect bite force project/Morphology/Odonata/Odonata.NWK",
  Orthoptera = "D:/Akash Ajay/insect bite force project/Morphology/Orthoptera/Orthoptera.NWK",
  Phasmatodea = "D:/Akash Ajay/insect bite force project/Morphology/Phasmatodea/Phasmatodea.NWK"
)
selected_traits <- c("Bite_F", "Head_W", "Head_H", "Head_L", "Thorax_W", "Body_L")
# Define mapping for renaming columns in data frame
rename_mapping <- c(
  "head.w" = "Head_W",
  "head.h" = "Head_H",
  "head.l" = "Head_L",
  "th.w" = "Thorax_W",
  "wing.l" = "Wing_L",
  "body.l" = "Body_L",
  "iBite" = "Bite_F"
)
data_traits <- selected_traits  # Expected column names after renaming
predictors <- selected_traits[-1]  # Predictors: Head_W, Head_H, Head_L, Thorax_W, Body_L

# Initialize warning log
warning_log <- c()

# Start capturing output to a text file
sink("phyr_log_bite_force.txt", append = TRUE, split = FALSE)

# Function to clean tree edge lengths
clean_tree <- function(tree, order) {
  invalid_edges <- is.na(tree$edge.length) | is.nan(tree$edge.length) | !is.finite(tree$edge.length) | tree$edge.length <= 0
  if (any(invalid_edges)) {
    warning_log <<- c(warning_log, paste("Cleaning", sum(invalid_edges), "invalid edges in", order))
    tree$edge.length[invalid_edges] <- 1e-6
  }
  return(tree)
}

# Function to process one order for a single predictor
process_order_phyr <- function(order, predictor) {
  cat("\n=== Processing", order, "with predictor", predictor, "for Bite_F ===\n")
  
  # Load data
  tryCatch({
    data <- as.data.frame(read_excel(file_path, sheet = order))
  }, error = function(e) {
    warning_log <<- c(warning_log, paste("Error loading data for", order, ":", e$message))
    cat("Error loading data:", e$message, "\n")
    return(NULL)
  })
  
  # Check for required columns
  if (!all(c("genus", "species") %in% colnames(data))) {
    warning_log <<- c(warning_log, paste("Error: 'genus' or 'species' column missing in", order))
    cat("Error: 'genus' or 'species' missing\n")
    return(NULL)
  }
  
  # Rename columns in data frame
  for (old_name in names(rename_mapping)) {
    if (old_name %in% colnames(data)) {
      colnames(data)[colnames(data) == old_name] <- rename_mapping[old_name]
    } else {
      cat("Warning: Column", old_name, "not found in", order, "\n")
    }
  }
  
  # Verify trait columns
  required_traits <- c("Bite_F", predictor)
  available_data_traits <- intersect(required_traits, colnames(data))
  cat("Available traits in", order, ":", paste(available_data_traits, collapse = ", "), "\n")
  if (length(available_data_traits) < 2) {
    warning_log <<- c(warning_log, paste("Error: Missing required traits (Bite_F or", predictor, ") in", order))
    cat("Error: Missing required traits\n")
    return(NULL)
  }
  
  # Create row names
  data$genus <- trimws(as.character(data$genus))
  data$species <- trimws(as.character(data$species))
  data <- data[!is.na(data$genus) & !is.na(data$species), , drop = FALSE]
  if (nrow(data) == 0) {
    warning_log <<- c(warning_log, paste("Error: No valid genus/species data in", order))
    cat("Error: No valid genus/species\n")
    return(NULL)
  }
  rownames(data) <- paste(data$genus, data$species, sep = "_")
  
  # Select traits and convert to numeric
  data <- data[, available_data_traits, drop = FALSE]
  cat("Data summary before conversion for", order, ":\n")
  print(summary(data))
  for (trait in available_data_traits) {
    original <- data[[trait]]
    data[[trait]] <- as.numeric(as.character(data[[trait]]))
    if (any(is.na(data[[trait]]) & !is.na(original))) {
      warning_log <<- c(warning_log, paste("NAs introduced for", trait, "in", order, ":", 
                                           paste(original[is.na(data[[trait]]) & !is.na(original)], collapse = ", ")))
      cat("Warning: NAs introduced for", trait, "\n")
    }
  }
  cat("Data dimensions before NA removal for", order, ":", nrow(data), "rows,", ncol(data), "columns\n")
  data <- na.omit(data)  # Remove rows with any NAs
  cat("Data dimensions after NA removal for", order, ":", nrow(data), "rows,", ncol(data), "columns\n")
  if (nrow(data) < 5) {
    warning_log <<- c(warning_log, paste("Error: Insufficient data after NA removal for", order, "(", nrow(data), "rows)"))
    cat("Error: Insufficient data (", nrow(data), "rows)\n")
    return(NULL)
  }
  
  # Load tree
  tryCatch({
    tree <- read.tree(tree_paths[[order]])
    tree <- multi2di(tree, random = TRUE)
    tree$node.label <- NULL
    tree <- clean_tree(tree, order)
  }, error = function(e) {
    warning_log <<- c(warning_log, paste("Error loading tree for", order, ":", e$message))
    cat("Error loading tree:", e$message, "\n")
    return(NULL)
  })
  
  # Match species
  common_species <- intersect(rownames(data), tree$tip.label)
  cat("Number of common species:", length(common_species), "\n")
  if (length(common_species) < 5) {
    warning_log <<- c(warning_log, paste("Error: Insufficient common species for", order, "(", length(common_species), ")"))
    cat("Error: Insufficient common species (", length(common_species), ")\n")
    return(NULL)
  }
  
  data <- data[common_species, , drop = FALSE]
  tree <- keep.tip(tree, common_species)
  
  # Prepare data for phyr
  response_col <- "Bite_F"
  predictor_cols <- predictor
  phyr_data <- data[, c(response_col, predictor_cols), drop = FALSE]
  phyr_data[, predictor_cols] <- scale(phyr_data[, predictor_cols])  # Scale predictor
  phyr_data$response <- phyr_data[[response_col]]  # Response variable
  phyr_data$sp <- rownames(phyr_data)  # Species ID for random effect
  
  # Fit phyr model
  formula_str <- paste("response ~", predictor_cols, "+ (1 | sp__)")
  cat("Model formula for", order, ":", formula_str, "\n")
  phyr_model <- tryCatch({
    pglmm(as.formula(formula_str), data = phyr_data, family = "gaussian", cov_ranef = list(sp = tree))
  }, error = function(e) {
    warning_log <<- c(warning_log, paste("Error fitting phyr for", order, "with", predictor, ":", e$message))
    cat("Error fitting phyr:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(phyr_model)) {
    # Extract fixed effects
    coef_names <- colnames(phyr_model$X)
    fixed_effects <- data.frame(
      Value = phyr_model$B,
      Std.Error = phyr_model$B.se,
      Zscore = phyr_model$B.zscore,
      Pvalue = phyr_model$B.pvalue
    )
    rownames(fixed_effects) <- coef_names
    
    # Extract random effects
    std_dev <- phyr_model$ss
    variance <- std_dev^2
    random_effects <- data.frame(
      Variance = variance,
      Std.Dev = std_dev
    )
    rownames(random_effects) <- names(phyr_model$ss)
    
    # Print model details
    cat("Dependent variable: Bite_F\n")
    cat("Independent variable:", predictor, "\n")
    cat("Random variables:", paste(rownames(random_effects), collapse = ", "), "\n")
    cat("logLik    AIC    BIC \n")
    cat(phyr_model$logLik, phyr_model$AIC, phyr_model$BIC, "\n")
    cat("Random effects:\n")
    print(random_effects)
    cat("Fixed effects:\n")
    print(fixed_effects)
    
    # Check if predictor is present
    if (!(predictor %in% rownames(fixed_effects))) {
      warning_log <<- c(warning_log, paste("Error: Predictor", predictor, "not found in fixed effects for", order))
      cat("Error: Predictor", predictor, "not found in fixed effects for", order, "\n")
      return(NULL)
    }
    
    phyr_importance <- data.frame(
      Order = order,
      Dependent = "Bite_F",
      Predictor = predictor,
      Coefficient = fixed_effects[predictor, "Value"],
      P_value = fixed_effects[predictor, "Pvalue"],
      LogLik = phyr_model$logLik,
      AIC = phyr_model$AIC,
      BIC = phyr_model$BIC,
      Random_effects = paste(rownames(random_effects), collapse = "; "),
      stringsAsFactors = FALSE
    )
    
    # Generate plot for coefficient importance
    phyr_plot <- ggplot(phyr_importance, aes(x = Predictor, y = abs(Coefficient))) +
      geom_col(fill = ifelse(phyr_importance$P_value < 0.05, "red", "gray")) +  # Red for significant (p < 0.05), gray otherwise
      theme_minimal() +
      labs(title = paste("Phyr Coefficient Importance for Bite Force in", order, "with", predictor), 
           x = "Predictor Trait", 
           y = "Absolute Coefficient (Red = p < 0.05)") +
      theme(
        axis.text.y = element_text(size = 12),  # Larger text for predictor names
        axis.text.x = element_text(size = 10),  # Text for coefficient values
        axis.title = element_text(size = 12, face = "bold"),  # Bold axis titles
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)  # Centered, bold title
      )
    # Save plot as PNG and PDF at 800 DPI
    ggsave(paste0(predictor, "_", order, "_phyr_importance_bite_force.png"), phyr_plot, width = 8, height = 6, dpi = 800)
    ggsave(paste0(predictor, "_", order, "_phyr_importance_bite_force.pdf"), phyr_plot, width = 8, height = 6, dpi = 800)
    
    cat("Result for", order, ":", 
        phyr_importance$Predictor, 
        "(p =", round(phyr_importance$P_value, 4), ", Coefficient =", 
        round(phyr_importance$Coefficient, 3), ")\n")
    
    return(phyr_importance)
  }
  return(NULL)
}

# Process each predictor
for (predictor in predictors) {
  cat("\n=== Processing models for predictor:", predictor, "for Bite_F ===\n")
  # Initialize results data frame for this predictor
  phyr_results_df <- data.frame(
    Order = character(),
    Dependent = character(),
    Predictor = character(),
    Coefficient = numeric(),
    P_value = numeric(),
    LogLik = numeric(),
    AIC = numeric(),
    BIC = numeric(),
    Random_effects = character(),
    stringsAsFactors = FALSE
  )
  
  # Process each order for this predictor
  for (order in orders) {
    result <- process_order_phyr(order, predictor)
    if (!is.null(result)) {
      phyr_results_df <- rbind(phyr_results_df, result)
    }
  }
  
  # Save results for this predictor
  write.csv(phyr_results_df, paste0(predictor, "_phyr_importance_per_order_bite_force.csv"), row.names = FALSE)
  cat("Results saved to:", paste0(predictor, "_phyr_importance_per_order_bite_force.csv"), "\n")
}

# Merge individual predictor tables
merged_results <- data.frame()
for (predictor in predictors) {
  file_name <- paste0(predictor, "_phyr_importance_per_order_bite_force.csv")
  if (file.exists(file_name)) {
    # Read the CSV file
    df <- read.csv(file_name, stringsAsFactors = FALSE)
    # Append to merged_results
    merged_results <- rbind(merged_results, df)
  } else {
    cat("Warning: File", file_name, "not found in working directory\n")
  }
}

# Save the merged results to a new CSV file
write.csv(merged_results, "merged_phyr_importance_bite_force.csv", row.names = FALSE)
cat("Merged results saved to: merged_phyr_importance_bite_force.csv\n")

# Print warnings
if (length(warning_log) > 0) {
  cat("\n=== Warnings ===\n")
  for (w in warning_log) cat(w, "\n")
}

# Stop capturing output
sink()

cat("\n=== Outputs ===\n")
cat("Log saved to: phyr_log_bite_force.txt\n")
cat("Results saved to: <predictor>_phyr_importance_per_order_bite_force.csv\n")
cat("Merged results saved to: merged_phyr_importance_bite_force.csv\n")
cat("Plots saved to: <predictor>_<order>_phyr_importance_bite_force.png/pdf (800 DPI)\n")

# Define predictors
predictors <- c("Head_W", "Head_H", "Head_L", "Thorax_W", "Body_L")

# Initialize an empty data frame to store the merged results
merged_results <- data.frame()

# Loop through each predictor to read and merge CSV files
for (predictor in predictors) {
  file_name <- paste0(predictor, "_phyr_importance_per_order_bite_force.csv")
  if (file.exists(file_name)) {
    # Read the CSV file
    df <- read.csv(file_name, stringsAsFactors = FALSE)
    # Append to merged_results
    merged_results <- rbind(merged_results, df)
  } else {
    cat("Warning: File", file_name, "not found in working directory\n")
  }
}

# Save the merged results to a new CSV file
write.csv(merged_results, "merged_phyr_importance_bite_force.csv", row.names = FALSE)
cat("Merged results saved to: merged_phyr_importance_bite_force.csv\n")