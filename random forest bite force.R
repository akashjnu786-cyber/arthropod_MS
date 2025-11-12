##Arthropod morphology code
## random forest method to asses variable importance

library(ape)
library(readxl)
library(dplyr)
library(randomForest)
library(ggplot2)

# Define paths and parameters
file_path <- "D:/Akash Ajay/insect bite force project/Morphology/Insect_Orders_Averaged.xlsx"
orders <- c("Blattodea", "Coleoptera", "Hymenoptera", "Mantodea", "Odonata", "Orthoptera")
tree_paths <- list(
  Blattodea = "D:/Akash Ajay/insect bite force project/Morphology/Blattodea/Blattodea.NWK",
  Coleoptera = "D:/Akash Ajay/insect bite force project/Morphology/Coleoptera/Coleoptera.NWK",
  Hymenoptera = "D:/Akash Ajay/insect bite force project/Morphology/Hymenoptera/Hymenoptera.NWK",
  Mantodea = "D:/Akash Ajay/insect bite force project/Morphology/Mantodea/Mantodea.NWK",
  Odonata = "D:/Akash Ajay/insect bite force project/Morphology/Odonata/Odonata.NWK",
  Orthoptera = "D:/Akash Ajay/insect bite force project/Morphology/Orthoptera/Orthoptera.NWK"
)
selected_traits <- c("Bite_F", "Head_W", "Head_H", "Head_L", "Thorax_W", "Wing_L", "Body_L")
trait_mapping <- c(
  "iBite" = "Bite_F",
  "head.w" = "Head_W",
  "head.h" = "Head_H",
  "head.l" = "Head_L",
  "th.w" = "Thorax_W",
  "wing.l" = "Wing_L",
  "body.l" = "Body_L"
)
data_traits <- names(trait_mapping)  # Data column names (e.g., "iBite")
predictors <- selected_traits[-1]  # Predictors: all except Bite_F
data_predictors <- names(trait_mapping)[names(trait_mapping) != "iBite"]  # Data predictor columns

# Initialize results data frame for RF importance
rf_importance_df <- data.frame(
  Order = character(),
  Trait = character(),  # Predictor trait
  Importance = numeric(),  # Mean decrease in accuracy
  stringsAsFactors = FALSE
)

# Function to clean tree edge lengths
clean_tree <- function(tree) {
  invalid_edges <- is.na(tree$edge.length) | is.nan(tree$edge.length) | !is.finite(tree$edge.length) | tree$edge.length <= 0
  if (any(invalid_edges)) {
    cat("Warning: Cleaning", sum(invalid_edges), "invalid edges in", order, "\n")
    tree$edge.length[invalid_edges] <- 1e-6
  }
  return(tree)
}

# Function to process one order
process_order_rf <- function(order) {
  cat("\n=== Processing", order, "===\n")
  
  # Load data
  tryCatch({
    data <- as.data.frame(read_excel(file_path, sheet = order))
  }, error = function(e) {
    cat("Error loading data for", order, ":", e$message, "\n")
    return(NULL)
  })
  
  # Check for required columns
  if (!all(c("genus", "species") %in% colnames(data))) {
    cat("Error: 'genus' or 'species' column missing in", order, "\n")
    return(NULL)
  }
  
  # Create row names
  data$genus <- trimws(as.character(data$genus))
  data$species <- trimws(as.character(data$species))
  data <- data[!is.na(data$genus) & !is.na(data$species), , drop = FALSE]
  if (nrow(data) == 0) {
    cat("Error: No valid genus/species data in", order, "\n")
    return(NULL)
  }
  rownames(data) <- paste(data$genus, data$species, sep = "_")
  
  # Select traits and convert to numeric
  available_data_traits <- intersect(data_traits, colnames(data))
  if (length(available_data_traits) == 0) {
    cat("Error: No valid traits available for", order, "\n")
    return(NULL)
  }
  data <- data[, available_data_traits, drop = FALSE]
  for (trait in available_data_traits) {
    original <- data[[trait]]
    data[[trait]] <- as.numeric(as.character(data[[trait]]))
    if (any(is.na(data[[trait]]) & !is.na(original))) {
      cat("Warning: NAs introduced for", trait_mapping[trait], "in", order, ":", 
          paste(original[is.na(data[[trait]]) & !is.na(original)], collapse = ", "), "\n")
    }
  }
  data <- na.omit(data)  # Remove rows with any NAs
  if (nrow(data) < 5) {
    cat("Error: Insufficient data after NA removal for", order, "(", nrow(data), "rows)\n")
    return(NULL)
  }
  
  # Load tree
  tryCatch({
    tree <- read.tree(tree_paths[[order]])
    tree <- multi2di(tree, random = TRUE)
    tree$node.label <- NULL
    tree <- clean_tree(tree)
  }, error = function(e) {
    cat("Error loading tree for", order, ":", e$message, "\n")
    return(NULL)
  })
  
  # Match species
  common_species <- intersect(rownames(data), tree$tip.label)
  cat("Number of common species:", length(common_species), "\n")
  if (length(common_species) < 5) {
    cat("Error: Insufficient common species for", order, "(", length(common_species), ")\n")
    return(NULL)
  }
  
  data <- data[common_species, , drop = FALSE]
  tree <- keep.tip(tree, common_species)
  
  # Prepare data for RF (Bite_F as response, others as predictors)
  if ("iBite" %in% colnames(data)) {
    response_col <- "iBite"
    predictor_cols <- setdiff(available_data_traits, "iBite")
    if (length(predictor_cols) == 0) {
      cat("Error: No predictors available for", order, "\n")
      return(NULL)
    }
    rf_data <- data[, c(response_col, predictor_cols), drop = FALSE]
    rf_data[, predictor_cols] <- scale(rf_data[, predictor_cols])  # Scale predictors
    rf_data$response <- rf_data[[response_col]]  # Response variable
    rf_data <- rf_data[, c("response", predictor_cols)]
    
    # Fit random forest
    rf_model <- tryCatch({
      randomForest(response ~ ., data = rf_data, ntree = 500, importance = TRUE)
    }, error = function(e) {
      cat("Error fitting RF for", order, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(rf_model)) {
      importance_scores <- importance(rf_model, type = 1)  # Mean decrease in accuracy
      importance_df <- data.frame(
        Order = order,
        Trait = trait_mapping[rownames(importance_scores)],  # Map to display names (e.g., Head_W)
        Importance = importance_scores[, 1],
        stringsAsFactors = FALSE
      )
      rf_importance_df <<- rbind(rf_importance_df, importance_df)
      
      # Plot importance
      importance_plot <- ggplot(importance_df, aes(x = reorder(Trait, Importance), y = Importance)) +
        geom_col(fill = "steelblue") +
        coord_flip() +
        theme_minimal() +
        labs(title = paste("RF Variable Importance for Bite_F in", order), 
             x = "Predictor Trait", y = "Mean Decrease in Accuracy") +
        theme(axis.text.y = element_text(size = 10))
      ggsave(paste0(order, "_rf_importance.png"), importance_plot, width = 8, height = 6, dpi = 300)
      ggsave(paste0(order, "_rf_importance.pdf"), importance_plot, width = 8, height = 6)
      
      cat("Top predictor for", order, ":", importance_df$Trait[which.max(importance_df$Importance)], 
          "(Importance =", round(max(importance_df$Importance), 2), ")\n")
    }
    return(rf_model)
  } else {
    cat("Error: Response variable 'iBite' missing for", order, "\n")
    return(NULL)
  }
}

# Process each order
for (order in orders) {
  result <- process_order_rf(order)
}

# Save results
write.csv(rf_importance_df, "rf_importance_per_order.csv", row.names = FALSE)

# Summary of top predictors
cat("\n=== RF Top Predictors ===\n")
for (order in orders) {
  order_importance <- rf_importance_df %>% filter(Order == order)
  if (nrow(order_importance) > 0) {
    top_predictor <- order_importance$Trait[which.max(order_importance$Importance)]
    cat(order, ":", top_predictor, "(Importance =", 
        round(max(order_importance$Importance), 2), ")\n")
  } else {
    cat(order, ": No results\n")
  }
}

cat("\n=== Outputs ===\n")
cat("Results saved to: rf_importance_per_order.csv\n")
cat("Plots saved to: <order>_rf_importance.png/pdf\n")
