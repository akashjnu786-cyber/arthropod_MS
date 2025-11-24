library(ape)
library(geiger)
library(phytools)
library(caper)
library(readxl)
library(dplyr)
library(tidyr)
library(openxlsx)

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

# Function to perform PGLS for an order
run_pgls <- function(order, data, tree, traits) {
  cat("\n=== Running PGLS for", order, "===\n")
  
  # Create species column and set row names
  data$species <- paste(data$genus, data$species, sep = "_")
  rownames(data) <- data$species
  
  # Match data and tree
  common_species <- intersect(data$species, tree$tip.label)
  cat("Number of species in analysis:", length(common_species), "\n")
  
  if (length(common_species) < 5) {
    warning(paste("Only", length(common_species), "species matched for", order, "! Skipping analysis."))
    return(NULL)
  }
  
  # Subset data and tree to common species
  data <- data[common_species, , drop = FALSE]
  tree <- keep.tip(tree, common_species)
  
  # Scale branch lengths to avoid singular matrix
  min_branch <- 1e-5 * max(node.depth.edgelength(tree))
  tree$edge.length[tree$edge.length < min_branch] <- min_branch
  
  # Ensure tree is ultrametric
  if (!is.ultrametric(tree)) {
    cat("Tree is not ultrametric, forcing ultrametric\n")
    tree <- force.ultrametric(tree, method = "extend")
  }
  
  # Create comparative data object
  comp_data <- comparative.data(phy = tree, data = data, names.col = species, vcv = TRUE)
  
  # Initialize results data frame
  pgls_results <- data.frame(
    Order = character(),
    Response_Trait = character(),
    Predictor_Trait = character(),
    Slope = numeric(),
    P_value = numeric(),
    R_squared = numeric(),
    F_value = numeric(),
    df_num = numeric(),
    df_den = numeric(),
    logLik = numeric(),
    AIC = numeric(),
    lambda = numeric(),
    kappa = numeric(),
    delta = numeric(),
    Best_Model = character(),
    Sample_Size = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Run PGLS for each response-predictor pair
  for (response in traits) {
    cat("  Analyzing response:", response, "\n")
    predictors <- setdiff(traits, response)
    
    for (predictor in predictors) {
      # Subset data for complete cases
      pair_data <- data[, c("species", response, predictor), drop = FALSE]
      pair_data <- pair_data[complete.cases(pair_data), ]
      sample_size <- nrow(pair_data)
      
      cat("    Sample size for", response, "~", predictor, ":", sample_size, "\n")
      
      if (sample_size < 3) {
        cat("    Skipping", response, "~", predictor, ": insufficient non-NA data\n")
        pgls_results <- rbind(pgls_results, data.frame(
          Order = order,
          Response_Trait = response,
          Predictor_Trait = predictor,
          Slope = NA,
          P_value = NA,
          R_squared = NA,
          F_value = NA,
          df_num = NA,
          df_den = NA,
          logLik = NA,
          AIC = NA,
          lambda = NA,
          kappa = NA,
          delta = NA,
          Best_Model = NA,
          Sample_Size = sample_size
        ))
        next
      }
      
      # Create temporary comparative data for the pair
      comp_data_pair <- comparative.data(
        phy = tree,
        data = pair_data,
        names.col = species,
        vcv = TRUE
      )
      
      # Run PGLS with fallback for singular matrix and kappa errors
      tryCatch({
        formula <- as.formula(paste(response, "~", predictor))
        pgls_model <- pgls(formula, data = comp_data_pair, lambda = "ML", kappa = "ML", delta = "ML")
        summary_pgls <- summary(pgls_model)
        
        # Extract results
        slope <- coef(summary_pgls)[2, "Estimate"]
        p_value <- coef(summary_pgls)[2, "Pr(>|t|)"]
        r_squared <- summary_pgls$r.squared
        f_statistic <- summary_pgls$fstatistic
        f_value <- f_statistic["value"]
        df_num <- summary_pgls$df[1]
        df_den <- summary_pgls$df[2]
        log_lik <- logLik(pgls_model)
        aic <- AIC(pgls_model)
        lambda <- pgls_model$param["lambda"]
        kappa <- pgls_model$param["kappa"]
        delta <- pgls_model$param["delta"]
        best_model <- "Full"
        
        # Append to results
        pgls_results <- rbind(pgls_results, data.frame(
          Order = order,
          Response_Trait = response,
          Predictor_Trait = predictor,
          Slope = round(slope, 3),
          P_value = round(p_value, 4),
          R_squared = round(r_squared, 3),
          F_value = round(f_value, 3),
          df_num = df_num,
          df_den = df_den,
          logLik = round(log_lik, 3),
          AIC = round(aic, 3),
          lambda = round(lambda, 3),
          kappa = round(kappa, 3),
          delta = round(delta, 3),
          Best_Model = best_model,
          Sample_Size = sample_size
        ))
      }, error = function(e) {
        if (grepl("system is singular|3D VCV.array needed for kappa transformation", e$message)) {
          cat("    Singular matrix or kappa error for", response, "~", predictor, ", trying alternative models\n")
          
          # Try alternative models: lambda and delta only
          aic_values <- numeric(2)
          models <- list()
          
          # Lambda only
          tryCatch({
            pgls_lambda <- pgls(formula, data = comp_data_pair, lambda = "ML", kappa = 1, delta = 1)
            aic_values[1] <- AIC(pgls_lambda)
            models[[1]] <- pgls_lambda
          }, error = function(e) { aic_values[1] <<- Inf })
          
          # Delta only
          tryCatch({
            pgls_delta <- pgls(formula, data = comp_data_pair, lambda = 1, kappa = 1, delta = "ML")
            aic_values[2] <- AIC(pgls_delta)
            models[[2]] <- pgls_delta
          }, error = function(e) { aic_values[2] <<- Inf })
          
          # Select best model
          if (all(is.infinite(aic_values))) {
            cat("    All alternative models failed for", response, "~", predictor, "\n")
            pgls_results <<- rbind(pgls_results, data.frame(
              Order = order,
              Response_Trait = response,
              Predictor_Trait = predictor,
              Slope = NA,
              P_value = NA,
              R_squared = NA,
              F_value = NA,
              df_num = NA,
              df_den = NA,
              logLik = NA,
              AIC = NA,
              lambda = NA,
              kappa = NA,
              delta = NA,
              Best_Model = NA,
              Sample_Size = sample_size
            ))
            return()
          }
          
          best_idx <- which.min(aic_values)
          best_model_name <- c("Lambda", "Delta")[best_idx]
          pgls_model <- models[[best_idx]]
          summary_pgls <- summary(pgls_model)
          
          # Extract results from best model
          slope <- coef(summary_pgls)[2, "Estimate"]
          p_value <- coef(summary_pgls)[2, "Pr(>|t|)"]
          r_squared <- summary_pgls$r.squared
          f_statistic <- summary_pgls$fstatistic
          f_value <- f_statistic["value"]
          df_num <- summary_pgls$df[1]
          df_den <- summary_pgls$df[2]
          log_lik <- logLik(pgls_model)
          aic <- AIC(pgls_model)
          lambda <- pgls_model$param["lambda"]
          kappa <- pgls_model$param["kappa"]
          delta <- pgls_model$param["delta"]
          
          # Append to results
          pgls_results <<- rbind(pgls_results, data.frame(
            Order = order,
            Response_Trait = response,
            Predictor_Trait = predictor,
            Slope = round(slope, 3),
            P_value = round(p_value, 4),
            R_squared = round(r_squared, 3),
            F_value = round(f_value, 3),
            df_num = df_num,
            df_den = df_den,
            logLik = round(log_lik, 3),
            AIC = round(aic, 3),
            lambda = round(lambda, 3),
            kappa = round(kappa, 3),
            delta = round(delta, 3),
            Best_Model = best_model_name,
            Sample_Size = sample_size
          ))
          cat("    Best model for", response, "~", predictor, ":", best_model_name, "\n")
        } else {
          cat("    Error in PGLS for", response, "~", predictor, ":", e$message, "\n")
          pgls_results <<- rbind(pgls_results, data.frame(
            Order = order,
            Response_Trait = response,
            Predictor_Trait = predictor,
            Slope = NA,
            P_value = NA,
            R_squared = NA,
            F_value = NA,
            df_num = NA,
            df_den = NA,
            logLik = NA,
            AIC = NA,
            lambda = NA,
            kappa = NA,
            delta = NA,
            Best_Model = NA,
            Sample_Size = sample_size
          ))
        }
      })
    }
  }
  
  return(pgls_results)
}

# Initialize list to store all PGLS results
all_pgls_results <- data.frame(
  Order = character(),
  Response_Trait = character(),
  Predictor_Trait = character(),
  Slope = numeric(),
  P_value = numeric(),
  R_squared = numeric(),
  F_value = numeric(),
  df_num = numeric(),
  df_den = numeric(),
  logLik = numeric(),
  AIC = numeric(),
  lambda = numeric(),
  kappa = numeric(),
  delta = numeric(),
  Best_Model = character(),
  Sample_Size = numeric(),
  stringsAsFactors = FALSE
)

# Process each order
for (order in orders) {
  cat("\n=== Processing", order, "===\n")
  
  # Load morphological data
  data <- as.data.frame(read_excel(file_path, sheet = order))
  
  # Rename columns using trait_mapping
  colnames(data) <- ifelse(colnames(data) %in% names(trait_mapping), 
                           trait_mapping[colnames(data)], 
                           colnames(data))
  
  # Select only the specified traits
  data <- data[, c("genus", "species", selected_traits), drop = FALSE]
  
  # Load phylogenetic tree
  tree <- read.tree(tree_paths[[order]])
  tree <- multi2di(tree)  # Convert polytomies to dichotomies
  tree$node.label <- NULL  # Remove node labels
  
  # Run PGLS
  order_results <- run_pgls(order, data, tree, selected_traits)
  
  if (!is.null(order_results)) {
    all_pgls_results <- rbind(all_pgls_results, order_results)
    # Save order-specific results
    write.xlsx(order_results, paste0(order, "_pgls_results.xlsx"))
  }
}

# Filter significant correlations (P < 0.05)
significant_pgls_results <- all_pgls_results %>%
  filter(P_value < 0.05 & !is.na(P_value))

# Save all and significant results
write.xlsx(
  x = list(All_PGLS_Results = all_pgls_results, Significant_PGLS_Results = significant_pgls_results),
  file = "all_pgls_results.xlsx"
)

# Print summary
cat("\n\n===== PGLS Correlation Summary =====\n")
for (order in orders) {
  order_data <- all_pgls_results %>% filter(Order == order)
  if (nrow(order_data) > 0) {
    cat("\n", order, ":\n")
    cat("  Response-Predictor pairs analyzed:", nrow(order_data), "\n")
    cat("  Significant correlations (P < 0.05):", sum(order_data$P_value < 0.05, na.rm = TRUE), "\n")
    cat("  Pairs with sample size < 3:", sum(order_data$Sample_Size < 3), "\n")
    cat("  Minimum sample size:", min(order_data$Sample_Size), "\n")
    cat("  Maximum sample size:", max(order_data$Sample_Size), "\n")
  } else {
    cat("\n", order, ": No results (insufficient data)\n")
  }
}

cat("\nResults saved to: all_pgls_results.xlsx (sheets: All_PGLS_Results, Significant_PGLS_Results)\n")
cat("Order-specific results saved to: <order>_pgls_results.xlsx\n")