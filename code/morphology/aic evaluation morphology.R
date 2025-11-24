### aic computation morphological values

library(ape)
library(geiger)
library(phytools)
library(readxl)
library(openxlsx)
library(dplyr)

# Create directory for outputs
bayou_dir <- "D:/Akash Ajay/insect bite force project/Morphology/bayou"
dir.create(bayou_dir, showWarnings = FALSE)

# Start capturing console output
sink(file.path(bayou_dir, "all_orders_bayou_output.txt"), type = "output")

# Define orders and their file paths
orders <- c("Odonata", "Orthoptera", "Blattodea", "Coleoptera", "Mantodea", "Hymenoptera")
tree_paths <- c(
  "D:/Akash Ajay/insect bite force project/Morphology/Odonata/Odonata.NWK",
  "D:/Akash Ajay/insect bite force project/Morphology/Orthoptera/Orthoptera.NWK", 
  "D:/Akash Ajay/insect bite force project/Morphology/Blattodea/Blattodea.NWK",
  "D:/Akash Ajay/insect bite force project/Morphology/Coleoptera/Coleoptera.NWK",
  "D:/Akash Ajay/insect bite force project/Morphology/Mantodea/Mantodea.NWK",
  "D:/Akash Ajay/insect bite force project/Morphology/Hymenoptera/Hymenoptera.NWK"
)

# Data file path
file_path <- "D:/Akash Ajay/insect bite force project/Morphology/Insect_Orders_Averaged.xlsx"

# Traits to analyze
traits <- c("iBite", "head.w", "head.h", "head.l", "th.w", "wing.l", "body.l")

# Initialize results data frame
bayou_results <- data.frame(
  Order = character(),
  Trait = character(),
  Sample_Size = numeric(),
  Dropped_Species = character(),
  AIC_PE1 = numeric(),
  AIC_PE1_BM = numeric(),
  AIC_PE2 = numeric(),
  AIC_PE2_BM = numeric(),
  AIC_PE3 = numeric(),
  AIC_PE3_BM = numeric(),
  AIC_BM = numeric(),
  Best_Model = character(),
  stringsAsFactors = FALSE
)

# Function to merge data and coerce traits to numeric
merge_data <- function(df) {
  # Coerce traits to numeric
  for (tr in traits) {
    if (tr %in% colnames(df)) {
      df[[tr]] <- suppressWarnings(as.numeric(df[[tr]]))
      if (sum(is.na(df[[tr]])) > 0) {
        cat("Warning: NAs introduced by coercion for trait", tr, "(count:", sum(is.na(df[[tr]])), ")\n")
      }
    }
  }
  df$species <- gsub(" ", "_", df$ID)
  df
}

# Function to fit BM using geiger with error handling
fit_standard <- function(tree, data, trait) {
  trait_data <- data[[trait]]
  names(trait_data) <- data$species
  valid <- !is.na(trait_data) & names(trait_data) %in% tree$tip.label
  dropped_na <- names(trait_data)[is.na(trait_data)]
  
  if (length(dropped_na) > 0) {
    cat("Dropped species with NA for", trait, ":", paste(dropped_na, collapse = ", "), "\n")
  }
  
  if (sum(valid) < 5) {
    cat("Skipping", trait, ": insufficient non-NA data (", sum(valid), "species)\n")
    return(NULL)
  }
  
  trait_data <- trait_data[valid]
  tree_sub <- keep.tip(tree, names(trait_data))
  
  # Additional tree stability check
  if (any(tree_sub$edge.length <= 0)) {
    tree_sub$edge.length[tree_sub$edge.length <= 0] <- replacement_value
  }
  
  fit_BM <- tryCatch({
    fitContinuous(tree_sub, trait_data, model = "BM")
  }, error = function(e) {
    cat("Warning: BM model failed for", trait, ":", e$message, "\n")
    return(NULL)
  })
  
  list(BM = fit_BM, Dropped_NA = dropped_na)
}

# Function for PE density
dpe <- function(x, t, lambda, sigma_jump, weights = NULL, sigma_bm = 0, max_k = 50, tol = 1e-10) {
  k <- 0:max_k
  probs <- dpois(k, lambda * t)
  vars <- sigma_bm^2 * t + k * sigma_jump^2
  vars[vars == 0] <- 1e-10
  dens <- probs * dnorm(x, mean = 0, sd = sqrt(vars))
  sum_dens <- sum(dens)
  if (sum(probs) < 1 - tol) warning("Poisson sum not converged; increase max_k")
  sum_dens
}

# Log-likelihood for PE or PE + BM
loglik_pe <- function(par, changes, lengths, n_rates) {
  lambda <- exp(par[1])
  sigma_jump <- exp(par[2:(1+n_rates)])
  weights <- if (n_rates > 1) exp(par[(2+n_rates):(1+2*n_rates-1)]) / sum(exp(par[(2+n_rates):(1+2*n_rates-1)])) else NULL
  sigma_bm <- if (length(par) > 1+2*n_rates-1) exp(par[length(par)]) else 0
  
  # Calculate log-likelihood for each edge
  log_liks <- sapply(1:length(changes), function(i) {
    log(dpe(changes[i], lengths[i], lambda, sigma_jump, weights, sigma_bm) + 1e-10)
  })
  
  sum(log_liks)
}

# Function to fit PE models with enhanced stability
fit_pulsed <- function(tree, data, trait) {
  trait_data <- data[[trait]]
  names(trait_data) <- data$species
  valid <- !is.na(trait_data) & names(trait_data) %in% tree$tip.label
  dropped_na <- names(trait_data)[is.na(trait_data)]
  
  if (length(dropped_na) > 0) {
    cat("Dropped species with NA for", trait, ":", paste(dropped_na, collapse = ", "), "\n")
  }
  
  if (sum(valid) < 5) {
    cat("Skipping", trait, ": insufficient non-NA data (", sum(valid), "species)\n")
    return(NULL)
  }
  
  trait_data <- trait_data[valid]
  tree_sub <- keep.tip(tree, names(trait_data))
  
  # Ensure tree stability
  if (any(tree_sub$edge.length <= 0)) {
    tree_sub$edge.length[tree_sub$edge.length <= 0] <- replacement_value
  }
  
  # Get ancestral states with error handling
  anc <- tryCatch({
    fastAnc(tree_sub, trait_data)
  }, error = function(e) {
    cat("Warning: fastAnc failed for", trait, ":", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(anc)) return(NULL)
  
  # Calculate edge changes
  edge_changes <- rep(NA, nrow(tree_sub$edge))
  for (i in 1:nrow(tree_sub$edge)) {
    parent_node <- as.character(tree_sub$edge[i, 1])
    child_node <- tree_sub$edge[i, 2]
    parent_anc <- anc[parent_node]
    
    if (child_node <= Ntip(tree_sub)) {
      child_val <- trait_data[tree_sub$tip.label[child_node]]
    } else {
      child_val <- anc[as.character(child_node)]
    }
    edge_changes[i] <- child_val - parent_anc
  }
  
  lengths <- tree_sub$edge.length
  valid_edges <- !is.na(edge_changes) & !is.na(lengths) & is.finite(edge_changes) & is.finite(lengths)
  
  if (sum(valid_edges) < 5) {
    cat("Skipping", trait, ": insufficient valid edges (", sum(valid_edges), ")\n")
    return(NULL)
  }
  
  edge_changes <- edge_changes[valid_edges]
  lengths <- lengths[valid_edges]
  
  # Check for variation
  if (sd(edge_changes, na.rm = TRUE) == 0 || sd(lengths, na.rm = TRUE) == 0) {
    cat("Skipping", trait, ": no variation in changes or lengths\n")
    return(NULL)
  }
  
  # Conservative parameter initialization
  init_lambda <- log(max(0.01, 1 / mean(lengths, na.rm = TRUE)))
  init_sigma_jump <- log(max(0.01, sd(edge_changes, na.rm = TRUE)))
  init_sigma_bm <- log(max(0.01, sd(edge_changes, na.rm = TRUE) / sqrt(mean(lengths, na.rm = TRUE))))
  init_weight <- log(1)
  
  # DEBUG: Print initial parameters
  cat("Initial parameters for", trait, ":\n")
  cat("  lambda:", exp(init_lambda), "\n")
  cat("  sigma_jump:", exp(init_sigma_jump), "\n")
  cat("  sigma_bm:", exp(init_sigma_bm), "\n")
  
  # PE1 (1 rate class) - 2 parameters
  opt_pe1 <- tryCatch({
    result <- optim(c(init_lambda, init_sigma_jump), 
                    function(par) loglik_pe(par, changes = edge_changes, lengths = lengths, n_rates = 1),
                    control = list(fnscale = -1, maxit = 1000), method = "BFGS")
    cat("PE1 optimization - convergence:", result$convergence, "value:", result$value, "\n")
    result
  }, error = function(e) {
    cat("PE1 model failed for", trait, ":", e$message, "\n")
    return(NULL)
  })
  
  # BM + PE1 - 3 parameters
  opt_pe1_bm <- tryCatch({
    result <- optim(c(init_lambda, init_sigma_jump, init_sigma_bm), 
                    function(par) loglik_pe(par, changes = edge_changes, lengths = lengths, n_rates = 1),
                    control = list(fnscale = -1, maxit = 1000), method = "BFGS")
    cat("PE1+BM optimization - convergence:", result$convergence, "value:", result$value, "\n")
    result
  }, error = function(e) {
    cat("BM+PE1 model failed for", trait, ":", e$message, "\n")
    return(NULL)
  })
  
  # PE2 (2 rate classes) - 4 parameters
  opt_pe2 <- tryCatch({
    result <- optim(c(init_lambda, rep(init_sigma_jump, 2), init_weight), 
                    function(par) loglik_pe(par, changes = edge_changes, lengths = lengths, n_rates = 2),
                    control = list(fnscale = -1, maxit = 2000), method = "BFGS")
    cat("PE2 optimization - convergence:", result$convergence, "value:", result$value, "\n")
    result
  }, error = function(e) {
    cat("PE2 model failed for", trait, ":", e$message, "\n")
    return(NULL)
  })
  
  # BM + PE2 - 5 parameters
  opt_pe2_bm <- tryCatch({
    result <- optim(c(init_lambda, rep(init_sigma_jump, 2), init_weight, init_sigma_bm), 
                    function(par) loglik_pe(par, changes = edge_changes, lengths = lengths, n_rates = 2),
                    control = list(fnscale = -1, maxit = 2000), method = "BFGS")
    cat("PE2+BM optimization - convergence:", result$convergence, "value:", result$value, "\n")
    result
  }, error = function(e) {
    cat("BM+PE2 model failed for", trait, ":", e$message, "\n")
    return(NULL)
  })
  
  # PE3 (3 rate classes) - 6 parameters
  opt_pe3 <- tryCatch({
    result <- optim(c(init_lambda, rep(init_sigma_jump, 3), rep(init_weight, 2)), 
                    function(par) loglik_pe(par, changes = edge_changes, lengths = lengths, n_rates = 3),
                    control = list(fnscale = -1, maxit = 3000), method = "BFGS")
    cat("PE3 optimization - convergence:", result$convergence, "value:", result$value, "\n")
    result
  }, error = function(e) {
    cat("PE3 model failed for", trait, ":", e$message, "\n")
    return(NULL)
  })
  
  # BM + PE3 - 7 parameters
  opt_pe3_bm <- tryCatch({
    result <- optim(c(init_lambda, rep(init_sigma_jump, 3), rep(init_weight, 2), init_sigma_bm), 
                    function(par) loglik_pe(par, changes = edge_changes, lengths = lengths, n_rates = 3),
                    control = list(fnscale = -1, maxit = 3000), method = "BFGS")
    cat("PE3+BM optimization - convergence:", result$convergence, "value:", result$value, "\n")
    result
  }, error = function(e) {
    cat("BM+PE3 model failed for", trait, ":", e$message, "\n")
    return(NULL)
  })
  
  list(PE1 = opt_pe1, PE1_BM = opt_pe1_bm, PE2 = opt_pe2, PE2_BM = opt_pe2_bm, 
       PE3 = opt_pe3, PE3_BM = opt_pe3_bm, Dropped_NA = dropped_na)
}

# Global zero-length branch replacement value
replacement_value <- 1e-6

# Process all orders
for (i in 1:length(orders)) {
  order_name <- orders[i]
  tree_path <- tree_paths[i]
  
  cat("\n")
  cat(rep("=", 60), "\n")
  cat("Processing order:", order_name, "\n")
  cat(rep("=", 60), "\n")
  
  # Load tree with robust zero-length correction
  tree <- read.tree(tree_path)
  
  # ROBUST ZERO-LENGTH CORRECTION
  cat("Original tree - Tips:", length(tree$tip.label), "Edges:", nrow(tree$edge), "\n")
  cat("Zero-length branches before correction:", sum(tree$edge.length == 0), "\n")
  cat("NA branches before correction:", sum(is.na(tree$edge.length)), "\n")
  
  # Convert to dichotomous tree first
  tree <- multi2di(tree)
  tree$node.label <- NULL
  
  # Multiple strategies for zero-length correction
  non_zero_lengths <- tree$edge.length[tree$edge.length > 0]
  if (length(non_zero_lengths) > 0) {
    min_non_zero <- min(non_zero_lengths)
    replacement_value <- min_non_zero * 1e-6
  } else {
    replacement_value <- 1e-6
  }
  
  # Apply correction
  tree$edge.length[tree$edge.length <= 0 | is.na(tree$edge.length)] <- replacement_value
  
  cat("After correction - Tips:", length(tree$tip.label), "Edges:", nrow(tree$edge), "\n")
  cat("Zero-length branches after correction:", sum(tree$edge.length == 0), "\n")
  
  # Load data for this order
  data <- read_excel(file_path, sheet = order_name)
  data_processed <- merge_data(data)
  
  if (is.null(data_processed)) {
    cat("Skipping", order_name, ": data processing failed\n")
    next
  }
  
  for (trait in traits) {
    cat("\n")
    cat(rep("-", 50), "\n")
    cat("Processing", order_name, "for", trait, "\n")
    cat(rep("-", 50), "\n")
    
    # Match species
    common_species <- intersect(data_processed$species, tree$tip.label)
    dropped_unmatched <- setdiff(tree$tip.label, data_processed$species)
    
    if (length(dropped_unmatched) > 0) {
      cat("Dropped tree tips for", order_name, ":", paste(dropped_unmatched, collapse = ", "), "\n")
    }
    
    sample_size <- length(common_species)
    if (sample_size < 5) {
      cat("Skipping", order_name, ":", trait, ": insufficient common species (", sample_size, ")\n")
      next
    }
    
    data_subset <- data_processed[data_processed$species %in% common_species, ]
    
    # Fit models
    standard <- fit_standard(tree, data_subset, trait)
    pulsed <- fit_pulsed(tree, data_subset, trait)
    
    # Collect results with DEBUG output
    dropped_na_standard <- if (!is.null(standard)) standard$Dropped_NA else character()
    dropped_na_pulsed <- if (!is.null(pulsed)) pulsed$Dropped_NA else character()
    dropped_all <- unique(c(dropped_na_standard, dropped_na_pulsed, dropped_unmatched))
    
    aic_bm <- if (!is.null(standard$BM)) standard$BM$opt$aic else NA
    
    # DEBUG: Print likelihood values before AIC calculation
    cat("Likelihood values for", trait, ":\n")
    cat("  PE1 likelihood:", if (!is.null(pulsed$PE1)) pulsed$PE1$value else "NULL", "\n")
    cat("  PE1_BM likelihood:", if (!is.null(pulsed$PE1_BM)) pulsed$PE1_BM$value else "NULL", "\n")
    cat("  PE2 likelihood:", if (!is.null(pulsed$PE2)) pulsed$PE2$value else "NULL", "\n")
    cat("  PE3 likelihood:", if (!is.null(pulsed$PE3)) pulsed$PE3$value else "NULL", "\n")
    
    # Calculate AICs with proper parameter counts
    aic_pe1 <- if (!is.null(pulsed$PE1) && is.finite(pulsed$PE1$value)) -2 * pulsed$PE1$value + 2 * 2 else NA
    aic_pe1_bm <- if (!is.null(pulsed$PE1_BM) && is.finite(pulsed$PE1_BM$value)) -2 * pulsed$PE1_BM$value + 2 * 3 else NA
    aic_pe2 <- if (!is.null(pulsed$PE2) && is.finite(pulsed$PE2$value)) -2 * pulsed$PE2$value + 2 * 4 else NA
    aic_pe2_bm <- if (!is.null(pulsed$PE2_BM) && is.finite(pulsed$PE2_BM$value)) -2 * pulsed$PE2_BM$value + 2 * 5 else NA
    aic_pe3 <- if (!is.null(pulsed$PE3) && is.finite(pulsed$PE3$value)) -2 * pulsed$PE3$value + 2 * 6 else NA
    aic_pe3_bm <- if (!is.null(pulsed$PE3_BM) && is.finite(pulsed$PE3_BM$value)) -2 * pulsed$PE3_BM$value + 2 * 7 else NA
    
    # DEBUG: Print AIC values
    cat("AIC values for", trait, ":\n")
    cat("  BM:", aic_bm, "\n")
    cat("  PE1:", aic_pe1, "\n")
    cat("  PE1_BM:", aic_pe1_bm, "\n")
    cat("  PE2:", aic_pe2, "\n")
    cat("  PE2_BM:", aic_pe2_bm, "\n")
    cat("  PE3:", aic_pe3, "\n")
    cat("  PE3_BM:", aic_pe3_bm, "\n")
    
    aics <- c(BM = aic_bm, PE1 = aic_pe1, PE1_BM = aic_pe1_bm, PE2 = aic_pe2, 
              PE2_BM = aic_pe2_bm, PE3 = aic_pe3, PE3_BM = aic_pe3_bm)
    best_model <- if (all(is.na(aics))) NA_character_ else names(aics)[which.min(aics)]
    
    cat("Best model for", trait, ":", best_model, "\n")
    
    bayou_results <- rbind(bayou_results, data.frame(
      Order = order_name,
      Trait = trait,
      Sample_Size = sample_size,
      Dropped_Species = paste(dropped_all, collapse = ", "),
      AIC_BM = aic_bm,
      AIC_PE1 = aic_pe1,
      AIC_PE1_BM = aic_pe1_bm,
      AIC_PE2 = aic_pe2,
      AIC_PE2_BM = aic_pe2_bm,
      AIC_PE3 = aic_pe3,
      AIC_PE3_BM = aic_pe3_bm,
      Best_Model = best_model
    ))
  }
}

# Save results
write.xlsx(bayou_results, file.path(bayou_dir, "all_orders_bayou_results.xlsx"))
cat("\nResults saved to:", file.path(bayou_dir, "all_orders_bayou_results.xlsx"), "\n")

# Print summary
cat("\n===== Bayou Analysis Summary =====\n")
if (nrow(bayou_results) > 0) {
  for (order in unique(bayou_results$Order)) {
    order_data <- bayou_results[bayou_results$Order == order, ]
    cat("\n", order, ":\n")
    cat("  Traits analyzed:", nrow(order_data), "\n")
    cat("  Best model is BM:", sum(order_data$Best_Model == "BM", na.rm = TRUE), "\n")
    cat("  Best model is PE1:", sum(order_data$Best_Model == "PE1", na.rm = TRUE), "\n")
    cat("  Best model is PE1_BM:", sum(order_data$Best_Model == "PE1_BM", na.rm = TRUE), "\n")
    cat("  Best model is PE2:", sum(order_data$Best_Model == "PE2", na.rm = TRUE), "\n")
    cat("  Best model is PE2_BM:", sum(order_data$Best_Model == "PE2_BM", na.rm = TRUE), "\n")
    cat("  Best model is PE3:", sum(order_data$Best_Model == "PE3", na.rm = TRUE), "\n")
    cat("  Best model is PE3_BM:", sum(order_data$Best_Model == "PE3_BM", na.rm = TRUE), "\n")
  }
} else {
  cat("No results generated\n")
}

# Stop capturing output
sink()
cat("Console output saved to:", file.path(bayou_dir, "all_orders_bayou_output.txt"), "\n")

# After calculating AICs and best_model, add parameter extraction:

# Extract parameters for the best model
best_params <- list()
if (!is.na(best_model)) {
  if (best_model == "BM" && !is.null(standard$BM)) {
    best_params <- list(
      sigma_bm = standard$BM$opt$sigsq,
      lambda = NA,
      sigma_jump1 = NA,
      sigma_jump2 = NA,
      sigma_jump3 = NA,
      weight1 = NA,
      weight2 = NA
    )
  } else {
    pulsed_model <- pulsed[[best_model]]
    if (!is.null(pulsed_model)) {
      par <- pulsed_model$par
      if (best_model == "PE1") {
        best_params <- list(
          sigma_bm = 0,
          lambda = exp(par[1]),
          sigma_jump1 = exp(par[2]),
          sigma_jump2 = NA,
          sigma_jump3 = NA,
          weight1 = NA,
          weight2 = NA
        )
      } else if (best_model == "PE1_BM") {
        best_params <- list(
          sigma_bm = exp(par[3]),
          lambda = exp(par[1]),
          sigma_jump1 = exp(par[2]),
          sigma_jump2 = NA,
          sigma_jump3 = NA,
          weight1 = NA,
          weight2 = NA
        )
      } else if (best_model == "PE2") {
        weights <- exp(par[4:5])
        weights <- weights / sum(weights)
        best_params <- list(
          sigma_bm = 0,
          lambda = exp(par[1]),
          sigma_jump1 = exp(par[2]),
          sigma_jump2 = exp(par[3]),
          sigma_jump3 = NA,
          weight1 = weights[1],
          weight2 = weights[2]
        )
      } else if (best_model == "PE2_BM") {
        weights <- exp(par[4:5])
        weights <- weights / sum(weights)
        best_params <- list(
          sigma_bm = exp(par[6]),
          lambda = exp(par[1]),
          sigma_jump1 = exp(par[2]),
          sigma_jump2 = exp(par[3]),
          sigma_jump3 = NA,
          weight1 = weights[1],
          weight2 = weights[2]
        )
      } else if (best_model == "PE3") {
        weights <- exp(par[5:7])
        weights <- weights / sum(weights)
        best_params <- list(
          sigma_bm = 0,
          lambda = exp(par[1]),
          sigma_jump1 = exp(par[2]),
          sigma_jump2 = exp(par[3]),
          sigma_jump3 = exp(par[4]),
          weight1 = weights[1],
          weight2 = weights[2]
        )
      } else if (best_model == "PE3_BM") {
        weights <- exp(par[5:7])
        weights <- weights / sum(weights)
        best_params <- list(
          sigma_bm = exp(par[8]),
          lambda = exp(par[1]),
          sigma_jump1 = exp(par[2]),
          sigma_jump2 = exp(par[3]),
          sigma_jump3 = exp(par[4]),
          weight1 = weights[1],
          weight2 = weights[2]
        )
      }
    }
  }
}

# Then add these to your results data frame:
bayou_results <- rbind(bayou_results, data.frame(
  Order = order_name,
  Trait = trait,
  Sample_Size = sample_size,
  Dropped_Species = paste(dropped_all, collapse = ", "),
  AIC_BM = aic_bm,
  AIC_PE1 = aic_pe1,
  AIC_PE1_BM = aic_pe1_bm,
  AIC_PE2 = aic_pe2,
  AIC_PE2_BM = aic_pe2_bm,
  AIC_PE3 = aic_pe3,
  AIC_PE3_BM = aic_pe3_bm,
  Best_Model = best_model,
  # Add parameter columns:
  Lambda = ifelse(!is.null(best_params$lambda), best_params$lambda, NA),
  Sigma_BM = ifelse(!is.null(best_params$sigma_bm), best_params$sigma_bm, NA),
  Sigma_Jump1 = ifelse(!is.null(best_params$sigma_jump1), best_params$sigma_jump1, NA),
  Sigma_Jump2 = ifelse(!is.null(best_params$sigma_jump2), best_params$sigma_jump2, NA),
  Sigma_Jump3 = ifelse(!is.null(best_params$sigma_jump3), best_params$sigma_jump3, NA),
  Weight1 = ifelse(!is.null(best_params$weight1), best_params$weight1, NA),
  Weight2 = ifelse(!is.null(best_params$weight2), best_params$weight2, NA
  ))