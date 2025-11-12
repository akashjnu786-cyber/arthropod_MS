library(ape)
library(geiger)
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(nortest)
library(purrr)

# Create directory for plots
plot_dir <- "D:/Akash Ajay/insect bite force project/sravya work/GC skew/distribution_plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Start capturing console output
sink("D:/Akash Ajay/insect bite force project/sravya work/GC skew/ks_ad_test_output.txt", type = "output")

# Define paths and parameters
file_path <- "D:/Akash Ajay/insect bite force project/sravya work/GC skew/genomic_features_with_gc_skew1.xlsx"
tree_paths <- list(
  Coleoptera = "D:/Akash Ajay/insect bite force project/genome features/Coleoptera/Coleoptera.NWK",
  Orthoptera = "D:/Akash Ajay/insect bite force project/genome features/Orthoptera/Orthoptera.NWK"
)

# PROPER TRAIT NAMES FOR GENOMIC TRAITS
traits <- c("Genome_Size", "Genome_GC", "Chromosome Number", "Coding genes", "Non - coding genes", "GC_skew")
trait_display <- list(
  "Genome_Size" = "Genome Size",
  "Genome_GC" = "Genome GC Content",
  "Chromosome Number" = "Chromosome Number",
  "Coding genes" = "Coding Genes",
  "Non - coding genes" = "Non-coding Genes",
  "GC_skew" = "GC Skew"
)

orders <- names(tree_paths)

# Load trees with zero branch length handling
trees <- list()
for (order in orders) {
  cat("Loading tree for", order, "...\n")
  tree <- tryCatch({
    t <- read.tree(tree_paths[[order]])
    
    # Handle zero branch lengths
    if (any(t$edge.length == 0)) {
      cat("  Found", sum(t$edge.length == 0), "zero branch lengths. Adding epsilon = 1e-8\n")
      t$edge.length[t$edge.length == 0] <- 1e-8
    }
    
    # Ensure tree is bifurcating
    t <- multi2di(t, random = FALSE)
    t$node.label <- NULL
    
    # Clean tip labels
    t$tip.label <- gsub("[[:space:]]+", "_", t$tip.label)
    t$tip.label <- gsub("[^[:alnum:]_]", "", t$tip.label)
    
    cat("  Tree loaded with", length(t$tip.label), "tips\n")
    t
  }, error = function(e) {
    cat("  Error loading tree for", order, ":", e$message, "\n")
    NULL
  })
  trees[[order]] <- tree
}

# Initialize results data frame
test_results <- data.frame(
  Order = character(),
  Trait = character(),
  Sample_Size = numeric(),
  AD_Statistic_BM = numeric(),
  AD_P_Value_BM = numeric(),
  AD_Statistic_OU = numeric(),
  AD_P_Value_OU = numeric(),
  KS_Statistic_BM = numeric(),
  KS_P_Value_BM = numeric(),
  KS_Statistic_OU = numeric(),
  KS_P_Value_OU = numeric(),
  stringsAsFactors = FALSE
)

# IMPROVED simulation functions with zero branch length protection
sim_pics_bm <- function(tree, nsim = 500, sigsq) {
  # Ensure no zero branch lengths
  tree_sim <- tree
  if (any(tree_sim$edge.length <= 0)) {
    tree_sim$edge.length[tree_sim$edge.length <= 0] <- 1e-8
  }
  
  sim_pics <- numeric(0)
  successful_sims <- 0
  max_attempts <- nsim * 2
  
  for(i in 1:max_attempts) {
    if(successful_sims >= nsim) break
    
    sim_data <- tryCatch({
      # Use fastBM for more robust simulation
      data <- fastBM(tree_sim, sig2 = sigsq)
      # Add small noise to avoid perfect correlations
      data + rnorm(length(data), 0, sd(data) * 0.001)
    }, error = function(e) NULL)
    
    if(!is.null(sim_data) && length(sim_data) == length(tree_sim$tip.label)) {
      pics <- tryCatch({
        result <- pic(sim_data, tree_sim)
        # Remove infinite or NA values
        result <- result[is.finite(result) & !is.na(result)]
        if(length(result) > 0) result else NULL
      }, error = function(e) NULL)
      
      if(!is.null(pics) && length(pics) > 5) {
        sim_pics <- c(sim_pics, pics)
        successful_sims <- successful_sims + 1
      }
    }
  }
  
  cat("    BM: Generated", successful_sims, "simulations (", length(sim_pics), "PICs)\n")
  return(sim_pics)
}

sim_pics_ou <- function(tree, nsim = 500, alpha, sigsq, theta) {
  # Ensure no zero branch lengths
  tree_sim <- tree
  if (any(tree_sim$edge.length <= 0)) {
    tree_sim$edge.length[tree_sim$edge.length <= 0] <- 1e-8
  }
  
  # Constrain alpha to reasonable bounds
  alpha_constrained <- max(min(alpha, 5), 1e-8)
  
  sim_pics <- numeric(0)
  successful_sims <- 0
  max_attempts <- nsim * 2
  
  for(i in 1:max_attempts) {
    if(successful_sims >= nsim) break
    
    sim_data <- tryCatch({
      data <- rTraitCont(tree_sim, model = "OU", 
                         alpha = alpha_constrained, 
                         sigma = sqrt(sigsq),
                         theta = theta)
      # Add small noise to avoid perfect correlations
      data + rnorm(length(data), 0, sd(data) * 0.001)
    }, error = function(e) NULL)
    
    if(!is.null(sim_data) && length(sim_data) == length(tree_sim$tip.label)) {
      pics <- tryCatch({
        result <- pic(sim_data, tree_sim)
        # Remove infinite or NA values
        result <- result[is.finite(result) & !is.na(result)]
        if(length(result) > 0) result else NULL
      }, error = function(e) NULL)
      
      if(!is.null(pics) && length(pics) > 5) {
        sim_pics <- c(sim_pics, pics)
        successful_sims <- successful_sims + 1
      }
    }
  }
  
  cat("    OU: Generated", successful_sims, "simulations (", length(sim_pics), "PICs)\n")
  return(sim_pics)
}

# IMPROVED test function with better parameter handling
robust_distribution_test <- function(observed_pics, tree, trait_data, nsim = 500) {
  cat("    Performing distribution tests with", length(observed_pics), "PICs\n")
  
  if(length(observed_pics) < 8) {
    cat("    Warning: Too few PICs for reliable analysis\n")
    return(list(
      bm_dist = numeric(0),
      ou_dist = numeric(0),
      ad_statistic_bm = NA, ad_p.value_bm = NA,
      ad_statistic_ou = NA, ad_p.value_ou = NA,
      ks_statistic_bm = NA, ks_p.value_bm = NA,
      ks_statistic_ou = NA, ks_p.value_ou = NA
    ))
  }
  
  # Initialize results
  results <- list(
    bm_dist = numeric(0),
    ou_dist = numeric(0),
    ad_statistic_bm = NA, ad_p.value_bm = NA,
    ad_statistic_ou = NA, ad_p.value_ou = NA,
    ks_statistic_bm = NA, ks_p.value_bm = NA,
    ks_statistic_ou = NA, ks_p.value_ou = NA
  )
  
  # Ensure tree has no zero branch lengths for model fitting
  tree_fit <- tree
  if (any(tree_fit$edge.length <= 0)) {
    tree_fit$edge.length[tree_fit$edge.length <= 0] <- 1e-8
  }
  
  # Fit BM model with error handling
  bm_fit <- tryCatch({
    cat("    Fitting BM model...\n")
    fit <- fitContinuous(tree_fit, trait_data, model = "BM")
    # Validate BM parameters
    if (fit$opt$sigsq <= 0) {
      cat("    Warning: BM sigsq parameter is non-positive\n")
      NULL
    } else {
      cat("    BM fit successful: sigsq =", fit$opt$sigsq, "\n")
      fit
    }
  }, error = function(e) {
    cat("    BM fit failed:", e$message, "\n")
    NULL
  })
  
  # Fit OU model with constrained parameters
  ou_fit <- tryCatch({
    cat("    Fitting OU model...\n")
    fit <- fitContinuous(tree_fit, trait_data, model = "OU", 
                         bounds = list(alpha = c(1e-8, 5)))  # Constrain alpha
    
    # Validate OU parameters
    if (fit$opt$sigsq <= 0 || fit$opt$alpha <= 0) {
      cat("    Warning: OU parameters are non-positive\n")
      NULL
    } else {
      cat("    OU fit successful: alpha =", fit$opt$alpha, "sigsq =", fit$opt$sigsq, "\n")
      fit
    }
  }, error = function(e) {
    cat("    OU fit failed:", e$message, "\n")
    NULL
  })
  
  # Simulate under BM if fit successful
  if(!is.null(bm_fit)) {
    results$bm_dist <- tryCatch({
      sim_pics_bm(tree, nsim, bm_fit$opt$sigsq)
    }, error = function(e) {
      cat("    BM simulation failed:", e$message, "\n")
      numeric(0)
    })
  }
  
  # Simulate under OU if fit successful
  if(!is.null(ou_fit)) {
    results$ou_dist <- tryCatch({
      sim_pics_ou(tree, nsim, ou_fit$opt$alpha, ou_fit$opt$sigsq, ou_fit$opt$z0)
    }, error = function(e) {
      cat("    OU simulation failed:", e$message, "\n")
      numeric(0)
    })
  }
  
  # Perform Anderson-Darling normality tests
  if(length(observed_pics) >= 8) {
    tryCatch({
      cat("    Performing Anderson-Darling test...\n")
      ad_test <- ad.test(observed_pics)
      results$ad_statistic_bm <- ad_test$statistic
      results$ad_p.value_bm <- ad_test$p.value
      results$ad_statistic_ou <- ad_test$statistic
      results$ad_p.value_ou <- ad_test$p.value
      cat("    AD test p-value:", ad_test$p.value, "\n")
    }, error = function(e) {
      cat("    AD test failed:", e$message, "\n")
    })
  }
  
  # KS tests against simulated distributions with tie handling
  if(length(results$bm_dist) > 30) {
    tryCatch({
      cat("    Performing KS test for BM...\n")
      # Add small jitter to handle ties
      obs_jitter <- observed_pics + rnorm(length(observed_pics), 0, sd(observed_pics) * 0.0001)
      bm_jitter <- results$bm_dist + rnorm(length(results$bm_dist), 0, sd(results$bm_dist) * 0.0001)
      
      ks_test_bm <- ks.test(obs_jitter, bm_jitter)
      results$ks_statistic_bm <- ks_test_bm$statistic
      results$ks_p.value_bm <- ks_test_bm$p.value
      cat("    BM KS test p-value:", ks_test_bm$p.value, "\n")
    }, error = function(e) {
      cat("    KS test for BM failed:", e$message, "\n")
    })
  }
  
  if(length(results$ou_dist) > 30) {
    tryCatch({
      cat("    Performing KS test for OU...\n")
      # Add small jitter to handle ties
      obs_jitter <- observed_pics + rnorm(length(observed_pics), 0, sd(observed_pics) * 0.0001)
      ou_jitter <- results$ou_dist + rnorm(length(results$ou_dist), 0, sd(results$ou_dist) * 0.0001)
      
      ks_test_ou <- ks.test(obs_jitter, ou_jitter)
      results$ks_statistic_ou <- ks_test_ou$statistic
      results$ks_p.value_ou <- ks_test_ou$p.value
      cat("    OU KS test p-value:", ks_test_ou$p.value, "\n")
    }, error = function(e) {
      cat("    KS test for OU failed:", e$message, "\n")
    })
  }
  
  return(results)
}

# Enhanced PIC computation with zero branch length handling
safe_pic_computation <- function(trait_data, tree) {
  # Ensure no zero branch lengths
  tree_safe <- tree
  if (any(tree_safe$edge.length <= 0)) {
    zero_count <- sum(tree_safe$edge.length <= 0)
    cat("    Found", zero_count, "zero/negative branch lengths. Applying correction...\n")
    tree_safe$edge.length[tree_safe$edge.length <= 0] <- 1e-8
  }
  
  # Compute PICs with error handling
  pic_values <- tryCatch({
    pics <- pic(trait_data, tree_safe)
    # Remove any infinite or NA values
    pics <- pics[is.finite(pics) & !is.na(pics)]
    
    if(length(pics) == 0) {
      stop("No valid PICs computed")
    }
    
    cat("    Computed", length(pics), "valid PICs\n")
    pics
  }, error = function(e) {
    cat("    PIC computation error:", e$message, "\n")
    numeric(0)
  })
  
  return(pic_values)
}

# Enhanced plotting function
create_distribution_plot <- function(order, trait, observed_pics, test_results, bm_dist, ou_dist) {
  if(length(observed_pics) < 3) return(NULL)
  
  # Create data frame for plotting
  plot_data <- data.frame(
    Value = c(observed_pics, bm_dist, ou_dist),
    Type = factor(c(
      rep("Observed PICs", length(observed_pics)),
      rep("BM Simulated", length(bm_dist)),
      rep("OU Simulated", length(ou_dist))
    ), levels = c("Observed PICs", "BM Simulated", "OU Simulated"))
  )
  
  # Get proper trait name for display
  trait_name <- trait_display[[trait]]
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = Value, fill = Type, color = Type)) +
    # Density curves for all types
    geom_density(data = subset(plot_data, Type == "Observed PICs"), alpha = 0.7, size = 1.2) +
    geom_density(data = subset(plot_data, Type == "BM Simulated"), alpha = 0.5, size = 1, linetype = "dashed") +
    geom_density(data = subset(plot_data, Type == "OU Simulated"), alpha = 0.5, size = 1, linetype = "dotted") +
    
    # Custom colors
    scale_fill_manual(values = c("Observed PICs" = "#6BAED6", "BM Simulated" = "#E74C3C", "OU Simulated" = "#27AE60")) +
    scale_color_manual(values = c("Observed PICs" = "#3182BD", "BM Simulated" = "#C0392B", "OU Simulated" = "#229954")) +
    
    # Labels and title with test results
    labs(
      title = paste(order, ":", trait_name, "Distribution Comparison"),
      subtitle = paste(
        "AD Test p =", sprintf("%.3f", test_results$ad_p.value_bm),
        "| KS-BM p =", sprintf("%.3f", test_results$ks_p.value_bm),
        "| KS-OU p =", sprintf("%.3f", test_results$ks_p.value_ou)
      ),
      x = "Standardized PIC Values",
      y = "Density",
      fill = "Distribution",
      color = "Distribution"
    ) +
    
    # Theme customization
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "darkgray"),
      legend.position = "top",
      panel.grid.major = element_line(color = "gray90", size = 0.2),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    ) +
    
    # Add sample size annotation
    annotate("text", x = Inf, y = Inf, 
             label = paste("n =", length(observed_pics)), 
             hjust = 1.1, vjust = 1.1, size = 4, color = "black")
  
  return(p)
}

# Main function to compute PICs and create enhanced plots
perform_distribution_analysis <- function(order, tree, data, trait) {
  cat("  Analyzing trait:", trait, "\n")
  
  # Data validation and processing
  data[[trait]] <- suppressWarnings(as.numeric(data[[trait]]))
  if (all(is.na(data[[trait]]))) {
    cat("  Warning: All data for", trait, "is NA or non-numeric\n")
    return(NULL)
  }
  
  trait_data <- data[[trait]]
  names(trait_data) <- data$Species
  trait_data <- trait_data[tree$tip.label]
  
  # Remove NAs
  valid <- !is.na(trait_data)
  valid_count <- sum(valid)
  
  cat("  Valid data points:", valid_count, "/", length(valid), "\n")
  
  if (valid_count < 15) {
    cat("  Skipping: insufficient non-NA data\n")
    return(NULL)
  }
  
  # Subset tree and data
  trait_data <- trait_data[valid]
  tree_sub <- drop.tip(tree, tree$tip.label[!valid])
  
  cat("  Subset tree has", length(tree_sub$tip.label), "tips\n")
  
  # Check tree quality
  if (any(tree_sub$edge.length <= 0)) {
    zero_count <- sum(tree_sub$edge.length <= 0)
    cat("  Found", zero_count, "zero/negative branch lengths. Applying correction...\n")
    tree_sub$edge.length[tree_sub$edge.length <= 0] <- 1e-8
  }
  
  # Compute PICs safely
  pic_values <- safe_pic_computation(trait_data, tree_sub)
  
  if (length(pic_values) < 8) {
    cat("  Warning: Too few valid PICs for analysis (", length(pic_values), ")\n")
    return(NULL)
  }
  
  # Standardize PICs
  pic_values_std <- scale(pic_values)[,1]
  
  # Perform distribution tests
  test_output <- robust_distribution_test(pic_values_std, tree_sub, trait_data, nsim = 300)
  
  # Standardize simulated distributions
  bm_dist_std <- if(length(test_output$bm_dist) > 0) scale(test_output$bm_dist)[,1] else numeric(0)
  ou_dist_std <- if(length(test_output$ou_dist) > 0) scale(test_output$ou_dist)[,1] else numeric(0)
  
  # Create results row
  result_row <- data.frame(
    Order = order,
    Trait = trait,
    Sample_Size = length(pic_values),
    AD_Statistic_BM = test_output$ad_statistic_bm,
    AD_P_Value_BM = test_output$ad_p.value_bm,
    AD_Statistic_OU = test_output$ad_statistic_ou,
    AD_P_Value_OU = test_output$ad_p.value_ou,
    KS_Statistic_BM = test_output$ks_statistic_bm,
    KS_P_Value_BM = test_output$ks_p.value_bm,
    KS_Statistic_OU = test_output$ks_statistic_ou,
    KS_P_Value_OU = test_output$ks_p.value_ou,
    stringsAsFactors = FALSE
  )
  
  # Create and save plot if we have enough data
  if(length(bm_dist_std) > 30 && length(ou_dist_std) > 30) {
    plot_data <- list(
      ad_p.value_bm = result_row$AD_P_Value_BM,
      ad_p.value_ou = result_row$AD_P_Value_OU,
      ks_p.value_bm = result_row$KS_P_Value_BM,
      ks_p.value_ou = result_row$KS_P_Value_OU
    )
    
    p <- create_distribution_plot(order, trait, pic_values_std, plot_data, bm_dist_std, ou_dist_std)
    
    if(!is.null(p)) {
      trait_name_clean <- gsub(" ", "_", trait_display[[trait]])
      plot_filename <- file.path(plot_dir, paste0(order, "_", trait_name_clean, "_distribution_comparison.png"))
      ggsave(plot_filename, plot = p, width = 10, height = 8, dpi = 300)
      cat("  Plot saved:", plot_filename, "\n")
    }
  } else {
    cat("  Insufficient simulated data for plotting (BM:", length(bm_dist_std), "OU:", length(ou_dist_std), ")\n")
  }
  
  return(result_row)
}

# Process each order
for (order in orders) {
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("PROCESSING ORDER:", order, "\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  # Load data
  data <- tryCatch({
    cat("Loading data from Excel...\n")
    as.data.frame(read_excel(file_path, sheet = order))
  }, error = function(e) {
    cat("Error loading data for", order, ":", e$message, "\n")
    # List available sheets for debugging
    sheets <- excel_sheets(file_path)
    cat("Available sheets:", paste(sheets, collapse = ", "), "\n")
    next
  })
  
  if (is.null(data) || nrow(data) == 0) {
    cat("No data loaded for", order, "\n")
    next
  }
  
  cat("Data loaded with", nrow(data), "rows and", ncol(data), "columns\n")
  cat("Column names:", paste(colnames(data), collapse = ", "), "\n")
  
  # Data processing
  if ("Species" %in% colnames(data)) {
    data$Species <- trimws(as.character(data$Species))
    data$Species <- gsub("[[:space:]]+", "_", data$Species)
    data$Species <- gsub("[^[:alnum:]_]", "", data$Species)
    valid_rows <- !is.na(data$Species) & data$Species != ""
    cat("Valid species rows after cleaning:", sum(valid_rows), "\n")
    
    if (sum(valid_rows) == 0) {
      cat("Error: No valid Species entries after cleaning\n")
      next
    }
    data <- data[valid_rows, ]
    data$species_lower <- tolower(data$Species)
  } else {
    cat("Error: Missing 'Species' column in data\n")
    next
  }
  
  # Get tree
  tree <- trees[[order]]
  if (is.null(tree)) {
    cat("Error: Tree is NULL or failed to load\n")
    next
  }
  
  cat("Tree has", length(tree$tip.label), "tips\n")
  
  # Match species between data and tree
  tree_tip_lower <- tolower(tree$tip.label)
  data_tip_lower <- data$species_lower
  match_indices <- match(data_tip_lower, tree_tip_lower)
  common_indices_data <- !is.na(match_indices)
  common_species <- data$Species[common_indices_data]
  common_indices_tree <- tree$tip.label %in% common_species
  common_species <- intersect(common_species, tree$tip.label[common_indices_tree])
  
  cat("Common species between data and tree:", length(common_species), "\n")
  
  if (length(common_species) < 15) {
    cat("Skipping: insufficient common species\n")
    next
  }
  
  # Subset data and tree to common species
  data <- data[common_indices_data, ]
  tree <- keep.tip(tree, common_species)
  
  cat("Final data has", nrow(data), "rows, tree has", length(tree$tip.label), "tips\n")
  
  if (nrow(data) != length(tree$tip.label) || !all(data$Species %in% tree$tip.label)) {
    cat("Error: Mismatch in data and tree species after subsetting\n")
    next
  }
  
  # Analyze each trait
  for (trait in traits) {
    if (trait %in% colnames(data)) {
      cat("\n--- Analyzing", trait, "---\n")
      result <- perform_distribution_analysis(order, tree, data, trait)
      if(!is.null(result)) {
        test_results <- rbind(test_results, result)
        cat("Successfully analyzed", trait, "\n")
      } else {
        cat("Analysis failed for", trait, "\n")
      }
    } else {
      cat("Trait", trait, "not found in data columns\n")
    }
  }
}

# Save results to Excel
if (nrow(test_results) > 0) {
  excel_path <- "D:/Akash Ajay/insect bite force project/sravya work/GC skew/genomic_traits_distribution_test_results.xlsx"
  write.xlsx(test_results, excel_path, rowNames = FALSE)
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("ANALYSIS COMPLETED SUCCESSFULLY!\n")
  cat("Results saved to:", excel_path, "\n")
  cat("Total results:", nrow(test_results), "rows\n")
  cat("All plots saved to:", plot_dir, "\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  # Print summary
  cat("\n===== DISTRIBUTION TEST SUMMARY =====\n")
  for (order in orders) {
    order_data <- test_results %>% filter(Order == order)
    if (nrow(order_data) > 0) {
      cat("\n", order, ":\n")
      cat("  Traits analyzed:", nrow(order_data), "\n")
      cat("  Significant deviation from normality (AD p < 0.05):", sum(order_data$AD_P_Value_BM < 0.05, na.rm = TRUE), "\n")
      cat("  Significant difference from BM (KS p < 0.05):", sum(order_data$KS_P_Value_BM < 0.05, na.rm = TRUE), "\n")
      cat("  Significant difference from OU (KS p < 0.05):", sum(order_data$KS_P_Value_OU < 0.05, na.rm = TRUE), "\n")
    } else {
      cat("\n", order, ": No results\n")
    }
  }
} else {
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("ANALYSIS COMPLETED BUT NO RESULTS GENERATED!\n")
  cat("Please check the error messages above.\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
}

sink()

cat("\nR SCRIPT EXECUTION COMPLETED!\n")
cat("Check the output file for detailed logs: ks_ad_test_output.txt\n")
cat("Results Excel file: genomic_traits_distribution_test_results.xlsx\n")
cat("Plots directory: distribution_plots/\n")

# Show any warnings
if(length(warnings()) > 0) {
  cat("\nWarnings encountered:\n")
  print(warnings())
} else {
  cat("\nNo warnings encountered.\n")
}